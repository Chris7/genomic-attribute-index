//! GNI, the GFF Name Index.
//!
//! GNI is a small, explicit binary index for configured GFF3 attribute values.
//! It is intentionally not a feature-identity index: `ID` is just another
//! attribute, and only names supplied in [`NameIndexOptions`] are searchable.
//! Coordinates in the on-disk format are zero-based, half-open `start +
//! length` tuples.  TBI/CSI remains responsible for locating source records;
//! GNI stores no BGZF virtual offsets.

#[cfg(test)]
use std::collections::BTreeSet;
use std::{
    cmp::Reverse,
    collections::{BTreeMap, BinaryHeap, HashMap, HashSet},
    convert::TryFrom,
    fs::{self, File},
    io::{self, BufRead, BufReader, Cursor, Read, Seek, Write},
    num::NonZeroUsize,
    path::{Path, PathBuf},
    sync::{
        Arc,
        atomic::{AtomicU64, Ordering},
    },
    time::{Duration, Instant},
};

use bgzf::io::Seek as _;
use memmap2::Mmap;
use noodles::{
    bgzf,
    core::{Position, Region},
    csi::{
        self, BinningIndex,
        binning_index::index::{
            header::{Format, format::CoordinateSystem},
            reference_sequence::bin::Chunk,
        },
    },
    gff, tabix,
};
use sha2::{Digest, Sha256};

const MAGIC: [u8; 4] = *b"GNI\x01";
const MAJOR_VERSION: u16 = 1;
const MINOR_VERSION: u16 = 0;
const BYTE_ORDER_LITTLE: u8 = 1;
const COORDINATE_ZERO_BASED_HALF_OPEN: u8 = 1;
const NORMALIZATION_ASCII_LOWER: u8 = 0;
const NORMALIZATION_CASE_SENSITIVE: u8 = 1;
const HEADER_SIZE: usize = 256;
const DIRECTORY_ENTRY_SIZE: usize = 40;
const POSTINGS_DIRECTORY_ENTRY_SIZE: usize = 32;
const SPAN_DIRECTORY_ENTRY_SIZE: usize = 52;
const MAX_SECTION_BYTES: u64 = 1 << 40;
const MAX_BLOCK_BYTES: u64 = 256 << 20;
const DEFAULT_POSTINGS_BLOCK_TARGET: usize = 64 * 1024;
const DEFAULT_SPANS_PER_BLOCK: usize = 4096;
const MAX_RUN_FANIN: usize = 64;
const PROGRESS_RECORD_INTERVAL: u64 = 250_000;
const RUN_MAGIC: [u8; 8] = *b"GNIR\x01\x00\x00\x00";
static RUN_COUNTER: AtomicU64 = AtomicU64::new(0);

/// Errors returned by GNI construction, reading, and indexed querying.
#[derive(Debug, thiserror::Error)]
pub enum Error {
    /// A source or index file could not be read.
    #[error("I/O error: {0}")]
    Io(#[from] io::Error),
    /// The supplied input is malformed or inconsistent.
    #[error("invalid input: {0}")]
    InvalidInput(String),
    /// The GNI byte stream is malformed or unsafe to decode.
    #[error("corrupt GNI: {0}")]
    Corrupt(String),
    /// The source GFF or coordinate index does not match the GNI metadata.
    #[error("stale GNI: {0}")]
    Stale(String),
    /// A compressed block could not be decoded.
    #[error("compression error: {0}")]
    Compression(String),
    /// A value cannot be represented in the requested coordinate type.
    #[error("invalid coordinate")]
    InvalidCoordinate,
}

/// Options controlling which GFF3 attributes are indexed.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct NameIndexOptions {
    /// Attribute tags whose values should become searchable.  The list must
    /// be nonempty; exact duplicate tags are retained only once.
    pub attributes: Vec<String>,
    /// If true, trim values without applying ASCII lowercasing.
    pub case_sensitive: bool,
}

impl NameIndexOptions {
    /// Creates options and validates the explicit attribute list.
    pub fn new<I, S>(attributes: I, case_sensitive: bool) -> Result<Self>
    where
        I: IntoIterator<Item = S>,
        S: Into<String>,
    {
        let mut unique = Vec::new();
        let mut seen = HashSet::new();
        for attribute in attributes {
            let attribute = attribute.into();
            if attribute.is_empty() {
                return Err(Error::InvalidInput(
                    "attribute names must not be empty".to_string(),
                ));
            }
            if seen.insert(attribute.clone()) {
                unique.push(attribute);
            }
        }
        if unique.is_empty() {
            return Err(Error::InvalidInput(
                "at least one --attribute is required".to_string(),
            ));
        }
        Ok(Self {
            attributes: unique,
            case_sensitive,
        })
    }
}

/// Statistics collected while building a GNI.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct IndexStats {
    /// Number of non-comment GFF records read.
    pub records_processed: u64,
    /// Number of records containing at least one configured nonempty value.
    pub records_indexed: u64,
    /// Number of distinct normalized terms.
    pub distinct_terms: u64,
    /// Number of globally unique `(reference, start, length)` tuples.
    pub unique_spans: u64,
    /// Number of unique term-to-span relationships.
    pub postings: u64,
    /// Duplicate term-to-span relationships removed.
    pub duplicate_postings_removed: u64,
    /// Duplicate span observations removed.
    pub duplicate_spans_removed: u64,
    /// Final index size in bytes.
    pub index_bytes: u64,
    /// Postings bytes before block compression.
    pub postings_bytes_before_compression: u64,
    /// Postings bytes after block compression.
    pub postings_bytes_after_compression: u64,
    /// Span bytes represented by fixed-width triples.
    pub span_bytes_fixed_width: u64,
    /// Span bytes after structural integer encoding.
    pub span_bytes_structural: u64,
    /// Span bytes after optional block compression.
    pub span_bytes_after_compression: u64,
    /// Final index bytes divided by the number of indexed terms.
    pub bytes_per_term: f64,
    /// Final index bytes divided by the number of unique postings.
    pub bytes_per_posting: f64,
    /// Final index bytes divided by the number of unique spans.
    pub bytes_per_unique_span: f64,
    /// Timing for each major build phase.
    pub timings: BuildTimings,
    /// Peak size of the bounded scan working set, excluding output sections.
    pub peak_working_set_bytes: u64,
}

/// A major phase reported by a GNI build.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum BuildPhase {
    /// Streaming source parse and extraction.
    Scan,
    /// Sorting and writing a bounded-memory run.
    Spill,
    /// Deterministic k-way merge of sorted runs.
    Merge,
    /// Posting dictionary and block encoding.
    EncodePostings,
    /// Reference-local span block encoding.
    EncodeSpans,
    /// Final section serialization and atomic replacement.
    Serialize,
    /// Build completed successfully.
    Complete,
}

/// A monotonic progress snapshot emitted by a build callback.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct BuildProgress {
    /// Current build phase.
    pub phase: BuildPhase,
    /// Number of source records processed so far.
    pub records_processed: u64,
    /// Exact source bytes consumed by the hashing reader so far.
    pub bytes_read: u64,
    /// Elapsed wall time since the build began.
    pub elapsed: Duration,
}

/// Wall-clock durations for major build phases.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct BuildTimings {
    /// Streaming source scan duration.
    pub scan: Duration,
    /// Spill sorting and run-writing duration.
    pub spill: Duration,
    /// Run merge duration.
    pub merge: Duration,
    /// Posting encoding and compression duration.
    pub encode_postings: Duration,
    /// Span encoding and compression duration.
    pub encode_spans: Duration,
    /// Final serialization and atomic write duration.
    pub serialize: Duration,
    /// Total build duration.
    pub total: Duration,
}

/// Callback invoked with monotonic build progress snapshots.
pub type ProgressCallback = Arc<dyn Fn(BuildProgress) + Send + Sync + 'static>;

/// Resource and observability controls for a GNI build.
#[derive(Clone)]
pub struct BuildOptions {
    /// Approximate maximum scan/run working set before a sorted run is spilled.
    pub memory_budget_bytes: usize,
    /// Number of workers used for independent block compression.
    pub compression_threads: usize,
    /// Number of BGZF decompression workers. One keeps the reader single-threaded.
    pub bgzf_threads: usize,
    /// Optional callback; the library never writes progress to stderr itself.
    pub progress: Option<ProgressCallback>,
}

impl Default for BuildOptions {
    fn default() -> Self {
        let workers = std::thread::available_parallelism()
            .map(NonZeroUsize::get)
            .unwrap_or(1)
            .min(8);
        Self {
            memory_budget_bytes: 64 * 1024 * 1024,
            compression_threads: workers,
            bgzf_threads: workers,
            progress: None,
        }
    }
}

impl BuildOptions {
    /// Sets the approximate bounded scan memory budget.
    pub fn with_memory_budget(mut self, bytes: usize) -> Self {
        self.memory_budget_bytes = bytes;
        self
    }

    /// Sets the number of deterministic block-compression workers.
    pub fn with_compression_threads(mut self, threads: usize) -> Self {
        self.compression_threads = threads;
        self
    }

    /// Sets the BGZF decompression worker count.
    pub fn with_bgzf_threads(mut self, threads: usize) -> Self {
        self.bgzf_threads = threads;
        self
    }

    /// Sets a progress callback shared by the build operation.
    pub fn with_progress<F>(mut self, callback: F) -> Self
    where
        F: Fn(BuildProgress) + Send + Sync + 'static,
    {
        self.progress = Some(Arc::new(callback));
        self
    }
}

/// A parsed GFF3 record returned by [`IndexedGff::query_name`].
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct GffRecord {
    /// Reference sequence name.
    pub reference_sequence_name: String,
    /// Source column.
    pub source: String,
    /// Type column.
    pub ty: String,
    /// One-based inclusive start as written in GFF3.
    pub start: u64,
    /// One-based inclusive end as written in GFF3.
    pub end: u64,
    /// Score column.
    pub score: String,
    /// Strand column.
    pub strand: String,
    /// Phase column.
    pub phase: String,
    /// Decoded attribute pairs in source order.  Repeated tags are retained.
    pub attributes: Vec<(String, Vec<String>)>,
    /// Original record text without the trailing line terminator.
    pub raw_line: String,
}

impl GffRecord {
    /// Returns all decoded values for an attribute tag.
    pub fn attribute_values(&self, tag: &str) -> impl Iterator<Item = &str> {
        self.attributes
            .iter()
            .filter(move |(name, _)| name == tag)
            .flat_map(|(_, values)| values.iter().map(String::as_str))
    }
}

/// A zero-based, half-open GNI span.
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub struct Span {
    /// Reference dictionary identifier.
    pub reference_id: u32,
    /// Zero-based start.
    pub start: u64,
    /// Number of bases; the end is `start + length`.
    pub length: u64,
}

impl Span {
    /// Returns the checked zero-based half-open end.
    pub fn end(self) -> Result<u64> {
        self.start
            .checked_add(self.length)
            .ok_or(Error::InvalidCoordinate)
    }

    /// Converts one-based inclusive GFF coordinates at the parser boundary.
    pub fn from_gff(reference_id: u32, gff_start: u64, gff_end: u64) -> Result<Self> {
        let (start, length) = gff_to_span(gff_start, gff_end)?;
        Ok(Self {
            reference_id,
            start,
            length,
        })
    }
}

/// Converts GFF3's one-based inclusive interval to a zero-based half-open
/// `start + length` pair.
pub fn gff_to_span(gff_start: u64, gff_end: u64) -> Result<(u64, u64)> {
    if gff_start == 0 || gff_end < gff_start {
        return Err(Error::InvalidCoordinate);
    }
    let start = gff_start.checked_sub(1).ok_or(Error::InvalidCoordinate)?;
    let length = gff_end.checked_sub(start).ok_or(Error::InvalidCoordinate)?;
    if length == 0 {
        return Err(Error::InvalidCoordinate);
    }
    start.checked_add(length).ok_or(Error::InvalidCoordinate)?;
    Ok((start, length))
}

/// Metadata reported by [`NameIndexReader::inspect`].
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct IndexMetadata {
    /// Major format version.
    pub major_version: u16,
    /// Minor format version.
    pub minor_version: u16,
    /// Whether values use case-sensitive normalization.
    pub case_sensitive: bool,
    /// Explicit configured tags.
    pub attributes: Vec<String>,
    /// SHA-256 of the source GFF bytes.
    pub gff_fingerprint: [u8; 32],
    /// SHA-256 of the TBI/CSI bytes.
    pub coordinate_index_fingerprint: [u8; 32],
    /// SHA-256 of the reference dictionary.
    pub reference_dictionary_fingerprint: [u8; 32],
    /// Number of terms.
    pub term_count: u64,
    /// Number of spans.
    pub unique_span_count: u64,
    /// Number of term-to-span postings.
    pub posting_count: u64,
    /// Number of postings blocks.
    pub postings_block_count: u64,
    /// Number of span blocks.
    pub span_block_count: u64,
    /// Number of reference sequences in the coordinate-index dictionary.
    pub reference_count: u32,
    /// Configured target number of rows per reference-specific span block.
    pub span_block_size: u32,
    /// On-disk byte length.
    pub file_size: u64,
    /// Attribute section size.
    pub attribute_section_bytes: u64,
    /// FST term dictionary size.
    pub term_dictionary_bytes: u64,
    /// Fixed postings-directory size.
    pub postings_directory_bytes: u64,
    /// Sum of uncompressed postings block lengths.
    pub postings_uncompressed_bytes: u64,
    /// Compressed postings-data size.
    pub postings_data_bytes: u64,
    /// Fixed span-directory size.
    pub span_directory_bytes: u64,
    /// Sum of uncompressed span block lengths.
    pub span_uncompressed_bytes: u64,
    /// Compressed span-data size.
    pub span_data_bytes: u64,
    /// Number of postings blocks using zstd.
    pub compressed_postings_blocks: u64,
    /// Number of span blocks using zstd.
    pub compressed_span_blocks: u64,
    /// Number of span blocks using delta-varint starts.
    pub delta_start_blocks: u64,
    /// Number of span blocks using frame-of-reference starts.
    pub for_start_blocks: u64,
    /// Number of span blocks using ordinary-varint lengths.
    pub varint_length_blocks: u64,
    /// Number of span blocks using frame-of-reference lengths.
    pub for_length_blocks: u64,
}

/// Bytes decompressed while resolving one term. This is useful for measuring
/// query amplification without exposing block offsets as part of the format.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct LookupStats {
    /// Uncompressed postings bytes read.
    pub postings_bytes_decompressed: u64,
    /// Unique uncompressed span-block bytes read.
    pub span_bytes_decompressed: u64,
}

/// Instrumentation for one indexed name query.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct QueryStats {
    /// Number of exact GNI spans requested by the term posting.
    pub requested_spans: u64,
    /// Number of distinct GNI span blocks decoded while resolving the posting.
    pub distinct_span_blocks_decoded: u64,
    /// Number of exact coordinate-index interval queries issued.
    pub exact_interval_queries: u64,
    /// Number of raw BGZF chunks returned by those interval queries.
    pub raw_chunks: u64,
    /// Number of merged, non-overlapping BGZF chunks read from the source.
    pub merged_chunks: u64,
    /// Number of unique source records parsed, keyed by BGZF virtual position.
    pub unique_candidate_records: u64,
    /// Number of source records that matched the exact span and configured term.
    pub matching_records: u64,
    /// Uncompressed bytes returned while reading merged chunks.
    pub bytes_read: u64,
}

/// A result alias using the crate's corruption-safe error type.
pub type Result<T> = std::result::Result<T, Error>;

#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
struct SpanKey {
    reference_id: u32,
    start: u64,
    length: u64,
}

impl From<SpanKey> for Span {
    fn from(value: SpanKey) -> Self {
        Self {
            reference_id: value.reference_id,
            start: value.start,
            length: value.length,
        }
    }
}

#[derive(Clone, Debug)]
struct ParsedRecord {
    record: GffRecord,
}

#[derive(Clone, Debug)]
struct CoordinateDictionary {
    names: Vec<String>,
    fingerprint: [u8; 32],
}

enum CoordinateIndex {
    Tabix(tabix::Index),
    Csi(csi::Index),
}

impl CoordinateIndex {
    fn dictionary(&self) -> Result<CoordinateDictionary> {
        let format = match self {
            Self::Tabix(index) => index
                .header()
                .ok_or_else(|| Error::InvalidInput("coordinate index has no header".into()))?
                .format(),
            Self::Csi(index) => index
                .header()
                .ok_or_else(|| Error::InvalidInput("coordinate index has no header".into()))?
                .format(),
        };
        validate_coordinate_index_format(format)?;
        let names = match self {
            Self::Tabix(index) => index
                .header()
                .ok_or_else(|| Error::InvalidInput("coordinate index has no header".into()))?
                .reference_sequence_names()
                .iter()
                .map(|name| {
                    String::from_utf8(<_ as AsRef<[u8]>>::as_ref(name).to_vec()).map_err(|_| {
                        Error::InvalidInput("reference names must be valid UTF-8".into())
                    })
                })
                .collect::<Result<Vec<_>>>()?,
            Self::Csi(index) => index
                .header()
                .ok_or_else(|| Error::InvalidInput("coordinate index has no header".into()))?
                .reference_sequence_names()
                .iter()
                .map(|name| {
                    String::from_utf8(<_ as AsRef<[u8]>>::as_ref(name).to_vec()).map_err(|_| {
                        Error::InvalidInput("reference names must be valid UTF-8".into())
                    })
                })
                .collect::<Result<Vec<_>>>()?,
        };
        if names.iter().any(String::is_empty) {
            return Err(Error::InvalidInput(
                "coordinate index contains an empty reference name".into(),
            ));
        }
        Ok(CoordinateDictionary {
            fingerprint: fingerprint_reference_dictionary(&names),
            names,
        })
    }
}

fn fingerprint_reference_dictionary(names: &[String]) -> [u8; 32] {
    let mut hasher = Sha256::new();
    for name in names {
        hasher.update((name.len() as u64).to_le_bytes());
        hasher.update(name.as_bytes());
    }
    hasher.finalize().into()
}

fn fingerprint_file(path: &Path) -> Result<[u8; 32]> {
    let mut reader = File::open(path)?;
    let mut hasher = Sha256::new();
    let mut buffer = [0_u8; 128 * 1024];
    loop {
        let length = reader.read(&mut buffer)?;
        if length == 0 {
            break;
        }
        hasher.update(&buffer[..length]);
    }
    Ok(hasher.finalize().into())
}

/// A source reader that hashes exactly the bytes returned by `Read` while
/// exposing a shared byte counter for progress reporting. For BGZF input the
/// bytes are compressed source bytes; the BGZF decoder consumes this reader
/// without changing the fingerprint semantics.
struct HashingReader<R> {
    inner: R,
    hasher: Sha256,
    bytes_read: Arc<AtomicU64>,
}

impl<R> HashingReader<R> {
    fn with_counter(inner: R, bytes_read: Arc<AtomicU64>) -> Self {
        Self {
            inner,
            hasher: Sha256::new(),
            bytes_read,
        }
    }

    fn finish(self) -> ([u8; 32], u64) {
        (
            self.hasher.finalize().into(),
            self.bytes_read.load(Ordering::Relaxed),
        )
    }
}

impl<R> Read for HashingReader<R>
where
    R: Read,
{
    fn read(&mut self, buffer: &mut [u8]) -> io::Result<usize> {
        let length = self.inner.read(buffer)?;
        if length != 0 {
            self.hasher.update(&buffer[..length]);
            self.bytes_read.fetch_add(length as u64, Ordering::Relaxed);
        }
        Ok(length)
    }
}

trait HashingInput: BufRead {
    fn drain_and_finish(self) -> io::Result<([u8; 32], u64)>;
}

impl HashingInput for BufReader<HashingReader<File>> {
    fn drain_and_finish(mut self) -> io::Result<([u8; 32], u64)> {
        io::copy(&mut self, &mut io::sink())?;
        Ok(self.into_inner().finish())
    }
}

impl HashingInput for bgzf::io::Reader<HashingReader<File>> {
    fn drain_and_finish(mut self) -> io::Result<([u8; 32], u64)> {
        io::copy(&mut self, &mut io::sink())?;
        Ok(self.into_inner().finish())
    }
}

impl HashingInput for bgzf::io::MultithreadedReader<HashingReader<File>> {
    fn drain_and_finish(mut self) -> io::Result<([u8; 32], u64)> {
        let drain_result = io::copy(&mut self, &mut io::sink());
        let inner = self.finish()?;
        drain_result?;
        Ok(inner.finish())
    }
}

fn scan_gff_reader<R, F>(mut reader: gff::io::Reader<R>, mut process: F) -> Result<([u8; 32], u64)>
where
    R: HashingInput,
    F: FnMut(gff::feature::RecordBuf) -> Result<()>,
{
    let mut records = reader.record_bufs();
    for result in records.by_ref() {
        let record =
            result.map_err(|error| Error::InvalidInput(format!("invalid GFF record: {error}")))?;
        process(record)?;
    }
    drop(records);
    reader.into_inner().drain_and_finish().map_err(Error::from)
}

fn read_coordinate_index_with_fingerprint(path: &Path) -> Result<(CoordinateIndex, [u8; 32])> {
    // Parse from one in-memory byte copy so the exact coordinate-index
    // fingerprint and the noodles reader share a single file read. The index
    // itself is retained by noodles in owned structures, not by this byte
    // buffer after parsing.
    let bytes = fs::read(path)?;
    let fingerprint = Sha256::digest(&bytes).into();
    // Explicit --coordinate-index paths may have arbitrary names. Try both
    // compressed readers by content and retain useful diagnostics if neither
    // accepts the bytes.
    let tabix_error = match tabix::io::Reader::new(Cursor::new(bytes.clone())).read_index() {
        Ok(index) => return Ok((CoordinateIndex::Tabix(index), fingerprint)),
        Err(error) => error.to_string(),
    };
    let csi_error = match csi::io::Reader::new(Cursor::new(bytes)).read_index() {
        Ok(index) => return Ok((CoordinateIndex::Csi(index), fingerprint)),
        Err(error) => error.to_string(),
    };
    Err(Error::InvalidInput(format!(
        "coordinate index {} is neither TBI nor CSI (TBI: {tabix_error}; CSI: {csi_error})",
        path.display()
    )))
}

fn validate_coordinate_index_format(format: Format) -> Result<()> {
    if format != Format::Generic(CoordinateSystem::Gff) {
        return Err(Error::InvalidInput(
            "coordinate index is not a generic GFF coordinate index".into(),
        ));
    }
    Ok(())
}

/// Applies the GNI normalization policy: trim Unicode whitespace, then
/// lowercase ASCII letters unless `case_sensitive` is true.
pub fn normalize_value(value: &str, case_sensitive: bool) -> String {
    let trimmed = value.trim_matches(char::is_whitespace);
    if case_sensitive {
        trimmed.to_owned()
    } else {
        trimmed
            .chars()
            .map(|character| character.to_ascii_lowercase())
            .collect()
    }
}

fn parsed_record_from_record_buf(
    record: &gff::feature::RecordBuf,
    raw_line: String,
    reference_ids: &HashMap<String, u32>,
) -> Result<ParsedRecord> {
    let reference_sequence_name = String::from_utf8(record.reference_sequence_name().to_vec())
        .map_err(|_| Error::InvalidInput("reference sequence name is not UTF-8".into()))?;
    if !reference_ids.contains_key(&reference_sequence_name) {
        return Err(Error::InvalidInput(format!(
            "GFF reference sequence {reference_sequence_name:?} is absent from coordinate index"
        )));
    }
    let start = record.start().get() as u64;
    let end = record.end().get() as u64;
    gff_to_span(start, end)?;
    let source = String::from_utf8(record.source().to_vec())
        .map_err(|_| Error::InvalidInput("GFF source is not UTF-8".into()))?;
    let ty = String::from_utf8(record.ty().to_vec())
        .map_err(|_| Error::InvalidInput("GFF type is not UTF-8".into()))?;
    let attributes = record
        .attributes()
        .as_ref()
        .iter()
        .map(|(tag, value)| {
            let tag = String::from_utf8(<_ as AsRef<[u8]>>::as_ref(tag).to_vec())
                .map_err(|_| Error::InvalidInput("GFF attribute tag is not UTF-8".into()))?;
            let values = value
                .iter()
                .map(|value| {
                    String::from_utf8(<_ as AsRef<[u8]>>::as_ref(value).to_vec())
                        .map_err(|_| Error::InvalidInput("GFF attribute value is not UTF-8".into()))
                })
                .collect::<Result<Vec<_>>>()?;
            Ok((tag, values))
        })
        .collect::<Result<Vec<_>>>()?;
    Ok(ParsedRecord {
        record: GffRecord {
            reference_sequence_name,
            source,
            ty,
            start,
            end,
            score: record
                .score()
                .map(|score| score.to_string())
                .unwrap_or_else(|| ".".to_string()),
            strand: match record.strand() {
                gff::feature::record::Strand::None => ".".to_string(),
                gff::feature::record::Strand::Forward => "+".to_string(),
                gff::feature::record::Strand::Reverse => "-".to_string(),
                gff::feature::record::Strand::Unknown => "?".to_string(),
            },
            phase: record
                .phase()
                .map(|phase| match phase {
                    gff::feature::record::Phase::Zero => "0".to_string(),
                    gff::feature::record::Phase::One => "1".to_string(),
                    gff::feature::record::Phase::Two => "2".to_string(),
                })
                .unwrap_or_else(|| ".".to_string()),
            attributes,
            raw_line,
        },
    })
}

fn feature_record_span<R>(record: &R, reference_ids: &HashMap<String, u32>) -> Result<SpanKey>
where
    R: gff::feature::Record + ?Sized,
{
    let reference_sequence_name = String::from_utf8(record.reference_sequence_name().to_vec())
        .map_err(|_| Error::InvalidInput("GFF reference sequence name is not UTF-8".into()))?;
    let reference_id = *reference_ids.get(&reference_sequence_name).ok_or_else(|| {
        Error::InvalidInput(format!(
            "GFF reference sequence {reference_sequence_name:?} is absent from coordinate index"
        ))
    })?;
    let start = record
        .feature_start()
        .map_err(|error| Error::InvalidInput(format!("invalid GFF start: {error}")))?
        .get() as u64;
    let end = record
        .feature_end()
        .map_err(|error| Error::InvalidInput(format!("invalid GFF end: {error}")))?
        .get() as u64;
    let (start, length) = gff_to_span(start, end)?;
    Ok(SpanKey {
        reference_id,
        start,
        length,
    })
}

fn feature_record_matches_term<R>(
    record: &R,
    configured_attributes: &HashSet<String>,
    case_sensitive: bool,
    normalized_query: &str,
) -> Result<bool>
where
    R: gff::feature::Record + ?Sized,
{
    for result in record.attributes().iter() {
        let (tag, value) = result
            .map_err(|error| Error::InvalidInput(format!("invalid GFF attributes: {error}")))?;
        let tag = std::str::from_utf8(tag.as_ref())
            .map_err(|_| Error::InvalidInput("GFF attribute tag is not UTF-8".into()))?;
        if !configured_attributes.contains(tag) {
            continue;
        }
        for value in value.iter() {
            let value = value.map_err(|error| {
                Error::InvalidInput(format!("invalid GFF attribute value: {error}"))
            })?;
            let value = std::str::from_utf8(value.as_ref())
                .map_err(|_| Error::InvalidInput("GFF attribute value is not UTF-8".into()))?;
            if normalize_value(value, case_sensitive) == normalized_query {
                return Ok(true);
            }
        }
    }
    Ok(false)
}

fn parsed_record_from_feature_record<R>(
    record: &R,
    raw_line: String,
    reference_ids: &HashMap<String, u32>,
) -> Result<ParsedRecord>
where
    R: gff::feature::Record + ?Sized,
{
    let record = gff::feature::RecordBuf::try_from_feature_record(record)
        .map_err(|error| Error::InvalidInput(format!("invalid GFF record: {error}")))?;
    parsed_record_from_record_buf(&record, raw_line, reference_ids)
}

#[derive(Clone, Debug)]
struct ExtractedRecord {
    span: SpanKey,
    terms: Vec<String>,
}

/// Extracts only the coordinate and configured values needed by the builder.
/// The lossless [`GffRecord`] conversion above remains query-only; this path
/// intentionally does not clone source, score, strand, phase, or raw text.
fn extract_record(
    record: &gff::feature::RecordBuf,
    reference_ids: &HashMap<String, u32>,
    configured: &HashSet<&str>,
    case_sensitive: bool,
) -> Result<ExtractedRecord> {
    let reference_sequence_name = String::from_utf8(record.reference_sequence_name().to_vec())
        .map_err(|_| Error::InvalidInput("reference sequence name is not UTF-8".into()))?;
    let reference_id = *reference_ids.get(&reference_sequence_name).ok_or_else(|| {
        Error::InvalidInput(format!(
            "GFF reference sequence {reference_sequence_name:?} is absent from coordinate index"
        ))
    })?;
    let (start, length) = gff_to_span(record.start().get() as u64, record.end().get() as u64)?;
    let mut terms = Vec::new();
    for (tag, value) in record.attributes().as_ref() {
        let tag = std::str::from_utf8(tag.as_ref())
            .map_err(|_| Error::InvalidInput("GFF attribute tag is not UTF-8".into()))?;
        if !configured.contains(tag) {
            continue;
        }
        for value in value.iter() {
            let value = String::from_utf8(<_ as AsRef<[u8]>>::as_ref(value).to_vec())
                .map_err(|_| Error::InvalidInput("GFF attribute value is not UTF-8".into()))?;
            let normalized = normalize_value(&value, case_sensitive);
            if !normalized.is_empty() {
                terms.push(normalized);
            }
        }
    }
    Ok(ExtractedRecord {
        span: SpanKey {
            reference_id,
            start,
            length,
        },
        terms,
    })
}

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
enum SectionKind {
    Attributes = 1,
    Terms = 2,
    PostingsDirectory = 3,
    PostingsData = 4,
    SpansDirectory = 5,
    SpansData = 6,
}

impl TryFrom<u32> for SectionKind {
    type Error = Error;

    fn try_from(value: u32) -> Result<Self> {
        match value {
            1 => Ok(Self::Attributes),
            2 => Ok(Self::Terms),
            3 => Ok(Self::PostingsDirectory),
            4 => Ok(Self::PostingsData),
            5 => Ok(Self::SpansDirectory),
            6 => Ok(Self::SpansData),
            _ => Err(Error::Corrupt(format!("unknown section kind {value}"))),
        }
    }
}

#[derive(Clone, Debug)]
struct SectionDirectoryEntry {
    kind: SectionKind,
    flags: u32,
    offset: u64,
    length: u64,
    item_count: u64,
    checksum: u32,
}

#[derive(Clone, Debug)]
struct EncodedBlock {
    compressed: Vec<u8>,
    uncompressed_length: u32,
    checksum: u32,
    compression: u8,
}

#[derive(Clone, Debug, Eq, PartialEq)]
struct PostingDirectoryEntry {
    compressed_offset: u64,
    compressed_length: u32,
    uncompressed_length: u32,
    checksum: u32,
    compression: u8,
}

#[derive(Clone, Debug, Eq, PartialEq)]
struct SpanDirectoryEntry {
    first_span_id: u64,
    span_count: u32,
    reference_id: u32,
    first_start: u64,
    compressed_offset: u64,
    compressed_length: u32,
    uncompressed_length: u32,
    checksum: u32,
    start_encoding: u8,
    length_encoding: u8,
    compression: u8,
}

fn put_u32(buffer: &mut Vec<u8>, value: u32) {
    buffer.extend_from_slice(&value.to_le_bytes());
}

fn put_u64(buffer: &mut Vec<u8>, value: u64) {
    buffer.extend_from_slice(&value.to_le_bytes());
}

fn read_u8(bytes: &[u8], offset: &mut usize, context: &str) -> Result<u8> {
    let value = *bytes
        .get(*offset)
        .ok_or_else(|| Error::Corrupt(format!("truncated {context}")))?;
    *offset += 1;
    Ok(value)
}

fn read_u16(bytes: &[u8], offset: &mut usize, context: &str) -> Result<u16> {
    let end = offset
        .checked_add(2)
        .ok_or_else(|| Error::Corrupt(format!("{context} offset overflow")))?;
    let value = u16::from_le_bytes(
        bytes
            .get(*offset..end)
            .ok_or_else(|| Error::Corrupt(format!("truncated {context}")))?
            .try_into()
            .map_err(|_| Error::Corrupt(format!("invalid {context}")))?,
    );
    *offset = end;
    Ok(value)
}

fn read_u32(bytes: &[u8], offset: &mut usize, context: &str) -> Result<u32> {
    let end = offset
        .checked_add(4)
        .ok_or_else(|| Error::Corrupt(format!("{context} offset overflow")))?;
    let value = u32::from_le_bytes(
        bytes
            .get(*offset..end)
            .ok_or_else(|| Error::Corrupt(format!("truncated {context}")))?
            .try_into()
            .map_err(|_| Error::Corrupt(format!("invalid {context}")))?,
    );
    *offset = end;
    Ok(value)
}

fn read_u64(bytes: &[u8], offset: &mut usize, context: &str) -> Result<u64> {
    let end = offset
        .checked_add(8)
        .ok_or_else(|| Error::Corrupt(format!("{context} offset overflow")))?;
    let value = u64::from_le_bytes(
        bytes
            .get(*offset..end)
            .ok_or_else(|| Error::Corrupt(format!("truncated {context}")))?
            .try_into()
            .map_err(|_| Error::Corrupt(format!("invalid {context}")))?,
    );
    *offset = end;
    Ok(value)
}

fn read_array<const N: usize>(bytes: &[u8], offset: &mut usize, context: &str) -> Result<[u8; N]> {
    let end = offset
        .checked_add(N)
        .ok_or_else(|| Error::Corrupt(format!("{context} offset overflow")))?;
    let value = bytes
        .get(*offset..end)
        .ok_or_else(|| Error::Corrupt(format!("truncated {context}")))?
        .try_into()
        .map_err(|_| Error::Corrupt(format!("invalid {context}")))?;
    *offset = end;
    Ok(value)
}

fn write_varint(buffer: &mut Vec<u8>, mut value: u64) {
    while value >= 0x80 {
        buffer.push((value as u8) | 0x80);
        value >>= 7;
    }
    buffer.push(value as u8);
}

fn read_varint(bytes: &[u8], offset: &mut usize, context: &str) -> Result<u64> {
    let mut value = 0_u64;
    for shift in (0..64).step_by(7) {
        let byte = read_u8(bytes, offset, context)?;
        let payload = (byte & 0x7f) as u64;
        if shift == 63 && payload > 1 {
            return Err(Error::Corrupt(format!("overflowing varint in {context}")));
        }
        value |= payload
            .checked_shl(shift)
            .ok_or_else(|| Error::Corrupt(format!("overflowing varint in {context}")))?;
        if byte & 0x80 == 0 {
            return Ok(value);
        }
    }
    Err(Error::Corrupt(format!("unterminated varint in {context}")))
}

fn checksum(bytes: &[u8]) -> u32 {
    crc32fast::hash(bytes)
}

fn compress_block(bytes: &[u8]) -> Result<EncodedBlock> {
    let checksum = checksum(bytes);
    let compressed =
        zstd::bulk::compress(bytes, 3).map_err(|error| Error::Compression(error.to_string()))?;
    if compressed.len() < bytes.len() {
        Ok(EncodedBlock {
            compressed,
            uncompressed_length: u32::try_from(bytes.len())
                .map_err(|_| Error::InvalidInput("block exceeds 4 GiB".into()))?,
            checksum,
            compression: 1,
        })
    } else {
        Ok(EncodedBlock {
            compressed: bytes.to_vec(),
            uncompressed_length: u32::try_from(bytes.len())
                .map_err(|_| Error::InvalidInput("block exceeds 4 GiB".into()))?,
            checksum,
            compression: 0,
        })
    }
}

fn decompress_block(
    compressed: &[u8],
    compression: u8,
    uncompressed_length: u32,
    expected_checksum: u32,
    context: &str,
) -> Result<Vec<u8>> {
    if u64::from(uncompressed_length) > MAX_BLOCK_BYTES {
        return Err(Error::Corrupt(format!(
            "{context} uncompressed length exceeds safety limit"
        )));
    }
    let bytes = match compression {
        0 => {
            if compressed.len() != usize::try_from(uncompressed_length).unwrap_or(usize::MAX) {
                return Err(Error::Corrupt(format!(
                    "{context} uncompressed-size mismatch"
                )));
            }
            compressed.to_vec()
        }
        1 => zstd::bulk::decompress(compressed, uncompressed_length as usize)
            .map_err(|error| Error::Compression(format!("{context}: {error}")))?,
        other => {
            return Err(Error::Corrupt(format!(
                "{context} unknown compression {other}"
            )));
        }
    };
    if bytes.len() != uncompressed_length as usize {
        return Err(Error::Corrupt(format!(
            "{context} decompressed-size mismatch"
        )));
    }
    if checksum(&bytes) != expected_checksum {
        return Err(Error::Corrupt(format!("{context} checksum mismatch")));
    }
    Ok(bytes)
}

fn encode_attributes(attributes: &[String]) -> Result<Vec<u8>> {
    let mut bytes = Vec::new();
    put_u32(
        &mut bytes,
        u32::try_from(attributes.len())
            .map_err(|_| Error::InvalidInput("too many attributes".into()))?,
    );
    for attribute in attributes {
        let raw = attribute.as_bytes();
        put_u32(
            &mut bytes,
            u32::try_from(raw.len())
                .map_err(|_| Error::InvalidInput("attribute name exceeds 4 GiB".into()))?,
        );
        bytes.extend_from_slice(raw);
    }
    Ok(bytes)
}

fn decode_attributes(bytes: &[u8]) -> Result<Vec<String>> {
    let mut offset = 0;
    let count = read_u32(bytes, &mut offset, "attribute count")? as usize;
    if count > 1_000_000 {
        return Err(Error::Corrupt("excessive attribute count".into()));
    }
    let mut attributes = Vec::with_capacity(count);
    let mut seen = HashSet::with_capacity(count);
    for _ in 0..count {
        let length = read_u32(bytes, &mut offset, "attribute length")? as usize;
        let end = offset
            .checked_add(length)
            .ok_or_else(|| Error::Corrupt("attribute length overflow".into()))?;
        let name = String::from_utf8(
            bytes
                .get(offset..end)
                .ok_or_else(|| Error::Corrupt("truncated attribute name".into()))?
                .to_vec(),
        )
        .map_err(|_| Error::Corrupt("attribute name is not UTF-8".into()))?;
        if name.is_empty() || !seen.insert(name.clone()) {
            return Err(Error::Corrupt("invalid or duplicate attribute name".into()));
        }
        attributes.push(name);
        offset = end;
    }
    if offset != bytes.len() {
        return Err(Error::Corrupt("trailing bytes in attribute section".into()));
    }
    Ok(attributes)
}

type PostingEncoding = (Vec<u8>, Vec<PostingDirectoryEntry>, Vec<u8>, u64, u64, u64);

struct PostingEncoder {
    raw_blocks: Vec<Vec<u8>>,
    current: Vec<u8>,
    block_id: u32,
    locators: Vec<(String, u32, u32)>,
    postings_bytes_before_compression: u64,
    posting_count: u64,
}

impl PostingEncoder {
    fn new() -> Self {
        Self {
            raw_blocks: Vec::new(),
            current: Vec::new(),
            block_id: 0,
            locators: Vec::new(),
            postings_bytes_before_compression: 0,
            posting_count: 0,
        }
    }

    fn add(&mut self, term: String, spans: &[u64]) -> Result<()> {
        if spans.is_empty() {
            return Ok(());
        }
        let mut record = Vec::new();
        write_varint(
            &mut record,
            u64::try_from(spans.len())
                .map_err(|_| Error::InvalidInput("too many spans for a term".into()))?,
        );
        let mut previous = 0_u64;
        for (index, span_id) in spans.iter().copied().enumerate() {
            if index == 0 {
                write_varint(&mut record, span_id);
            } else {
                let delta = span_id
                    .checked_sub(previous)
                    .ok_or_else(|| Error::InvalidInput("span IDs are not sorted".into()))?;
                if delta == 0 {
                    return Err(Error::InvalidInput(
                        "posting spans are not deduplicated".into(),
                    ));
                }
                write_varint(&mut record, delta);
            }
            previous = span_id;
        }
        if !self.current.is_empty()
            && self
                .current
                .len()
                .checked_add(record.len())
                .is_some_and(|length| length > DEFAULT_POSTINGS_BLOCK_TARGET)
        {
            self.raw_blocks.push(std::mem::take(&mut self.current));
            self.block_id = self
                .block_id
                .checked_add(1)
                .ok_or(Error::InvalidCoordinate)?;
        }
        let record_offset = u32::try_from(self.current.len())
            .map_err(|_| Error::InvalidInput("posting offset exceeds 4 GiB".into()))?;
        self.current.extend_from_slice(&record);
        self.locators.push((term, self.block_id, record_offset));
        self.postings_bytes_before_compression = self
            .postings_bytes_before_compression
            .checked_add(record.len() as u64)
            .ok_or(Error::InvalidCoordinate)?;
        self.posting_count = self
            .posting_count
            .checked_add(spans.len() as u64)
            .ok_or(Error::InvalidCoordinate)?;
        if self.current.len() >= DEFAULT_POSTINGS_BLOCK_TARGET {
            self.raw_blocks.push(std::mem::take(&mut self.current));
            self.block_id = self
                .block_id
                .checked_add(1)
                .ok_or(Error::InvalidCoordinate)?;
        }
        Ok(())
    }

    fn finish(mut self, compression_threads: usize) -> Result<PostingEncoding> {
        if !self.current.is_empty() || self.locators.is_empty() {
            self.raw_blocks.push(self.current);
        }
        let encoded_blocks = compress_blocks_parallel(&self.raw_blocks, compression_threads)?;
        let mut posting_directory = Vec::with_capacity(encoded_blocks.len());
        let mut posting_data = Vec::new();
        for block in &encoded_blocks {
            let compressed_offset = u64::try_from(posting_data.len())
                .map_err(|_| Error::InvalidInput("postings section exceeds 4 GiB".into()))?;
            posting_data.extend_from_slice(&block.compressed);
            posting_directory.push(PostingDirectoryEntry {
                compressed_offset,
                compressed_length: u32::try_from(block.compressed.len())
                    .map_err(|_| Error::InvalidInput("compressed block exceeds 4 GiB".into()))?,
                uncompressed_length: block.uncompressed_length,
                checksum: block.checksum,
                compression: block.compression,
            });
        }
        let mut fst_bytes = Vec::new();
        let term_count = self.locators.len() as u64;
        {
            let mut builder = fst::MapBuilder::new(&mut fst_bytes).map_err(|error| {
                Error::InvalidInput(format!("could not create term FST: {error}"))
            })?;
            for (term, block_id, record_offset) in self.locators {
                let value = (u64::from(block_id) << 32) | u64::from(record_offset);
                builder.insert(term, value).map_err(|error| {
                    Error::InvalidInput(format!("could not build term FST: {error}"))
                })?;
            }
            builder.finish().map_err(|error| {
                Error::InvalidInput(format!("could not finish term FST: {error}"))
            })?;
        }
        Ok((
            fst_bytes,
            posting_directory,
            posting_data,
            self.postings_bytes_before_compression,
            self.posting_count,
            term_count,
        ))
    }
}

fn compress_blocks_parallel(blocks: &[Vec<u8>], thread_count: usize) -> Result<Vec<EncodedBlock>> {
    use rayon::prelude::*;

    let thread_count = thread_count.max(1);
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(thread_count)
        .build()
        .map_err(|error| {
            Error::InvalidInput(format!("invalid compression thread count: {error}"))
        })?;
    pool.install(|| {
        blocks
            .par_iter()
            .map(|block| compress_block(block))
            .collect()
    })
}

#[cfg(test)]
fn encode_postings_groups<I>(groups: I, compression_threads: usize) -> Result<PostingEncoding>
where
    I: IntoIterator<Item = (String, Vec<u64>)>,
{
    let mut encoder = PostingEncoder::new();
    for (term, spans) in groups {
        encoder.add(term, &spans)?;
    }
    encoder.finish(compression_threads)
}

#[allow(clippy::type_complexity)]
#[cfg(test)]
fn encode_postings(
    term_spans: &BTreeMap<String, BTreeSet<u64>>,
) -> Result<(Vec<u8>, Vec<PostingDirectoryEntry>, Vec<u8>, u64)> {
    let groups = term_spans
        .iter()
        .map(|(term, spans)| (term.clone(), spans.iter().copied().collect::<Vec<_>>()));
    let (terms, directory, data, before, _, _) = encode_postings_groups(groups, 1)?;
    Ok((terms, directory, data, before))
}

fn encode_for_values(values: &[u64]) -> Result<(Vec<u8>, u64, u8)> {
    if values.is_empty() {
        return Ok((Vec::new(), 0, 0));
    }
    let base = *values
        .iter()
        .min()
        .ok_or_else(|| Error::InvalidInput("empty integer stream".into()))?;
    encode_for_values_with_base(values, base)
}

fn encode_for_values_with_base(values: &[u64], base: u64) -> Result<(Vec<u8>, u64, u8)> {
    if values.is_empty() {
        return Ok((Vec::new(), base, 0));
    }
    let maximum = values
        .iter()
        .map(|value| value.checked_sub(base))
        .collect::<Option<Vec<_>>>()
        .ok_or_else(|| Error::InvalidInput("frame value below base".into()))?
        .into_iter()
        .max()
        .unwrap_or(0);
    let bit_width = if maximum == 0 {
        0
    } else {
        64 - maximum.leading_zeros() as u8
    };
    if bit_width > 63 {
        return Ok((Vec::new(), base, 64));
    }
    let bit_count = values
        .len()
        .checked_mul(bit_width as usize)
        .ok_or(Error::InvalidCoordinate)?;
    let byte_count = bit_count.checked_add(7).ok_or(Error::InvalidCoordinate)? / 8;
    let mut packed = vec![0_u8; byte_count];
    if bit_width > 0 {
        for (index, value) in values.iter().enumerate() {
            let relative = value
                .checked_sub(base)
                .ok_or_else(|| Error::InvalidInput("frame value below base".into()))?;
            let bit_offset = index
                .checked_mul(bit_width as usize)
                .ok_or(Error::InvalidCoordinate)?;
            for bit in 0..bit_width as usize {
                if relative & (1_u64 << bit) != 0 {
                    let target = bit_offset + bit;
                    packed[target / 8] |= 1 << (target % 8);
                }
            }
        }
    }
    Ok((packed, base, bit_width))
}

fn encode_span_block(spans: &[SpanKey]) -> Result<(Vec<u8>, u8, u8, u64)> {
    if spans.is_empty() {
        return Err(Error::InvalidInput(
            "cannot encode an empty span block".into(),
        ));
    }
    let first_start = spans[0].start;
    let mut start_deltas = Vec::with_capacity(spans.len().saturating_sub(1));
    let mut previous_start = first_start;
    for span in spans.iter().skip(1) {
        let delta = span
            .start
            .checked_sub(previous_start)
            .ok_or_else(|| Error::InvalidInput("span block is not sorted".into()))?;
        start_deltas.push(delta);
        previous_start = span.start;
    }
    let lengths = spans.iter().map(|span| span.length).collect::<Vec<_>>();
    let mut start_varints = Vec::new();
    for delta in &start_deltas {
        write_varint(&mut start_varints, *delta);
    }
    let mut length_varints = Vec::new();
    for length in &lengths {
        write_varint(&mut length_varints, *length);
    }
    // The first start is already stored in the fixed directory entry, so both
    // candidate streams encode only the remaining rows.  Absolute start
    // varints (encoding 0) are dominated by these sorted deltas: every delta
    // is no larger than the corresponding absolute start, and the first row
    // is omitted from both representations.
    let relative_starts = spans
        .iter()
        .skip(1)
        .map(|span| span.start - first_start)
        .collect::<Vec<_>>();
    let (start_for, _, start_width) = encode_for_values_with_base(&relative_starts, 0)?;
    let (length_for, length_base, length_width) = encode_for_values(&lengths)?;
    let (start_encoding, start_stream, start_bit_width) =
        if start_width <= 63 && start_for.len() < start_varints.len() {
            (2_u8, start_for, start_width)
        } else {
            (1_u8, start_varints, 0_u8)
        };
    let (length_encoding, length_stream, length_base, length_bit_width) =
        if length_width <= 63 && length_for.len() < length_varints.len() {
            (2_u8, length_for, length_base, length_width)
        } else {
            (0_u8, length_varints, 0_u64, 0_u8)
        };
    let mut payload = Vec::new();
    put_u32(
        &mut payload,
        u32::try_from(start_stream.len())
            .map_err(|_| Error::InvalidInput("span start stream exceeds 4 GiB".into()))?,
    );
    put_u32(
        &mut payload,
        u32::try_from(length_stream.len())
            .map_err(|_| Error::InvalidInput("span length stream exceeds 4 GiB".into()))?,
    );
    put_u64(&mut payload, length_base);
    payload.push(start_encoding);
    payload.push(length_encoding);
    payload.push(start_bit_width);
    payload.push(length_bit_width);
    payload.extend_from_slice(&start_stream);
    payload.extend_from_slice(&length_stream);
    Ok((payload, start_encoding, length_encoding, first_start))
}

#[allow(clippy::type_complexity)]
#[cfg(test)]
fn encode_span_blocks(
    spans: &[SpanKey],
    spans_per_block: usize,
) -> Result<(Vec<SpanDirectoryEntry>, Vec<u8>, u64, u64, u64)> {
    encode_span_blocks_with_threads(spans, spans_per_block, 1)
}

#[allow(clippy::type_complexity)]
fn encode_span_blocks_with_threads(
    spans: &[SpanKey],
    spans_per_block: usize,
    compression_threads: usize,
) -> Result<(Vec<SpanDirectoryEntry>, Vec<u8>, u64, u64, u64)> {
    use rayon::prelude::*;

    if spans_per_block == 0 {
        return Err(Error::InvalidInput(
            "span block size must be nonzero".into(),
        ));
    }
    let mut ranges = Vec::new();
    let mut block_start = 0;
    while block_start < spans.len() {
        let reference_id = spans[block_start].reference_id;
        let mut block_end = block_start;
        while block_end < spans.len()
            && spans[block_end].reference_id == reference_id
            && block_end - block_start < spans_per_block
        {
            block_end += 1;
        }
        ranges.push((block_start, block_end));
        block_start = block_end;
    }

    let thread_count = compression_threads.max(1);
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(thread_count)
        .build()
        .map_err(|error| {
            Error::InvalidInput(format!("invalid compression thread count: {error}"))
        })?;
    let encoded_blocks = pool.install(|| {
        ranges
            .par_iter()
            .map(|(block_start, block_end)| {
                let block = &spans[*block_start..*block_end];
                let (payload, start_encoding, length_encoding, first_start) =
                    encode_span_block(block)?;
                let encoded = compress_block(&payload)?;
                Ok::<_, Error>((
                    block.len(),
                    block[0].reference_id,
                    first_start,
                    start_encoding,
                    length_encoding,
                    payload.len() as u64,
                    encoded,
                ))
            })
            .collect::<Result<Vec<_>>>()
    })?;

    let mut directory = Vec::with_capacity(encoded_blocks.len());
    let mut data = Vec::new();
    let mut fixed_width_bytes = 0_u64;
    let mut structural_bytes = 0_u64;
    let mut compressed_bytes = 0_u64;
    let mut first_span_id = 0_u64;
    for (
        span_count,
        reference_id,
        first_start,
        start_encoding,
        length_encoding,
        payload_length,
        encoded,
    ) in encoded_blocks
    {
        let compressed_offset = u64::try_from(data.len())
            .map_err(|_| Error::InvalidInput("span data exceeds 4 GiB".into()))?;
        data.extend_from_slice(&encoded.compressed);
        let span_count = u32::try_from(span_count)
            .map_err(|_| Error::InvalidInput("span block has too many rows".into()))?;
        directory.push(SpanDirectoryEntry {
            first_span_id,
            span_count,
            reference_id,
            first_start,
            compressed_offset,
            compressed_length: u32::try_from(encoded.compressed.len())
                .map_err(|_| Error::InvalidInput("compressed span block exceeds 4 GiB".into()))?,
            uncompressed_length: encoded.uncompressed_length,
            checksum: encoded.checksum,
            start_encoding,
            length_encoding,
            compression: encoded.compression,
        });
        fixed_width_bytes = fixed_width_bytes
            .checked_add(u64::from(span_count).saturating_mul(20))
            .ok_or(Error::InvalidCoordinate)?;
        structural_bytes = structural_bytes
            .checked_add(payload_length)
            .ok_or(Error::InvalidCoordinate)?;
        compressed_bytes = compressed_bytes
            .checked_add(encoded.compressed.len() as u64)
            .ok_or(Error::InvalidCoordinate)?;
        first_span_id = first_span_id
            .checked_add(u64::from(span_count))
            .ok_or(Error::InvalidCoordinate)?;
    }
    Ok((
        directory,
        data,
        fixed_width_bytes,
        structural_bytes,
        compressed_bytes,
    ))
}

fn encode_span_directory(entries: &[SpanDirectoryEntry]) -> Result<Vec<u8>> {
    let mut bytes = Vec::with_capacity(
        entries
            .len()
            .checked_mul(SPAN_DIRECTORY_ENTRY_SIZE)
            .ok_or(Error::InvalidCoordinate)?,
    );
    for entry in entries {
        put_u64(&mut bytes, entry.first_span_id);
        put_u32(&mut bytes, entry.span_count);
        put_u32(&mut bytes, entry.reference_id);
        put_u64(&mut bytes, entry.first_start);
        put_u64(&mut bytes, entry.compressed_offset);
        put_u32(&mut bytes, entry.compressed_length);
        put_u32(&mut bytes, entry.uncompressed_length);
        put_u32(&mut bytes, entry.checksum);
        bytes.push(entry.start_encoding);
        bytes.push(entry.length_encoding);
        bytes.push(entry.compression);
        bytes.push(0);
        put_u32(&mut bytes, 0);
    }
    Ok(bytes)
}

fn encode_posting_directory(entries: &[PostingDirectoryEntry]) -> Result<Vec<u8>> {
    let mut bytes = Vec::with_capacity(
        entries
            .len()
            .checked_mul(POSTINGS_DIRECTORY_ENTRY_SIZE)
            .ok_or(Error::InvalidCoordinate)?,
    );
    for entry in entries {
        put_u64(&mut bytes, entry.compressed_offset);
        put_u32(&mut bytes, entry.compressed_length);
        put_u32(&mut bytes, entry.uncompressed_length);
        put_u32(&mut bytes, entry.checksum);
        bytes.push(entry.compression);
        bytes.extend_from_slice(&[0; 3]);
        put_u64(&mut bytes, 0);
    }
    Ok(bytes)
}

fn set_u16(bytes: &mut [u8], offset: usize, value: u16) {
    bytes[offset..offset + 2].copy_from_slice(&value.to_le_bytes());
}

fn set_u32(bytes: &mut [u8], offset: usize, value: u32) {
    bytes[offset..offset + 4].copy_from_slice(&value.to_le_bytes());
}

fn set_u64(bytes: &mut [u8], offset: usize, value: u64) {
    bytes[offset..offset + 8].copy_from_slice(&value.to_le_bytes());
}

fn encode_directory(entries: &[SectionDirectoryEntry]) -> Result<Vec<u8>> {
    let mut bytes = Vec::with_capacity(
        entries
            .len()
            .checked_mul(DIRECTORY_ENTRY_SIZE)
            .ok_or(Error::InvalidCoordinate)?,
    );
    for entry in entries {
        put_u32(&mut bytes, entry.kind as u32);
        put_u32(&mut bytes, entry.flags);
        put_u64(&mut bytes, entry.offset);
        put_u64(&mut bytes, entry.length);
        put_u64(&mut bytes, entry.item_count);
        put_u32(&mut bytes, entry.checksum);
        put_u32(&mut bytes, 0);
    }
    Ok(bytes)
}

#[allow(clippy::too_many_arguments)]
fn serialize_index(
    attributes: &[String],
    case_sensitive: bool,
    gff_fingerprint: [u8; 32],
    coordinate_index_fingerprint: [u8; 32],
    reference_dictionary_fingerprint: [u8; 32],
    term_count: u64,
    unique_span_count: u64,
    posting_count: u64,
    postings_block_count: u64,
    span_block_count: u64,
    reference_count: u32,
    terms: Vec<u8>,
    postings_directory: &[PostingDirectoryEntry],
    postings_data: Vec<u8>,
    spans_directory: Vec<u8>,
    spans_data: Vec<u8>,
    spans_per_block: usize,
) -> Result<Vec<u8>> {
    let attributes_section = encode_attributes(attributes)?;
    let postings_data_length = postings_data.len() as u64;
    let spans_data_length = spans_data.len() as u64;
    let sections = [
        (
            SectionKind::Attributes,
            attributes_section,
            attributes.len() as u64,
        ),
        (SectionKind::Terms, terms, term_count),
        (
            SectionKind::PostingsDirectory,
            encode_posting_directory(postings_directory)?,
            postings_block_count,
        ),
        (
            SectionKind::PostingsData,
            postings_data,
            postings_data_length,
        ),
        (
            SectionKind::SpansDirectory,
            spans_directory,
            span_block_count,
        ),
        (SectionKind::SpansData, spans_data, spans_data_length),
    ];
    let section_directory_offset = HEADER_SIZE as u64;
    let section_directory_length = u64::try_from(
        sections
            .len()
            .checked_mul(DIRECTORY_ENTRY_SIZE)
            .ok_or(Error::InvalidCoordinate)?,
    )
    .map_err(|_| Error::InvalidInput("section directory exceeds 4 GiB".into()))?;
    let mut output = vec![0_u8; HEADER_SIZE + section_directory_length as usize];
    let mut directory = Vec::with_capacity(sections.len());
    let section_count = sections.len();
    for (kind, section, item_count) in &sections {
        let offset = u64::try_from(output.len())
            .map_err(|_| Error::InvalidInput("GNI exceeds 4 GiB".into()))?;
        if section.len() as u64 > MAX_SECTION_BYTES {
            return Err(Error::InvalidInput(
                "GNI section exceeds safety limit".into(),
            ));
        }
        let section_checksum = checksum(section);
        output.extend_from_slice(section);
        directory.push(SectionDirectoryEntry {
            kind: *kind,
            flags: 0,
            offset,
            length: section.len() as u64,
            item_count: *item_count,
            checksum: section_checksum,
        });
    }
    let directory_bytes = encode_directory(&directory)?;
    let directory_start = HEADER_SIZE;
    output[directory_start..directory_start + directory_bytes.len()]
        .copy_from_slice(&directory_bytes);

    output[..4].copy_from_slice(&MAGIC);
    set_u16(&mut output, 4, MAJOR_VERSION);
    set_u16(&mut output, 6, MINOR_VERSION);
    set_u32(&mut output, 8, 1); // zstd blocks are selected independently.
    output[12] = BYTE_ORDER_LITTLE;
    output[13] = COORDINATE_ZERO_BASED_HALF_OPEN;
    output[14] = if case_sensitive {
        NORMALIZATION_CASE_SENSITIVE
    } else {
        NORMALIZATION_ASCII_LOWER
    };
    output[15] = 0;
    set_u32(&mut output, 16, HEADER_SIZE as u32);
    set_u32(&mut output, 20, DIRECTORY_ENTRY_SIZE as u32);
    set_u32(&mut output, 24, section_count as u32);
    set_u32(&mut output, 28, 0);
    set_u64(&mut output, 32, term_count);
    set_u64(&mut output, 40, unique_span_count);
    set_u64(&mut output, 48, posting_count);
    set_u64(&mut output, 56, postings_block_count);
    set_u64(&mut output, 64, span_block_count);
    output[72..104].copy_from_slice(&gff_fingerprint);
    output[104..136].copy_from_slice(&coordinate_index_fingerprint);
    output[136..168].copy_from_slice(&reference_dictionary_fingerprint);
    set_u32(&mut output, 168, attributes.len() as u32);
    set_u32(
        &mut output,
        172,
        u32::try_from(spans_per_block)
            .map_err(|_| Error::InvalidInput("span block size exceeds 4 GiB".into()))?,
    );
    set_u32(&mut output, 200, reference_count);
    set_u64(&mut output, 176, section_directory_offset);
    set_u64(&mut output, 184, section_directory_length);
    let output_len =
        u64::try_from(output.len()).map_err(|_| Error::InvalidInput("GNI exceeds 4 GiB".into()))?;
    set_u64(&mut output, 192, output_len);
    // Bytes 204..256 are reserved and remain zero for forward compatibility.
    Ok(output)
}

#[derive(Clone, Debug)]
struct TermSpanObservation {
    term: String,
    span: SpanKey,
}

struct SpillCollector {
    memory_budget_bytes: usize,
    parent: PathBuf,
    prefix: String,
    spans: Vec<SpanKey>,
    pairs: Vec<TermSpanObservation>,
    estimated_bytes: usize,
    peak_working_set_bytes: usize,
    run_paths: Vec<PathBuf>,
}

fn create_spill_file(parent: &Path, prefix: &str) -> Result<(File, PathBuf)> {
    let counter = RUN_COUNTER.fetch_add(1, Ordering::Relaxed);
    let path = parent.join(format!(
        ".{}-run-{}-{}",
        prefix,
        std::process::id(),
        counter
    ));
    let file = fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&path)?;
    Ok((file, path))
}

struct SpillRunGuard {
    path: Option<PathBuf>,
}

impl Drop for SpillRunGuard {
    fn drop(&mut self) {
        if let Some(path) = self.path.take() {
            let _ = fs::remove_file(path);
        }
    }
}

impl SpillCollector {
    fn new(memory_budget_bytes: usize, parent: &Path, prefix: &str) -> Self {
        Self {
            memory_budget_bytes,
            parent: parent.to_path_buf(),
            prefix: prefix.to_owned(),
            spans: Vec::new(),
            pairs: Vec::new(),
            estimated_bytes: 0,
            peak_working_set_bytes: 0,
            run_paths: Vec::new(),
        }
    }

    fn add(&mut self, span: SpanKey, terms: &[String]) -> Result<bool> {
        self.spans.push(span);
        self.estimated_bytes = self
            .estimated_bytes
            .saturating_add(std::mem::size_of::<SpanKey>());
        for term in terms {
            self.estimated_bytes = self
                .estimated_bytes
                .saturating_add(std::mem::size_of::<TermSpanObservation>() + term.len());
            self.pairs.push(TermSpanObservation {
                term: term.clone(),
                span,
            });
        }
        self.peak_working_set_bytes = self.peak_working_set_bytes.max(self.estimated_bytes);
        Ok(self.estimated_bytes >= self.memory_budget_bytes)
    }

    fn spill(&mut self) -> Result<Option<PathBuf>> {
        if self.spans.is_empty() && self.pairs.is_empty() {
            return Ok(None);
        }
        self.spans.sort_unstable();
        self.spans.dedup();
        self.pairs.sort_unstable_by(|left, right| {
            left.term
                .cmp(&right.term)
                .then_with(|| left.span.cmp(&right.span))
        });
        self.pairs
            .dedup_by(|left, right| left.term == right.term && left.span == right.span);

        let (file, path) = create_spill_file(&self.parent, &self.prefix)?;
        let mut guard = SpillRunGuard {
            path: Some(path.clone()),
        };
        let mut writer = io::BufWriter::new(file);
        writer.write_all(&RUN_MAGIC)?;
        writer.write_all(&(self.spans.len() as u64).to_le_bytes())?;
        for span in &self.spans {
            writer.write_all(&span.reference_id.to_le_bytes())?;
            writer.write_all(&span.start.to_le_bytes())?;
            writer.write_all(&span.length.to_le_bytes())?;
        }
        writer.write_all(&(self.pairs.len() as u64).to_le_bytes())?;
        for pair in &self.pairs {
            writer.write_all(
                &u32::try_from(pair.term.len())
                    .map_err(|_| Error::InvalidInput("term exceeds 4 GiB".into()))?
                    .to_le_bytes(),
            )?;
            writer.write_all(pair.term.as_bytes())?;
            writer.write_all(&pair.span.reference_id.to_le_bytes())?;
            writer.write_all(&pair.span.start.to_le_bytes())?;
            writer.write_all(&pair.span.length.to_le_bytes())?;
        }
        writer.flush()?;
        self.run_paths.push(path.clone());
        guard.path = None;
        self.compact_runs_if_needed()?;
        self.spans.clear();
        self.pairs.clear();
        self.estimated_bytes = 0;
        Ok(Some(path))
    }

    fn finish(&mut self) -> Result<()> {
        self.spill().map(|_| ())
    }

    fn run_paths(&self) -> &[PathBuf] {
        &self.run_paths
    }

    fn compact_runs_if_needed(&mut self) -> Result<()> {
        while self.run_paths.len() > MAX_RUN_FANIN {
            let group = self.run_paths.drain(..MAX_RUN_FANIN).collect::<Vec<_>>();
            let merged = match merge_run_group(&group, &self.parent, &self.prefix) {
                Ok(path) => path,
                Err(error) => {
                    self.run_paths.splice(0..0, group);
                    return Err(error);
                }
            };
            for path in group {
                let _ = fs::remove_file(path);
            }
            self.run_paths.push(merged);
        }
        Ok(())
    }

    fn peak_working_set_bytes(&self) -> u64 {
        self.peak_working_set_bytes as u64
    }
}

impl Drop for SpillCollector {
    fn drop(&mut self) {
        for path in &self.run_paths {
            let _ = fs::remove_file(path);
        }
    }
}

fn read_run_u32(reader: &mut impl Read, context: &str) -> Result<u32> {
    let mut bytes = [0_u8; 4];
    reader
        .read_exact(&mut bytes)
        .map_err(|error| Error::InvalidInput(format!("truncated {context}: {error}")))?;
    Ok(u32::from_le_bytes(bytes))
}

fn read_run_u64(reader: &mut impl Read, context: &str) -> Result<u64> {
    let mut bytes = [0_u8; 8];
    reader
        .read_exact(&mut bytes)
        .map_err(|error| Error::InvalidInput(format!("truncated {context}: {error}")))?;
    Ok(u64::from_le_bytes(bytes))
}

struct SpanRunReader {
    reader: io::BufReader<File>,
    remaining: u64,
}

impl SpanRunReader {
    fn open(path: &Path) -> Result<Self> {
        let mut reader = io::BufReader::new(File::open(path)?);
        let mut magic = [0_u8; 8];
        reader.read_exact(&mut magic)?;
        if magic != RUN_MAGIC {
            return Err(Error::InvalidInput("invalid GNI spill run magic".into()));
        }
        let remaining = read_run_u64(&mut reader, "spill span count")?;
        Ok(Self { reader, remaining })
    }

    fn next(&mut self) -> Result<Option<SpanKey>> {
        if self.remaining == 0 {
            return Ok(None);
        }
        let span = SpanKey {
            reference_id: read_run_u32(&mut self.reader, "spill reference ID")?,
            start: read_run_u64(&mut self.reader, "spill span start")?,
            length: read_run_u64(&mut self.reader, "spill span length")?,
        };
        self.remaining -= 1;
        Ok(Some(span))
    }
}

struct TermRunReader {
    reader: io::BufReader<File>,
    remaining: u64,
}

impl TermRunReader {
    fn open(path: &Path) -> Result<Self> {
        let mut reader = io::BufReader::new(File::open(path)?);
        let mut magic = [0_u8; 8];
        reader.read_exact(&mut magic)?;
        if magic != RUN_MAGIC {
            return Err(Error::InvalidInput("invalid GNI spill run magic".into()));
        }
        let span_count = read_run_u64(&mut reader, "spill span count")?;
        let span_bytes = span_count
            .checked_mul(20)
            .ok_or_else(|| Error::InvalidInput("spill span section is too large".into()))?;
        reader.seek(io::SeekFrom::Current(i64::try_from(span_bytes).map_err(
            |_| Error::InvalidInput("spill span section exceeds seek range".into()),
        )?))?;
        let remaining = read_run_u64(&mut reader, "spill posting pair count")?;
        Ok(Self { reader, remaining })
    }

    fn next(&mut self) -> Result<Option<TermSpanObservation>> {
        if self.remaining == 0 {
            return Ok(None);
        }
        let term_length = read_run_u32(&mut self.reader, "spill term length")? as usize;
        if term_length > 1 << 30 {
            return Err(Error::InvalidInput("spill term is too large".into()));
        }
        let mut term_bytes = vec![0_u8; term_length];
        self.reader.read_exact(&mut term_bytes)?;
        let term = String::from_utf8(term_bytes)
            .map_err(|_| Error::InvalidInput("spill term is not UTF-8".into()))?;
        let span = SpanKey {
            reference_id: read_run_u32(&mut self.reader, "spill reference ID")?,
            start: read_run_u64(&mut self.reader, "spill span start")?,
            length: read_run_u64(&mut self.reader, "spill span length")?,
        };
        self.remaining -= 1;
        Ok(Some(TermSpanObservation { term, span }))
    }
}

fn merge_run_group(group: &[PathBuf], parent: &Path, prefix: &str) -> Result<PathBuf> {
    let (file, path) = create_spill_file(parent, prefix)?;
    let mut guard = SpillRunGuard {
        path: Some(path.clone()),
    };
    let mut writer = io::BufWriter::new(file);
    writer.write_all(&RUN_MAGIC)?;
    writer.write_all(&0_u64.to_le_bytes())?;

    let mut span_readers = group
        .iter()
        .map(|run| SpanRunReader::open(run))
        .collect::<Result<Vec<_>>>()?;
    let mut span_heap = BinaryHeap::<Reverse<(SpanKey, usize)>>::new();
    for (index, reader) in span_readers.iter_mut().enumerate() {
        if let Some(span) = reader.next()? {
            span_heap.push(Reverse((span, index)));
        }
    }
    let mut span_count = 0_u64;
    let mut previous_span = None;
    while let Some(Reverse((span, index))) = span_heap.pop() {
        if previous_span != Some(span) {
            writer.write_all(&span.reference_id.to_le_bytes())?;
            writer.write_all(&span.start.to_le_bytes())?;
            writer.write_all(&span.length.to_le_bytes())?;
            span_count += 1;
            previous_span = Some(span);
        }
        if let Some(next) = span_readers[index].next()? {
            span_heap.push(Reverse((next, index)));
        }
    }
    writer.flush()?;
    let mut file = writer
        .into_inner()
        .map_err(|error| Error::Io(error.into_error()))?;
    file.seek(io::SeekFrom::Start(8))?;
    file.write_all(&span_count.to_le_bytes())?;
    file.seek(io::SeekFrom::End(0))?;
    let pair_count_offset = file.stream_position()?;
    file.write_all(&0_u64.to_le_bytes())?;
    let mut writer = io::BufWriter::new(file);

    let mut term_readers = group
        .iter()
        .map(|run| TermRunReader::open(run))
        .collect::<Result<Vec<_>>>()?;
    let mut term_heap = BinaryHeap::<Reverse<(String, SpanKey, usize)>>::new();
    for (index, reader) in term_readers.iter_mut().enumerate() {
        if let Some(pair) = reader.next()? {
            term_heap.push(Reverse((pair.term, pair.span, index)));
        }
    }
    let mut pair_count = 0_u64;
    let mut previous_pair: Option<(String, SpanKey)> = None;
    while let Some(Reverse((term, span, index))) = term_heap.pop() {
        if previous_pair
            .as_ref()
            .is_none_or(|previous| previous.0 != term || previous.1 != span)
        {
            writer.write_all(
                &u32::try_from(term.len())
                    .map_err(|_| Error::InvalidInput("term exceeds 4 GiB".into()))?
                    .to_le_bytes(),
            )?;
            writer.write_all(term.as_bytes())?;
            writer.write_all(&span.reference_id.to_le_bytes())?;
            writer.write_all(&span.start.to_le_bytes())?;
            writer.write_all(&span.length.to_le_bytes())?;
            pair_count += 1;
            previous_pair = Some((term, span));
        }
        if let Some(next) = term_readers[index].next()? {
            term_heap.push(Reverse((next.term, next.span, index)));
        }
    }
    writer.flush()?;
    let mut file = writer
        .into_inner()
        .map_err(|error| Error::Io(error.into_error()))?;
    file.seek(io::SeekFrom::Start(pair_count_offset))?;
    file.write_all(&pair_count.to_le_bytes())?;
    file.flush()?;
    guard.path = None;
    Ok(path)
}

fn merge_unique_spans(run_paths: &[PathBuf]) -> Result<Vec<SpanKey>> {
    let mut readers = run_paths
        .iter()
        .map(|path| SpanRunReader::open(path))
        .collect::<Result<Vec<_>>>()?;
    let mut heap = BinaryHeap::<Reverse<(SpanKey, usize)>>::new();
    for (index, reader) in readers.iter_mut().enumerate() {
        if let Some(span) = reader.next()? {
            heap.push(Reverse((span, index)));
        }
    }
    let mut spans = Vec::new();
    let mut previous = None;
    while let Some(Reverse((span, index))) = heap.pop() {
        if previous != Some(span) {
            spans.push(span);
            previous = Some(span);
        }
        if let Some(next) = readers[index].next()? {
            heap.push(Reverse((next, index)));
        }
    }
    Ok(spans)
}

fn encode_postings_from_runs(
    run_paths: &[PathBuf],
    coordinate_spans: &[SpanKey],
    compression_threads: usize,
) -> Result<PostingEncoding> {
    let mut readers = run_paths
        .iter()
        .map(|path| TermRunReader::open(path))
        .collect::<Result<Vec<_>>>()?;
    let mut heap = BinaryHeap::<Reverse<(String, SpanKey, usize)>>::new();
    for (index, reader) in readers.iter_mut().enumerate() {
        if let Some(pair) = reader.next()? {
            heap.push(Reverse((pair.term, pair.span, index)));
        }
    }

    let mut encoder = PostingEncoder::new();
    while let Some(Reverse((term, span, index))) = heap.pop() {
        let mut spans = Vec::new();
        let mut previous_span = None;
        let mut push_span = |span: SpanKey, spans: &mut Vec<u64>| -> Result<()> {
            if previous_span == Some(span) {
                return Ok(());
            }
            let span_id = coordinate_spans
                .binary_search(&span)
                .map_err(|_| Error::InvalidInput("term references an unknown span".into()))?;
            spans.push(span_id as u64);
            previous_span = Some(span);
            Ok(())
        };
        push_span(span, &mut spans)?;
        if let Some(next) = readers[index].next()? {
            heap.push(Reverse((next.term, next.span, index)));
        }
        while heap.peek().is_some_and(|entry| entry.0.0 == term) {
            let Reverse((_, next_span, next_index)) = heap.pop().unwrap();
            push_span(next_span, &mut spans)?;
            if let Some(next) = readers[next_index].next()? {
                heap.push(Reverse((next.term, next.span, next_index)));
            }
        }
        encoder.add(term, &spans)?;
    }
    encoder.finish(compression_threads)
}

fn report_progress(
    options: &BuildOptions,
    phase: BuildPhase,
    records_processed: u64,
    bytes_read: &Arc<AtomicU64>,
    started: Instant,
) {
    if let Some(callback) = &options.progress {
        callback(BuildProgress {
            phase,
            records_processed,
            bytes_read: bytes_read.load(Ordering::Relaxed),
            elapsed: started.elapsed(),
        });
    }
}

fn scan_gff_path<F>(
    path: &Path,
    bgzf_threads: usize,
    bytes_read: Arc<AtomicU64>,
    process: F,
) -> Result<([u8; 32], u64)>
where
    F: FnMut(gff::feature::RecordBuf) -> Result<()>,
{
    let mut source = File::open(path)?;
    let mut magic = [0_u8; 2];
    let magic_length = source.read(&mut magic)?;
    source.rewind()?;
    let hashing = HashingReader::with_counter(source, bytes_read);
    if magic_length == magic.len() && magic == [0x1f, 0x8b] {
        if bgzf_threads > 1 {
            let workers = NonZeroUsize::new(bgzf_threads)
                .ok_or_else(|| Error::InvalidInput("BGZF worker count must be positive".into()))?;
            scan_gff_reader(
                gff::io::Reader::new(bgzf::io::MultithreadedReader::with_worker_count(
                    workers, hashing,
                )),
                process,
            )
        } else {
            scan_gff_reader(
                gff::io::Reader::new(bgzf::io::Reader::new(hashing)),
                process,
            )
        }
    } else {
        scan_gff_reader(gff::io::Reader::new(BufReader::new(hashing)), process)
    }
}

/// Builds a deterministic GNI beside a BGZF or plain GFF3 source.
pub fn build_name_index(
    gff_path: impl AsRef<Path>,
    coordinate_index_path: impl AsRef<Path>,
    destination: impl AsRef<Path>,
    options: &NameIndexOptions,
) -> Result<IndexStats> {
    build_name_index_with_options(
        gff_path,
        coordinate_index_path,
        destination,
        options,
        &BuildOptions::default(),
    )
}

/// Builds a deterministic GNI with explicit resource and progress controls.
pub fn build_name_index_with_options(
    gff_path: impl AsRef<Path>,
    coordinate_index_path: impl AsRef<Path>,
    destination: impl AsRef<Path>,
    options: &NameIndexOptions,
    build_options: &BuildOptions,
) -> Result<IndexStats> {
    build_name_index_with_options_and_span_block_size(
        gff_path,
        coordinate_index_path,
        destination,
        options,
        DEFAULT_SPANS_PER_BLOCK,
        build_options,
    )
}

/// Builds a GNI with an explicit reference-specific span-block row target.
/// The default [`build_name_index`] value is 4,096; this variant is provided
/// for reproducible compression experiments and deployment tuning.
pub fn build_name_index_with_span_block_size(
    gff_path: impl AsRef<Path>,
    coordinate_index_path: impl AsRef<Path>,
    destination: impl AsRef<Path>,
    options: &NameIndexOptions,
    spans_per_block: usize,
) -> Result<IndexStats> {
    build_name_index_with_options_and_span_block_size(
        gff_path,
        coordinate_index_path,
        destination,
        options,
        spans_per_block,
        &BuildOptions::default(),
    )
}

/// Builds a GNI with explicit span-block and resource controls.
pub fn build_name_index_with_options_and_span_block_size(
    gff_path: impl AsRef<Path>,
    coordinate_index_path: impl AsRef<Path>,
    destination: impl AsRef<Path>,
    options: &NameIndexOptions,
    spans_per_block: usize,
    build_options: &BuildOptions,
) -> Result<IndexStats> {
    if !(1..=1_000_000).contains(&spans_per_block) {
        return Err(Error::InvalidInput(
            "span block size must be between 1 and 1,000,000".into(),
        ));
    }
    if build_options.memory_budget_bytes == 0 {
        return Err(Error::InvalidInput("memory budget must be positive".into()));
    }
    if build_options.compression_threads == 0 {
        return Err(Error::InvalidInput(
            "compression worker count must be positive".into(),
        ));
    }
    if build_options.bgzf_threads == 0 {
        return Err(Error::InvalidInput(
            "BGZF worker count must be positive".into(),
        ));
    }
    let started = Instant::now();
    let gff_path = gff_path.as_ref();
    let coordinate_index_path = coordinate_index_path.as_ref();
    let destination = destination.as_ref();
    let name_options = NameIndexOptions::new(options.attributes.clone(), options.case_sensitive)?;
    let destination_parent = destination.parent().unwrap_or_else(|| Path::new("."));
    fs::create_dir_all(destination_parent)?;
    let (coordinate_index, coordinate_index_fingerprint) =
        read_coordinate_index_with_fingerprint(coordinate_index_path)?;
    let dictionary = coordinate_index.dictionary()?;
    let reference_ids = dictionary
        .names
        .iter()
        .enumerate()
        .map(|(index, name)| {
            Ok((
                name.clone(),
                u32::try_from(index)
                    .map_err(|_| Error::InvalidInput("too many reference sequences".into()))?,
            ))
        })
        .collect::<Result<HashMap<_, _>>>()?;

    let configured: HashSet<&str> = name_options.attributes.iter().map(String::as_str).collect();
    let mut collector = SpillCollector::new(
        build_options.memory_budget_bytes,
        destination_parent,
        &destination
            .file_name()
            .map(|name| name.to_string_lossy())
            .unwrap_or_else(|| "gni".into()),
    );
    let bytes_read = Arc::new(AtomicU64::new(0));
    let mut records_processed = 0_u64;
    let mut records_indexed = 0_u64;
    let mut span_observations = 0_u64;
    let mut pair_observations = 0_u64;
    let mut spill_duration = Duration::ZERO;
    let mut scan_spill_duration = Duration::ZERO;
    let scan_started = Instant::now();
    let mut process_record = |record: gff::feature::RecordBuf| -> Result<()> {
        records_processed = records_processed
            .checked_add(1)
            .ok_or(Error::InvalidCoordinate)?;
        let extracted = extract_record(
            &record,
            &reference_ids,
            &configured,
            name_options.case_sensitive,
        )?;
        if extracted.terms.is_empty() {
            if records_processed.is_multiple_of(PROGRESS_RECORD_INTERVAL) {
                report_progress(
                    build_options,
                    BuildPhase::Scan,
                    records_processed,
                    &bytes_read,
                    started,
                );
            }
            return Ok(());
        }
        records_indexed = records_indexed
            .checked_add(1)
            .ok_or(Error::InvalidCoordinate)?;
        span_observations = span_observations
            .checked_add(1)
            .ok_or(Error::InvalidCoordinate)?;
        pair_observations = pair_observations
            .checked_add(extracted.terms.len() as u64)
            .ok_or(Error::InvalidCoordinate)?;
        let should_spill = collector.add(extracted.span, &extracted.terms)?;
        if should_spill {
            let spill_started = Instant::now();
            collector.spill()?;
            let elapsed = spill_started.elapsed();
            spill_duration = spill_duration.saturating_add(elapsed);
            scan_spill_duration = scan_spill_duration.saturating_add(elapsed);
            report_progress(
                build_options,
                BuildPhase::Spill,
                records_processed,
                &bytes_read,
                started,
            );
        } else if records_processed.is_multiple_of(PROGRESS_RECORD_INTERVAL) {
            report_progress(
                build_options,
                BuildPhase::Scan,
                records_processed,
                &bytes_read,
                started,
            );
        }
        Ok(())
    };
    let (gff_fingerprint, source_bytes) = scan_gff_path(
        gff_path,
        build_options.bgzf_threads,
        Arc::clone(&bytes_read),
        &mut process_record,
    )?;
    let scan_duration = scan_started.elapsed();
    let final_spill_started = Instant::now();
    collector.finish()?;
    spill_duration = spill_duration.saturating_add(final_spill_started.elapsed());
    report_progress(
        build_options,
        BuildPhase::Scan,
        records_processed,
        &bytes_read,
        started,
    );

    let merge_started = Instant::now();
    let coordinate_spans = merge_unique_spans(collector.run_paths())?;
    let merge_duration = merge_started.elapsed();
    report_progress(
        build_options,
        BuildPhase::Merge,
        records_processed,
        &bytes_read,
        started,
    );

    let postings_started = Instant::now();
    let (terms, posting_directory, postings_data, postings_before, unique_postings, term_count) =
        encode_postings_from_runs(
            collector.run_paths(),
            &coordinate_spans,
            build_options.compression_threads,
        )?;
    let postings_duration = postings_started.elapsed();
    report_progress(
        build_options,
        BuildPhase::EncodePostings,
        records_processed,
        &bytes_read,
        started,
    );

    let spans_started = Instant::now();
    let (span_directory, spans_data, spans_fixed, spans_structural, spans_compressed) =
        encode_span_blocks_with_threads(
            &coordinate_spans,
            spans_per_block,
            build_options.compression_threads,
        )?;
    let spans_duration = spans_started.elapsed();
    let span_directory_bytes = encode_span_directory(&span_directory)?;
    report_progress(
        build_options,
        BuildPhase::EncodeSpans,
        records_processed,
        &bytes_read,
        started,
    );

    let serialize_started = Instant::now();
    let serialized = serialize_index(
        &name_options.attributes,
        name_options.case_sensitive,
        gff_fingerprint,
        coordinate_index_fingerprint,
        dictionary.fingerprint,
        term_count,
        coordinate_spans.len() as u64,
        unique_postings,
        posting_directory.len() as u64,
        span_directory.len() as u64,
        u32::try_from(dictionary.names.len())
            .map_err(|_| Error::InvalidInput("too many reference sequences".into()))?,
        terms,
        &posting_directory,
        postings_data,
        span_directory_bytes,
        spans_data,
        spans_per_block,
    )?;
    let index_bytes = serialized.len() as u64;
    atomic_write(destination, &serialized)?;
    let serialize_duration = serialize_started.elapsed();
    report_progress(
        build_options,
        BuildPhase::Serialize,
        records_processed,
        &bytes_read,
        started,
    );
    let timings = BuildTimings {
        scan: scan_duration.saturating_sub(scan_spill_duration),
        spill: spill_duration,
        merge: merge_duration,
        encode_postings: postings_duration,
        encode_spans: spans_duration,
        serialize: serialize_duration,
        total: started.elapsed(),
    };
    report_progress(
        build_options,
        BuildPhase::Complete,
        records_processed,
        &bytes_read,
        started,
    );
    let postings_after = posting_directory
        .iter()
        .map(|entry| u64::from(entry.compressed_length))
        .sum();
    let stats = IndexStats {
        records_processed,
        records_indexed,
        distinct_terms: term_count,
        unique_spans: coordinate_spans.len() as u64,
        postings: unique_postings,
        duplicate_postings_removed: pair_observations.saturating_sub(unique_postings),
        duplicate_spans_removed: span_observations.saturating_sub(coordinate_spans.len() as u64),
        index_bytes,
        postings_bytes_before_compression: postings_before,
        postings_bytes_after_compression: postings_after,
        span_bytes_fixed_width: spans_fixed,
        span_bytes_structural: spans_structural,
        span_bytes_after_compression: spans_compressed,
        bytes_per_term: index_bytes as f64 / (term_count.max(1) as f64),
        bytes_per_posting: index_bytes as f64 / (unique_postings.max(1) as f64),
        bytes_per_unique_span: index_bytes as f64 / (coordinate_spans.len().max(1) as f64),
        timings,
        peak_working_set_bytes: collector.peak_working_set_bytes(),
    };
    debug_assert_eq!(source_bytes, bytes_read.load(Ordering::Relaxed));
    Ok(stats)
}

fn atomic_write(destination: &Path, bytes: &[u8]) -> Result<()> {
    let parent = destination.parent().unwrap_or_else(|| Path::new("."));
    fs::create_dir_all(parent)?;
    let mut temporary = tempfile_path(parent, destination.file_name().unwrap_or_default())?;
    let result = (|| {
        temporary.file.write_all(bytes)?;
        temporary.file.flush()?;
        temporary.file.sync_all()?;
        Ok::<_, io::Error>(())
    })();
    if let Err(error) = result {
        let _ = fs::remove_file(&temporary.path);
        return Err(error.into());
    }
    let temporary_path = temporary.path.clone();
    drop(temporary);
    fs::rename(&temporary_path, destination).map_err(|error| {
        let _ = fs::remove_file(&temporary_path);
        Error::Io(error)
    })
}

fn tempfile_path(parent: &Path, destination_name: &std::ffi::OsStr) -> Result<TempFile> {
    let mut path = parent.to_path_buf();
    let mut name = destination_name.to_os_string();
    name.push(format!(".tmp-{}", std::process::id()));
    path.push(name);
    // A stale temporary from an interrupted build must never be reused.
    let file = fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&path)?;
    Ok(TempFile { file, path })
}

struct TempFile {
    file: File,
    path: PathBuf,
}

enum IndexStorage {
    Owned(Vec<u8>),
    Mapped(Mmap),
}

impl AsRef<[u8]> for IndexStorage {
    fn as_ref(&self) -> &[u8] {
        match self {
            Self::Owned(bytes) => bytes,
            Self::Mapped(bytes) => bytes,
        }
    }
}

/// A bounds-checked GNI reader. It only decompresses the posting and span
/// blocks needed by an individual lookup. [`Self::open_mmap`] keeps the FST
/// and fixed directories in an OS memory map instead of copying the file.
pub struct NameIndexReader {
    bytes: IndexStorage,
    metadata: IndexMetadata,
    sections: BTreeMap<SectionKind, SectionDirectoryEntry>,
    posting_directory: Vec<PostingDirectoryEntry>,
    span_directory: Vec<SpanDirectoryEntry>,
}

impl NameIndexReader {
    /// Opens and validates a GNI file without opening its source files.
    pub fn open(path: impl AsRef<Path>) -> Result<Self> {
        let bytes = fs::read(path)?;
        Self::from_bytes(bytes)
    }

    /// Opens and validates a GNI using a read-only memory map.
    pub fn open_mmap(path: impl AsRef<Path>) -> Result<Self> {
        let file = File::open(path)?;
        // SAFETY: the file descriptor remains valid while the map is created;
        // the resulting Mmap owns the mapping for the reader's lifetime.
        let bytes = unsafe { Mmap::map(&file)? };
        Self::from_storage(IndexStorage::Mapped(bytes))
    }

    /// Parses a GNI byte stream.  This is useful for corruption tests and for
    /// callers that memory-map the file themselves before handing it to GNI.
    pub fn from_bytes(bytes: Vec<u8>) -> Result<Self> {
        Self::from_storage(IndexStorage::Owned(bytes))
    }

    fn from_storage(storage: IndexStorage) -> Result<Self> {
        let bytes = storage.as_ref();
        if bytes.len() < HEADER_SIZE {
            return Err(Error::Corrupt("truncated GNI header".into()));
        }
        if bytes[..4] != MAGIC {
            return Err(Error::Corrupt("invalid GNI magic".into()));
        }
        let mut offset = 4;
        let major_version = read_u16(bytes, &mut offset, "major version")?;
        let minor_version = read_u16(bytes, &mut offset, "minor version")?;
        if major_version != MAJOR_VERSION {
            return Err(Error::Corrupt(format!(
                "unsupported GNI major version {major_version}"
            )));
        }
        if minor_version > MINOR_VERSION {
            return Err(Error::Corrupt(format!(
                "unsupported GNI minor version {minor_version}"
            )));
        }
        let flags = read_u32(bytes, &mut offset, "flags")?;
        if flags & !1 != 0 {
            return Err(Error::Corrupt("unknown GNI flags".into()));
        }
        let byte_order = read_u8(bytes, &mut offset, "byte order")?;
        let coordinate_convention = read_u8(bytes, &mut offset, "coordinate convention")?;
        let normalization = read_u8(bytes, &mut offset, "normalization policy")?;
        let reserved = read_u8(bytes, &mut offset, "header reserved byte")?;
        if byte_order != BYTE_ORDER_LITTLE
            || coordinate_convention != COORDINATE_ZERO_BASED_HALF_OPEN
            || reserved != 0
        {
            return Err(Error::Corrupt("unsupported GNI header conventions".into()));
        }
        if normalization > NORMALIZATION_CASE_SENSITIVE {
            return Err(Error::Corrupt("unknown normalization policy".into()));
        }
        let header_size = read_u32(bytes, &mut offset, "header size")? as usize;
        let directory_entry_size = read_u32(bytes, &mut offset, "directory entry size")? as usize;
        let section_count = read_u32(bytes, &mut offset, "section count")? as usize;
        let reserved = read_u32(bytes, &mut offset, "header reserved field")?;
        if header_size != HEADER_SIZE
            || directory_entry_size != DIRECTORY_ENTRY_SIZE
            || reserved != 0
        {
            return Err(Error::Corrupt("unsupported GNI header layout".into()));
        }
        if !(1..=32).contains(&section_count) {
            return Err(Error::Corrupt("invalid section count".into()));
        }
        let term_count = read_u64(bytes, &mut offset, "term count")?;
        let unique_span_count = read_u64(bytes, &mut offset, "span count")?;
        let posting_count = read_u64(bytes, &mut offset, "posting count")?;
        let postings_block_count = read_u64(bytes, &mut offset, "postings block count")?;
        let span_block_count = read_u64(bytes, &mut offset, "span block count")?;
        let max_posting_count = term_count
            .checked_mul(unique_span_count)
            .ok_or_else(|| Error::Corrupt("posting count overflow".into()))?;
        if posting_count > max_posting_count
            || (term_count == 0 && posting_count != 0)
            || (term_count > 0 && posting_count < term_count)
            || (term_count > 0 && postings_block_count == 0)
            || (unique_span_count > 0 && span_block_count == 0)
        {
            return Err(Error::Corrupt("inconsistent GNI counts".into()));
        }
        let gff_fingerprint = read_array::<32>(bytes, &mut offset, "GFF fingerprint")?;
        let coordinate_index_fingerprint =
            read_array::<32>(bytes, &mut offset, "coordinate-index fingerprint")?;
        let reference_dictionary_fingerprint =
            read_array::<32>(bytes, &mut offset, "reference-dictionary fingerprint")?;
        let attribute_count = read_u32(bytes, &mut offset, "attribute count")? as usize;
        let spans_per_block = read_u32(bytes, &mut offset, "span block size")?;
        if spans_per_block == 0 || spans_per_block > 1_000_000 {
            return Err(Error::Corrupt("invalid span block size".into()));
        }
        let section_directory_offset = read_u64(bytes, &mut offset, "section directory offset")?;
        let section_directory_length = read_u64(bytes, &mut offset, "section directory length")?;
        let file_size = read_u64(bytes, &mut offset, "file size")?;
        let reference_count = read_u32(bytes, &mut offset, "reference count")?;
        if reference_count > 1_000_000 {
            return Err(Error::Corrupt("excessive reference count".into()));
        }
        if bytes[204..HEADER_SIZE].iter().any(|byte| *byte != 0) {
            return Err(Error::Corrupt("nonzero reserved header bytes".into()));
        }
        if file_size != bytes.len() as u64 {
            return Err(Error::Corrupt("GNI file size mismatch".into()));
        }
        let directory_end = section_directory_offset
            .checked_add(section_directory_length)
            .ok_or_else(|| Error::Corrupt("section directory overflow".into()))?;
        if section_directory_offset < HEADER_SIZE as u64
            || directory_end > bytes.len() as u64
            || section_directory_length
                != u64::from(section_count as u32)
                    .checked_mul(DIRECTORY_ENTRY_SIZE as u64)
                    .ok_or_else(|| Error::Corrupt("section directory size overflow".into()))?
        {
            return Err(Error::Corrupt("invalid section directory range".into()));
        }
        let mut sections = BTreeMap::new();
        let mut section_ranges = Vec::with_capacity(section_count);
        let mut directory_offset = section_directory_offset as usize;
        for _ in 0..section_count {
            let kind =
                SectionKind::try_from(read_u32(bytes, &mut directory_offset, "section kind")?)?;
            let flags = read_u32(bytes, &mut directory_offset, "section flags")?;
            if flags != 0 {
                return Err(Error::Corrupt("unknown section flags".into()));
            }
            let section_offset = read_u64(bytes, &mut directory_offset, "section offset")?;
            let section_length = read_u64(bytes, &mut directory_offset, "section length")?;
            let item_count = read_u64(bytes, &mut directory_offset, "section item count")?;
            let section_checksum = read_u32(bytes, &mut directory_offset, "section checksum")?;
            let reserved = read_u32(bytes, &mut directory_offset, "section reserved")?;
            if reserved != 0 || section_length > MAX_SECTION_BYTES {
                return Err(Error::Corrupt("invalid section metadata".into()));
            }
            let section_end = section_offset
                .checked_add(section_length)
                .ok_or_else(|| Error::Corrupt("section range overflow".into()))?;
            if section_offset < directory_end || section_end > bytes.len() as u64 {
                return Err(Error::Corrupt("section is outside file".into()));
            }
            section_ranges.push((section_offset, section_end));
            let section = bytes
                .get(section_offset as usize..section_end as usize)
                .ok_or_else(|| Error::Corrupt("section range is not addressable".into()))?;
            if checksum(section) != section_checksum {
                return Err(Error::Corrupt(format!(
                    "{kind:?} section checksum mismatch"
                )));
            }
            if sections
                .insert(
                    kind,
                    SectionDirectoryEntry {
                        kind,
                        flags,
                        offset: section_offset,
                        length: section_length,
                        item_count,
                        checksum: section_checksum,
                    },
                )
                .is_some()
            {
                return Err(Error::Corrupt("duplicate section directory entry".into()));
            }
        }
        section_ranges.sort_unstable();
        if section_ranges
            .windows(2)
            .any(|ranges| ranges[0].1 > ranges[1].0)
        {
            return Err(Error::Corrupt("overlapping GNI sections".into()));
        }
        let required = [
            SectionKind::Attributes,
            SectionKind::Terms,
            SectionKind::PostingsDirectory,
            SectionKind::PostingsData,
            SectionKind::SpansDirectory,
            SectionKind::SpansData,
        ];
        if required.iter().any(|kind| !sections.contains_key(kind)) {
            return Err(Error::Corrupt("missing required GNI section".into()));
        }
        let postings_data_length = sections
            .get(&SectionKind::PostingsData)
            .map(|section| section.length)
            .ok_or_else(|| Error::Corrupt("missing postings data section".into()))?;
        let spans_data_length = sections
            .get(&SectionKind::SpansData)
            .map(|section| section.length)
            .ok_or_else(|| Error::Corrupt("missing spans data section".into()))?;
        let expected_items = [
            (SectionKind::Attributes, attribute_count as u64),
            (SectionKind::Terms, term_count),
            (SectionKind::PostingsDirectory, postings_block_count),
            (SectionKind::PostingsData, postings_data_length),
            (SectionKind::SpansDirectory, span_block_count),
            (SectionKind::SpansData, spans_data_length),
        ];
        for (kind, expected) in expected_items {
            if sections
                .get(&kind)
                .is_some_and(|section| section.item_count != expected)
            {
                return Err(Error::Corrupt(format!("{kind:?} item count mismatch")));
            }
        }
        let section_bytes = |kind: SectionKind| -> Result<&[u8]> {
            let section = sections
                .get(&kind)
                .ok_or_else(|| Error::Corrupt("missing required section".into()))?;
            let end = section
                .offset
                .checked_add(section.length)
                .ok_or_else(|| Error::Corrupt("section range overflow".into()))?;
            bytes
                .get(section.offset as usize..end as usize)
                .ok_or_else(|| Error::Corrupt("section range out of bounds".into()))
        };
        let attributes = decode_attributes(section_bytes(SectionKind::Attributes)?)?;
        if attributes.len() != attribute_count || attributes.is_empty() {
            return Err(Error::Corrupt("attribute count mismatch".into()));
        }
        let terms_bytes = section_bytes(SectionKind::Terms)?;
        let term_map = fst::Map::new(terms_bytes)
            .map_err(|error| Error::Corrupt(format!("invalid term FST: {error}")))?;
        if term_map.len() as u64 != term_count {
            return Err(Error::Corrupt("term count mismatch".into()));
        }
        let posting_directory = decode_posting_directory(
            section_bytes(SectionKind::PostingsDirectory)?,
            section_bytes(SectionKind::PostingsData)?.len() as u64,
            postings_block_count,
        )?;
        let span_directory = decode_span_directory(
            section_bytes(SectionKind::SpansDirectory)?,
            section_bytes(SectionKind::SpansData)?.len() as u64,
            span_block_count,
            unique_span_count,
            spans_per_block,
            reference_count,
        )?;
        let section_length =
            |kind: SectionKind| -> u64 { sections.get(&kind).map_or(0, |section| section.length) };
        let postings_uncompressed_bytes = posting_directory
            .iter()
            .try_fold(0_u64, |total, entry| {
                total.checked_add(u64::from(entry.uncompressed_length))
            })
            .ok_or_else(|| Error::Corrupt("postings uncompressed size overflow".into()))?;
        let span_uncompressed_bytes = span_directory
            .iter()
            .try_fold(0_u64, |total, entry| {
                total.checked_add(u64::from(entry.uncompressed_length))
            })
            .ok_or_else(|| Error::Corrupt("span uncompressed size overflow".into()))?;
        let file_size = bytes.len() as u64;
        Ok(Self {
            bytes: storage,
            metadata: IndexMetadata {
                major_version,
                minor_version,
                case_sensitive: normalization == NORMALIZATION_CASE_SENSITIVE,
                attributes,
                gff_fingerprint,
                coordinate_index_fingerprint,
                reference_dictionary_fingerprint,
                term_count,
                unique_span_count,
                posting_count,
                postings_block_count,
                span_block_count,
                reference_count,
                span_block_size: spans_per_block,
                file_size,
                attribute_section_bytes: section_length(SectionKind::Attributes),
                term_dictionary_bytes: section_length(SectionKind::Terms),
                postings_directory_bytes: section_length(SectionKind::PostingsDirectory),
                postings_uncompressed_bytes,
                postings_data_bytes: section_length(SectionKind::PostingsData),
                span_directory_bytes: section_length(SectionKind::SpansDirectory),
                span_uncompressed_bytes,
                span_data_bytes: section_length(SectionKind::SpansData),
                compressed_postings_blocks: posting_directory
                    .iter()
                    .filter(|entry| entry.compression == 1)
                    .count() as u64,
                compressed_span_blocks: span_directory
                    .iter()
                    .filter(|entry| entry.compression == 1)
                    .count() as u64,
                delta_start_blocks: span_directory
                    .iter()
                    .filter(|entry| entry.start_encoding == 1)
                    .count() as u64,
                for_start_blocks: span_directory
                    .iter()
                    .filter(|entry| entry.start_encoding == 2)
                    .count() as u64,
                varint_length_blocks: span_directory
                    .iter()
                    .filter(|entry| entry.length_encoding == 0)
                    .count() as u64,
                for_length_blocks: span_directory
                    .iter()
                    .filter(|entry| entry.length_encoding == 2)
                    .count() as u64,
            },
            sections,
            posting_directory,
            span_directory,
        })
    }

    /// Returns format metadata and configured attributes.
    pub fn metadata(&self) -> &IndexMetadata {
        &self.metadata
    }

    /// Returns a copy of metadata for inspection interfaces.
    pub fn inspect(&self) -> IndexMetadata {
        self.metadata.clone()
    }

    /// Looks up a normalized term and returns its exact spans.
    pub fn lookup_spans(&self, term: &str) -> Result<Vec<Span>> {
        let normalized = normalize_value(term, self.metadata.case_sensitive);
        if normalized.is_empty() {
            return Ok(Vec::new());
        }
        let terms_section = self.section(SectionKind::Terms)?;
        let term_map = fst::Map::new(terms_section)
            .map_err(|error| Error::Corrupt(format!("invalid term FST: {error}")))?;
        let Some(locator) = term_map.get(normalized) else {
            return Ok(Vec::new());
        };
        let block_id = usize::try_from(locator >> 32)
            .map_err(|_| Error::Corrupt("posting block ID overflows usize".into()))?;
        let record_offset = usize::try_from(locator & u64::from(u32::MAX))
            .map_err(|_| Error::Corrupt("posting offset overflows usize".into()))?;
        let posting_bytes = self.decode_posting_block(block_id)?;
        let span_ids = decode_posting_record(
            &posting_bytes,
            record_offset,
            self.metadata.unique_span_count,
        )?;
        self.resolve_span_ids(&span_ids)
    }

    /// Looks up spans and reports the independently decoded bytes used by the
    /// lookup. Each span block is counted once even when several span IDs in
    /// the posting share it.
    pub fn lookup_spans_with_stats(&self, term: &str) -> Result<(Vec<Span>, LookupStats)> {
        let normalized = normalize_value(term, self.metadata.case_sensitive);
        if normalized.is_empty() {
            return Ok((Vec::new(), LookupStats::default()));
        }
        let terms_section = self.section(SectionKind::Terms)?;
        let term_map = fst::Map::new(terms_section)
            .map_err(|error| Error::Corrupt(format!("invalid term FST: {error}")))?;
        let Some(locator) = term_map.get(normalized) else {
            return Ok((Vec::new(), LookupStats::default()));
        };
        let block_id = usize::try_from(locator >> 32)
            .map_err(|_| Error::Corrupt("posting block ID overflows usize".into()))?;
        let record_offset = usize::try_from(locator & u64::from(u32::MAX))
            .map_err(|_| Error::Corrupt("posting offset overflows usize".into()))?;
        let posting_bytes = self.decode_posting_block(block_id)?;
        let span_ids = decode_posting_record(
            &posting_bytes,
            record_offset,
            self.metadata.unique_span_count,
        )?;
        let (spans, _, span_bytes_decompressed) = self.resolve_span_ids_with_stats(&span_ids)?;
        Ok((
            spans,
            LookupStats {
                postings_bytes_decompressed: posting_bytes.len() as u64,
                span_bytes_decompressed,
            },
        ))
    }

    /// Returns the sorted span IDs for an exact normalized term.
    pub fn lookup_span_ids(&self, term: &str) -> Result<Vec<u64>> {
        let normalized = normalize_value(term, self.metadata.case_sensitive);
        if normalized.is_empty() {
            return Ok(Vec::new());
        }
        let terms_section = self.section(SectionKind::Terms)?;
        let term_map = fst::Map::new(terms_section)
            .map_err(|error| Error::Corrupt(format!("invalid term FST: {error}")))?;
        let Some(locator) = term_map.get(normalized) else {
            return Ok(Vec::new());
        };
        let block_id = usize::try_from(locator >> 32)
            .map_err(|_| Error::Corrupt("posting block ID overflows usize".into()))?;
        let record_offset = usize::try_from(locator & u64::from(u32::MAX))
            .map_err(|_| Error::Corrupt("posting offset overflows usize".into()))?;
        let posting_bytes = self.decode_posting_block(block_id)?;
        decode_posting_record(
            &posting_bytes,
            record_offset,
            self.metadata.unique_span_count,
        )
    }

    /// Resolves one span ID using a binary search over fixed-width block entries.
    pub fn resolve_span_id(&self, span_id: u64) -> Result<Span> {
        self.resolve_span_ids(&[span_id])?
            .into_iter()
            .next()
            .ok_or_else(|| Error::Corrupt("missing span ID".into()))
    }

    /// Resolves a batch of span IDs while decoding each referenced span block
    /// at most once. The input and output retain the same order.
    pub fn resolve_span_ids(&self, span_ids: &[u64]) -> Result<Vec<Span>> {
        self.resolve_span_ids_with_stats(span_ids)
            .map(|(spans, _, _)| spans)
    }

    fn resolve_span_ids_with_stats(&self, span_ids: &[u64]) -> Result<(Vec<Span>, u64, u64)> {
        let mut requests = BTreeMap::<usize, Vec<(usize, usize)>>::new();
        for (output_index, &span_id) in span_ids.iter().enumerate() {
            if span_id >= self.metadata.unique_span_count {
                return Err(Error::Corrupt("invalid span ID".into()));
            }
            let block_index = self
                .span_directory
                .partition_point(|entry| entry.first_span_id <= span_id)
                .checked_sub(1)
                .ok_or_else(|| Error::Corrupt("span ID is outside span directory".into()))?;
            let entry = &self.span_directory[block_index];
            let row = usize::try_from(span_id - entry.first_span_id)
                .map_err(|_| Error::Corrupt("span row overflows usize".into()))?;
            if row >= entry.span_count as usize {
                return Err(Error::Corrupt("span ID is outside span block".into()));
            }
            requests
                .entry(block_index)
                .or_default()
                .push((output_index, row));
        }

        let mut spans = vec![None; span_ids.len()];
        let mut span_bytes_decompressed = 0_u64;
        for (block_index, rows) in requests.iter() {
            let entry = &self.span_directory[*block_index];
            let payload = self.decode_span_block(*block_index)?;
            span_bytes_decompressed = span_bytes_decompressed
                .checked_add(payload.len() as u64)
                .ok_or(Error::InvalidCoordinate)?;
            let decoded = decode_span_rows(&payload, entry)?;
            for &(output_index, row) in rows {
                let span = *decoded
                    .get(row)
                    .ok_or_else(|| Error::Corrupt("span row is out of bounds".into()))?;
                spans[output_index] = Some(span);
            }
        }

        let spans = spans
            .into_iter()
            .map(|span| span.ok_or_else(|| Error::Corrupt("missing span ID".into())))
            .collect::<Result<Vec<_>>>()?;
        Ok((spans, requests.len() as u64, span_bytes_decompressed))
    }

    fn section(&self, kind: SectionKind) -> Result<&[u8]> {
        let section = self
            .sections
            .get(&kind)
            .ok_or_else(|| Error::Corrupt("missing GNI section".into()))?;
        let end = section
            .offset
            .checked_add(section.length)
            .ok_or_else(|| Error::Corrupt("section range overflow".into()))?;
        self.bytes
            .as_ref()
            .get(section.offset as usize..end as usize)
            .ok_or_else(|| Error::Corrupt("section range out of bounds".into()))
    }

    fn decode_posting_block(&self, block_id: usize) -> Result<Vec<u8>> {
        let entry = self
            .posting_directory
            .get(block_id)
            .ok_or_else(|| Error::Corrupt("invalid posting block ID".into()))?;
        let data = self.section(SectionKind::PostingsData)?;
        let end = entry
            .compressed_offset
            .checked_add(u64::from(entry.compressed_length))
            .ok_or_else(|| Error::Corrupt("posting block range overflow".into()))?;
        let compressed = data
            .get(entry.compressed_offset as usize..end as usize)
            .ok_or_else(|| Error::Corrupt("posting block range out of bounds".into()))?;
        decompress_block(
            compressed,
            entry.compression,
            entry.uncompressed_length,
            entry.checksum,
            "posting block",
        )
    }

    fn decode_span_block(&self, block_id: usize) -> Result<Vec<u8>> {
        let entry = self
            .span_directory
            .get(block_id)
            .ok_or_else(|| Error::Corrupt("invalid span block ID".into()))?;
        let data = self.section(SectionKind::SpansData)?;
        let end = entry
            .compressed_offset
            .checked_add(u64::from(entry.compressed_length))
            .ok_or_else(|| Error::Corrupt("span block range overflow".into()))?;
        let compressed = data
            .get(entry.compressed_offset as usize..end as usize)
            .ok_or_else(|| Error::Corrupt("span block range out of bounds".into()))?;
        decompress_block(
            compressed,
            entry.compression,
            entry.uncompressed_length,
            entry.checksum,
            "span block",
        )
    }
}

fn decode_posting_record(
    bytes: &[u8],
    record_offset: usize,
    span_count_limit: u64,
) -> Result<Vec<u64>> {
    let mut offset = record_offset;
    if offset >= bytes.len() {
        return Err(Error::Corrupt(
            "posting record offset is outside block".into(),
        ));
    }
    let count = read_varint(bytes, &mut offset, "posting span count")?;
    if count > span_count_limit || count > 100_000_000 {
        return Err(Error::Corrupt("excessive posting span count".into()));
    }
    let mut ids = Vec::with_capacity(count as usize);
    let mut previous = 0_u64;
    for index in 0..count {
        let value = read_varint(bytes, &mut offset, "posting span ID")?;
        let id = if index == 0 {
            value
        } else {
            if value == 0 {
                return Err(Error::Corrupt("posting span IDs are not sorted".into()));
            }
            previous
                .checked_add(value)
                .ok_or_else(|| Error::Corrupt("posting span ID overflow".into()))?
        };
        if id >= span_count_limit {
            return Err(Error::Corrupt("posting references invalid span ID".into()));
        }
        ids.push(id);
        previous = id;
    }
    Ok(ids)
}

fn decode_posting_directory(
    bytes: &[u8],
    data_length: u64,
    expected_count: u64,
) -> Result<Vec<PostingDirectoryEntry>> {
    if expected_count > 100_000_000 {
        return Err(Error::Corrupt("excessive postings block count".into()));
    }
    if !bytes.len().is_multiple_of(POSTINGS_DIRECTORY_ENTRY_SIZE)
        || bytes.len() / POSTINGS_DIRECTORY_ENTRY_SIZE != expected_count as usize
    {
        return Err(Error::Corrupt("postings directory size mismatch".into()));
    }
    let mut offset = 0;
    let mut entries = Vec::with_capacity(expected_count as usize);
    let mut previous_end = 0_u64;
    for _ in 0..expected_count {
        let compressed_offset = read_u64(bytes, &mut offset, "posting compressed offset")?;
        let compressed_length = read_u32(bytes, &mut offset, "posting compressed length")?;
        let uncompressed_length = read_u32(bytes, &mut offset, "posting uncompressed length")?;
        let checksum = read_u32(bytes, &mut offset, "posting checksum")?;
        let compression = read_u8(bytes, &mut offset, "posting compression")?;
        let reserved = read_array::<3>(bytes, &mut offset, "posting reserved")?;
        let reserved_tail = read_u64(bytes, &mut offset, "posting reserved tail")?;
        if reserved != [0; 3] || reserved_tail != 0 || compression > 1 {
            return Err(Error::Corrupt("invalid postings directory entry".into()));
        }
        if u64::from(uncompressed_length) > MAX_BLOCK_BYTES {
            return Err(Error::Corrupt("postings block is too large".into()));
        }
        let end = compressed_offset
            .checked_add(u64::from(compressed_length))
            .ok_or_else(|| Error::Corrupt("postings block range overflow".into()))?;
        if compressed_offset < previous_end || end > data_length {
            return Err(Error::Corrupt(
                "postings block is outside data section".into(),
            ));
        }
        previous_end = end;
        entries.push(PostingDirectoryEntry {
            compressed_offset,
            compressed_length,
            uncompressed_length,
            checksum,
            compression,
        });
    }
    Ok(entries)
}

fn decode_span_directory(
    bytes: &[u8],
    data_length: u64,
    expected_count: u64,
    expected_span_count: u64,
    spans_per_block: u32,
    reference_count: u32,
) -> Result<Vec<SpanDirectoryEntry>> {
    if expected_count > 100_000_000 || expected_span_count > 100_000_000 {
        return Err(Error::Corrupt("excessive span directory count".into()));
    }
    if !bytes.len().is_multiple_of(SPAN_DIRECTORY_ENTRY_SIZE)
        || bytes.len() / SPAN_DIRECTORY_ENTRY_SIZE != expected_count as usize
    {
        return Err(Error::Corrupt("span directory size mismatch".into()));
    }
    let mut offset = 0;
    let mut entries = Vec::with_capacity(expected_count as usize);
    let mut previous_span_id = 0_u64;
    let mut previous_data_end = 0_u64;
    let mut previous_reference = None;
    let mut previous_start = 0_u64;
    let mut total_rows = 0_u64;
    for index in 0..expected_count {
        let first_span_id = read_u64(bytes, &mut offset, "span first ID")?;
        let span_count = read_u32(bytes, &mut offset, "span count")?;
        let reference_id = read_u32(bytes, &mut offset, "span reference ID")?;
        let first_start = read_u64(bytes, &mut offset, "span first start")?;
        let compressed_offset = read_u64(bytes, &mut offset, "span compressed offset")?;
        let compressed_length = read_u32(bytes, &mut offset, "span compressed length")?;
        let uncompressed_length = read_u32(bytes, &mut offset, "span uncompressed length")?;
        let checksum = read_u32(bytes, &mut offset, "span checksum")?;
        let start_encoding = read_u8(bytes, &mut offset, "span start encoding")?;
        let length_encoding = read_u8(bytes, &mut offset, "span length encoding")?;
        let compression = read_u8(bytes, &mut offset, "span compression")?;
        let reserved = read_u8(bytes, &mut offset, "span reserved")?;
        let reserved_tail = read_u32(bytes, &mut offset, "span reserved tail")?;
        if span_count == 0
            || span_count > spans_per_block
            || reference_id >= reference_count
            || start_encoding > 2
            || length_encoding > 2
            || compression > 1
            || reserved != 0
            || reserved_tail != 0
        {
            return Err(Error::Corrupt("invalid span directory entry".into()));
        }
        let expected_first_id = if index == 0 { 0 } else { previous_span_id };
        if first_span_id != expected_first_id {
            return Err(Error::Corrupt(
                "span directory IDs are not contiguous".into(),
            ));
        }
        if let Some(previous_reference) = previous_reference
            && (reference_id < previous_reference
                || (reference_id == previous_reference && first_start < previous_start))
        {
            return Err(Error::Corrupt(
                "span directory is not coordinate ordered".into(),
            ));
        }
        let end = compressed_offset
            .checked_add(u64::from(compressed_length))
            .ok_or_else(|| Error::Corrupt("span block range overflow".into()))?;
        if compressed_offset < previous_data_end || end > data_length {
            return Err(Error::Corrupt("span block is outside data section".into()));
        }
        if u64::from(uncompressed_length) > MAX_BLOCK_BYTES {
            return Err(Error::Corrupt("span block is too large".into()));
        }
        total_rows = total_rows
            .checked_add(u64::from(span_count))
            .ok_or_else(|| Error::Corrupt("span row count overflow".into()))?;
        previous_span_id = first_span_id
            .checked_add(u64::from(span_count))
            .ok_or_else(|| Error::Corrupt("span ID range overflow".into()))?;
        previous_data_end = end;
        previous_reference = Some(reference_id);
        previous_start = first_start;
        entries.push(SpanDirectoryEntry {
            first_span_id,
            span_count,
            reference_id,
            first_start,
            compressed_offset,
            compressed_length,
            uncompressed_length,
            checksum,
            start_encoding,
            length_encoding,
            compression,
        });
    }
    if total_rows != expected_span_count {
        return Err(Error::Corrupt("span row count mismatch".into()));
    }
    Ok(entries)
}

fn decode_for_stream(
    bytes: &[u8],
    offset: &mut usize,
    length: usize,
    count: usize,
    base: u64,
    bit_width: u8,
    context: &str,
) -> Result<Vec<u64>> {
    if bit_width > 63 {
        return Err(Error::Corrupt(format!("invalid {context} bit width")));
    }
    let bit_count = count
        .checked_mul(bit_width as usize)
        .ok_or_else(|| Error::Corrupt(format!("{context} bit count overflow")))?;
    let expected_length = bit_count
        .checked_add(7)
        .ok_or_else(|| Error::Corrupt(format!("{context} length overflow")))?
        / 8;
    if expected_length != length {
        return Err(Error::Corrupt(format!("{context} byte length mismatch")));
    }
    let end = offset
        .checked_add(length)
        .ok_or_else(|| Error::Corrupt(format!("{context} range overflow")))?;
    let packed = bytes
        .get(*offset..end)
        .ok_or_else(|| Error::Corrupt(format!("truncated {context}")))?;
    *offset = end;
    let mut values = Vec::with_capacity(count);
    for index in 0..count {
        let bit_offset = index
            .checked_mul(bit_width as usize)
            .ok_or_else(|| Error::Corrupt(format!("{context} bit offset overflow")))?;
        let mut relative = 0_u64;
        for bit in 0..bit_width as usize {
            let target = bit_offset + bit;
            if packed[target / 8] & (1 << (target % 8)) != 0 {
                relative |= 1_u64 << bit;
            }
        }
        values.push(
            base.checked_add(relative)
                .ok_or_else(|| Error::Corrupt(format!("{context} value overflow")))?,
        );
    }
    Ok(values)
}

fn decode_span_rows(payload: &[u8], entry: &SpanDirectoryEntry) -> Result<Vec<Span>> {
    let mut offset = 0;
    let start_length = read_u32(payload, &mut offset, "span start stream length")? as usize;
    let length_length = read_u32(payload, &mut offset, "span length stream length")? as usize;
    let length_base = read_u64(payload, &mut offset, "span length base")?;
    let start_encoding = read_u8(payload, &mut offset, "span payload start encoding")?;
    let length_encoding = read_u8(payload, &mut offset, "span payload length encoding")?;
    let start_bit_width = read_u8(payload, &mut offset, "span start bit width")?;
    let length_bit_width = read_u8(payload, &mut offset, "span length bit width")?;
    let start_base = entry.first_start;
    if start_encoding != entry.start_encoding || length_encoding != entry.length_encoding {
        return Err(Error::Corrupt("span payload metadata mismatch".into()));
    }
    if (start_encoding != 2 && start_bit_width != 0)
        || (length_encoding != 2 && length_bit_width != 0)
        || (length_encoding != 2 && length_base != 0)
    {
        return Err(Error::Corrupt("invalid span payload metadata".into()));
    }
    let count = entry.span_count as usize;
    let starts = match start_encoding {
        0 => {
            let end = offset
                .checked_add(start_length)
                .ok_or_else(|| Error::Corrupt("span start range overflow".into()))?;
            let stream = payload
                .get(offset..end)
                .ok_or_else(|| Error::Corrupt("truncated span start stream".into()))?;
            offset = end;
            let mut stream_offset = 0;
            let mut values = Vec::with_capacity(count);
            for _ in 0..count {
                values.push(read_varint(stream, &mut stream_offset, "span start")?);
            }
            if stream_offset != stream.len() {
                return Err(Error::Corrupt("trailing span start bytes".into()));
            }
            values
        }
        1 => {
            if count == 0 {
                return Err(Error::Corrupt("empty span block".into()));
            }
            let end = offset
                .checked_add(start_length)
                .ok_or_else(|| Error::Corrupt("span start range overflow".into()))?;
            let stream = payload
                .get(offset..end)
                .ok_or_else(|| Error::Corrupt("truncated span start stream".into()))?;
            offset = end;
            let mut stream_offset = 0;
            let mut values = Vec::with_capacity(count);
            values.push(start_base);
            let mut previous = start_base;
            for _ in 1..count {
                let delta = read_varint(stream, &mut stream_offset, "span start delta")?;
                let value = previous
                    .checked_add(delta)
                    .ok_or_else(|| Error::Corrupt("span start overflow".into()))?;
                values.push(value);
                previous = value;
            }
            if stream_offset != stream.len() {
                return Err(Error::Corrupt("trailing span start bytes".into()));
            }
            values
        }
        2 => {
            if count == 0 {
                return Err(Error::Corrupt("empty span block".into()));
            }
            let mut values = Vec::with_capacity(count);
            values.push(start_base);
            for relative in decode_for_stream(
                payload,
                &mut offset,
                start_length,
                count - 1,
                0,
                start_bit_width,
                "span start",
            )? {
                values.push(
                    start_base
                        .checked_add(relative)
                        .ok_or_else(|| Error::Corrupt("span start overflow".into()))?,
                );
            }
            values
        }
        _ => return Err(Error::Corrupt("unknown span start encoding".into())),
    };
    let lengths = match length_encoding {
        0 => {
            let end = offset
                .checked_add(length_length)
                .ok_or_else(|| Error::Corrupt("span length range overflow".into()))?;
            let stream = payload
                .get(offset..end)
                .ok_or_else(|| Error::Corrupt("truncated span length stream".into()))?;
            offset = end;
            let mut stream_offset = 0;
            let mut values = Vec::with_capacity(count);
            for _ in 0..count {
                values.push(read_varint(stream, &mut stream_offset, "span length")?);
            }
            if stream_offset != stream.len() {
                return Err(Error::Corrupt("trailing span length bytes".into()));
            }
            values
        }
        1 => {
            let end = offset
                .checked_add(length_length)
                .ok_or_else(|| Error::Corrupt("span length range overflow".into()))?;
            let stream = payload
                .get(offset..end)
                .ok_or_else(|| Error::Corrupt("truncated span length stream".into()))?;
            offset = end;
            let mut stream_offset = 0;
            let mut values = Vec::with_capacity(count);
            let mut previous = 0_u64;
            for index in 0..count {
                let delta = read_varint(stream, &mut stream_offset, "span length delta")?;
                let value = if index == 0 {
                    delta
                } else {
                    previous
                        .checked_add(delta)
                        .ok_or_else(|| Error::Corrupt("span length overflow".into()))?
                };
                values.push(value);
                previous = value;
            }
            if stream_offset != stream.len() {
                return Err(Error::Corrupt("trailing span length bytes".into()));
            }
            values
        }
        2 => decode_for_stream(
            payload,
            &mut offset,
            length_length,
            count,
            length_base,
            length_bit_width,
            "span length",
        )?,
        _ => return Err(Error::Corrupt("unknown span length encoding".into())),
    };
    if offset != payload.len() {
        return Err(Error::Corrupt("trailing span payload bytes".into()));
    }
    if starts.first().copied() != Some(entry.first_start)
        || starts.windows(2).any(|values| values[0] > values[1])
    {
        return Err(Error::Corrupt(
            "span starts are not coordinate ordered".into(),
        ));
    }
    for (start, length) in starts.iter().zip(&lengths) {
        if *length == 0 {
            return Err(Error::Corrupt("span length is zero".into()));
        }
        start.checked_add(*length).ok_or(Error::InvalidCoordinate)?;
    }
    Ok(starts
        .into_iter()
        .zip(lengths)
        .map(|(start, length)| Span {
            reference_id: entry.reference_id,
            start,
            length,
        })
        .collect())
}

#[cfg(test)]
fn decode_span_row(payload: &[u8], entry: &SpanDirectoryEntry, row: usize) -> Result<Span> {
    decode_span_rows(payload, entry)?
        .get(row)
        .copied()
        .ok_or_else(|| Error::Corrupt("span row is out of bounds".into()))
}

/// An indexed GFF3 source, its TBI/CSI coordinate index, and a GNI reader.
pub struct IndexedGff {
    gff_path: PathBuf,
    coordinate_index_path: PathBuf,
    coordinate_index: CoordinateIndex,
    dictionary: CoordinateDictionary,
    reference_ids: HashMap<String, u32>,
    configured_attributes: HashSet<String>,
    name_index: NameIndexReader,
}

impl IndexedGff {
    /// Opens a source, coordinate index, and GNI and rejects stale pairs by
    /// checking all source, index, and reference-dictionary fingerprints.
    pub fn open(
        gff_path: impl AsRef<Path>,
        coordinate_index_path: impl AsRef<Path>,
        gni_path: impl AsRef<Path>,
    ) -> Result<Self> {
        Self::open_inner(
            gff_path.as_ref(),
            coordinate_index_path.as_ref(),
            gni_path.as_ref(),
            None,
        )
    }

    /// Opens an indexed source while additionally checking the caller's
    /// attribute and normalization configuration against the stored header.
    pub fn open_with_options(
        gff_path: impl AsRef<Path>,
        coordinate_index_path: impl AsRef<Path>,
        gni_path: impl AsRef<Path>,
        options: &NameIndexOptions,
    ) -> Result<Self> {
        let options = NameIndexOptions::new(options.attributes.clone(), options.case_sensitive)?;
        Self::open_inner(
            gff_path.as_ref(),
            coordinate_index_path.as_ref(),
            gni_path.as_ref(),
            Some(&options),
        )
    }

    fn open_inner(
        gff_path: &Path,
        coordinate_index_path: &Path,
        gni_path: &Path,
        options: Option<&NameIndexOptions>,
    ) -> Result<Self> {
        let gff_path = gff_path.to_path_buf();
        let coordinate_index_path = coordinate_index_path.to_path_buf();
        let name_index = NameIndexReader::open_mmap(gni_path)?;
        let (coordinate_index, coordinate_index_fingerprint) =
            read_coordinate_index_with_fingerprint(&coordinate_index_path)?;
        let dictionary = coordinate_index.dictionary()?;
        if fingerprint_file(&gff_path)? != name_index.metadata.gff_fingerprint {
            return Err(Error::Stale("source GFF fingerprint does not match".into()));
        }
        if coordinate_index_fingerprint != name_index.metadata.coordinate_index_fingerprint {
            return Err(Error::Stale("TBI/CSI fingerprint does not match".into()));
        }
        if dictionary.fingerprint != name_index.metadata.reference_dictionary_fingerprint {
            return Err(Error::Stale(
                "reference dictionary fingerprint does not match".into(),
            ));
        }
        if name_index.metadata.reference_count
            != u32::try_from(dictionary.names.len())
                .map_err(|_| Error::Corrupt("reference dictionary is too large".into()))?
        {
            return Err(Error::Corrupt(
                "GNI reference count does not match coordinate index".into(),
            ));
        }
        if let Some(options) = options
            && (options.attributes != name_index.metadata.attributes
                || options.case_sensitive != name_index.metadata.case_sensitive)
        {
            return Err(Error::Stale(
                "configured attributes or normalization policy does not match GNI".into(),
            ));
        }
        for span_directory_entry in &name_index.span_directory {
            if usize::try_from(span_directory_entry.reference_id)
                .ok()
                .is_none_or(|reference_id| reference_id >= dictionary.names.len())
            {
                return Err(Error::Corrupt(
                    "span block has an invalid reference ID".into(),
                ));
            }
        }
        let reference_ids = dictionary
            .names
            .iter()
            .enumerate()
            .map(|(index, name)| (name.clone(), index as u32))
            .collect();
        let configured_attributes = name_index.metadata.attributes.iter().cloned().collect();
        Ok(Self {
            gff_path,
            coordinate_index_path,
            coordinate_index,
            dictionary,
            reference_ids,
            configured_attributes,
            name_index,
        })
    }

    /// Returns GNI metadata.
    pub fn metadata(&self) -> &IndexMetadata {
        self.name_index.metadata()
    }

    /// Queries a configured attribute value and returns matching records in
    /// source coordinate order. Each posted span is queried as an exact
    /// coordinate interval; returned chunks are then merged and read once.
    /// Unknown terms are successful empty queries.
    pub fn query_name(&mut self, term: &str) -> Result<Vec<GffRecord>> {
        self.query_name_with_stats(term).map(|(records, _)| records)
    }

    /// Queries a configured attribute value and returns records plus bounded
    /// query instrumentation. The instrumentation is useful for measuring
    /// index amplification without changing ordinary CLI output.
    pub fn query_name_with_stats(&mut self, term: &str) -> Result<(Vec<GffRecord>, QueryStats)> {
        let normalized_query = normalize_value(term, self.name_index.metadata.case_sensitive);
        if normalized_query.is_empty() {
            return Ok((Vec::new(), QueryStats::default()));
        }
        let span_ids = self.name_index.lookup_span_ids(&normalized_query)?;
        let mut stats = QueryStats {
            requested_spans: span_ids.len() as u64,
            ..QueryStats::default()
        };
        if span_ids.is_empty() {
            return Ok((Vec::new(), stats));
        }

        let (spans, distinct_span_blocks, _) =
            self.name_index.resolve_span_ids_with_stats(&span_ids)?;
        stats.distinct_span_blocks_decoded = distinct_span_blocks;
        let requested_spans = spans
            .iter()
            .map(|span| SpanKey {
                reference_id: span.reference_id,
                start: span.start,
                length: span.length,
            })
            .collect::<HashSet<_>>();
        let mut raw_chunks = Vec::new();
        for span in &spans {
            let reference_sequence_id = usize::try_from(span.reference_id)
                .map_err(|_| Error::Corrupt("span reference ID overflows usize".into()))?;
            let reference_name = self
                .dictionary
                .names
                .get(reference_sequence_id)
                .ok_or_else(|| Error::Corrupt("span references invalid reference ID".into()))?;
            let query_start = Position::try_from(
                usize::try_from(span.start.checked_add(1).ok_or(Error::InvalidCoordinate)?)
                    .map_err(|_| Error::InvalidCoordinate)?,
            )
            .map_err(|_| Error::InvalidCoordinate)?;
            let query_end = Position::try_from(
                usize::try_from(span.end()?).map_err(|_| Error::InvalidCoordinate)?,
            )
            .map_err(|_| Error::InvalidCoordinate)?;
            let region = Region::new(reference_name.as_str(), query_start..=query_end);
            let chunks = match &self.coordinate_index {
                CoordinateIndex::Tabix(index) => {
                    index.query(reference_sequence_id, region.interval())?
                }
                CoordinateIndex::Csi(index) => {
                    index.query(reference_sequence_id, region.interval())?
                }
            };
            stats.exact_interval_queries += 1;
            stats.raw_chunks += chunks.len() as u64;
            raw_chunks.extend(chunks);
        }
        let merged_chunks = merge_query_chunks(raw_chunks);
        stats.merged_chunks = merged_chunks.len() as u64;
        let context = QueryReadContext {
            requested_spans: &requested_spans,
            reference_ids: &self.reference_ids,
            configured_attributes: &self.configured_attributes,
            case_sensitive: self.name_index.metadata.case_sensitive,
            normalized_query: &normalized_query,
        };
        let records = read_query_chunks(&self.gff_path, &merged_chunks, &context, &mut stats)?;
        Ok((records, stats))
    }

    /// Returns the path of the coordinate index used for this reader.
    pub fn coordinate_index_path(&self) -> &Path {
        &self.coordinate_index_path
    }

    /// Returns the GNI reader for callers needing direct span lookup.
    pub fn name_index(&self) -> &NameIndexReader {
        &self.name_index
    }
}

fn merge_query_chunks(mut chunks: Vec<Chunk>) -> Vec<Chunk> {
    chunks.sort_unstable_by_key(|chunk| (chunk.start(), chunk.end()));
    let mut merged: Vec<Chunk> = Vec::with_capacity(chunks.len());
    for chunk in chunks {
        if let Some(last) = merged.last_mut()
            && chunk.start() <= last.end()
        {
            if chunk.end() > last.end() {
                *last = Chunk::new(last.start(), chunk.end());
            }
            continue;
        }
        merged.push(chunk);
    }
    merged
}

struct QueryReadContext<'a> {
    requested_spans: &'a HashSet<SpanKey>,
    reference_ids: &'a HashMap<String, u32>,
    configured_attributes: &'a HashSet<String>,
    case_sensitive: bool,
    normalized_query: &'a str,
}

fn read_query_chunks(
    gff_path: &Path,
    chunks: &[Chunk],
    context: &QueryReadContext<'_>,
    stats: &mut QueryStats,
) -> Result<Vec<GffRecord>> {
    let source = File::open(gff_path)?;
    let mut reader = gff::io::Reader::new(bgzf::io::Reader::new(source));
    let mut line = gff::Line::default();
    let mut positions = HashSet::new();
    let mut records = Vec::new();
    for chunk in chunks {
        reader.get_mut().seek_to_virtual_position(chunk.start())?;
        loop {
            let source_position = reader.get_ref().virtual_position();
            if source_position >= chunk.end() {
                break;
            }
            let length = reader.read_line(&mut line)?;
            if length == 0 {
                break;
            }
            stats.bytes_read = stats
                .bytes_read
                .checked_add(length as u64)
                .ok_or(Error::InvalidCoordinate)?;
            let Some(result) = line.as_record() else {
                continue;
            };
            if !positions.insert(u64::from(source_position)) {
                continue;
            }
            stats.unique_candidate_records += 1;
            let record = result
                .map_err(|error| Error::InvalidInput(format!("invalid GFF record: {error}")))?;
            let span = feature_record_span(&record, context.reference_ids)?;
            if !context.requested_spans.contains(&span) {
                continue;
            }
            if !feature_record_matches_term(
                &record,
                context.configured_attributes,
                context.case_sensitive,
                context.normalized_query,
            )? {
                continue;
            }
            let raw_line = String::from_utf8(line.as_ref().to_vec())
                .map_err(|_| Error::InvalidInput("GFF record is not UTF-8".into()))?;
            let parsed =
                parsed_record_from_feature_record(&record, raw_line, context.reference_ids)?;
            stats.matching_records += 1;
            records.push(parsed.record);
        }
    }
    Ok(records)
}

#[cfg(test)]
mod tests {
    use std::{
        collections::BTreeSet,
        sync::{Arc, Mutex},
    };

    use noodles::{
        bgzf,
        core::Position,
        csi::{
            self,
            binning_index::{
                Indexer,
                index::reference_sequence::{bin::Chunk, index::BinnedIndex},
            },
        },
        tabix,
    };
    use tempfile::tempdir;

    use super::*;

    fn fixture_records() -> [&'static str; 6] {
        [
            "##gff-version 3",
            "##sequence-region chr1 1 1000",
            "chr1\tsrc\tgene\t100\t150\t.\t+\t.\tName=BRCA1;Alias=BRCC1,RNF53;gene_name=BRCA1",
            "chr1\tsrc\tgene\t100\t150\t.\t+\t.\tAlias=BRCA1;ID=without_special_status",
            "chr1\tsrc\tgene\t200\t210\t.\t-\t.\tName=Other%20Gene",
            "chr2\tsrc\tgene\t10\t20\t.\t+\t.\tName=BRCA1",
        ]
    }

    fn write_fixture(directory: &Path) -> (PathBuf, PathBuf) {
        let source_path = directory.join("fixture.gff3.gz");
        let index_path = directory.join("fixture.gff3.gz.tbi");
        let records = fixture_records();
        let mut writer = File::create(&source_path)
            .map(bgzf::io::Writer::new)
            .expect("should create BGZF source");
        let mut indexer = tabix::index::Indexer::default();
        indexer.set_header(csi::binning_index::index::header::Builder::gff().build());
        for (record_index, line) in records.iter().enumerate() {
            if line.starts_with('#') {
                writeln!(writer, "{line}").expect("should write directive");
                continue;
            }
            let fields = line.split('\t').collect::<Vec<_>>();
            let start = Position::try_from(fields[3].parse::<usize>().unwrap()).unwrap();
            let end = Position::try_from(fields[4].parse::<usize>().unwrap()).unwrap();
            let start_position = writer.virtual_position();
            writeln!(writer, "{line}").expect("should write record");
            let end_position = writer.virtual_position();
            indexer
                .add_record(
                    fields[0],
                    start,
                    end,
                    Chunk::new(start_position, end_position),
                )
                .expect("should index record");
            assert!(record_index > 0);
        }
        writer.finish().expect("should finish BGZF");
        let index = indexer.build();
        let mut index_writer = File::create(&index_path)
            .map(tabix::io::Writer::new)
            .expect("should create TBI");
        index_writer.write_index(&index).expect("should write TBI");
        (source_path, index_path)
    }

    fn write_tbi_lines(directory: &Path, stem: &str, lines: &[&str]) -> (PathBuf, PathBuf) {
        let source_path = directory.join(format!("{stem}.gff3.gz"));
        let index_path = directory.join(format!("{stem}.gff3.gz.tbi"));
        let mut writer = File::create(&source_path)
            .map(bgzf::io::Writer::new)
            .expect("should create BGZF source");
        let mut indexer = tabix::index::Indexer::default();
        indexer.set_header(csi::binning_index::index::header::Builder::gff().build());
        for line in lines {
            if line.starts_with('#') {
                writeln!(writer, "{line}").expect("should write directive");
                continue;
            }
            let fields = line.split('\t').collect::<Vec<_>>();
            let start = Position::try_from(fields[3].parse::<usize>().unwrap()).unwrap();
            let end = Position::try_from(fields[4].parse::<usize>().unwrap()).unwrap();
            let start_position = writer.virtual_position();
            writeln!(writer, "{line}").expect("should write record");
            let end_position = writer.virtual_position();
            indexer
                .add_record(
                    fields[0],
                    start,
                    end,
                    Chunk::new(start_position, end_position),
                )
                .expect("should index record");
        }
        writer.finish().expect("should finish BGZF");
        let index = indexer.build();
        let mut index_writer = File::create(&index_path)
            .map(tabix::io::Writer::new)
            .expect("should create TBI");
        index_writer.write_index(&index).expect("should write TBI");
        drop(index_writer);
        (source_path, index_path)
    }

    fn write_csi_fixture(directory: &Path) -> (PathBuf, PathBuf) {
        let source_path = directory.join("fixture-csi.gff3.gz");
        let index_path = directory.join("fixture-csi.gff3.gz.csi");
        let records = fixture_records();
        let mut writer = File::create(&source_path)
            .map(bgzf::io::Writer::new)
            .expect("should create BGZF source");
        let mut names = csi::binning_index::index::header::ReferenceSequenceNames::new();
        names.insert("chr1".into());
        names.insert("chr2".into());
        let header = csi::binning_index::index::Header::builder()
            .set_format(Format::Generic(CoordinateSystem::Gff))
            .set_reference_sequence_names(names)
            .build();
        let mut indexer = Indexer::<BinnedIndex>::default().set_header(header);
        for line in records {
            if line.starts_with('#') {
                writeln!(writer, "{line}").expect("should write directive");
                continue;
            }
            let fields = line.split('\t').collect::<Vec<_>>();
            let reference_id = if fields[0] == "chr1" { 0 } else { 1 };
            let start = Position::try_from(fields[3].parse::<usize>().unwrap()).unwrap();
            let end = Position::try_from(fields[4].parse::<usize>().unwrap()).unwrap();
            let start_position = writer.virtual_position();
            writeln!(writer, "{line}").expect("should write record");
            let end_position = writer.virtual_position();
            indexer
                .add_record(
                    Some((reference_id, start, end, true)),
                    Chunk::new(start_position, end_position),
                )
                .expect("should index record");
        }
        writer.finish().expect("should finish BGZF");
        let index = indexer.build(2);
        let mut index_writer =
            csi::io::Writer::new(File::create(&index_path).expect("should create CSI file"));
        index_writer.write_index(&index).expect("should write CSI");
        (source_path, index_path)
    }

    fn mutate_section(bytes: &mut [u8], kind: SectionKind, mutate: impl FnOnce(&mut [u8])) {
        let directory_offset =
            usize::try_from(u64::from_le_bytes(bytes[176..184].try_into().unwrap())).unwrap();
        for index in 0..6 {
            let entry = directory_offset + index * DIRECTORY_ENTRY_SIZE;
            let entry_kind = u32::from_le_bytes(bytes[entry..entry + 4].try_into().unwrap());
            if entry_kind != kind as u32 {
                continue;
            }
            let section_offset = usize::try_from(u64::from_le_bytes(
                bytes[entry + 8..entry + 16].try_into().unwrap(),
            ))
            .unwrap();
            let section_length = usize::try_from(u64::from_le_bytes(
                bytes[entry + 16..entry + 24].try_into().unwrap(),
            ))
            .unwrap();
            mutate(&mut bytes[section_offset..section_offset + section_length]);
            let section_checksum =
                checksum(&bytes[section_offset..section_offset + section_length]);
            bytes[entry + 32..entry + 36].copy_from_slice(&section_checksum.to_le_bytes());
            return;
        }
        panic!("missing section {kind:?}");
    }

    fn mutate_section_directory_item_count(bytes: &mut [u8], kind: SectionKind, count: u64) {
        let directory_offset =
            usize::try_from(u64::from_le_bytes(bytes[176..184].try_into().unwrap())).unwrap();
        for index in 0..6 {
            let entry = directory_offset + index * DIRECTORY_ENTRY_SIZE;
            let entry_kind = u32::from_le_bytes(bytes[entry..entry + 4].try_into().unwrap());
            if entry_kind == kind as u32 {
                bytes[entry + 24..entry + 32].copy_from_slice(&count.to_le_bytes());
                return;
            }
        }
        panic!("missing section {kind:?}");
    }

    fn mutate_section_directory_offset(bytes: &mut [u8], kind: SectionKind, offset: u64) {
        let directory_offset =
            usize::try_from(u64::from_le_bytes(bytes[176..184].try_into().unwrap())).unwrap();
        for index in 0..6 {
            let entry = directory_offset + index * DIRECTORY_ENTRY_SIZE;
            let entry_kind = u32::from_le_bytes(bytes[entry..entry + 4].try_into().unwrap());
            if entry_kind == kind as u32 {
                bytes[entry + 8..entry + 16].copy_from_slice(&offset.to_le_bytes());
                return;
            }
        }
        panic!("missing section {kind:?}");
    }

    #[test]
    fn test_name_index_round_trip_and_query() {
        let directory = tempdir().expect("should create temp directory");
        let (source, coordinate_index) = write_fixture(directory.path());
        let destination = directory.path().join("fixture.gni");
        let options = NameIndexOptions::new(["Name", "Alias", "gene_name", "Name"], false)
            .expect("should validate attributes");
        let stats = build_name_index(&source, &coordinate_index, &destination, &options)
            .expect("should build GNI");
        let second_destination = directory.path().join("fixture-second.gni");
        build_name_index(&source, &coordinate_index, &second_destination, &options)
            .expect("should rebuild GNI deterministically");
        assert_eq!(
            fs::read(&destination).unwrap(),
            fs::read(&second_destination).unwrap()
        );
        assert_eq!(stats.records_processed, 4);
        assert_eq!(stats.records_indexed, 4);
        assert_eq!(stats.distinct_terms, 4);
        assert_eq!(stats.unique_spans, 3);

        let reader = NameIndexReader::open(&destination).expect("should open GNI");
        let mapped = NameIndexReader::open_mmap(&destination).expect("should mmap GNI");
        assert_eq!(mapped.lookup_span_ids("BRCA1").unwrap(), vec![0, 2]);
        assert_eq!(reader.lookup_span_ids("BRCA1").unwrap(), vec![0, 2]);
        assert_eq!(reader.lookup_span_ids("brcc1").unwrap(), vec![0]);
        assert_eq!(reader.lookup_span_ids("rnf53").unwrap(), vec![0]);
        assert_eq!(reader.lookup_span_ids("other gene").unwrap(), vec![1]);
        assert_eq!(reader.resolve_span_id(0).unwrap().start, 99);
        assert_eq!(reader.resolve_span_id(0).unwrap().length, 51);
        assert!(
            reader
                .lookup_span_ids("without_special_status")
                .unwrap()
                .is_empty()
        );

        let mut indexed = IndexedGff::open(&source, &coordinate_index, &destination)
            .expect("should open indexed source");
        let records = indexed.query_name(" BRCA1 ").expect("should query name");
        assert_eq!(records.len(), 3);
        assert_eq!(records[0].start, 100);
        assert_eq!(records[2].reference_sequence_name, "chr2");
        assert!(indexed.query_name("missing").unwrap().is_empty());
    }

    #[test]
    fn test_coordinate_and_normalization_boundaries() {
        assert_eq!(gff_to_span(100, 150).unwrap(), (99, 51));
        assert!(gff_to_span(0, 1).is_err());
        assert!(gff_to_span(2, 1).is_err());
        assert_eq!(normalize_value("  BRCA1\u{2003}", false), "brca1");
        assert_eq!(normalize_value("  BRCA1\u{2003}", true), "BRCA1");
        let directory = tempdir().expect("should create temp directory");
        let (source, coordinate_index) = write_fixture(directory.path());
        let destination = directory.path().join("case-sensitive.gni");
        let options = NameIndexOptions::new(["Name"], true).expect("should validate attributes");
        build_name_index(&source, &coordinate_index, &destination, &options)
            .expect("should build case-sensitive GNI");
        let reader = NameIndexReader::open(destination).expect("should open GNI");
        assert!(reader.lookup_span_ids("brca1").unwrap().is_empty());
        assert_eq!(reader.lookup_span_ids("BRCA1").unwrap(), vec![0, 2]);
    }

    #[test]
    fn test_corruption_is_an_error() {
        let directory = tempdir().expect("should create temp directory");
        let (source, coordinate_index) = write_fixture(directory.path());
        let destination = directory.path().join("fixture.gni");
        let options = NameIndexOptions::new(["Name"], false).expect("should validate attributes");
        build_name_index(&source, &coordinate_index, &destination, &options)
            .expect("should build GNI");
        let mut bytes = fs::read(&destination).expect("should read GNI");
        bytes[0] = b'X';
        assert!(NameIndexReader::from_bytes(bytes).is_err());
        let stale_source = directory.path().join("stale.gff3.gz");
        fs::copy(&source, &stale_source).expect("should copy source");
        let mut source_bytes = fs::read(&stale_source).expect("should read source");
        source_bytes.push(0);
        fs::write(&stale_source, source_bytes).expect("should alter source");
        assert!(matches!(
            IndexedGff::open(&stale_source, &coordinate_index, &destination),
            Err(Error::Stale(_))
        ));
        let stale_coordinate_index = directory.path().join("stale.gff3.gz.tbi");
        fs::copy(&coordinate_index, &stale_coordinate_index).expect("should copy TBI");
        let mut coordinate_bytes = fs::read(&stale_coordinate_index).unwrap();
        coordinate_bytes.push(0);
        fs::write(&stale_coordinate_index, coordinate_bytes).unwrap();
        assert!(IndexedGff::open(&source, &stale_coordinate_index, &destination).is_err());
    }

    #[test]
    fn test_varint_boundaries_and_malformed_values() {
        for value in [0, 1, 127, 128, 255, 16_384, u32::MAX as u64, u64::MAX] {
            let mut bytes = Vec::new();
            write_varint(&mut bytes, value);
            let mut offset = 0;
            assert_eq!(read_varint(&bytes, &mut offset, "test").unwrap(), value);
            assert_eq!(offset, bytes.len());
        }

        let mut offset = 0;
        assert!(read_varint(&[0x80; 10], &mut offset, "unterminated").is_err());
        let mut overflowing = vec![0xff; 9];
        overflowing.push(0x02);
        offset = 0;
        assert!(read_varint(&overflowing, &mut offset, "overflow").is_err());

        let mut term_spans = BTreeMap::new();
        term_spans.insert("term".to_string(), BTreeSet::from([1_u64, 4, 9]));
        let (fst_bytes, directory, data, _) = encode_postings(&term_spans).unwrap();
        assert_eq!(directory.len(), 1);
        assert_eq!(directory[0].compression, 0);
        assert_eq!(data, vec![3, 1, 3, 5]);
        let map = fst::Map::new(&fst_bytes).unwrap();
        let locator = map.get("term").unwrap();
        assert_eq!(
            decode_posting_record(&data, (locator & u64::from(u32::MAX)) as usize, 10).unwrap(),
            vec![1, 4, 9]
        );
    }

    #[test]
    fn test_adaptive_span_encodings_round_trip() {
        let packed_values = (0..64)
            .map(|index| SpanKey {
                reference_id: 7,
                start: 1_000_000 + index / 2,
                length: 100 + index % 3,
            })
            .collect::<Vec<_>>();
        let (packed_payload, start_encoding, length_encoding, first_start) =
            encode_span_block(&packed_values).expect("should encode packed spans");
        assert_eq!(start_encoding, 2, "clustered starts should use FOR");
        assert_eq!(length_encoding, 2, "constant lengths should use FOR");
        assert_eq!(
            u32::from_le_bytes(packed_payload[..4].try_into().unwrap()) as usize,
            ((packed_values.len() - 1) * 5).div_ceil(8),
            "the FOR start stream omits the directory-stored first row"
        );
        let packed_entry = SpanDirectoryEntry {
            first_span_id: 0,
            span_count: packed_values.len() as u32,
            reference_id: 7,
            first_start,
            compressed_offset: 0,
            compressed_length: packed_payload.len() as u32,
            uncompressed_length: packed_payload.len() as u32,
            checksum: checksum(&packed_payload),
            start_encoding,
            length_encoding,
            compression: 0,
        };
        for (row, expected) in packed_values.iter().enumerate() {
            assert_eq!(
                decode_span_row(&packed_payload, &packed_entry, row).unwrap(),
                Span::from(*expected)
            );
        }

        let varint_values = vec![
            SpanKey {
                reference_id: 7,
                start: 10,
                length: 1,
            },
            SpanKey {
                reference_id: 7,
                start: 1_u64 << 40,
                length: 1_u64 << 40,
            },
        ];
        let (varint_payload, start_encoding, length_encoding, first_start) =
            encode_span_block(&varint_values).expect("should encode sparse spans");
        assert_eq!(
            start_encoding, 2,
            "the exact FOR stream should win for this gap"
        );
        assert_eq!(length_encoding, 0, "sparse lengths should use varints");
        let varint_entry = SpanDirectoryEntry {
            first_span_id: 0,
            span_count: 2,
            reference_id: 7,
            first_start,
            compressed_offset: 0,
            compressed_length: varint_payload.len() as u32,
            uncompressed_length: varint_payload.len() as u32,
            checksum: checksum(&varint_payload),
            start_encoding,
            length_encoding,
            compression: 0,
        };
        for (row, expected) in varint_values.iter().enumerate() {
            assert_eq!(
                decode_span_row(&varint_payload, &varint_entry, row).unwrap(),
                Span::from(*expected)
            );
        }

        let values = [100_u64, 101, 103, 103];
        let (packed, base, width) = encode_for_values(&values).unwrap();
        assert_eq!(base, 100);
        assert_eq!(width, 2);
        let mut offset = 0;
        assert_eq!(
            decode_for_stream(
                &packed,
                &mut offset,
                packed.len(),
                values.len(),
                base,
                width,
                "test"
            )
            .unwrap(),
            values
        );
        assert_eq!(offset, packed.len());
        let (empty, base, width) = encode_for_values(&[0, u64::MAX]).unwrap();
        assert!(empty.is_empty());
        assert_eq!((base, width), (0, 64));
        assert!(decode_for_stream(&[], &mut 0, 0, 2, 0, 64, "overflow").is_err());
    }

    #[test]
    fn test_span_blocks_are_reference_local_and_deterministic() {
        let mut spans = Vec::new();
        for index in 0..5_000_u64 {
            spans.push(SpanKey {
                reference_id: 0,
                start: index * 10,
                length: 5,
            });
        }
        for index in 0..3_u64 {
            spans.push(SpanKey {
                reference_id: 1,
                start: index * 10,
                length: 7,
            });
        }
        let first = encode_span_blocks(&spans, 1_024).unwrap();
        let second = encode_span_blocks(&spans, 1_024).unwrap();
        assert_eq!(first.0.len(), 6);
        assert_eq!(first.0, second.0);
        assert_eq!(first.1, second.1);
        assert_eq!(
            first
                .0
                .iter()
                .map(|entry| entry.span_count)
                .collect::<Vec<_>>(),
            vec![1024, 1024, 1024, 1024, 904, 3]
        );
        assert!(first.0.windows(2).all(|entries| {
            entries[0].reference_id <= entries[1].reference_id
                && entries[0].first_span_id + u64::from(entries[0].span_count)
                    == entries[1].first_span_id
        }));
        assert_eq!(first.0[4].reference_id, 0);
        assert_eq!(first.0[5].reference_id, 1);
        assert!(encode_span_blocks(&spans, 0).is_err());
    }

    #[test]
    fn test_configured_attributes_and_disjoint_query_spans() {
        let directory = tempdir().expect("should create temp directory");
        let (source, coordinate_index) = write_fixture(directory.path());
        let destination = directory.path().join("id-only.gni");
        let options = NameIndexOptions::new(["ID"], false).unwrap();
        build_name_index(&source, &coordinate_index, &destination, &options).unwrap();
        let reader = NameIndexReader::open(&destination).unwrap();
        assert!(reader.lookup_span_ids("brca1").unwrap().is_empty());
        assert_eq!(
            reader.lookup_span_ids("without_special_status").unwrap(),
            vec![0]
        );
        let duplicate_options = NameIndexOptions::new(["Name", "Name", "Alias"], false).unwrap();
        assert_eq!(duplicate_options.attributes, vec!["Name", "Alias"]);
        assert!(NameIndexOptions::new([""], false).is_err());
        assert!(NameIndexOptions::new(std::iter::empty::<&str>(), false).is_err());

        // The two records below share a term but not a span.  This exercises
        // exact disjoint postings instead of a bounding interval.
        let source_path = directory.path().join("disjoint.gff3.gz");
        let index_path = directory.path().join("disjoint.gff3.gz.tbi");
        let lines = [
            "##gff-version 3",
            "chr1\tsrc\tgene\t10\t10\t.\t+\t.\tName=shared",
            "chr1\tsrc\tgene\t100\t100\t.\t+\t.\tName=shared",
        ];
        let mut writer = File::create(&source_path)
            .map(bgzf::io::Writer::new)
            .unwrap();
        let mut indexer = tabix::index::Indexer::default();
        indexer.set_header(csi::binning_index::index::header::Builder::gff().build());
        for line in lines {
            if line.starts_with('#') {
                writeln!(writer, "{line}").unwrap();
                continue;
            }
            let fields = line.split('\t').collect::<Vec<_>>();
            let position = Position::try_from(fields[3].parse::<usize>().unwrap()).unwrap();
            let start_position = writer.virtual_position();
            writeln!(writer, "{line}").unwrap();
            let end_position = writer.virtual_position();
            indexer
                .add_record(
                    fields[0],
                    position,
                    position,
                    Chunk::new(start_position, end_position),
                )
                .unwrap();
        }
        writer.finish().unwrap();
        let index = indexer.build();
        let mut index_writer = File::create(&index_path)
            .map(tabix::io::Writer::new)
            .unwrap();
        index_writer.write_index(&index).unwrap();
        drop(index_writer);
        let destination = directory.path().join("disjoint.gni");
        build_name_index(
            &source_path,
            &index_path,
            &destination,
            &NameIndexOptions::new(["Name"], false).unwrap(),
        )
        .unwrap();
        let reader = NameIndexReader::open(&destination).unwrap();
        assert_eq!(reader.lookup_span_ids("shared").unwrap(), vec![0, 1]);
        let mut indexed = IndexedGff::open(&source_path, &index_path, &destination).unwrap();
        let (records, stats) = indexed.query_name_with_stats("shared").unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(stats.requested_spans, 2);
        assert_eq!(stats.exact_interval_queries, 2);
        assert_eq!(stats.matching_records, 2);
        assert_eq!(stats.unique_candidate_records, 2);
    }

    #[test]
    fn test_query_exact_spans_preserves_identical_records_and_overlap_order() {
        let directory = tempdir().expect("should create temporary directory");
        let first = "chr1\tsrc\tgene\t10\t20\t.\t+\t.\tName=overlap";
        let overlapping = "chr1\tsrc\tgene\t15\t25\t.\t+\t.\tName=overlap";
        let (source, coordinate_index) = write_tbi_lines(
            directory.path(),
            "overlap",
            &["##gff-version 3", first, first, overlapping],
        );
        let destination = directory.path().join("overlap.gni");
        build_name_index(
            &source,
            &coordinate_index,
            &destination,
            &NameIndexOptions::new(["Name"], false).unwrap(),
        )
        .expect("should build overlap GNI");

        let mut indexed = IndexedGff::open(&source, &coordinate_index, &destination)
            .expect("should open overlap GNI");
        let (records, stats) = indexed
            .query_name_with_stats("overlap")
            .expect("should query overlap term");
        assert_eq!(records.len(), 3);
        assert_eq!(records[0].raw_line, first);
        assert_eq!(records[1].raw_line, first);
        assert_eq!(records[2].raw_line, overlapping);
        assert_eq!(records[0].start, 10);
        assert_eq!(records[1].start, 10);
        assert_eq!(records[2].start, 15);
        assert_eq!(stats.requested_spans, 2);
        assert_eq!(stats.exact_interval_queries, 2);
        assert_eq!(stats.distinct_span_blocks_decoded, 1);
        assert_eq!(stats.unique_candidate_records, 3);
        assert_eq!(stats.matching_records, 3);
    }

    #[test]
    fn test_query_chunk_union_and_virtual_position_deduplication() {
        let merged = merge_query_chunks(vec![
            Chunk::new(
                bgzf::VirtualPosition::from(10),
                bgzf::VirtualPosition::from(20),
            ),
            Chunk::new(
                bgzf::VirtualPosition::from(15),
                bgzf::VirtualPosition::from(25),
            ),
            Chunk::new(
                bgzf::VirtualPosition::from(30),
                bgzf::VirtualPosition::from(35),
            ),
            Chunk::new(
                bgzf::VirtualPosition::from(31),
                bgzf::VirtualPosition::from(32),
            ),
        ]);
        assert_eq!(
            merged,
            vec![
                Chunk::new(
                    bgzf::VirtualPosition::from(10),
                    bgzf::VirtualPosition::from(25)
                ),
                Chunk::new(
                    bgzf::VirtualPosition::from(30),
                    bgzf::VirtualPosition::from(35)
                ),
            ]
        );

        let directory = tempdir().expect("should create temp directory");
        let first = "chr1\tsrc\tgene\t10\t20\t.\t+\t.\tName=dedup";
        let (source, coordinate_index) = write_tbi_lines(
            directory.path(),
            "virtual-position",
            &["##gff-version 3", first, first],
        );
        let (index, _) = read_coordinate_index_with_fingerprint(&coordinate_index).unwrap();
        let spans = [SpanKey {
            reference_id: 0,
            start: 9,
            length: 11,
        }];
        let reference_ids = HashMap::from([(String::from("chr1"), 0_u32)]);
        let configured_attributes = HashSet::from([String::from("Name")]);
        let reference_name = "chr1";
        let query_start = Position::try_from(10).unwrap();
        let query_end = Position::try_from(20).unwrap();
        let region = Region::new(reference_name, query_start..=query_end);
        let chunks = match &index {
            CoordinateIndex::Tabix(index) => index.query(0, region.interval()).unwrap(),
            CoordinateIndex::Csi(index) => index.query(0, region.interval()).unwrap(),
        };
        assert!(!chunks.is_empty());
        let duplicated_chunks = [chunks[0], chunks[0]];
        let mut stats = QueryStats::default();
        let requested_spans = spans.into_iter().collect::<HashSet<_>>();
        let context = QueryReadContext {
            requested_spans: &requested_spans,
            reference_ids: &reference_ids,
            configured_attributes: &configured_attributes,
            case_sensitive: false,
            normalized_query: "dedup",
        };
        let records = read_query_chunks(&source, &duplicated_chunks, &context, &mut stats).unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(stats.unique_candidate_records, 2);
        assert_eq!(stats.matching_records, 2);
    }

    #[test]
    fn test_query_same_start_preserves_source_length_order() {
        let directory = tempdir().expect("should create temporary directory");
        let longer = "chr1\tsrc\tgene\t10\t30\t.\t+\t.\tName=same-start";
        let shorter = "chr1\tsrc\tgene\t10\t20\t.\t+\t.\tName=same-start";
        let (source, coordinate_index) = write_tbi_lines(
            directory.path(),
            "same-start",
            &["##gff-version 3", longer, shorter],
        );
        let destination = directory.path().join("same-start.gni");
        build_name_index(
            &source,
            &coordinate_index,
            &destination,
            &NameIndexOptions::new(["Name"], false).unwrap(),
        )
        .expect("should build same-start GNI");
        let mut indexed = IndexedGff::open(&source, &coordinate_index, &destination)
            .expect("should open same-start GNI");
        let records = indexed.query_name("same-start").unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].raw_line, longer);
        assert_eq!(records[1].raw_line, shorter);
    }

    #[test]
    fn test_coordinate_index_path_and_gff_header_validation() {
        let directory = tempdir().expect("should create temporary directory");
        let (source, coordinate_index) = write_fixture(directory.path());
        let renamed_index = directory.path().join("coordinates.data");
        fs::copy(&coordinate_index, &renamed_index).expect("should copy coordinate index");
        let destination = directory.path().join("renamed.gni");
        build_name_index(
            &source,
            &renamed_index,
            &destination,
            &NameIndexOptions::new(["Name"], false).unwrap(),
        )
        .expect("should parse a valid index with a nonstandard filename");
        let mut indexed = IndexedGff::open(&source, &renamed_index, &destination)
            .expect("should open renamed coordinate index");
        assert_eq!(indexed.query_name("brca1").unwrap().len(), 2);

        let non_gff_source = directory.path().join("bed-index.gff3.gz");
        let non_gff_index = directory.path().join("bed-index.data");
        let line = "chr1\tsrc\tgene\t10\t20\t.\t+\t.\tName=wrong-header";
        let mut writer = File::create(&non_gff_source)
            .map(bgzf::io::Writer::new)
            .expect("should create non-GFF source");
        let mut indexer = tabix::index::Indexer::default();
        indexer.set_header(csi::binning_index::index::header::Builder::bed().build());
        let fields = line.split('\t').collect::<Vec<_>>();
        let start = Position::try_from(10).unwrap();
        let end = Position::try_from(20).unwrap();
        let start_position = writer.virtual_position();
        writeln!(writer, "{line}").unwrap();
        let end_position = writer.virtual_position();
        indexer
            .add_record(
                fields[0],
                start,
                end,
                Chunk::new(start_position, end_position),
            )
            .unwrap();
        writer.finish().unwrap();
        let index = indexer.build();
        let mut index_writer = File::create(&non_gff_index)
            .map(tabix::io::Writer::new)
            .unwrap();
        index_writer.write_index(&index).unwrap();
        drop(index_writer);
        let non_gff_destination = directory.path().join("non-gff.gni");
        assert!(matches!(
            build_name_index(
                &non_gff_source,
                &non_gff_index,
                &non_gff_destination,
                &NameIndexOptions::new(["Name"], false).unwrap(),
            ),
            Err(Error::InvalidInput(message)) if message.contains("generic GFF")
        ));
    }

    #[test]
    fn test_attribute_and_data_section_item_counts_are_validated() {
        let directory = tempdir().expect("should create temp directory");
        let (source, coordinate_index) = write_fixture(directory.path());
        let destination = directory.path().join("counts.gni");
        build_name_index(
            &source,
            &coordinate_index,
            &destination,
            &NameIndexOptions::new(["Name"], false).unwrap(),
        )
        .unwrap();
        let original = fs::read(&destination).unwrap();

        let mut empty_attribute = original.clone();
        mutate_section(&mut empty_attribute, SectionKind::Attributes, |section| {
            section[4..8].copy_from_slice(&0_u32.to_le_bytes());
        });
        assert!(matches!(
            NameIndexReader::from_bytes(empty_attribute),
            Err(Error::Corrupt(message)) if message.contains("attribute")
        ));

        let two_attribute_destination = directory.path().join("two-attributes.gni");
        build_name_index(
            &source,
            &coordinate_index,
            &two_attribute_destination,
            &NameIndexOptions::new(["Name", "Type"], false).unwrap(),
        )
        .unwrap();
        let mut duplicate_attribute = fs::read(&two_attribute_destination).unwrap();
        mutate_section(
            &mut duplicate_attribute,
            SectionKind::Attributes,
            |section| {
                let first_length = u32::from_le_bytes(section[4..8].try_into().unwrap()) as usize;
                let second_offset = 8 + first_length;
                let second_length = u32::from_le_bytes(
                    section[second_offset..second_offset + 4]
                        .try_into()
                        .unwrap(),
                ) as usize;
                assert_eq!(first_length, second_length);
                let first = section[8..8 + first_length].to_vec();
                section[second_offset + 4..second_offset + 4 + second_length]
                    .copy_from_slice(&first);
            },
        );
        assert!(matches!(
            NameIndexReader::from_bytes(duplicate_attribute),
            Err(Error::Corrupt(message)) if message.contains("attribute")
        ));

        let mut bad_postings_items = original.clone();
        mutate_section_directory_item_count(
            &mut bad_postings_items,
            SectionKind::PostingsData,
            u64::MAX,
        );
        assert!(matches!(
            NameIndexReader::from_bytes(bad_postings_items),
            Err(Error::Corrupt(message)) if message.contains("PostingsData")
        ));

        let mut bad_spans_items = original;
        mutate_section_directory_item_count(&mut bad_spans_items, SectionKind::SpansData, u64::MAX);
        assert!(matches!(
            NameIndexReader::from_bytes(bad_spans_items),
            Err(Error::Corrupt(message)) if message.contains("SpansData")
        ));
    }

    #[test]
    fn test_feature_crossing_bgzf_block_boundary_is_queryable() {
        let directory = tempdir().expect("should create temporary directory");
        let source = directory.path().join("large.gff3.gz");
        let coordinate_index = directory.path().join("large.gff3.gz.tbi");
        let large_note = "x".repeat(70_000);
        let line = format!("chr1\tsrc\tgene\t100\t200\t.\t+\t.\tName=large;Note={large_note}");
        let mut writer = File::create(&source)
            .map(bgzf::io::Writer::new)
            .expect("should create large BGZF source");
        writeln!(writer, "##gff-version 3").unwrap();
        let start_position = writer.virtual_position();
        writeln!(writer, "{line}").unwrap();
        let end_position = writer.virtual_position();
        assert!(
            end_position.compressed() > start_position.compressed(),
            "large feature should cross BGZF blocks"
        );
        let mut indexer = tabix::index::Indexer::default();
        indexer.set_header(csi::binning_index::index::header::Builder::gff().build());
        indexer
            .add_record(
                "chr1",
                Position::try_from(100).unwrap(),
                Position::try_from(200).unwrap(),
                Chunk::new(start_position, end_position),
            )
            .unwrap();
        writer.finish().unwrap();
        let index = indexer.build();
        let mut index_writer = File::create(&coordinate_index)
            .map(tabix::io::Writer::new)
            .unwrap();
        index_writer.write_index(&index).unwrap();
        drop(index_writer);

        let destination = directory.path().join("large.gni");
        build_name_index(
            &source,
            &coordinate_index,
            &destination,
            &NameIndexOptions::new(["Name"], false).unwrap(),
        )
        .expect("should build large-feature GNI");
        let mut indexed = IndexedGff::open(&source, &coordinate_index, &destination)
            .expect("should open large-feature GNI");
        let records = indexed.query_name("large").unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].start, 100);
        assert_eq!(records[0].end, 200);
        assert_eq!(
            records[0].attribute_values("Note").next(),
            Some(large_note.as_str())
        );
    }

    #[test]
    fn test_spill_budgets_and_compression_threads_are_byte_deterministic() {
        let directory = tempdir().expect("should create temp directory");
        let (source, coordinate_index) = write_fixture(directory.path());
        let options = NameIndexOptions::new(["Name", "Alias", "gene_name"], false).unwrap();
        let first_destination = directory.path().join("deterministic-first.gni");
        let second_destination = directory.path().join("deterministic-second.gni");
        let first_options = BuildOptions::default()
            .with_memory_budget(1)
            .with_compression_threads(1)
            .with_bgzf_threads(1);
        let second_options = BuildOptions::default()
            .with_memory_budget(1024)
            .with_compression_threads(4)
            .with_bgzf_threads(2);
        let first_stats = build_name_index_with_options(
            &source,
            &coordinate_index,
            &first_destination,
            &options,
            &first_options,
        )
        .unwrap();
        let second_stats = build_name_index_with_options(
            &source,
            &coordinate_index,
            &second_destination,
            &options,
            &second_options,
        )
        .unwrap();
        assert_eq!(
            fs::read(&first_destination).unwrap(),
            fs::read(&second_destination).unwrap()
        );
        assert_eq!(first_stats.distinct_terms, second_stats.distinct_terms);
        assert_eq!(first_stats.unique_spans, second_stats.unique_spans);
        assert!(first_stats.peak_working_set_bytes >= first_options.memory_budget_bytes as u64);
        let first_phase_total = first_stats.timings.scan
            + first_stats.timings.spill
            + first_stats.timings.merge
            + first_stats.timings.encode_postings
            + first_stats.timings.encode_spans
            + first_stats.timings.serialize;
        let second_phase_total = second_stats.timings.scan
            + second_stats.timings.spill
            + second_stats.timings.merge
            + second_stats.timings.encode_postings
            + second_stats.timings.encode_spans
            + second_stats.timings.serialize;
        assert!(first_phase_total <= first_stats.timings.total);
        assert!(second_phase_total <= second_stats.timings.total);
        assert!(
            fs::read_dir(directory.path())
                .unwrap()
                .filter_map(|entry| entry.ok())
                .all(|entry| !entry.file_name().to_string_lossy().contains("-run-"))
        );
    }

    #[test]
    fn test_streaming_progress_and_fasta_fingerprint() {
        let directory = tempdir().expect("should create temp directory");
        let source = directory.path().join("fasta.gff3.gz");
        let coordinate_index = directory.path().join("fasta.gff3.gz.tbi");
        let destination = directory.path().join("fasta.gni");
        let record = "chr1\tsrc\tgene\t10\t20\t.\t+\t.\tName=fasta-term";

        let write_source = |sequence: &str| {
            let mut writer = File::create(&source)
                .map(bgzf::io::Writer::new)
                .expect("should create FASTA source");
            let mut indexer = tabix::index::Indexer::default();
            indexer.set_header(csi::binning_index::index::header::Builder::gff().build());
            writeln!(writer, "##gff-version 3").unwrap();
            let start_position = writer.virtual_position();
            writeln!(writer, "{record}").unwrap();
            let end_position = writer.virtual_position();
            indexer
                .add_record(
                    "chr1",
                    Position::try_from(10).unwrap(),
                    Position::try_from(20).unwrap(),
                    Chunk::new(start_position, end_position),
                )
                .unwrap();
            writeln!(writer, "##FASTA").unwrap();
            writeln!(writer, ">chr1").unwrap();
            writeln!(writer, "{sequence}").unwrap();
            writer.finish().unwrap();
            let index = indexer.build();
            let mut index_writer = File::create(&coordinate_index)
                .map(tabix::io::Writer::new)
                .unwrap();
            index_writer.write_index(&index).unwrap();
        };
        write_source("ACGTACGT");

        let progress = Arc::new(Mutex::new(Vec::<BuildProgress>::new()));
        let progress_for_callback = Arc::clone(&progress);
        let build_options = BuildOptions::default()
            .with_memory_budget(1)
            .with_bgzf_threads(2)
            .with_progress(move |event| {
                progress_for_callback.lock().unwrap().push(event);
            });
        let stats = build_name_index_with_options(
            &source,
            &coordinate_index,
            &destination,
            &NameIndexOptions::new(["Name"], false).unwrap(),
            &build_options,
        )
        .unwrap();
        assert_eq!(stats.records_processed, 1);
        assert_eq!(
            NameIndexReader::open(&destination)
                .unwrap()
                .lookup_span_ids("fasta-term")
                .unwrap()
                .len(),
            1
        );
        let events = progress.lock().unwrap().clone();
        assert!(events.iter().any(|event| event.phase == BuildPhase::Scan));
        assert!(
            events
                .iter()
                .any(|event| event.phase == BuildPhase::Complete)
        );
        assert!(events.windows(2).all(|window| {
            window[0].records_processed <= window[1].records_processed
                && window[0].bytes_read <= window[1].bytes_read
                && window[0].elapsed <= window[1].elapsed
        }));
        let source_bytes = fs::metadata(&source).unwrap().len();
        assert_eq!(events.last().unwrap().bytes_read, source_bytes);

        write_source("TTTTCCCC");
        assert!(matches!(
            IndexedGff::open(&source, &coordinate_index, &destination),
            Err(Error::Stale(message)) if message.contains("source GFF")
        ));
    }

    #[test]
    fn test_scan_progress_is_rate_limited() {
        let directory = tempdir().expect("should create temp directory");
        let mut storage = Vec::with_capacity(20_001);
        storage.push("##gff-version 3".to_string());
        for index in 0..20_000 {
            let position = index + 1;
            storage.push(format!(
                "chr1\tsrc\tgene\t{position}\t{position}\t.\t+\t.\tName=progress-{index}"
            ));
        }
        let lines = storage.iter().map(String::as_str).collect::<Vec<_>>();
        let (source, coordinate_index) = write_tbi_lines(directory.path(), "progress", &lines);
        let destination = directory.path().join("progress.gni");
        let progress = Arc::new(Mutex::new(Vec::<BuildProgress>::new()));
        let progress_for_callback = Arc::clone(&progress);
        let options = BuildOptions::default()
            .with_memory_budget(1 << 30)
            .with_bgzf_threads(1)
            .with_progress(move |event| {
                progress_for_callback.lock().unwrap().push(event);
            });
        let stats = build_name_index_with_options(
            &source,
            &coordinate_index,
            &destination,
            &NameIndexOptions::new(["Name"], false).unwrap(),
            &options,
        )
        .unwrap();
        assert_eq!(stats.records_processed, 20_000);
        let events = progress.lock().unwrap().clone();
        assert_eq!(
            events
                .iter()
                .filter(|event| event.phase == BuildPhase::Scan)
                .count(),
            1,
            "the final scan event should be the only scan callback below the threshold"
        );
        assert_eq!(
            events
                .iter()
                .filter(|event| event.phase == BuildPhase::Complete)
                .count(),
            1
        );
        assert!(
            events.len() <= 8,
            "unexpected progress event count: {}",
            events.len()
        );
    }

    #[test]
    fn test_spill_run_fan_in_compaction_is_equivalent() {
        let directory = tempdir().expect("should create temp directory");
        let mut storage = Vec::new();
        for index in 0..200 {
            storage.push(format!(
                "chr1\tsrc\tgene\t{}\t{}\t.\t+\t.\tName=run-{index:03}",
                index + 1,
                index + 1
            ));
        }
        let mut lines = vec!["##gff-version 3"];
        lines.extend(storage.iter().map(String::as_str));
        let (source, coordinate_index) = write_tbi_lines(directory.path(), "fan-in", &lines);
        let options = NameIndexOptions::new(["Name"], false).unwrap();
        let compacted = directory.path().join("compacted.gni");
        let baseline = directory.path().join("fan-in-baseline.gni");
        build_name_index_with_options(
            &source,
            &coordinate_index,
            &compacted,
            &options,
            &BuildOptions::default()
                .with_memory_budget(1)
                .with_compression_threads(1)
                .with_bgzf_threads(1),
        )
        .unwrap();
        build_name_index(&source, &coordinate_index, &baseline, &options).unwrap();
        assert_eq!(fs::read(&compacted).unwrap(), fs::read(&baseline).unwrap());
        assert!(
            fs::read_dir(directory.path())
                .unwrap()
                .filter_map(|entry| entry.ok())
                .all(|entry| !entry.file_name().to_string_lossy().contains("-run-"))
        );
    }

    #[test]
    fn test_plain_gff_uses_the_same_streaming_parser() {
        let directory = tempdir().expect("should create temp directory");
        let (_bgzf_source, coordinate_index) = write_fixture(directory.path());
        let source = directory.path().join("plain.gff3");
        let text = fixture_records().join("\n") + "\n";
        fs::write(&source, text).unwrap();
        let destination = directory.path().join("plain.gni");
        let stats = build_name_index_with_options(
            &source,
            &coordinate_index,
            &destination,
            &NameIndexOptions::new(["Name", "Alias"], false).unwrap(),
            &BuildOptions::default().with_bgzf_threads(1),
        )
        .unwrap();
        assert_eq!(stats.records_processed, 4);
        assert_eq!(stats.records_indexed, 4);
        assert_eq!(stats.distinct_terms, 4);
    }

    #[test]
    fn test_valid_empty_index_has_no_terms_or_spans() {
        let directory = tempdir().unwrap();
        let (source, coordinate_index) = write_fixture(directory.path());
        let destination = directory.path().join("empty.gni");
        let options = NameIndexOptions::new(["not_present"], false).unwrap();
        let stats = build_name_index(&source, &coordinate_index, &destination, &options).unwrap();
        assert_eq!(stats.distinct_terms, 0);
        assert_eq!(stats.unique_spans, 0);
        let reader = NameIndexReader::open(&destination).unwrap();
        assert_eq!(reader.metadata().term_count, 0);
        assert!(reader.lookup_span_ids("anything").unwrap().is_empty());
        let mut indexed = IndexedGff::open(&source, &coordinate_index, &destination).unwrap();
        let (records, stats) = indexed.query_name_with_stats("anything").unwrap();
        assert!(records.is_empty());
        assert_eq!(stats, QueryStats::default());
    }

    #[test]
    fn test_csi_build_query_and_index_fingerprints() {
        let directory = tempdir().expect("should create temp directory");
        let (source, coordinate_index) = write_csi_fixture(directory.path());
        let destination = directory.path().join("fixture-csi.gni");
        build_name_index(
            &source,
            &coordinate_index,
            &destination,
            &NameIndexOptions::new(["Name", "Alias"], false).unwrap(),
        )
        .unwrap();
        let mut indexed = IndexedGff::open(&source, &coordinate_index, &destination).unwrap();
        let (records, stats) = indexed.query_name_with_stats("brca1").unwrap();
        assert_eq!(records.len(), 3);
        assert_eq!(records[2].reference_sequence_name, "chr2");
        assert_eq!(stats.requested_spans, 2);
        assert_eq!(stats.exact_interval_queries, 2);
        assert_eq!(stats.matching_records, 3);
        let metadata = indexed.metadata();
        assert_eq!(
            metadata.coordinate_index_fingerprint,
            fingerprint_file(&coordinate_index).unwrap()
        );
        assert!(metadata.reference_dictionary_fingerprint != [0; 32]);
    }

    #[test]
    fn test_corruption_classes_are_rejected_or_bounded() {
        let directory = tempdir().unwrap();
        let (source, coordinate_index) = write_fixture(directory.path());
        let destination = directory.path().join("corruptions.gni");
        build_name_index(
            &source,
            &coordinate_index,
            &destination,
            &NameIndexOptions::new(["Name"], false).unwrap(),
        )
        .unwrap();
        let original = fs::read(&destination).unwrap();

        let mut truncated = original.clone();
        truncated.truncate(truncated.len() - 1);
        assert!(NameIndexReader::from_bytes(truncated).is_err());

        let mut bad_offset = original.clone();
        mutate_section_directory_offset(&mut bad_offset, SectionKind::Terms, u64::MAX);
        assert!(NameIndexReader::from_bytes(bad_offset).is_err());

        let mut bad_version = original.clone();
        bad_version[4..6].copy_from_slice(&2_u16.to_le_bytes());
        assert!(NameIndexReader::from_bytes(bad_version).is_err());

        let mut bad_count = original.clone();
        mutate_section_directory_item_count(&mut bad_count, SectionKind::Terms, u64::MAX);
        assert!(NameIndexReader::from_bytes(bad_count).is_err());

        let mut bad_checksum = original.clone();
        mutate_section(&mut bad_checksum, SectionKind::Attributes, |section| {
            section[section.len() - 1] ^= 1;
        });
        // Restoring the section checksum above models a valid directory with
        // damaged payload; an un-restored checksum is checked at open time.
        let mut untrusted_checksum = original.clone();
        mutate_section(
            &mut untrusted_checksum,
            SectionKind::Attributes,
            |section| {
                section[section.len() - 1] ^= 1;
            },
        );
        // mutate_section refreshes the checksum, so explicitly damage it for
        // the section-level checksum case.
        let directory_offset = usize::try_from(u64::from_le_bytes(
            untrusted_checksum[176..184].try_into().unwrap(),
        ))
        .unwrap();
        untrusted_checksum[directory_offset + 32] ^= 1;
        assert!(NameIndexReader::from_bytes(untrusted_checksum).is_err());
        assert!(NameIndexReader::from_bytes(bad_checksum).is_ok());

        let mut bad_postings_payload = original.clone();
        mutate_section(
            &mut bad_postings_payload,
            SectionKind::PostingsData,
            |section| {
                if !section.is_empty() {
                    section[0] ^= 0xff;
                }
            },
        );
        let posting_reader = NameIndexReader::from_bytes(bad_postings_payload).unwrap();
        assert!(posting_reader.lookup_span_ids("brca1").is_err());

        let mut bad_decompressed_size = original.clone();
        mutate_section(
            &mut bad_decompressed_size,
            SectionKind::PostingsDirectory,
            |section| {
                let current = u32::from_le_bytes(section[12..16].try_into().unwrap());
                section[12..16].copy_from_slice(&current.saturating_add(1).to_le_bytes());
            },
        );
        let size_reader = NameIndexReader::from_bytes(bad_decompressed_size).unwrap();
        assert!(size_reader.lookup_span_ids("brca1").is_err());

        let mut bad_span_ref = original.clone();
        mutate_section(&mut bad_span_ref, SectionKind::SpansDirectory, |section| {
            section[12..16].copy_from_slice(&u32::MAX.to_le_bytes());
        });
        let bad_span_destination = directory.path().join("bad-span-ref.gni");
        fs::write(&bad_span_destination, &bad_span_ref).unwrap();
        assert!(matches!(
            IndexedGff::open(&source, &coordinate_index, &bad_span_destination),
            Err(Error::Corrupt(_))
        ));

        let mut bad_reference_count = original.clone();
        bad_reference_count[200..204].copy_from_slice(&3_u32.to_le_bytes());
        let bad_reference_count_destination = directory.path().join("bad-reference-count.gni");
        fs::write(&bad_reference_count_destination, &bad_reference_count).unwrap();
        assert!(matches!(
            IndexedGff::open(&source, &coordinate_index, &bad_reference_count_destination),
            Err(Error::Corrupt(_))
        ));

        assert!(decode_posting_record(&[0x80; 10], 0, 1).is_err());
        let mut excessive_count = Vec::new();
        write_varint(&mut excessive_count, 100_000_001);
        assert!(decode_posting_record(&excessive_count, 0, 1).is_err());
        let mut invalid_span_id = Vec::new();
        write_varint(&mut invalid_span_id, 1);
        write_varint(&mut invalid_span_id, 1);
        assert!(decode_posting_record(&invalid_span_id, 0, 1).is_err());
        assert!(decode_for_stream(&[0], &mut 0, 1, 1, 0, 64, "bit width").is_err());
        assert!(
            Span {
                reference_id: 0,
                start: u64::MAX,
                length: 1,
            }
            .end()
            .is_err()
        );

        let bad_source = directory.path().join("bad.gff3");
        fs::write(&bad_source, b"not a GFF record\n").unwrap();
        let atomic_destination = directory.path().join("atomic.gni");
        fs::write(&atomic_destination, b"previous valid output").unwrap();
        assert!(
            build_name_index(
                &bad_source,
                &coordinate_index,
                &atomic_destination,
                &NameIndexOptions::new(["Name"], false).unwrap(),
            )
            .is_err()
        );
        assert_eq!(
            fs::read(&atomic_destination).unwrap(),
            b"previous valid output"
        );

        for length in 0..256_usize {
            let bytes = (0..length)
                .map(|index| (index as u8).wrapping_mul(37))
                .collect::<Vec<_>>();
            assert!(
                std::panic::catch_unwind(|| NameIndexReader::from_bytes(bytes)).is_ok(),
                "random corruption length {length} panicked"
            );
        }
    }
}
