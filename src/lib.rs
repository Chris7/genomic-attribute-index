//! GAI, the Genomic Attribute Index.
//!
//! GAI is a small, explicit binary index for configured GFF3 attribute values.
//! It is intentionally not a feature-identity index: `ID` is just another
//! attribute, and only names supplied in [`NameIndexOptions`] are searchable.
//! Coordinates in the on-disk format are zero-based, half-open `start +
//! length` tuples.  TBI/CSI remains responsible for locating source records;
//! GAI stores no BGZF virtual offsets.

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
use fst::{Automaton, Streamer};
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
use regex::{Regex, RegexBuilder};
use sha2::{Digest, Sha256};

mod index;
mod query;
mod sort;

pub use index::*;
pub use query::{IndexedSource, NameIndexReader};
#[cfg(test)]
pub(crate) use query::{
    QueryReadContext, decode_delta_start_payload, decode_for_stream, decode_length_payload,
    decode_posting_record, decode_span_row, decode_span_rows, merge_query_chunks,
    read_query_chunks,
};
pub use sort::{SortFormat, sort_bed, sort_file, sort_gff};

const MAGIC: [u8; 4] = *b"GAI\x01";
const MAJOR_VERSION: u16 = 1;
const MINOR_VERSION: u16 = 0;
const BYTE_ORDER_LITTLE: u8 = 1;
const COORDINATE_ZERO_BASED_HALF_OPEN: u8 = 1;
const NORMALIZATION_ASCII_LOWER: u8 = 0;
const NORMALIZATION_CASE_SENSITIVE: u8 = 1;
const HEADER_SIZE: usize = 256;
const DIRECTORY_ENTRY_SIZE: usize = 40;
const POSTINGS_DIRECTORY_ENTRY_SIZE: usize = 32;
const SPAN_DIRECTORY_ENTRY_SIZE: usize = 72;
const LENGTH_PAYLOAD_HEADER_SIZE: usize = 20;
const START_ENCODING_DELTA: u8 = 1;
const MAX_SECTION_BYTES: u64 = 1 << 40;
const MAX_BLOCK_BYTES: u64 = 256 << 20;
const DEFAULT_POSTINGS_BLOCK_TARGET: usize = 64 * 1024;
const DEFAULT_SPANS_PER_BLOCK: usize = 4096;
const MAX_RUN_FANIN: usize = 64;
const PROGRESS_RECORD_INTERVAL: u64 = 250_000;
const RUN_MAGIC: [u8; 8] = *b"GAIR\x01\x00\x00\x00";
static RUN_COUNTER: AtomicU64 = AtomicU64::new(0);

/// Errors returned by GAI construction, reading, and indexed querying.
#[derive(Debug, thiserror::Error)]
pub enum Error {
    /// A source or index file could not be read.
    #[error("I/O error: {0}")]
    Io(#[from] io::Error),
    /// The supplied input is malformed or inconsistent.
    #[error("invalid input: {0}")]
    InvalidInput(String),
    /// The GAI byte stream is malformed or unsafe to decode.
    #[error("corrupt GAI: {0}")]
    Corrupt(String),
    /// The source GFF or coordinate index does not match the GAI metadata.
    #[error("stale GAI: {0}")]
    Stale(String),
    /// A compressed block could not be decoded.
    #[error("compression error: {0}")]
    Compression(String),
    /// A value cannot be represented in the requested coordinate type.
    #[error("invalid coordinate")]
    InvalidCoordinate,
}

/// Controls how a normalized annotation query is matched against indexed
/// terms.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum MatchMode {
    /// Match one complete normalized configured value.
    Exact,
    /// Match every normalized value beginning with the normalized query.
    Prefix,
    /// Match normalized values containing the query as a literal substring.
    Contains,
    /// Search normalized values with a regular expression. The regex engine
    /// uses Unicode-aware case folding when the index is case-insensitive.
    Regex,
}

impl MatchMode {
    /// Parses the explicit mode names used by the CLI and Python bindings.
    pub fn parse(value: &str) -> Result<Self> {
        match value {
            "exact" => Ok(Self::Exact),
            "prefix" => Ok(Self::Prefix),
            "contains" => Ok(Self::Contains),
            "regex" => Ok(Self::Regex),
            _ => Err(Error::InvalidInput(
                "match must be one of 'exact', 'prefix', 'contains', or 'regex'".into(),
            )),
        }
    }
}

enum CompiledMatch {
    Exact(String),
    Prefix(String),
    Contains(String),
    Regex { pattern: String, regex: Regex },
}

impl CompiledMatch {
    fn new(term: &str, mode: MatchMode, case_sensitive: bool) -> Result<Self> {
        match mode {
            MatchMode::Exact => Ok(Self::Exact(normalize_value(term, case_sensitive))),
            MatchMode::Prefix => Ok(Self::Prefix(normalize_value(term, case_sensitive))),
            MatchMode::Contains => Ok(Self::Contains(normalize_value(term, case_sensitive))),
            MatchMode::Regex => {
                let pattern = term.trim_matches(char::is_whitespace).to_owned();
                let regex = RegexBuilder::new(&pattern)
                    .case_insensitive(!case_sensitive)
                    .build()
                    .map_err(|error| {
                        Error::InvalidInput(format!("invalid regex query: {error}"))
                    })?;
                Ok(Self::Regex { pattern, regex })
            }
        }
    }

    fn is_empty(&self) -> bool {
        match self {
            Self::Exact(query) | Self::Prefix(query) | Self::Contains(query) => query.is_empty(),
            Self::Regex { pattern, .. } => pattern.is_empty(),
        }
    }

    fn matches_normalized(&self, value: &str) -> bool {
        match self {
            Self::Exact(query) => value == query,
            Self::Prefix(query) => value.starts_with(query),
            Self::Contains(query) => value.contains(query),
            Self::Regex { regex, .. } => regex.is_match(value),
        }
    }
}

/// Options controlling which annotation values are indexed.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct NameIndexOptions {
    /// GFF3 attribute tags whose values should become searchable. BED input
    /// always indexes the `name` field (column 4) instead.
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

    /// Creates options for BED input, whose searchable value is always the
    /// `name` field in column 4.
    pub fn bed(case_sensitive: bool) -> Self {
        Self {
            attributes: vec!["name".to_string()],
            case_sensitive,
        }
    }
}

/// Statistics collected while building a GAI.
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
    /// Span bytes represented by fixed-width row metadata before structural encoding.
    pub span_bytes_fixed_width: u64,
    /// Span bytes after structural integer encoding.
    pub span_bytes_structural: u64,
    /// Span bytes after optional block compression.
    pub span_bytes_after_compression: u64,
    /// Delta-varint start-coordinate bytes before block compression.
    pub span_starts_bytes_before_compression: u64,
    /// Delta-varint start-coordinate bytes after block compression.
    pub span_starts_bytes_after_compression: u64,
    /// Length-coordinate bytes before block compression.
    pub span_lengths_bytes_before_compression: u64,
    /// Length-coordinate bytes after block compression.
    pub span_lengths_bytes_after_compression: u64,
    /// Number of span blocks using delta-varint starts.
    pub delta_start_blocks: u64,
    /// Number of span blocks using varint lengths.
    pub length_varint_blocks: u64,
    /// Number of span blocks using FOR lengths.
    pub length_for_blocks: u64,
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

/// A major phase reported by a GAI build.
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

/// Resource and observability controls for a GAI build.
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

/// A parsed source record returned by [`IndexedSource::query_name`].
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

/// A zero-based, half-open GAI span.
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
    pub source_fingerprint: [u8; 32],
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
    /// Starts-data section size.
    pub starts_data_bytes: u64,
    /// Lengths-data section size.
    pub lengths_data_bytes: u64,
    /// Sum of uncompressed starts block lengths.
    pub starts_uncompressed_bytes: u64,
    /// Sum of uncompressed lengths block lengths.
    pub lengths_uncompressed_bytes: u64,
    /// Number of postings blocks using zstd.
    pub compressed_postings_blocks: u64,
    /// Number of span blocks using delta-varint starts.
    pub delta_start_blocks: u64,
    /// Number of span blocks using ordinary-varint lengths.
    pub varint_length_blocks: u64,
    /// Number of span blocks using frame-of-reference lengths.
    pub for_length_blocks: u64,
    /// Number of starts blocks using zstd.
    pub compressed_start_blocks: u64,
    /// Number of lengths blocks using zstd.
    pub compressed_length_blocks: u64,
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
    /// Number of exact GAI spans requested by the term posting.
    pub requested_spans: u64,
    /// Number of distinct GAI span blocks decoded while resolving the posting.
    pub distinct_span_blocks_decoded: u64,
    /// Number of exact coordinate-index interval queries issued.
    pub exact_interval_queries: u64,
    /// Number of raw BGZF chunks returned by those interval queries.
    pub raw_chunks: u64,
    /// Number of merged, non-overlapping BGZF chunks read from the source.
    pub merged_chunks: u64,
    /// Number of unique source records parsed, keyed by BGZF virtual position.
    pub unique_candidate_records: u64,
    /// Number of source records that matched a requested span and configured
    /// value under the selected query mode.
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
    fn dictionary(&self, source_format: SortFormat) -> Result<CoordinateDictionary> {
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
        validate_coordinate_index_format(format, source_format)?;
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

fn validate_coordinate_index_format(format: Format, source_format: SortFormat) -> Result<()> {
    let expected = match source_format {
        SortFormat::Gff => CoordinateSystem::Gff,
        SortFormat::Bed => CoordinateSystem::Bed,
    };
    if format != Format::Generic(expected) {
        let source_name = match source_format {
            SortFormat::Gff => "GFF",
            SortFormat::Bed => "BED",
        };
        return Err(Error::InvalidInput(format!(
            "coordinate index is not a generic {source_name} coordinate index"
        )));
    }
    Ok(())
}

/// Applies the GAI normalization policy: trim Unicode whitespace, then
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

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
enum SectionKind {
    Attributes = 1,
    Terms = 2,
    PostingsDirectory = 3,
    PostingsData = 4,
    SpansDirectory = 5,
    StartsData = 6,
    LengthsData = 7,
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
            6 => Ok(Self::StartsData),
            7 => Ok(Self::LengthsData),
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
    starts_compressed_offset: u64,
    starts_compressed_length: u32,
    starts_uncompressed_length: u32,
    starts_checksum: u32,
    lengths_compressed_offset: u64,
    lengths_compressed_length: u32,
    lengths_uncompressed_length: u32,
    lengths_checksum: u32,
    start_encoding: u8,
    length_encoding: u8,
    starts_compression: u8,
    lengths_compression: u8,
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
        let record = encode_delta_posting_record(spans)?;
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

fn encode_delta_posting_record(spans: &[u64]) -> Result<Vec<u8>> {
    let count = u64::try_from(spans.len())
        .map_err(|_| Error::InvalidInput("too many spans for a term".into()))?;
    if count == 0 {
        return Err(Error::InvalidInput(
            "posting spans must not be empty".into(),
        ));
    }
    let mut record = Vec::new();
    write_varint(&mut record, count);
    let mut previous = 0_u64;
    for (index, span_id) in spans.iter().copied().enumerate() {
        if index == 0 {
            write_varint(&mut record, span_id);
        } else {
            let delta = span_id
                .checked_sub(previous)
                .ok_or_else(|| Error::InvalidInput("posting spans are not sorted".into()))?;
            if delta == 0 {
                return Err(Error::InvalidInput(
                    "posting spans are not deduplicated".into(),
                ));
            }
            write_varint(&mut record, delta);
        }
        previous = span_id;
    }
    Ok(record)
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

/// Encodes reference-local, nondecreasing starts as canonical unsigned
/// LEB128 deltas. The first absolute start is stored in the block directory,
/// so this payload contains exactly one delta for every row after the first.
/// Equal starts therefore encode as a zero delta, and a single-row block has
/// an empty payload.
fn encode_delta_start_payload(spans: &[SpanKey]) -> Result<Vec<u8>> {
    if spans.is_empty() {
        return Err(Error::InvalidInput(
            "cannot encode an empty span block".into(),
        ));
    }
    let mut payload = Vec::new();
    let mut previous = spans[0].start;
    for span in spans.iter().skip(1) {
        let delta = span
            .start
            .checked_sub(previous)
            .ok_or_else(|| Error::InvalidInput("span starts are not nondecreasing".into()))?;
        write_varint(&mut payload, delta);
        previous = span.start;
    }
    if payload.len() as u64 > MAX_BLOCK_BYTES {
        return Err(Error::InvalidInput(
            "delta start payload is too large".into(),
        ));
    }
    Ok(payload)
}

fn canonical_varint_length(value: u64) -> usize {
    if value == 0 {
        1
    } else {
        (64 - value.leading_zeros()).div_ceil(7) as usize
    }
}

fn read_canonical_varint(bytes: &[u8], offset: &mut usize, context: &str) -> Result<u64> {
    let start = *offset;
    let value = read_varint(bytes, offset, context)?;
    if offset.saturating_sub(start) != canonical_varint_length(value) {
        return Err(Error::Corrupt(format!("noncanonical varint in {context}")));
    }
    Ok(value)
}

fn encode_length_payload(spans: &[SpanKey]) -> Result<(Vec<u8>, u8)> {
    let lengths = spans.iter().map(|span| span.length).collect::<Vec<_>>();
    let mut varints = Vec::new();
    for length in &lengths {
        write_varint(&mut varints, *length);
    }
    let (for_stream, base, bit_width) = encode_for_values(&lengths)?;
    let (encoding, stream, base, bit_width) = if bit_width <= 63 && for_stream.len() < varints.len()
    {
        (2_u8, for_stream, base, bit_width)
    } else {
        (0_u8, varints, 0_u64, 0_u8)
    };
    let mut payload = Vec::with_capacity(
        LENGTH_PAYLOAD_HEADER_SIZE
            .checked_add(stream.len())
            .ok_or(Error::InvalidCoordinate)?,
    );
    put_u32(
        &mut payload,
        u32::try_from(spans.len()).map_err(|_| Error::InvalidInput("too many span rows".into()))?,
    );
    payload.push(encoding);
    payload.push(bit_width);
    payload.extend_from_slice(&[0; 2]);
    put_u64(&mut payload, base);
    put_u32(
        &mut payload,
        u32::try_from(stream.len())
            .map_err(|_| Error::InvalidInput("length stream is too large".into()))?,
    );
    payload.extend_from_slice(&stream);
    if payload.len() as u64 > MAX_BLOCK_BYTES {
        return Err(Error::InvalidInput("length payload is too large".into()));
    }
    Ok((payload, encoding))
}

#[derive(Clone, Copy, Debug, Default)]
struct SpanEncodingStats {
    fixed_width_bytes: u64,
    starts_structural_bytes: u64,
    lengths_structural_bytes: u64,
    starts_compressed_bytes: u64,
    lengths_compressed_bytes: u64,
    delta_start_blocks: u64,
    length_varint_blocks: u64,
    length_for_blocks: u64,
}

#[allow(clippy::type_complexity)]
fn encode_span_blocks_with_threads(
    spans: &[SpanKey],
    spans_per_block: usize,
    compression_threads: usize,
) -> Result<(Vec<SpanDirectoryEntry>, Vec<u8>, Vec<u8>, SpanEncodingStats)> {
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
                let starts_payload = encode_delta_start_payload(block)?;
                let (lengths_payload, length_encoding) = encode_length_payload(block)?;
                let starts_encoded = compress_block(&starts_payload)?;
                let lengths_encoded = compress_block(&lengths_payload)?;
                Ok::<_, Error>((
                    block.len(),
                    block[0].reference_id,
                    block[0].start,
                    length_encoding,
                    starts_payload.len() as u64,
                    lengths_payload.len() as u64,
                    starts_encoded,
                    lengths_encoded,
                ))
            })
            .collect::<Result<Vec<_>>>()
    })?;

    let mut directory = Vec::with_capacity(encoded_blocks.len());
    let mut starts_data = Vec::new();
    let mut lengths_data = Vec::new();
    let mut stats = SpanEncodingStats::default();
    let mut first_span_id = 0_u64;
    for (
        span_count,
        reference_id,
        first_start,
        length_encoding,
        starts_payload_length,
        lengths_payload_length,
        starts_encoded,
        lengths_encoded,
    ) in encoded_blocks
    {
        let starts_compressed_offset = u64::try_from(starts_data.len())
            .map_err(|_| Error::InvalidInput("starts data exceeds 4 GiB".into()))?;
        starts_data.extend_from_slice(&starts_encoded.compressed);
        let lengths_compressed_offset = u64::try_from(lengths_data.len())
            .map_err(|_| Error::InvalidInput("lengths data exceeds 4 GiB".into()))?;
        lengths_data.extend_from_slice(&lengths_encoded.compressed);
        let span_count = u32::try_from(span_count)
            .map_err(|_| Error::InvalidInput("span block has too many rows".into()))?;
        directory.push(SpanDirectoryEntry {
            first_span_id,
            span_count,
            reference_id,
            first_start,
            starts_compressed_offset,
            starts_compressed_length: u32::try_from(starts_encoded.compressed.len())
                .map_err(|_| Error::InvalidInput("compressed starts block exceeds 4 GiB".into()))?,
            starts_uncompressed_length: starts_encoded.uncompressed_length,
            starts_checksum: starts_encoded.checksum,
            lengths_compressed_offset,
            lengths_compressed_length: u32::try_from(lengths_encoded.compressed.len()).map_err(
                |_| Error::InvalidInput("compressed lengths block exceeds 4 GiB".into()),
            )?,
            lengths_uncompressed_length: lengths_encoded.uncompressed_length,
            lengths_checksum: lengths_encoded.checksum,
            start_encoding: START_ENCODING_DELTA,
            length_encoding,
            starts_compression: starts_encoded.compression,
            lengths_compression: lengths_encoded.compression,
        });
        stats.fixed_width_bytes = stats
            .fixed_width_bytes
            .checked_add(u64::from(span_count).saturating_mul(20))
            .ok_or(Error::InvalidCoordinate)?;
        stats.starts_structural_bytes = stats
            .starts_structural_bytes
            .checked_add(starts_payload_length)
            .ok_or(Error::InvalidCoordinate)?;
        stats.lengths_structural_bytes = stats
            .lengths_structural_bytes
            .checked_add(lengths_payload_length)
            .ok_or(Error::InvalidCoordinate)?;
        stats.starts_compressed_bytes = stats
            .starts_compressed_bytes
            .checked_add(starts_encoded.compressed.len() as u64)
            .ok_or(Error::InvalidCoordinate)?;
        stats.lengths_compressed_bytes = stats
            .lengths_compressed_bytes
            .checked_add(lengths_encoded.compressed.len() as u64)
            .ok_or(Error::InvalidCoordinate)?;
        stats.delta_start_blocks += 1;
        if length_encoding == 0 {
            stats.length_varint_blocks += 1;
        } else if length_encoding == 2 {
            stats.length_for_blocks += 1;
        }
        first_span_id = first_span_id
            .checked_add(u64::from(span_count))
            .ok_or(Error::InvalidCoordinate)?;
    }
    Ok((directory, starts_data, lengths_data, stats))
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
        put_u64(&mut bytes, entry.starts_compressed_offset);
        put_u32(&mut bytes, entry.starts_compressed_length);
        put_u32(&mut bytes, entry.starts_uncompressed_length);
        put_u32(&mut bytes, entry.starts_checksum);
        put_u64(&mut bytes, entry.lengths_compressed_offset);
        put_u32(&mut bytes, entry.lengths_compressed_length);
        put_u32(&mut bytes, entry.lengths_uncompressed_length);
        put_u32(&mut bytes, entry.lengths_checksum);
        bytes.push(entry.start_encoding);
        bytes.push(entry.length_encoding);
        bytes.push(entry.starts_compression);
        bytes.push(entry.lengths_compression);
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
    source_fingerprint: [u8; 32],
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
    starts_data: Vec<u8>,
    lengths_data: Vec<u8>,
    spans_per_block: usize,
) -> Result<Vec<u8>> {
    let attributes_section = encode_attributes(attributes)?;
    let postings_data_length = postings_data.len() as u64;
    let starts_data_length = starts_data.len() as u64;
    let lengths_data_length = lengths_data.len() as u64;
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
        (SectionKind::StartsData, starts_data, starts_data_length),
        (SectionKind::LengthsData, lengths_data, lengths_data_length),
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
            .map_err(|_| Error::InvalidInput("GAI exceeds 4 GiB".into()))?;
        if section.len() as u64 > MAX_SECTION_BYTES {
            return Err(Error::InvalidInput(
                "GAI section exceeds safety limit".into(),
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
    output[72..104].copy_from_slice(&source_fingerprint);
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
        u64::try_from(output.len()).map_err(|_| Error::InvalidInput("GAI exceeds 4 GiB".into()))?;
    set_u64(&mut output, 192, output_len);
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
            return Err(Error::InvalidInput("invalid GAI spill run magic".into()));
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
            return Err(Error::InvalidInput("invalid GAI spill run magic".into()));
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

    fn write_bed_fixture(directory: &Path) -> (PathBuf, PathBuf) {
        let source_path = directory.join("fixture.bed.gz");
        let index_path = directory.join("fixture.bed.gz.tbi");
        let records = [
            "chr1\t0\t10\tAlpha",
            "chr1\t20\t25\tBeta",
            "chr2\t5\t7\tAlpha",
            "chr2\t8\t9",
        ];
        let mut writer = File::create(&source_path)
            .map(bgzf::io::Writer::new)
            .expect("should create BGZF BED source");
        let mut indexer = tabix::index::Indexer::default();
        indexer.set_header(csi::binning_index::index::header::Builder::bed().build());
        for line in records {
            let fields = line.split('\t').collect::<Vec<_>>();
            let bed_start = fields[1].parse::<usize>().unwrap();
            let bed_end = fields[2].parse::<usize>().unwrap();
            let start = Position::try_from(bed_start + 1).unwrap();
            let end = Position::try_from(bed_end).unwrap();
            let start_position = writer.virtual_position();
            writeln!(writer, "{line}").expect("should write BED record");
            let end_position = writer.virtual_position();
            indexer
                .add_record(
                    fields[0],
                    start,
                    end,
                    Chunk::new(start_position, end_position),
                )
                .expect("should index BED record");
        }
        writer.finish().expect("should finish BGZF BED source");
        let index = indexer.build();
        let mut index_writer = File::create(&index_path)
            .map(tabix::io::Writer::new)
            .expect("should create BED TBI");
        index_writer
            .write_index(&index)
            .expect("should write BED TBI");
        (source_path, index_path)
    }

    fn mutate_section(bytes: &mut [u8], kind: SectionKind, mutate: impl FnOnce(&mut [u8])) {
        let directory_offset =
            usize::try_from(u64::from_le_bytes(bytes[176..184].try_into().unwrap())).unwrap();
        let section_count = u32::from_le_bytes(bytes[24..28].try_into().unwrap()) as usize;
        for index in 0..section_count {
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
        let section_count = u32::from_le_bytes(bytes[24..28].try_into().unwrap()) as usize;
        for index in 0..section_count {
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
        let section_count = u32::from_le_bytes(bytes[24..28].try_into().unwrap()) as usize;
        for index in 0..section_count {
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
        let destination = directory.path().join("fixture.gai");
        let options = NameIndexOptions::new(["Name", "Alias", "gene_name", "Name"], false)
            .expect("should validate attributes");
        let stats = build_name_index(&source, &coordinate_index, &destination, &options)
            .expect("should build GAI");
        let second_destination = directory.path().join("fixture-second.gai");
        build_name_index(&source, &coordinate_index, &second_destination, &options)
            .expect("should rebuild GAI deterministically");
        assert_eq!(
            fs::read(&destination).unwrap(),
            fs::read(&second_destination).unwrap()
        );
        assert_eq!(stats.records_processed, 4);
        assert_eq!(stats.records_indexed, 4);
        assert_eq!(stats.distinct_terms, 4);
        assert_eq!(stats.unique_spans, 3);
        let reader = NameIndexReader::open(&destination).expect("should open GAI");
        assert_eq!(reader.metadata().minor_version, 0);
        let mapped = NameIndexReader::open_mmap(&destination).expect("should mmap GAI");
        assert_eq!(mapped.lookup_span_ids("BRCA1").unwrap(), vec![0, 2]);
        assert_eq!(reader.lookup_span_ids("BRCA1").unwrap(), vec![0, 2]);
        assert_eq!(reader.lookup_span_ids("brcc1").unwrap(), vec![0]);
        assert_eq!(reader.lookup_span_ids("rnf53").unwrap(), vec![0]);
        assert_eq!(reader.lookup_span_ids("other gene").unwrap(), vec![1]);
        assert_eq!(
            reader
                .lookup_span_ids_with_mode("br", MatchMode::Prefix)
                .unwrap(),
            vec![0, 2]
        );
        assert_eq!(
            reader
                .lookup_span_ids_with_mode("brcc", MatchMode::Prefix)
                .unwrap(),
            vec![0]
        );
        assert_eq!(
            reader
                .lookup_span_ids_with_mode(" RCA ", MatchMode::Contains)
                .unwrap(),
            vec![0, 2]
        );
        assert_eq!(
            reader
                .lookup_span_ids_with_mode("^brc(a1|c1)$", MatchMode::Regex)
                .unwrap(),
            vec![0, 2]
        );
        assert_eq!(
            reader
                .lookup_span_ids_with_mode(r"^brc[a-z][0-9]$", MatchMode::Regex)
                .unwrap(),
            vec![0, 2]
        );
        assert_eq!(
            reader
                .lookup_span_ids_with_mode(r"\S", MatchMode::Regex)
                .unwrap(),
            vec![0, 1, 2]
        );
        assert!(
            reader
                .lookup_span_ids_with_mode("", MatchMode::Contains)
                .unwrap()
                .is_empty()
        );
        assert!(matches!(
            reader.lookup_span_ids_with_mode("missing[", MatchMode::Regex),
            Err(Error::InvalidInput(message)) if message.contains("regex")
        ));
        assert_eq!(reader.resolve_span_id(0).unwrap().start, 99);
        assert_eq!(reader.resolve_span_id(0).unwrap().length, 51);
        assert!(
            reader
                .lookup_span_ids("without_special_status")
                .unwrap()
                .is_empty()
        );

        let mut indexed = IndexedSource::open(&source, &coordinate_index, &destination)
            .expect("should open indexed source");
        let records = indexed.query_name(" BRCA1 ").expect("should query name");
        assert_eq!(records.len(), 3);
        assert_eq!(records[0].start, 100);
        assert_eq!(records[2].reference_sequence_name, "chr2");
        let (prefix_records, prefix_stats) = indexed
            .query_name_with_mode_and_stats(" br ", MatchMode::Prefix)
            .expect("should query prefix");
        assert_eq!(prefix_records.len(), 3);
        assert_eq!(prefix_stats.requested_spans, 2);
        assert_eq!(prefix_stats.exact_interval_queries, 2);
        assert_eq!(prefix_stats.matching_records, 3);
        assert_eq!(
            prefix_records
                .iter()
                .map(|record| record.raw_line.as_str())
                .collect::<Vec<_>>(),
            vec![
                fixture_records()[2],
                fixture_records()[3],
                fixture_records()[5]
            ]
        );
        let (contains_records, contains_stats) = indexed
            .query_name_with_mode_and_stats(" RCA ", MatchMode::Contains)
            .expect("should query literal substring");
        assert_eq!(contains_records.len(), 3);
        assert_eq!(contains_stats.requested_spans, 2);
        assert_eq!(contains_stats.exact_interval_queries, 2);
        assert_eq!(contains_stats.matching_records, 3);
        let (regex_records, regex_stats) = indexed
            .query_name_with_mode_and_stats("^BR C[A-Z]$", MatchMode::Regex)
            .expect("should query regex");
        assert!(regex_records.is_empty());
        assert_eq!(regex_stats.requested_spans, 0);
        let (regex_records, regex_stats) = indexed
            .query_name_with_mode_and_stats("^BRC(A1|C1)$", MatchMode::Regex)
            .expect("should query anchored regex alternation");
        assert_eq!(regex_records.len(), 3);
        assert_eq!(regex_stats.requested_spans, 2);
        assert_eq!(regex_stats.matching_records, 3);
        let (escape_records, _) = indexed
            .query_name_with_mode_and_stats(r"\S", MatchMode::Regex)
            .expect("should preserve uppercase regex escape");
        assert_eq!(escape_records.len(), 4);
        assert!(indexed.query_name("missing").unwrap().is_empty());
        assert!(
            indexed
                .query_name_with_mode("missing", MatchMode::Prefix)
                .unwrap()
                .is_empty()
        );
    }

    #[test]
    fn test_bed_name_index_builds_from_column_four() {
        let directory = tempdir().expect("should create temp directory");
        let (source, coordinate_index) = write_bed_fixture(directory.path());
        let destination = directory.path().join("fixture-bed.gai");

        // BED ignores configured GFF attributes and always indexes column 4 (`name`).
        let options = NameIndexOptions::new(["DefinitelyIgnored"], false).unwrap();
        let stats = build_name_index(&source, &coordinate_index, &destination, &options)
            .expect("should build BED GAI");
        assert_eq!(stats.records_processed, 4);
        assert_eq!(stats.records_indexed, 3);
        assert_eq!(stats.distinct_terms, 2);
        assert_eq!(stats.unique_spans, 3);

        let reader = NameIndexReader::open(&destination).expect("should open BED GAI");
        assert_eq!(reader.metadata().attributes, vec!["name"]);
        assert_eq!(reader.lookup_span_ids("alpha").unwrap(), vec![0, 2]);
        assert_eq!(reader.lookup_span_ids("beta").unwrap(), vec![1]);
        assert!(
            reader
                .lookup_span_ids("DefinitelyIgnored")
                .unwrap()
                .is_empty()
        );
        assert_eq!(reader.resolve_span_id(0).unwrap().start, 0);
        assert_eq!(reader.resolve_span_id(0).unwrap().length, 10);
        assert_eq!(reader.resolve_span_id(1).unwrap().start, 20);
        assert_eq!(reader.resolve_span_id(1).unwrap().length, 5);

        let case_destination = directory.path().join("fixture-bed-case.gai");
        let case_options = NameIndexOptions::new(["Ignored"], true).unwrap();
        build_name_index(&source, &coordinate_index, &case_destination, &case_options)
            .expect("should build case-sensitive BED GAI");
        let case_reader = NameIndexReader::open(&case_destination).unwrap();
        assert_eq!(case_reader.lookup_span_ids("Alpha").unwrap(), vec![0, 2]);
        assert!(case_reader.lookup_span_ids("alpha").unwrap().is_empty());

        let mut indexed = IndexedSource::open(&source, &coordinate_index, &destination)
            .expect("should open BED indexed source");
        let records = indexed.query_name("alpha").expect("should query BED name");
        assert_eq!(
            records
                .iter()
                .map(|record| record.raw_line.as_str())
                .collect::<Vec<_>>(),
            vec!["chr1\t0\t10\tAlpha", "chr2\t5\t7\tAlpha"]
        );
        assert_eq!(records[0].start, 1);
        assert_eq!(records[0].end, 10);
        assert_eq!(
            records[0].attribute_values("name").collect::<Vec<_>>(),
            vec!["Alpha"]
        );
        let prefix = indexed
            .query_name_with_mode("al", MatchMode::Prefix)
            .expect("should query BED name prefix");
        assert_eq!(prefix.len(), 2);

        let (contains, contains_stats) = indexed
            .query_name_with_mode_and_stats(" ph ", MatchMode::Contains)
            .expect("should query BED name substring");
        assert_eq!(contains.len(), 2);
        assert_eq!(contains_stats.requested_spans, 2);
        assert_eq!(contains_stats.matching_records, 2);

        let (regex, regex_stats) = indexed
            .query_name_with_mode_and_stats("^(ALPHA|BETA)$", MatchMode::Regex)
            .expect("should query BED name regex");
        assert_eq!(regex.len(), 3);
        assert_eq!(regex_stats.matching_records, 3);
        let (class_regex, _) = indexed
            .query_name_with_mode_and_stats("^ALP[A-Z]+$", MatchMode::Regex)
            .expect("should query BED regex character class");
        assert_eq!(class_regex.len(), 2);
        let (escape_regex, _) = indexed
            .query_name_with_mode_and_stats(r"\S", MatchMode::Regex)
            .expect("should preserve BED regex uppercase escape");
        assert_eq!(escape_regex.len(), 3);
        assert!(
            indexed
                .query_name_with_mode("^missing$", MatchMode::Regex)
                .unwrap()
                .is_empty()
        );
        assert!(matches!(
            indexed.query_name_with_mode("missing[", MatchMode::Regex),
            Err(Error::InvalidInput(message)) if message.contains("regex")
        ));

        let mut case_indexed = IndexedSource::open(&source, &coordinate_index, &case_destination)
            .expect("should open case-sensitive BED GAI");
        assert_eq!(
            case_indexed
                .query_name_with_mode("^Alpha$", MatchMode::Regex)
                .unwrap()
                .len(),
            2
        );
        assert!(
            case_indexed
                .query_name_with_mode("^alpha$", MatchMode::Regex)
                .unwrap()
                .is_empty()
        );
        assert!(
            case_indexed
                .query_name_with_mode("PH", MatchMode::Contains)
                .unwrap()
                .is_empty()
        );
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
        let destination = directory.path().join("case-sensitive.gai");
        let options = NameIndexOptions::new(["Name"], true).expect("should validate attributes");
        build_name_index(&source, &coordinate_index, &destination, &options)
            .expect("should build case-sensitive GAI");
        let reader = NameIndexReader::open(destination).expect("should open GAI");
        assert!(reader.lookup_span_ids("brca1").unwrap().is_empty());
        assert_eq!(reader.lookup_span_ids("BRCA1").unwrap(), vec![0, 2]);
        assert!(
            reader
                .lookup_span_ids_with_mode("br", MatchMode::Prefix)
                .unwrap()
                .is_empty()
        );
        assert_eq!(
            reader
                .lookup_span_ids_with_mode("BR", MatchMode::Prefix)
                .unwrap(),
            vec![0, 2]
        );
        assert!(
            reader
                .lookup_span_ids_with_mode("^brca1$", MatchMode::Regex)
                .unwrap()
                .is_empty()
        );
        assert_eq!(
            reader
                .lookup_span_ids_with_mode("^BRCA1$", MatchMode::Regex)
                .unwrap(),
            vec![0, 2]
        );
    }

    #[test]
    fn test_corruption_is_an_error() {
        let directory = tempdir().expect("should create temp directory");
        let (source, coordinate_index) = write_fixture(directory.path());
        let destination = directory.path().join("fixture.gai");
        let options = NameIndexOptions::new(["Name"], false).expect("should validate attributes");
        build_name_index(&source, &coordinate_index, &destination, &options)
            .expect("should build GAI");
        let mut bytes = fs::read(&destination).expect("should read GAI");
        assert_eq!(&bytes[..4], b"GAI\x01");
        let mut legacy_magic = bytes.clone();
        legacy_magic[..4].copy_from_slice(&[b'G', b'N', b'I', 1]);
        assert!(NameIndexReader::from_bytes(legacy_magic).is_err());
        bytes[0] = b'X';
        assert!(NameIndexReader::from_bytes(bytes).is_err());
        let stale_source = directory.path().join("stale.gff3.gz");
        fs::copy(&source, &stale_source).expect("should copy source");
        let mut source_bytes = fs::read(&stale_source).expect("should read source");
        source_bytes.push(0);
        fs::write(&stale_source, source_bytes).expect("should alter source");
        assert!(matches!(
            IndexedSource::open(&stale_source, &coordinate_index, &destination),
            Err(Error::Stale(_))
        ));
        let stale_coordinate_index = directory.path().join("stale.gff3.gz.tbi");
        fs::copy(&coordinate_index, &stale_coordinate_index).expect("should copy TBI");
        let mut coordinate_bytes = fs::read(&stale_coordinate_index).unwrap();
        coordinate_bytes.push(0);
        fs::write(&stale_coordinate_index, coordinate_bytes).unwrap();
        assert!(IndexedSource::open(&source, &stale_coordinate_index, &destination).is_err());
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
        let example = [10_u64, 10, 15, 20]
            .into_iter()
            .map(|start| SpanKey {
                reference_id: 0,
                start,
                length: 1,
            })
            .collect::<Vec<_>>();
        assert_eq!(example[0].start, 10);
        assert_eq!(encode_delta_start_payload(&example).unwrap(), vec![0, 5, 5]);

        let packed_values = (0..64)
            .map(|index| SpanKey {
                reference_id: 7,
                start: 1_000_000 + index / 2,
                length: 100 + index % 3,
            })
            .collect::<Vec<_>>();
        let packed_starts = encode_delta_start_payload(&packed_values).unwrap();
        let (packed_lengths, packed_length_encoding) =
            encode_length_payload(&packed_values).unwrap();
        assert_eq!(
            packed_length_encoding, 2,
            "clustered lengths should use FOR"
        );
        let packed_entry = SpanDirectoryEntry {
            first_span_id: 0,
            span_count: packed_values.len() as u32,
            reference_id: 7,
            first_start: packed_values[0].start,
            starts_compressed_offset: 0,
            starts_compressed_length: packed_starts.len() as u32,
            starts_uncompressed_length: packed_starts.len() as u32,
            starts_checksum: checksum(&packed_starts),
            lengths_compressed_offset: 0,
            lengths_compressed_length: packed_lengths.len() as u32,
            lengths_uncompressed_length: packed_lengths.len() as u32,
            lengths_checksum: checksum(&packed_lengths),
            start_encoding: START_ENCODING_DELTA,
            length_encoding: packed_length_encoding,
            starts_compression: 0,
            lengths_compression: 0,
        };
        for (row, expected) in packed_values.iter().enumerate() {
            assert_eq!(
                decode_span_row(&packed_starts, &packed_lengths, &packed_entry, row).unwrap(),
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
        let varint_starts = encode_delta_start_payload(&varint_values).unwrap();
        let (varint_lengths, varint_length_encoding) =
            encode_length_payload(&varint_values).unwrap();
        assert_eq!(
            varint_length_encoding, 0,
            "sparse lengths should use varints"
        );
        let varint_entry = SpanDirectoryEntry {
            first_span_id: 0,
            span_count: 2,
            reference_id: 7,
            first_start: varint_values[0].start,
            starts_compressed_offset: 0,
            starts_compressed_length: varint_starts.len() as u32,
            starts_uncompressed_length: varint_starts.len() as u32,
            starts_checksum: checksum(&varint_starts),
            lengths_compressed_offset: 0,
            lengths_compressed_length: varint_lengths.len() as u32,
            lengths_uncompressed_length: varint_lengths.len() as u32,
            lengths_checksum: checksum(&varint_lengths),
            start_encoding: START_ENCODING_DELTA,
            length_encoding: varint_length_encoding,
            starts_compression: 0,
            lengths_compression: 0,
        };
        for (row, expected) in varint_values.iter().enumerate() {
            assert_eq!(
                decode_span_row(&varint_starts, &varint_lengths, &varint_entry, row).unwrap(),
                Span::from(*expected)
            );
        }

        let equal_start_values = (0..4)
            .map(|index| SpanKey {
                reference_id: 7,
                start: 77,
                length: index + 1,
            })
            .collect::<Vec<_>>();
        let equal_starts = encode_delta_start_payload(&equal_start_values).unwrap();
        let equal_entry = SpanDirectoryEntry {
            first_span_id: 0,
            span_count: 4,
            reference_id: 7,
            first_start: 77,
            starts_compressed_offset: 0,
            starts_compressed_length: equal_starts.len() as u32,
            starts_uncompressed_length: equal_starts.len() as u32,
            starts_checksum: checksum(&equal_starts),
            lengths_compressed_offset: 0,
            lengths_compressed_length: 0,
            lengths_uncompressed_length: 0,
            lengths_checksum: 0,
            start_encoding: START_ENCODING_DELTA,
            length_encoding: 0,
            starts_compression: 0,
            lengths_compression: 0,
        };
        let (equal_lengths, equal_encoding) = encode_length_payload(&equal_start_values).unwrap();
        let equal_entry = SpanDirectoryEntry {
            lengths_compressed_length: equal_lengths.len() as u32,
            lengths_uncompressed_length: equal_lengths.len() as u32,
            lengths_checksum: checksum(&equal_lengths),
            length_encoding: equal_encoding,
            ..equal_entry
        };
        assert_eq!(
            decode_span_rows(&equal_starts, &equal_lengths, &equal_entry)
                .unwrap()
                .iter()
                .map(|span| span.start)
                .collect::<Vec<_>>(),
            vec![77; 4]
        );

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
    fn test_delta_start_payload_rejects_malformed_payloads_without_unbounded_work() {
        let values = [
            SpanKey {
                reference_id: 2,
                start: 100,
                length: 5,
            },
            SpanKey {
                reference_id: 2,
                start: 101,
                length: 6,
            },
            SpanKey {
                reference_id: 2,
                start: 103,
                length: 7,
            },
            SpanKey {
                reference_id: 2,
                start: 110,
                length: 8,
            },
        ];
        let starts = encode_delta_start_payload(&values).unwrap();
        let (lengths, length_encoding) = encode_length_payload(&values).unwrap();
        let entry = SpanDirectoryEntry {
            first_span_id: 0,
            span_count: values.len() as u32,
            reference_id: 2,
            first_start: 100,
            starts_compressed_offset: 0,
            starts_compressed_length: starts.len() as u32,
            starts_uncompressed_length: starts.len() as u32,
            starts_checksum: checksum(&starts),
            lengths_compressed_offset: 0,
            lengths_compressed_length: lengths.len() as u32,
            lengths_uncompressed_length: lengths.len() as u32,
            lengths_checksum: checksum(&lengths),
            start_encoding: START_ENCODING_DELTA,
            length_encoding,
            starts_compression: 0,
            lengths_compression: 0,
        };
        assert!(decode_delta_start_payload(&starts, &entry).is_ok());

        let mut truncated = starts.clone();
        truncated.pop();
        assert!(decode_delta_start_payload(&truncated, &entry).is_err());

        let mut trailing = starts.clone();
        trailing.push(0);
        assert!(decode_delta_start_payload(&trailing, &entry).is_err());

        // One-byte deltas are required to use their minimal LEB128 form.
        let noncanonical = vec![0x81, 0x00, 0x02, 0x07];
        assert!(decode_delta_start_payload(&noncanonical, &entry).is_err());
        assert!(decode_delta_start_payload(&[0x80], &entry).is_err());

        let mut count_mismatch = entry.clone();
        count_mismatch.span_count = 3;
        assert!(decode_delta_start_payload(&starts, &count_mismatch).is_err());

        let mut overflowing_start = entry.clone();
        overflowing_start.first_start = u64::MAX;
        assert!(decode_delta_start_payload(&starts, &overflowing_start).is_err());

        let mut bad_lengths = lengths.clone();
        bad_lengths[0..4].copy_from_slice(&u32::MAX.to_le_bytes());
        assert!(decode_length_payload(&bad_lengths, &entry).is_err());
        assert!(decode_span_rows(&starts, &lengths, &entry).is_ok());

        let single = [SpanKey {
            reference_id: 2,
            start: 42,
            length: 9,
        }];
        let single_starts = encode_delta_start_payload(&single).unwrap();
        assert!(single_starts.is_empty());
        let (single_lengths, single_encoding) = encode_length_payload(&single).unwrap();
        let single_entry = SpanDirectoryEntry {
            first_span_id: 0,
            span_count: 1,
            reference_id: 2,
            first_start: 42,
            starts_compressed_offset: 0,
            starts_compressed_length: single_starts.len() as u32,
            starts_uncompressed_length: single_starts.len() as u32,
            starts_checksum: checksum(&single_starts),
            lengths_compressed_offset: 0,
            lengths_compressed_length: single_lengths.len() as u32,
            lengths_uncompressed_length: single_lengths.len() as u32,
            lengths_checksum: checksum(&single_lengths),
            start_encoding: START_ENCODING_DELTA,
            length_encoding: single_encoding,
            starts_compression: 0,
            lengths_compression: 0,
        };
        assert_eq!(
            decode_span_rows(&single_starts, &single_lengths, &single_entry).unwrap(),
            vec![Span::from(single[0])]
        );
        assert!(decode_delta_start_payload(&[0], &single_entry).is_err());
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
        let first = encode_span_blocks_with_threads(&spans, 1_024, 1).unwrap();
        let second = encode_span_blocks_with_threads(&spans, 1_024, 1).unwrap();
        assert_eq!(first.0.len(), 6);
        assert_eq!(first.0, second.0);
        assert_eq!(first.1, second.1);
        assert_eq!(first.2, second.2);
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
        assert!(encode_span_blocks_with_threads(&spans, 0, 1).is_err());
    }

    #[test]
    fn test_configured_attributes_and_disjoint_query_spans() {
        let directory = tempdir().expect("should create temp directory");
        let (source, coordinate_index) = write_fixture(directory.path());
        let destination = directory.path().join("id-only.gai");
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
        let destination = directory.path().join("disjoint.gai");
        build_name_index(
            &source_path,
            &index_path,
            &destination,
            &NameIndexOptions::new(["Name"], false).unwrap(),
        )
        .unwrap();
        let reader = NameIndexReader::open(&destination).unwrap();
        assert_eq!(reader.lookup_span_ids("shared").unwrap(), vec![0, 1]);
        let mut indexed = IndexedSource::open(&source_path, &index_path, &destination).unwrap();
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
        let destination = directory.path().join("overlap.gai");
        build_name_index(
            &source,
            &coordinate_index,
            &destination,
            &NameIndexOptions::new(["Name"], false).unwrap(),
        )
        .expect("should build overlap GAI");

        let mut indexed = IndexedSource::open(&source, &coordinate_index, &destination)
            .expect("should open overlap GAI");
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
    fn test_contains_filters_nonmatching_records_at_a_shared_span() {
        let directory = tempdir().expect("should create temporary directory");
        let matching = "chr1\tsrc\tgene\t10\t20\t.\t+\t.\tName=Alpha";
        let nonmatching = "chr1\tsrc\tgene\t10\t20\t.\t+\t.\tName=Beta";
        let (source, coordinate_index) = write_tbi_lines(
            directory.path(),
            "contains-shared-span",
            &["##gff-version 3", matching, nonmatching],
        );
        let destination = directory.path().join("contains-shared-span.gai");
        build_name_index(
            &source,
            &coordinate_index,
            &destination,
            &NameIndexOptions::new(["Name"], false).unwrap(),
        )
        .expect("should build shared-span GAI");

        let mut indexed = IndexedSource::open(&source, &coordinate_index, &destination)
            .expect("should open shared-span GAI");
        let (records, stats) = indexed
            .query_name_with_mode_and_stats("ph", MatchMode::Contains)
            .expect("should query interior substring");
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].raw_line, matching);
        assert_eq!(stats.requested_spans, 1);
        assert_eq!(stats.distinct_span_blocks_decoded, 1);
        assert_eq!(stats.exact_interval_queries, 1);
        assert_eq!(stats.unique_candidate_records, 2);
        assert_eq!(stats.matching_records, 1);
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
        let matcher = CompiledMatch::new("dedup", MatchMode::Exact, false).unwrap();
        let context = QueryReadContext {
            requested_spans: &requested_spans,
            reference_ids: &reference_ids,
            configured_terms: &configured_attributes,
            case_sensitive: false,
            matcher: &matcher,
        };
        let records = read_query_chunks(
            &source,
            SortFormat::Gff,
            &duplicated_chunks,
            &context,
            &mut stats,
        )
        .unwrap();
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
        let destination = directory.path().join("same-start.gai");
        build_name_index(
            &source,
            &coordinate_index,
            &destination,
            &NameIndexOptions::new(["Name"], false).unwrap(),
        )
        .expect("should build same-start GAI");
        let mut indexed = IndexedSource::open(&source, &coordinate_index, &destination)
            .expect("should open same-start GAI");
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
        let destination = directory.path().join("renamed.gai");
        build_name_index(
            &source,
            &renamed_index,
            &destination,
            &NameIndexOptions::new(["Name"], false).unwrap(),
        )
        .expect("should parse a valid index with a nonstandard filename");
        let mut indexed = IndexedSource::open(&source, &renamed_index, &destination)
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
        let non_gff_destination = directory.path().join("non-gff.gai");
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
        let destination = directory.path().join("counts.gai");
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

        let two_attribute_destination = directory.path().join("two-attributes.gai");
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
        mutate_section_directory_item_count(
            &mut bad_spans_items,
            SectionKind::StartsData,
            u64::MAX,
        );
        assert!(matches!(
            NameIndexReader::from_bytes(bad_spans_items),
            Err(Error::Corrupt(message)) if message.contains("StartsData")
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

        let destination = directory.path().join("large.gai");
        build_name_index(
            &source,
            &coordinate_index,
            &destination,
            &NameIndexOptions::new(["Name"], false).unwrap(),
        )
        .expect("should build large-feature GAI");
        let mut indexed = IndexedSource::open(&source, &coordinate_index, &destination)
            .expect("should open large-feature GAI");
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
        let first_destination = directory.path().join("deterministic-first.gai");
        let second_destination = directory.path().join("deterministic-second.gai");
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
        let destination = directory.path().join("fasta.gai");
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
            IndexedSource::open(&source, &coordinate_index, &destination),
            Err(Error::Stale(message)) if message.contains("source fingerprint")
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
        let destination = directory.path().join("progress.gai");
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
        let compacted = directory.path().join("compacted.gai");
        let baseline = directory.path().join("fan-in-baseline.gai");
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
        let destination = directory.path().join("plain.gai");
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
        let destination = directory.path().join("empty.gai");
        let options = NameIndexOptions::new(["not_present"], false).unwrap();
        let stats = build_name_index(&source, &coordinate_index, &destination, &options).unwrap();
        assert_eq!(stats.distinct_terms, 0);
        assert_eq!(stats.unique_spans, 0);
        let reader = NameIndexReader::open(&destination).unwrap();
        assert_eq!(reader.metadata().term_count, 0);
        assert!(reader.lookup_span_ids("anything").unwrap().is_empty());
        let mut indexed = IndexedSource::open(&source, &coordinate_index, &destination).unwrap();
        let (records, stats) = indexed.query_name_with_stats("anything").unwrap();
        assert!(records.is_empty());
        assert_eq!(stats, QueryStats::default());
    }

    #[test]
    fn test_csi_build_query_and_index_fingerprints() {
        let directory = tempdir().expect("should create temp directory");
        let (source, coordinate_index) = write_csi_fixture(directory.path());
        let destination = directory.path().join("fixture-csi.gai");
        build_name_index(
            &source,
            &coordinate_index,
            &destination,
            &NameIndexOptions::new(["Name", "Alias"], false).unwrap(),
        )
        .unwrap();
        let mut indexed = IndexedSource::open(&source, &coordinate_index, &destination).unwrap();
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
        let destination = directory.path().join("corruptions.gai");
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

        for component in [SectionKind::StartsData, SectionKind::LengthsData] {
            let mut damaged_component = original.clone();
            mutate_section(&mut damaged_component, component, |section| {
                assert!(!section.is_empty());
                section[0] ^= 0xff;
            });
            let component_reader = NameIndexReader::from_bytes(damaged_component).unwrap();
            assert!(component_reader.resolve_span_id(0).is_err());
        }

        let mut bad_span_ref = original.clone();
        mutate_section(&mut bad_span_ref, SectionKind::SpansDirectory, |section| {
            section[12..16].copy_from_slice(&u32::MAX.to_le_bytes());
        });
        let bad_span_destination = directory.path().join("bad-span-ref.gai");
        fs::write(&bad_span_destination, &bad_span_ref).unwrap();
        assert!(matches!(
            IndexedSource::open(&source, &coordinate_index, &bad_span_destination),
            Err(Error::Corrupt(_))
        ));

        let mut bad_reference_count = original.clone();
        bad_reference_count[200..204].copy_from_slice(&3_u32.to_le_bytes());
        let bad_reference_count_destination = directory.path().join("bad-reference-count.gai");
        fs::write(&bad_reference_count_destination, &bad_reference_count).unwrap();
        assert!(matches!(
            IndexedSource::open(&source, &coordinate_index, &bad_reference_count_destination),
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
        let atomic_destination = directory.path().join("atomic.gai");
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
