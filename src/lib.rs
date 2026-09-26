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
    #[tracing::instrument(level = "trace", skip_all)]
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
    #[tracing::instrument(level = "trace", skip_all)]
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

    #[tracing::instrument(level = "trace", skip_all)]
    fn is_empty(&self) -> bool {
        match self {
            Self::Exact(query) | Self::Prefix(query) | Self::Contains(query) => query.is_empty(),
            Self::Regex { pattern, .. } => pattern.is_empty(),
        }
    }

    #[tracing::instrument(level = "trace", skip_all)]
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
    #[tracing::instrument(level = "trace", skip_all)]
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
    #[tracing::instrument(level = "trace", skip_all)]
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
    #[tracing::instrument(level = "trace", skip_all)]
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
    #[tracing::instrument(level = "trace", skip_all)]
    pub fn with_memory_budget(mut self, bytes: usize) -> Self {
        self.memory_budget_bytes = bytes;
        self
    }

    /// Sets the number of deterministic block-compression workers.
    #[tracing::instrument(level = "trace", skip_all)]
    pub fn with_compression_threads(mut self, threads: usize) -> Self {
        self.compression_threads = threads;
        self
    }

    /// Sets the BGZF decompression worker count.
    #[tracing::instrument(level = "trace", skip_all)]
    pub fn with_bgzf_threads(mut self, threads: usize) -> Self {
        self.bgzf_threads = threads;
        self
    }

    /// Sets a progress callback shared by the build operation.
    #[tracing::instrument(level = "trace", skip_all)]
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
    #[tracing::instrument(level = "trace", skip_all)]
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
    #[tracing::instrument(level = "trace", skip_all)]
    pub fn end(self) -> Result<u64> {
        self.start
            .checked_add(self.length)
            .ok_or(Error::InvalidCoordinate)
    }

    /// Converts one-based inclusive GFF coordinates at the parser boundary.
    #[tracing::instrument(level = "trace", skip_all)]
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
#[tracing::instrument(level = "trace", skip_all)]
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
    #[tracing::instrument(level = "trace", skip_all)]
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
    #[tracing::instrument(level = "trace", skip_all)]
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
#[tracing::instrument(level = "trace", skip_all)]
fn fingerprint_reference_dictionary(names: &[String]) -> [u8; 32] {
    let mut hasher = Sha256::new();
    for name in names {
        hasher.update((name.len() as u64).to_le_bytes());
        hasher.update(name.as_bytes());
    }
    hasher.finalize().into()
}
#[tracing::instrument(level = "trace", skip_all)]
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
    #[tracing::instrument(level = "trace", skip_all)]
    fn with_counter(inner: R, bytes_read: Arc<AtomicU64>) -> Self {
        Self {
            inner,
            hasher: Sha256::new(),
            bytes_read,
        }
    }

    #[tracing::instrument(level = "trace", skip_all)]
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
    #[tracing::instrument(level = "trace", skip_all)]
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
    #[tracing::instrument(level = "trace", skip_all)]
    fn drain_and_finish(mut self) -> io::Result<([u8; 32], u64)> {
        io::copy(&mut self, &mut io::sink())?;
        Ok(self.into_inner().finish())
    }
}

impl HashingInput for bgzf::io::Reader<HashingReader<File>> {
    #[tracing::instrument(level = "trace", skip_all)]
    fn drain_and_finish(mut self) -> io::Result<([u8; 32], u64)> {
        io::copy(&mut self, &mut io::sink())?;
        Ok(self.into_inner().finish())
    }
}

impl HashingInput for bgzf::io::MultithreadedReader<HashingReader<File>> {
    #[tracing::instrument(level = "trace", skip_all)]
    fn drain_and_finish(mut self) -> io::Result<([u8; 32], u64)> {
        let drain_result = io::copy(&mut self, &mut io::sink());
        let inner = self.finish()?;
        drain_result?;
        Ok(inner.finish())
    }
}
#[tracing::instrument(level = "trace", skip_all)]
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
#[tracing::instrument(level = "trace", skip_all)]
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
#[tracing::instrument(level = "trace", skip_all)]
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

    #[tracing::instrument(level = "trace", skip_all)]
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
#[tracing::instrument(level = "trace", skip_all)]
fn put_u32(buffer: &mut Vec<u8>, value: u32) {
    buffer.extend_from_slice(&value.to_le_bytes());
}
#[tracing::instrument(level = "trace", skip_all)]
fn put_u64(buffer: &mut Vec<u8>, value: u64) {
    buffer.extend_from_slice(&value.to_le_bytes());
}
#[tracing::instrument(level = "trace", skip_all)]
fn read_u8(bytes: &[u8], offset: &mut usize, context: &str) -> Result<u8> {
    let value = *bytes
        .get(*offset)
        .ok_or_else(|| Error::Corrupt(format!("truncated {context}")))?;
    *offset += 1;
    Ok(value)
}
#[tracing::instrument(level = "trace", skip_all)]
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
#[tracing::instrument(level = "trace", skip_all)]
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
#[tracing::instrument(level = "trace", skip_all)]
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
#[tracing::instrument(level = "trace", skip_all)]
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
#[tracing::instrument(level = "trace", skip_all)]
fn write_varint(buffer: &mut Vec<u8>, mut value: u64) {
    while value >= 0x80 {
        buffer.push((value as u8) | 0x80);
        value >>= 7;
    }
    buffer.push(value as u8);
}
#[tracing::instrument(level = "trace", skip_all)]
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
#[tracing::instrument(level = "trace", skip_all)]
fn checksum(bytes: &[u8]) -> u32 {
    crc32fast::hash(bytes)
}
#[tracing::instrument(level = "trace", skip_all)]
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
#[tracing::instrument(level = "trace", skip_all)]
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
#[tracing::instrument(level = "trace", skip_all)]
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
#[tracing::instrument(level = "trace", skip_all)]
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
    #[tracing::instrument(level = "trace", skip_all)]
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

    #[tracing::instrument(level = "trace", skip_all)]
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

    #[tracing::instrument(level = "trace", skip_all)]
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
#[tracing::instrument(level = "trace", skip_all)]
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
#[tracing::instrument(level = "trace", skip_all)]
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
#[tracing::instrument(level = "trace", skip_all)]
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
#[tracing::instrument(level = "trace", skip_all)]
fn encode_postings(
    term_spans: &BTreeMap<String, BTreeSet<u64>>,
) -> Result<(Vec<u8>, Vec<PostingDirectoryEntry>, Vec<u8>, u64)> {
    let groups = term_spans
        .iter()
        .map(|(term, spans)| (term.clone(), spans.iter().copied().collect::<Vec<_>>()));
    let (terms, directory, data, before, _, _) = encode_postings_groups(groups, 1)?;
    Ok((terms, directory, data, before))
}
#[tracing::instrument(level = "trace", skip_all)]
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
#[tracing::instrument(level = "trace", skip_all)]
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
#[tracing::instrument(level = "trace", skip_all)]
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
#[tracing::instrument(level = "trace", skip_all)]
fn canonical_varint_length(value: u64) -> usize {
    if value == 0 {
        1
    } else {
        (64 - value.leading_zeros()).div_ceil(7) as usize
    }
}
#[tracing::instrument(level = "trace", skip_all)]
fn read_canonical_varint(bytes: &[u8], offset: &mut usize, context: &str) -> Result<u64> {
    let start = *offset;
    let value = read_varint(bytes, offset, context)?;
    if offset.saturating_sub(start) != canonical_varint_length(value) {
        return Err(Error::Corrupt(format!("noncanonical varint in {context}")));
    }
    Ok(value)
}
#[tracing::instrument(level = "trace", skip_all)]
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
#[tracing::instrument(level = "trace", skip_all)]
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
#[tracing::instrument(level = "trace", skip_all)]
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
#[tracing::instrument(level = "trace", skip_all)]
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
#[tracing::instrument(level = "trace", skip_all)]
fn set_u16(bytes: &mut [u8], offset: usize, value: u16) {
    bytes[offset..offset + 2].copy_from_slice(&value.to_le_bytes());
}
#[tracing::instrument(level = "trace", skip_all)]
fn set_u32(bytes: &mut [u8], offset: usize, value: u32) {
    bytes[offset..offset + 4].copy_from_slice(&value.to_le_bytes());
}
#[tracing::instrument(level = "trace", skip_all)]
fn set_u64(bytes: &mut [u8], offset: usize, value: u64) {
    bytes[offset..offset + 8].copy_from_slice(&value.to_le_bytes());
}
#[tracing::instrument(level = "trace", skip_all)]
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
#[tracing::instrument(level = "trace", skip_all)]
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
#[tracing::instrument(level = "trace", skip_all)]
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
    #[tracing::instrument(level = "trace", skip_all)]
    fn drop(&mut self) {
        if let Some(path) = self.path.take() {
            let _ = fs::remove_file(path);
        }
    }
}

impl SpillCollector {
    #[tracing::instrument(level = "trace", skip_all)]
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

    #[tracing::instrument(level = "trace", skip_all)]
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

    #[tracing::instrument(level = "trace", skip_all)]
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

    #[tracing::instrument(level = "trace", skip_all)]
    fn finish(&mut self) -> Result<()> {
        self.spill().map(|_| ())
    }

    #[tracing::instrument(level = "trace", skip_all)]
    fn run_paths(&self) -> &[PathBuf] {
        &self.run_paths
    }

    #[tracing::instrument(level = "trace", skip_all)]
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

    #[tracing::instrument(level = "trace", skip_all)]
    fn peak_working_set_bytes(&self) -> u64 {
        self.peak_working_set_bytes as u64
    }
}

impl Drop for SpillCollector {
    #[tracing::instrument(level = "trace", skip_all)]
    fn drop(&mut self) {
        for path in &self.run_paths {
            let _ = fs::remove_file(path);
        }
    }
}
#[tracing::instrument(level = "trace", skip_all)]
fn read_run_u32(reader: &mut impl Read, context: &str) -> Result<u32> {
    let mut bytes = [0_u8; 4];
    reader
        .read_exact(&mut bytes)
        .map_err(|error| Error::InvalidInput(format!("truncated {context}: {error}")))?;
    Ok(u32::from_le_bytes(bytes))
}
#[tracing::instrument(level = "trace", skip_all)]
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
    #[tracing::instrument(level = "trace", skip_all)]
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

    #[tracing::instrument(level = "trace", skip_all)]
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
    #[tracing::instrument(level = "trace", skip_all)]
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

    #[tracing::instrument(level = "trace", skip_all)]
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
#[tracing::instrument(level = "trace", skip_all)]
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
#[tracing::instrument(level = "trace", skip_all)]
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
#[tracing::instrument(level = "trace", skip_all)]
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
#[tracing::instrument(level = "trace", skip_all)]
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
#[tracing::instrument(level = "trace", skip_all)]
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
#[tracing::instrument(level = "trace", skip_all)]
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
