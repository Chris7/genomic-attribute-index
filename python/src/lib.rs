use std::{
    fs::File,
    io::{BufWriter, Write},
    path::PathBuf,
    sync::Mutex,
};

use gai::{
    BuildOptions, Error, GffRecord, IndexMetadata, IndexStats, IndexedGff, MatchMode,
    NameIndexOptions, QueryStats, SortFormat, build_name_index_with_options, sort_file,
};
use pyo3::{
    Bound, PyErr, PyResult, Python, create_exception,
    exceptions::{PyException, PyUserWarning},
    prelude::PyModule,
    pyclass, pyfunction, pymethods, pymodule,
    types::PyModuleMethods,
    wrap_pyfunction,
};

create_exception!(
    _gai,
    GaiError,
    PyException,
    "Base error raised by the GAI bindings."
);
create_exception!(
    _gai,
    GaiInputError,
    GaiError,
    "Invalid GAI input or API argument."
);
create_exception!(_gai, GaiIoError, GaiError, "GAI input/output failure.");
create_exception!(_gai, GaiCorruptError, GaiError, "Corrupt GAI data.");
create_exception!(
    _gai,
    GaiStaleError,
    GaiError,
    "GAI/source/index fingerprints do not match."
);

fn to_py_error(error: Error) -> pyo3::PyErr {
    match error {
        Error::Io(error) => GaiIoError::new_err(error.to_string()),
        Error::InvalidInput(message) => GaiInputError::new_err(message),
        Error::InvalidCoordinate => GaiInputError::new_err("invalid coordinate"),
        Error::Corrupt(message) => GaiCorruptError::new_err(message),
        Error::Stale(message) => GaiStaleError::new_err(message),
        Error::Compression(message) => GaiError::new_err(format!("compression error: {message}")),
    }
}

fn hex(bytes: &[u8; 32]) -> String {
    bytes.iter().map(|byte| format!("{byte:02x}")).collect()
}

#[pyclass(name = "GffRecord", frozen)]
#[derive(Clone)]
struct PyGffRecord {
    #[pyo3(get)]
    reference_sequence_name: String,
    #[pyo3(get)]
    source: String,
    #[pyo3(get)]
    ty: String,
    #[pyo3(get)]
    start: u64,
    #[pyo3(get)]
    end: u64,
    #[pyo3(get)]
    score: String,
    #[pyo3(get)]
    strand: String,
    #[pyo3(get)]
    phase: String,
    #[pyo3(get)]
    attributes: Vec<(String, Vec<String>)>,
    #[pyo3(get)]
    raw_line: String,
}

impl From<GffRecord> for PyGffRecord {
    fn from(record: GffRecord) -> Self {
        Self {
            reference_sequence_name: record.reference_sequence_name,
            source: record.source,
            ty: record.ty,
            start: record.start,
            end: record.end,
            score: record.score,
            strand: record.strand,
            phase: record.phase,
            attributes: record.attributes,
            raw_line: record.raw_line,
        }
    }
}

#[pymethods]
impl PyGffRecord {
    fn attribute_values(&self, tag: &str) -> Vec<String> {
        self.attributes
            .iter()
            .filter(|(name, _)| name == tag)
            .flat_map(|(_, values)| values.iter().cloned())
            .collect()
    }

    fn __repr__(&self) -> String {
        format!(
            "GffRecord(reference_sequence_name={:?}, start={}, end={}, raw_line={:?})",
            self.reference_sequence_name, self.start, self.end, self.raw_line
        )
    }
}

#[pyclass(name = "IndexMetadata", frozen)]
#[derive(Clone)]
struct PyIndexMetadata {
    #[pyo3(get)]
    major_version: u16,
    #[pyo3(get)]
    minor_version: u16,
    #[pyo3(get)]
    case_sensitive: bool,
    #[pyo3(get)]
    attributes: Vec<String>,
    #[pyo3(get)]
    gff_fingerprint: String,
    #[pyo3(get)]
    coordinate_index_fingerprint: String,
    #[pyo3(get)]
    reference_dictionary_fingerprint: String,
    #[pyo3(get)]
    term_count: u64,
    #[pyo3(get)]
    unique_span_count: u64,
    #[pyo3(get)]
    posting_count: u64,
    #[pyo3(get)]
    postings_block_count: u64,
    #[pyo3(get)]
    span_block_count: u64,
    #[pyo3(get)]
    reference_count: u32,
    #[pyo3(get)]
    span_block_size: u32,
    #[pyo3(get)]
    file_size: u64,
    #[pyo3(get)]
    attribute_section_bytes: u64,
    #[pyo3(get)]
    term_dictionary_bytes: u64,
    #[pyo3(get)]
    postings_directory_bytes: u64,
    #[pyo3(get)]
    postings_uncompressed_bytes: u64,
    #[pyo3(get)]
    postings_data_bytes: u64,
    #[pyo3(get)]
    span_directory_bytes: u64,
    #[pyo3(get)]
    span_uncompressed_bytes: u64,
    #[pyo3(get)]
    starts_data_bytes: u64,
    #[pyo3(get)]
    lengths_data_bytes: u64,
    #[pyo3(get)]
    starts_uncompressed_bytes: u64,
    #[pyo3(get)]
    lengths_uncompressed_bytes: u64,
    #[pyo3(get)]
    compressed_postings_blocks: u64,
    #[pyo3(get)]
    delta_start_blocks: u64,
    #[pyo3(get)]
    varint_length_blocks: u64,
    #[pyo3(get)]
    for_length_blocks: u64,
    #[pyo3(get)]
    compressed_start_blocks: u64,
    #[pyo3(get)]
    compressed_length_blocks: u64,
}

impl From<IndexMetadata> for PyIndexMetadata {
    fn from(metadata: IndexMetadata) -> Self {
        Self {
            major_version: metadata.major_version,
            minor_version: metadata.minor_version,
            case_sensitive: metadata.case_sensitive,
            attributes: metadata.attributes,
            gff_fingerprint: hex(&metadata.gff_fingerprint),
            coordinate_index_fingerprint: hex(&metadata.coordinate_index_fingerprint),
            reference_dictionary_fingerprint: hex(&metadata.reference_dictionary_fingerprint),
            term_count: metadata.term_count,
            unique_span_count: metadata.unique_span_count,
            posting_count: metadata.posting_count,
            postings_block_count: metadata.postings_block_count,
            span_block_count: metadata.span_block_count,
            reference_count: metadata.reference_count,
            span_block_size: metadata.span_block_size,
            file_size: metadata.file_size,
            attribute_section_bytes: metadata.attribute_section_bytes,
            term_dictionary_bytes: metadata.term_dictionary_bytes,
            postings_directory_bytes: metadata.postings_directory_bytes,
            postings_uncompressed_bytes: metadata.postings_uncompressed_bytes,
            postings_data_bytes: metadata.postings_data_bytes,
            span_directory_bytes: metadata.span_directory_bytes,
            span_uncompressed_bytes: metadata.span_uncompressed_bytes,
            starts_data_bytes: metadata.starts_data_bytes,
            lengths_data_bytes: metadata.lengths_data_bytes,
            starts_uncompressed_bytes: metadata.starts_uncompressed_bytes,
            lengths_uncompressed_bytes: metadata.lengths_uncompressed_bytes,
            compressed_postings_blocks: metadata.compressed_postings_blocks,
            delta_start_blocks: metadata.delta_start_blocks,
            varint_length_blocks: metadata.varint_length_blocks,
            for_length_blocks: metadata.for_length_blocks,
            compressed_start_blocks: metadata.compressed_start_blocks,
            compressed_length_blocks: metadata.compressed_length_blocks,
        }
    }
}

#[pyclass(name = "BuildStats", frozen)]
#[derive(Clone)]
struct PyBuildStats {
    #[pyo3(get)]
    records_processed: u64,
    #[pyo3(get)]
    records_indexed: u64,
    #[pyo3(get)]
    distinct_terms: u64,
    #[pyo3(get)]
    unique_spans: u64,
    #[pyo3(get)]
    postings: u64,
    #[pyo3(get)]
    duplicate_postings_removed: u64,
    #[pyo3(get)]
    duplicate_spans_removed: u64,
    #[pyo3(get)]
    index_bytes: u64,
    #[pyo3(get)]
    postings_bytes_before_compression: u64,
    #[pyo3(get)]
    postings_bytes_after_compression: u64,
    #[pyo3(get)]
    span_bytes_fixed_width: u64,
    #[pyo3(get)]
    span_bytes_structural: u64,
    #[pyo3(get)]
    span_bytes_after_compression: u64,
    #[pyo3(get)]
    span_starts_bytes_before_compression: u64,
    #[pyo3(get)]
    span_starts_bytes_after_compression: u64,
    #[pyo3(get)]
    span_lengths_bytes_before_compression: u64,
    #[pyo3(get)]
    span_lengths_bytes_after_compression: u64,
    #[pyo3(get)]
    delta_start_blocks: u64,
    #[pyo3(get)]
    length_varint_blocks: u64,
    #[pyo3(get)]
    length_for_blocks: u64,
    #[pyo3(get)]
    bytes_per_term: f64,
    #[pyo3(get)]
    bytes_per_posting: f64,
    #[pyo3(get)]
    bytes_per_unique_span: f64,
    #[pyo3(get)]
    scan_seconds: f64,
    #[pyo3(get)]
    spill_seconds: f64,
    #[pyo3(get)]
    merge_seconds: f64,
    #[pyo3(get)]
    encode_postings_seconds: f64,
    #[pyo3(get)]
    encode_spans_seconds: f64,
    #[pyo3(get)]
    serialize_seconds: f64,
    #[pyo3(get)]
    total_seconds: f64,
    #[pyo3(get)]
    peak_working_set_bytes: u64,
}

impl From<IndexStats> for PyBuildStats {
    fn from(stats: IndexStats) -> Self {
        Self {
            records_processed: stats.records_processed,
            records_indexed: stats.records_indexed,
            distinct_terms: stats.distinct_terms,
            unique_spans: stats.unique_spans,
            postings: stats.postings,
            duplicate_postings_removed: stats.duplicate_postings_removed,
            duplicate_spans_removed: stats.duplicate_spans_removed,
            index_bytes: stats.index_bytes,
            postings_bytes_before_compression: stats.postings_bytes_before_compression,
            postings_bytes_after_compression: stats.postings_bytes_after_compression,
            span_bytes_fixed_width: stats.span_bytes_fixed_width,
            span_bytes_structural: stats.span_bytes_structural,
            span_bytes_after_compression: stats.span_bytes_after_compression,
            span_starts_bytes_before_compression: stats.span_starts_bytes_before_compression,
            span_starts_bytes_after_compression: stats.span_starts_bytes_after_compression,
            span_lengths_bytes_before_compression: stats.span_lengths_bytes_before_compression,
            span_lengths_bytes_after_compression: stats.span_lengths_bytes_after_compression,
            delta_start_blocks: stats.delta_start_blocks,
            length_varint_blocks: stats.length_varint_blocks,
            length_for_blocks: stats.length_for_blocks,
            bytes_per_term: stats.bytes_per_term,
            bytes_per_posting: stats.bytes_per_posting,
            bytes_per_unique_span: stats.bytes_per_unique_span,
            scan_seconds: stats.timings.scan.as_secs_f64(),
            spill_seconds: stats.timings.spill.as_secs_f64(),
            merge_seconds: stats.timings.merge.as_secs_f64(),
            encode_postings_seconds: stats.timings.encode_postings.as_secs_f64(),
            encode_spans_seconds: stats.timings.encode_spans.as_secs_f64(),
            serialize_seconds: stats.timings.serialize.as_secs_f64(),
            total_seconds: stats.timings.total.as_secs_f64(),
            peak_working_set_bytes: stats.peak_working_set_bytes,
        }
    }
}

#[pyclass(name = "QueryStats", frozen)]
#[derive(Clone)]
struct PyQueryStats {
    #[pyo3(get)]
    requested_spans: u64,
    #[pyo3(get)]
    distinct_span_blocks_decoded: u64,
    #[pyo3(get)]
    exact_interval_queries: u64,
    #[pyo3(get)]
    raw_chunks: u64,
    #[pyo3(get)]
    merged_chunks: u64,
    #[pyo3(get)]
    unique_candidate_records: u64,
    #[pyo3(get)]
    matching_records: u64,
    #[pyo3(get)]
    bytes_read: u64,
}

impl From<QueryStats> for PyQueryStats {
    fn from(stats: QueryStats) -> Self {
        Self {
            requested_spans: stats.requested_spans,
            distinct_span_blocks_decoded: stats.distinct_span_blocks_decoded,
            exact_interval_queries: stats.exact_interval_queries,
            raw_chunks: stats.raw_chunks,
            merged_chunks: stats.merged_chunks,
            unique_candidate_records: stats.unique_candidate_records,
            matching_records: stats.matching_records,
            bytes_read: stats.bytes_read,
        }
    }
}

#[pyclass(name = "IndexedGff")]
struct PyIndexedGff {
    inner: Mutex<IndexedGff>,
}

fn records(records: Vec<GffRecord>) -> Vec<PyGffRecord> {
    records.into_iter().map(PyGffRecord::from).collect()
}

fn open_inner(
    input: PathBuf,
    coordinate_index: PathBuf,
    gai: PathBuf,
) -> Result<IndexedGff, Error> {
    IndexedGff::open(input, coordinate_index, gai)
}

fn parse_match_mode(value: &str) -> Result<MatchMode, Error> {
    MatchMode::parse(value)
}

#[pymethods]
impl PyIndexedGff {
    /// Return metadata after the source/index/GAI fingerprints were checked.
    fn metadata(&self) -> PyResult<PyIndexMetadata> {
        let guard = self
            .inner
            .lock()
            .map_err(|_| GaiError::new_err("indexed reader lock is poisoned"))?;
        Ok(guard.metadata().clone().into())
    }

    /// Query one normalized configured attribute value.
    #[pyo3(signature = (term, *, r#match = "exact"))]
    fn query(&self, py: Python<'_>, term: &str, r#match: &str) -> PyResult<Vec<PyGffRecord>> {
        let term = term.to_owned();
        let match_mode = parse_match_mode(r#match).map_err(to_py_error)?;
        let result = py.allow_threads(|| {
            let mut guard = self
                .inner
                .lock()
                .map_err(|_| Error::InvalidInput("indexed reader lock is poisoned".into()))?;
            guard.query_name_with_mode(&term, match_mode)
        });
        result.map(records).map_err(to_py_error)
    }

    /// Query and return bounded I/O/candidate instrumentation.
    #[pyo3(signature = (term, *, r#match = "exact"))]
    fn query_with_stats(
        &self,
        py: Python<'_>,
        term: &str,
        r#match: &str,
    ) -> PyResult<(Vec<PyGffRecord>, PyQueryStats)> {
        let term = term.to_owned();
        let match_mode = parse_match_mode(r#match).map_err(to_py_error)?;
        let result = py.allow_threads(|| {
            let mut guard = self
                .inner
                .lock()
                .map_err(|_| Error::InvalidInput("indexed reader lock is poisoned".into()))?;
            guard.query_name_with_mode_and_stats(&term, match_mode)
        });
        result
            .map(|(values, stats)| (records(values), stats.into()))
            .map_err(to_py_error)
    }

    fn __repr__(&self) -> String {
        "IndexedGff(...)".to_string()
    }
}

#[pyfunction]
#[pyo3(signature = (input, output, *, disk_sort=false))]
fn sort(py: Python<'_>, input: PathBuf, output: PathBuf, disk_sort: bool) -> PyResult<()> {
    py.allow_threads(|| -> Result<(), Error> {
        let output = File::create(output)?;
        let mut writer = BufWriter::new(output);
        sort_file(input, disk_sort, &mut writer)?;
        writer.flush()?;
        Ok(())
    })
    .map_err(to_py_error)
}

#[pyfunction]
#[pyo3(signature = (input, coordinate_index, output, attributes, case_sensitive=false, memory_budget=67108864, compression_threads=None, bgzf_threads=None))]
#[allow(clippy::too_many_arguments)]
fn build_index(
    py: Python<'_>,
    input: PathBuf,
    coordinate_index: PathBuf,
    output: PathBuf,
    attributes: Vec<String>,
    case_sensitive: bool,
    memory_budget: usize,
    compression_threads: Option<usize>,
    bgzf_threads: Option<usize>,
) -> PyResult<PyBuildStats> {
    let format = SortFormat::from_path(&input).map_err(to_py_error)?;
    let options = match format {
        SortFormat::Gff => {
            NameIndexOptions::new(attributes, case_sensitive).map_err(to_py_error)?
        }
        SortFormat::Bed => {
            if !attributes.is_empty() {
                let warning = py.get_type::<PyUserWarning>();
                PyErr::warn(
                    py,
                    &warning,
                    c"attributes are ignored for BED input; BED indexes the name field (column 4)",
                    1,
                )?;
            }
            NameIndexOptions::bed(case_sensitive)
        }
    };
    let mut build_options = BuildOptions::default().with_memory_budget(memory_budget);
    if let Some(threads) = compression_threads {
        build_options = build_options.with_compression_threads(threads);
    }
    if let Some(threads) = bgzf_threads {
        build_options = build_options.with_bgzf_threads(threads);
    }
    py.allow_threads(|| {
        build_name_index_with_options(input, coordinate_index, output, &options, &build_options)
    })
    .map(PyBuildStats::from)
    .map_err(to_py_error)
}

#[pyfunction]
fn open_index(
    py: Python<'_>,
    input: PathBuf,
    coordinate_index: PathBuf,
    gai: PathBuf,
) -> PyResult<PyIndexedGff> {
    py.allow_threads(|| open_inner(input, coordinate_index, gai))
        .map(|inner| PyIndexedGff {
            inner: Mutex::new(inner),
        })
        .map_err(to_py_error)
}

#[pyfunction]
#[pyo3(signature = (input, coordinate_index, gai, term, *, r#match = "exact"))]
fn query_index(
    py: Python<'_>,
    input: PathBuf,
    coordinate_index: PathBuf,
    gai: PathBuf,
    term: String,
    r#match: &str,
) -> PyResult<Vec<PyGffRecord>> {
    let match_mode = parse_match_mode(r#match).map_err(to_py_error)?;
    py.allow_threads(|| {
        let mut indexed = open_inner(input, coordinate_index, gai)?;
        indexed.query_name_with_mode(&term, match_mode)
    })
    .map(records)
    .map_err(to_py_error)
}

#[pyfunction]
fn inspect_index(py: Python<'_>, gai: PathBuf) -> PyResult<PyIndexMetadata> {
    py.allow_threads(|| gai::NameIndexReader::open(gai).map(|reader| reader.inspect()))
        .map(Into::into)
        .map_err(to_py_error)
}

#[pymodule]
fn _gai(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add("GaiError", m.py().get_type::<GaiError>())?;
    m.add("GaiInputError", m.py().get_type::<GaiInputError>())?;
    m.add("GaiIoError", m.py().get_type::<GaiIoError>())?;
    m.add("GaiCorruptError", m.py().get_type::<GaiCorruptError>())?;
    m.add("GaiStaleError", m.py().get_type::<GaiStaleError>())?;
    m.add_class::<PyGffRecord>()?;
    m.add_class::<PyIndexMetadata>()?;
    m.add_class::<PyBuildStats>()?;
    m.add_class::<PyQueryStats>()?;
    m.add_class::<PyIndexedGff>()?;
    m.add_function(wrap_pyfunction!(sort, m)?)?;
    m.add_function(wrap_pyfunction!(build_index, m)?)?;
    m.add_function(wrap_pyfunction!(open_index, m)?)?;
    m.add_function(wrap_pyfunction!(query_index, m)?)?;
    m.add_function(wrap_pyfunction!(inspect_index, m)?)?;
    Ok(())
}
