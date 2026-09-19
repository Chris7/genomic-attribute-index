from pathlib import Path
from typing import Iterable, Optional, Union

PathLike = Union[str, Path]
__version__: str

class GniError(Exception): ...
class GniInputError(GniError): ...
class GniIoError(GniError): ...
class GniCorruptError(GniError): ...
class GniStaleError(GniError): ...

class GffRecord:
    reference_sequence_name: str
    source: str
    ty: str
    start: int
    end: int
    score: str
    strand: str
    phase: str
    attributes: list[tuple[str, list[str]]]
    raw_line: str
    def attribute_values(self, tag: str) -> Iterable[str]: ...

class IndexMetadata:
    major_version: int
    minor_version: int
    case_sensitive: bool
    attributes: list[str]
    gff_fingerprint: str
    coordinate_index_fingerprint: str
    reference_dictionary_fingerprint: str
    term_count: int
    unique_span_count: int
    posting_count: int
    postings_block_count: int
    span_block_count: int
    reference_count: int
    span_block_size: int
    file_size: int
    attribute_section_bytes: int
    term_dictionary_bytes: int
    postings_directory_bytes: int
    postings_uncompressed_bytes: int
    postings_data_bytes: int
    span_directory_bytes: int
    span_uncompressed_bytes: int
    span_data_bytes: int
    compressed_postings_blocks: int
    compressed_span_blocks: int
    delta_start_blocks: int
    for_start_blocks: int
    varint_length_blocks: int
    for_length_blocks: int

class BuildStats:
    records_processed: int
    records_indexed: int
    distinct_terms: int
    unique_spans: int
    postings: int
    duplicate_postings_removed: int
    duplicate_spans_removed: int
    index_bytes: int
    postings_bytes_before_compression: int
    postings_bytes_after_compression: int
    span_bytes_fixed_width: int
    span_bytes_structural: int
    span_bytes_after_compression: int
    bytes_per_term: float
    bytes_per_posting: float
    bytes_per_unique_span: float
    scan_seconds: float
    spill_seconds: float
    merge_seconds: float
    encode_postings_seconds: float
    encode_spans_seconds: float
    serialize_seconds: float
    total_seconds: float
    peak_working_set_bytes: int

class QueryStats:
    requested_spans: int
    distinct_span_blocks_decoded: int
    exact_interval_queries: int
    raw_chunks: int
    merged_chunks: int
    unique_candidate_records: int
    matching_records: int
    bytes_read: int

class IndexedGff:
    def metadata(self) -> IndexMetadata: ...
    def query(self, term: str) -> list[GffRecord]: ...
    def query_with_stats(self, term: str) -> tuple[list[GffRecord], QueryStats]: ...

def build_index(
    input: PathLike,
    coordinate_index: PathLike,
    output: PathLike,
    attributes: Iterable[str],
    *,
    case_sensitive: bool = False,
    memory_budget: int = 64 * 1024 * 1024,
    compression_threads: Optional[int] = None,
    bgzf_threads: Optional[int] = None,
) -> BuildStats: ...
def open_index(
    input: PathLike, coordinate_index: PathLike, gni: PathLike
) -> IndexedGff: ...
def query_index(
    input: PathLike, coordinate_index: PathLike, gni: PathLike, term: str
) -> list[GffRecord]: ...
def inspect_index(gni: PathLike) -> IndexMetadata: ...
