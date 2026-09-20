from pathlib import Path
from typing import Iterable, Literal, Optional, Union

PathLike = Union[str, Path]
MatchMode = Literal["exact", "prefix"]
__version__: str

class GaiError(Exception): ...
class GaiInputError(GaiError): ...
class GaiIoError(GaiError): ...
class GaiCorruptError(GaiError): ...
class GaiStaleError(GaiError): ...

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
    starts_data_bytes: int
    lengths_data_bytes: int
    starts_uncompressed_bytes: int
    lengths_uncompressed_bytes: int
    compressed_postings_blocks: int
    delta_start_blocks: int
    varint_length_blocks: int
    for_length_blocks: int
    compressed_start_blocks: int
    compressed_length_blocks: int

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
    span_starts_bytes_before_compression: int
    span_starts_bytes_after_compression: int
    span_lengths_bytes_before_compression: int
    span_lengths_bytes_after_compression: int
    delta_start_blocks: int
    length_varint_blocks: int
    length_for_blocks: int
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
    def query(self, term: str, *, match: MatchMode = "exact") -> list[GffRecord]: ...
    def query_with_stats(
        self, term: str, *, match: MatchMode = "exact"
    ) -> tuple[list[GffRecord], QueryStats]: ...

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
    input: PathLike, coordinate_index: PathLike, gai: PathLike
) -> IndexedGff: ...
def query_index(
    input: PathLike,
    coordinate_index: PathLike,
    gai: PathLike,
    term: str,
    *,
    match: MatchMode = "exact",
) -> list[GffRecord]: ...
def inspect_index(gai: PathLike) -> IndexMetadata: ...
