"""Python API for the deterministic GFF Name Index.

The heavy build and query work runs in Rust and releases the GIL. Paths may be
strings or :class:`pathlib.Path` instances.
"""

from importlib.metadata import PackageNotFoundError, version

from ._gni import (
    BuildStats,
    GffRecord,
    GniCorruptError,
    GniError,
    GniInputError,
    GniIoError,
    GniStaleError,
    IndexMetadata,
    IndexedGff,
    QueryStats,
    build_index as _build_index,
    inspect_index as _inspect_index,
    open_index as _open_index,
    query_index as _query_index,
)

__all__ = [
    "BuildStats",
    "GffRecord",
    "GniCorruptError",
    "GniError",
    "GniInputError",
    "GniIoError",
    "GniStaleError",
    "IndexMetadata",
    "IndexedGff",
    "QueryStats",
    "build_index",
    "inspect_index",
    "open_index",
    "query_index",
]
try:
    __version__ = version("genomic-attribute-index")
except PackageNotFoundError:
    # Source-tree imports before a wheel/editable install still expose a useful
    # version; release metadata is checked against Cargo and pyproject.
    __version__ = "0.1.0"


def build_index(
    input,
    coordinate_index,
    output,
    attributes,
    *,
    case_sensitive=False,
    memory_budget=64 * 1024 * 1024,
    compression_threads=None,
    bgzf_threads=None,
) -> BuildStats:
    """Build a GNI beside a BGZF or plain GFF3 source.

    ``attributes`` is a nonempty iterable of explicit GFF3 attribute tags.
    Duplicate tags are removed while preserving first occurrence. The result
    contains deterministic build counters and phase timings.
    """
    return _build_index(
        input,
        coordinate_index,
        output,
        list(attributes),
        case_sensitive,
        memory_budget,
        compression_threads,
        bgzf_threads,
    )


def open_index(input, coordinate_index, gni) -> IndexedGff:
    """Open a source, TBI/CSI index, and matching GNI with stale checks."""
    return _open_index(input, coordinate_index, gni)


def query_index(input, coordinate_index, gni, term: str) -> list[GffRecord]:
    """Open an indexed source and return records matching ``term``."""
    return _query_index(input, coordinate_index, gni, term)


def inspect_index(gni) -> IndexMetadata:
    """Read format and compression metadata from a GNI file."""
    return _inspect_index(gni)
