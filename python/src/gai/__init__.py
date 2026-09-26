"""Python API for the deterministic Genomic Attribute Index.

The heavy build and query work runs in Rust and releases the GIL. Paths may be
strings or :class:`pathlib.Path` instances.
"""

from importlib.metadata import PackageNotFoundError, version
from typing import Literal

from ._gai import (
    BuildStats,
    GffRecord,
    GaiCorruptError,
    GaiError,
    GaiInputError,
    GaiIoError,
    GaiStaleError,
    IndexMetadata,
    IndexedSource,
    QueryStats,
    __version__ as _extension_version,
    build_index as _build_index,
    inspect_index as _inspect_index,
    open_index as _open_index,
    query_index as _query_index,
    sort as _sort,
)

__all__ = [
    "BuildStats",
    "GffRecord",
    "GaiCorruptError",
    "GaiError",
    "GaiInputError",
    "GaiIoError",
    "GaiStaleError",
    "IndexMetadata",
    "IndexedSource",
    "MatchMode",
    "QueryStats",
    "build_index",
    "inspect_index",
    "open_index",
    "query_index",
    "sort",
]

MatchMode = Literal["exact", "prefix", "contains", "regex"]
try:
    __version__ = version("genomic-attribute-index")
except PackageNotFoundError:
    # Source-tree imports before a wheel/editable install still expose a useful
    # version from the compiled extension, which is compiled from Cargo's
    # package version.
    __version__ = _extension_version


def sort(input, output, *, disk_sort=False) -> None:
    """Sort a GFF/GFF3 or BED file into ``output``.

    The input format is inferred from its extension. Set ``disk_sort=True``
    for files that may not fit in memory.
    """
    return _sort(input, output, disk_sort=disk_sort)


def build_index(
    input,
    coordinate_index,
    output,
    attributes=None,
    *,
    case_sensitive=False,
    memory_budget=64 * 1024 * 1024,
    compression_threads=None,
    bgzf_threads=None,
) -> BuildStats:
    """Build a GAI beside a BGZF or plain GFF3 source, or a BGZF BED source.

    For GFF3, ``attributes`` is a nonempty iterable of explicit attribute
    tags. For BED, names from column 4 are indexed and ``attributes`` is
    ignored. Duplicate GFF3 tags are removed while preserving first
    occurrence. The result contains deterministic build counters and phase
    timings.
    """
    return _build_index(
        input,
        coordinate_index,
        output,
        [] if attributes is None else list(attributes),
        case_sensitive,
        memory_budget,
        compression_threads,
        bgzf_threads,
    )


def open_index(input, coordinate_index, gai) -> IndexedSource:
    """Open a source, TBI/CSI index, and matching GAI with stale checks."""
    return _open_index(input, coordinate_index, gai)


def query_index(
    input, coordinate_index, gai, term: str, *, match: MatchMode = "exact"
) -> list[GffRecord]:
    """Open an indexed GFF3 or BED source and return records matching ``term``.

    ``match`` accepts ``"exact"`` (the default), ``"prefix"``,
    ``"contains"`` for literal substring matching, or ``"regex"`` for a
    Unicode-aware regular expression search. GFF3 attribute values or BED
    column-4 names are normalized before matching; regex syntax is preserved
    and the pattern is trimmed only at its boundaries.
    """
    return _query_index(input, coordinate_index, gai, term, match=match)


def inspect_index(gai) -> IndexMetadata:
    """Read format and compression metadata from a GAI file."""
    return _inspect_index(gai)
