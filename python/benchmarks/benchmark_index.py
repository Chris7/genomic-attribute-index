#!/usr/bin/env python3
"""Measure GAI index builds and indexed GFF3/BED queries.

The script writes separate indexing and query tables as Markdown by default.
Opening the index (including source fingerprint checks) is outside each timed
query. Query time is the median over the configured number of calls.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import os
import statistics
import sys
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path
import gai


MATCH_MODES = ("exact", "prefix", "contains", "regex")


@dataclass(frozen=True)
class IndexResult:
    input_file: str
    file_type: str
    uncompressed_size_bytes: int
    attributes: tuple[str, ...]
    index_size_bytes: int
    build_seconds: float


@dataclass(frozen=True)
class QueryResult:
    input_file: str
    file_type: str
    attributes: tuple[str, ...]
    match_mode: str
    query: str
    records_returned: int
    query_seconds: float


def _positive_int(value: str) -> int:
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be at least 1")
    return parsed


def _nonnegative_int(value: str) -> int:
    parsed = int(value)
    if parsed < 0:
        raise argparse.ArgumentTypeError("must be zero or greater")
    return parsed


def _parse_query_mode(value: str) -> tuple[str, str]:
    mode, separator, query = value.partition("=")
    if not separator or not query:
        raise argparse.ArgumentTypeError("expected MODE=QUERY with a nonempty query")
    if mode not in MATCH_MODES:
        choices = ", ".join(MATCH_MODES)
        raise argparse.ArgumentTypeError(f"mode must be one of: {choices}")
    return mode, query


def _parse_attribute_set(value: str) -> tuple[str, ...]:
    attributes = tuple(part.strip() for part in value.split(",") if part.strip())
    if not attributes:
        raise argparse.ArgumentTypeError("attribute sets must contain at least one tag")
    if len(set(attributes)) != len(attributes):
        raise argparse.ArgumentTypeError("attribute tags must be unique within a set")
    return attributes


def _file_type(path: Path) -> str:
    suffixes = [suffix.lower() for suffix in path.suffixes]
    if suffixes and suffixes[-1] in {".gz", ".bgz", ".bgzf"}:
        suffixes.pop()
    if suffixes and suffixes[-1] in {".gff3", ".gff"}:
        return "GFF3" if suffixes[-1] == ".gff3" else "GFF"
    if suffixes and suffixes[-1] == ".bed":
        return "BED"
    raise ValueError(
        "input must have a .gff, .gff3, or .bed extension, optionally gzip-compressed"
    )


def _uncompressed_size(path: Path) -> int:
    """Return the decoded byte count without loading the whole input in memory."""
    with path.open("rb") as raw:
        compressed = raw.read(2) == b"\x1f\x8b"

    opener = gzip.open if compressed else Path.open
    count = 0
    with opener(path, "rb") as source:
        while chunk := source.read(1024 * 1024):
            count += len(chunk)
    return count


def _coordinate_index(input_path: Path, supplied: Path | None) -> Path:
    if supplied is not None:
        if not supplied.is_file():
            raise ValueError(f"coordinate index does not exist: {supplied}")
        return supplied

    candidates = [Path(f"{input_path}.tbi"), Path(f"{input_path}.csi")]
    existing = [candidate for candidate in candidates if candidate.is_file()]
    if len(existing) == 1:
        return existing[0]
    if len(existing) > 1:
        raise ValueError(
            f"both {candidates[0]} and {candidates[1]} exist; pass --coordinate-index"
        )
    raise ValueError(
        f"could not find {candidates[0]} or {candidates[1]}; pass --coordinate-index"
    )


def _validate_output_path(
    output: Path | None, input_path: Path, coordinate_index: Path
) -> Path | None:
    if output is None:
        return None
    resolved = output.expanduser().resolve()
    protected_paths = (input_path, coordinate_index.resolve())
    if resolved in protected_paths or (
        resolved.exists()
        and any(os.path.samefile(resolved, protected) for protected in protected_paths)
    ):
        raise ValueError(
            "report output must not overwrite the input or coordinate index"
        )
    if not resolved.parent.is_dir():
        raise ValueError(f"report output directory does not exist: {resolved.parent}")
    if not os.access(resolved.parent, os.W_OK | os.X_OK):
        raise ValueError(f"report output directory is not writable: {resolved.parent}")
    if resolved.exists() and not resolved.is_file():
        raise ValueError(f"report output is not a regular file: {resolved}")
    if resolved.exists() and not os.access(resolved, os.W_OK):
        raise ValueError(f"report output is not writable: {resolved}")
    return resolved


def _markdown_cell(value: object) -> str:
    return str(value).replace("|", "\\|").replace("\n", " ")


def _render_markdown(
    indexes: list[IndexResult],
    queries: list[QueryResult],
    *,
    warmups: int,
    repeats: int,
    build_repeats: int,
) -> str:
    lines = [
        f"Build times are medians of {build_repeats} build(s) per attribute set. "
        "Index opening and source fingerprint checks are excluded from query timings. "
        f"Query times are medians of {repeats} calls after {warmups} warmup call(s).",
        "",
        "### Indexing",
        "",
        "| File | Type | Uncompressed size (bytes) | Attributes indexed | Attribute count | Index size (bytes) | Build time (s) |",
        "| --- | --- | ---: | --- | ---: | ---: | ---: |",
    ]
    for result in indexes:
        values = (
            result.input_file,
            result.file_type,
            result.uncompressed_size_bytes,
            ", ".join(result.attributes),
            len(result.attributes),
            result.index_size_bytes,
            f"{result.build_seconds:.6f}",
        )
        lines.append(
            "| " + " | ".join(_markdown_cell(value) for value in values) + " |"
        )

    lines.extend(
        [
            "",
            "### Querying",
            "",
            "| File | Type | Attributes indexed | Query type | Query | Records returned | Query time (s) |",
            "| --- | --- | --- | --- | --- | ---: | ---: |",
        ]
    )
    for result in queries:
        values = (
            result.input_file,
            result.file_type,
            ", ".join(result.attributes),
            result.match_mode,
            result.query,
            result.records_returned,
            f"{result.query_seconds:.6f}",
        )
        lines.append(
            "| " + " | ".join(_markdown_cell(value) for value in values) + " |"
        )
    return "\n".join(lines) + "\n"


def _render_csv(
    indexes: list[IndexResult],
    queries: list[QueryResult],
    *,
    warmups: int,
    repeats: int,
    build_repeats: int,
) -> str:
    buffer = io.StringIO(newline="")
    writer = csv.writer(buffer)
    writer.writerow(
        [
            "section",
            "file",
            "file_type",
            "uncompressed_size_bytes",
            "attributes_indexed",
            "attribute_count",
            "index_size_bytes",
            "build_time_seconds",
            "build_repeats",
            "query_type",
            "query",
            "records_returned",
            "query_time_seconds",
            "query_repeats",
            "warmup_calls",
        ]
    )
    for result in indexes:
        writer.writerow(
            [
                "indexing",
                result.input_file,
                result.file_type,
                result.uncompressed_size_bytes,
                ", ".join(result.attributes),
                len(result.attributes),
                result.index_size_bytes,
                f"{result.build_seconds:.6f}",
                build_repeats,
                "",
                "",
                "",
                "",
                "",
                "",
            ]
        )
    for result in queries:
        writer.writerow(
            [
                "querying",
                result.input_file,
                result.file_type,
                "",
                ", ".join(result.attributes),
                len(result.attributes),
                "",
                "",
                "",
                result.match_mode,
                result.query,
                result.records_returned,
                f"{result.query_seconds:.6f}",
                repeats,
                warmups,
            ]
        )
    return buffer.getvalue()


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "input", type=Path, help="sorted, tabix-indexed GFF/GFF3 or BED file"
    )
    parser.add_argument(
        "--coordinate-index",
        type=Path,
        help="TBI/CSI path; otherwise an unambiguous sibling is discovered",
    )
    parser.add_argument(
        "--attribute-set",
        action="append",
        type=_parse_attribute_set,
        metavar="TAG[,TAG...]",
        help="GFF attributes for one index; repeat to compare index sizes and query speeds",
    )
    parser.add_argument(
        "--query",
        action="append",
        metavar="TEXT_OR_PATTERN",
        help="query text; each query is measured with all selected modes",
    )
    parser.add_argument(
        "--query-mode",
        action="append",
        type=_parse_query_mode,
        metavar="MODE=TEXT_OR_PATTERN",
        help="measure a query with one mode; repeat for mode-specific comparisons",
    )
    parser.add_argument(
        "--match",
        action="append",
        choices=MATCH_MODES,
        help="query mode to measure (repeatable; default: all four modes)",
    )
    parser.add_argument("--repeats", type=_positive_int, default=5)
    parser.add_argument("--warmups", type=_nonnegative_int, default=1)
    parser.add_argument("--build-repeats", type=_positive_int, default=1)
    parser.add_argument("--memory-budget", type=_positive_int, default=64 * 1024 * 1024)
    parser.add_argument(
        "--threads",
        type=_positive_int,
        help="set both compression and BGZF thread counts; default uses GAI settings",
    )
    parser.add_argument("--format", choices=("markdown", "csv"), default="markdown")
    parser.add_argument(
        "--output", type=Path, help="write the report to a file instead of stdout"
    )
    return parser


def run(args: argparse.Namespace) -> tuple[list[IndexResult], list[QueryResult]]:
    input_path = args.input.expanduser().resolve()
    if not input_path.is_file():
        raise ValueError(f"input file does not exist: {input_path}")
    file_type = _file_type(input_path)
    coordinate_index = _coordinate_index(
        input_path,
        args.coordinate_index.expanduser().resolve() if args.coordinate_index else None,
    )
    if file_type == "BED":
        if args.attribute_set:
            raise ValueError(
                "BED always indexes its name column as `name`; omit --attribute-set"
            )
        attribute_sets = [("name",)]
    else:
        if not args.attribute_set:
            raise ValueError("GFF/GFF3 input requires at least one --attribute-set")
        attribute_sets = args.attribute_set
    if len(set(attribute_sets)) != len(attribute_sets):
        raise ValueError("attribute sets must be distinct")
    modes = args.match if args.match else list(MATCH_MODES)
    if args.query_mode:
        if args.query or args.match:
            raise ValueError("use --query-mode by itself, without --query or --match")
        query_cases = args.query_mode
    else:
        if not args.query:
            raise ValueError("provide at least one --query or --query-mode")
        query_cases = [(mode, query) for query in args.query for mode in modes]
    _validate_output_path(args.output, input_path, coordinate_index)

    uncompressed_size_bytes = _uncompressed_size(input_path)
    indexes: list[IndexResult] = []
    queries: list[QueryResult] = []

    with tempfile.TemporaryDirectory(prefix="gai-benchmark-") as temporary_directory:
        temporary_path = Path(temporary_directory)
        for set_number, attributes in enumerate(attribute_sets, start=1):
            gai_path: Path | None = None
            build_times: list[float] = []
            for build_number in range(1, args.build_repeats + 1):
                next_gai_path = (
                    temporary_path / f"attributes-{set_number}-build-{build_number}.gai"
                )
                started = time.perf_counter()
                gai.build_index(
                    input_path,
                    coordinate_index,
                    next_gai_path,
                    None if file_type == "BED" else attributes,
                    memory_budget=args.memory_budget,
                    compression_threads=args.threads,
                    bgzf_threads=args.threads,
                )
                build_times.append(time.perf_counter() - started)
                if gai_path is not None:
                    gai_path.unlink()
                gai_path = next_gai_path

            assert gai_path is not None  # build_repeats is validated as positive

            indexes.append(
                IndexResult(
                    input_file=input_path.name,
                    file_type=file_type,
                    uncompressed_size_bytes=uncompressed_size_bytes,
                    attributes=attributes,
                    index_size_bytes=gai_path.stat().st_size,
                    build_seconds=statistics.median(build_times),
                )
            )

            # This validates fingerprints once, outside all timed query calls.
            indexed = gai.open_index(input_path, coordinate_index, gai_path)
            for match_mode, query in query_cases:
                for _ in range(args.warmups):
                    indexed.query(query, match=match_mode)

                query_times: list[float] = []
                records_returned = 0
                for _ in range(args.repeats):
                    started = time.perf_counter()
                    records = indexed.query(query, match=match_mode)
                    query_times.append(time.perf_counter() - started)
                    records_returned = len(records)
                    del records

                queries.append(
                    QueryResult(
                        input_file=input_path.name,
                        file_type=file_type,
                        attributes=attributes,
                        match_mode=match_mode,
                        query=query,
                        records_returned=records_returned,
                        query_seconds=statistics.median(query_times),
                    )
                )
    return indexes, queries


def main(argv: list[str] | None = None) -> int:
    parser = _parser()
    args = parser.parse_args(argv)
    try:
        indexes, queries = run(args)
    except (OSError, ValueError, gai.GaiError) as error:
        parser.error(str(error))

    render = _render_markdown if args.format == "markdown" else _render_csv
    report = render(
        indexes,
        queries,
        warmups=args.warmups,
        repeats=args.repeats,
        build_repeats=args.build_repeats,
    )
    if args.output:
        try:
            output_path = args.output.expanduser().resolve()
            output_path.write_text(report, encoding="utf-8")
        except OSError as error:
            parser.error(f"could not write report to {args.output}: {error}")
    else:
        sys.stdout.write(report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
