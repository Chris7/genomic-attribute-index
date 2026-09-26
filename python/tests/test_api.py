from __future__ import annotations

import gzip
from importlib.metadata import version as distribution_version
from pathlib import Path

import pytest

import gai


def test_distribution_version_is_exposed():
    expected = distribution_version("genomic-attribute-index")
    assert gai.__version__ == expected
    assert gai._gai.__version__ == expected


def _build(
    paths: dict[str, Path], output: Path, *, case_sensitive: bool = False
) -> gai.BuildStats:
    return gai.build_index(
        paths["source"],
        paths["tbi"],
        output,
        ["Name", "Alias", "Name"],
        case_sensitive=case_sensitive,
        memory_budget=1,
        compression_threads=1,
        bgzf_threads=1,
    )


def test_sort(fixture_paths, tmp_path):
    output = tmp_path / "sorted.gff3"
    gai.sort(fixture_paths["unsorted_gff"], output)
    with gzip.open(fixture_paths["source"], "rb") as handle:
        expected = handle.read()
    with open(output, "rb") as handle:
        results = handle.read()
    assert expected == results


def test_build_query_inspect_tbi_and_source_order(fixture_paths, tmp_path):
    destination = tmp_path / "fixture.gai"
    stats = _build(fixture_paths, destination)
    assert stats.records_processed == 5
    assert stats.records_indexed == 5
    assert stats.distinct_terms >= 4
    assert stats.total_seconds >= 0

    indexed = gai.open_index(fixture_paths["source"], fixture_paths["tbi"], destination)
    metadata = indexed.metadata()
    assert metadata.attributes == ["Name", "Alias"]
    assert metadata.major_version == 1
    assert metadata.minor_version == 0
    assert len(metadata.source_fingerprint) == 64

    records, query_stats = indexed.query_with_stats(" ALPHA ")
    assert [record.raw_line for record in records] == [
        "chr1\tsrc\tgene\t10\t20\t.\t+\t.\tName=Alpha;Alias=Beta,Gamma",
        "chr1\tsrc\tgene\t10\t20\t.\t+\t.\tName=Alpha;ID=identical",
        "chr2\tsrc\tgene\t5\t15\t.\t+\t.\tAlias=Alpha",
        "chr2\tsrc\tgene\t5\t25\t.\t+\t.\tName=Alpha",
    ]
    assert query_stats.matching_records == 4
    assert query_stats.exact_interval_queries == query_stats.requested_spans
    assert records[0].attribute_values("Alias") == ["Beta", "Gamma"]
    assert records[2].attributes == [("Alias", ["Alpha"])]
    assert (
        gai.query_index(
            fixture_paths["source"], fixture_paths["tbi"], destination, "other gene"
        )[0].start
        == 30
    )
    prefix_records = indexed.query("alp", match="prefix")
    assert len(prefix_records) == 4
    prefix_records, prefix_stats = indexed.query_with_stats("alp", match="prefix")
    assert len(prefix_records) == 4
    assert prefix_stats.matching_records == 4
    assert (
        len(
            gai.query_index(
                fixture_paths["source"],
                fixture_paths["tbi"],
                destination,
                "alp",
                match="prefix",
            )
        )
        == 4
    )
    with pytest.raises(gai.GaiInputError, match="exact.*prefix"):
        indexed.query("alp", match="substring")
    with pytest.raises(TypeError):
        indexed.query("alp", "prefix")
    with pytest.raises(TypeError):
        gai.query_index(
            fixture_paths["source"], fixture_paths["tbi"], destination, "alp", "prefix"
        )

    inspected = gai.inspect_index(destination)
    assert inspected.file_size == metadata.file_size
    assert inspected.starts_data_bytes >= 0
    assert inspected.lengths_data_bytes >= 0
    assert inspected.delta_start_blocks == metadata.span_block_count
    assert inspected.compressed_start_blocks >= 0
    assert inspected.compressed_length_blocks >= 0


def test_csi_case_sensitive_unknown_and_no_id(fixture_paths, tmp_path):
    destination = tmp_path / "fixture-csi.gai"
    gai.build_index(
        fixture_paths["source"],
        fixture_paths["csi"],
        destination,
        ["Name", "Alias"],
        case_sensitive=True,
        compression_threads=1,
        bgzf_threads=1,
    )
    indexed = gai.open_index(fixture_paths["source"], fixture_paths["csi"], destination)
    assert len(indexed.query("Alpha")) == 4
    assert indexed.query("alpha") == []
    assert len(indexed.query("Al", match="prefix")) == 4
    assert indexed.query("al", match="prefix") == []
    assert len(indexed.query("ph", match="contains")) == 4
    assert indexed.query("PH", match="contains") == []
    assert len(indexed.query(r"^Alpha$", match="regex")) == 4
    assert indexed.query(r"^alpha$", match="regex") == []
    assert indexed.query("missing") == []

    id_only = tmp_path / "id-only.gai"
    gai.build_index(fixture_paths["source"], fixture_paths["tbi"], id_only, ["ID"])
    assert (
        gai.open_index(fixture_paths["source"], fixture_paths["tbi"], id_only).query(
            "Alpha"
        )
        == []
    )


def test_contains_regex_modes_and_query_stats(fixture_paths, tmp_path):
    destination = tmp_path / "matching.gai"
    _build(fixture_paths, destination)
    indexed = gai.open_index(fixture_paths["source"], fixture_paths["tbi"], destination)

    contains_records, contains_stats = indexed.query_with_stats(
        " PH ", match="contains"
    )
    alpha_records = indexed.query("Alpha")
    assert [record.raw_line for record in contains_records] == [
        record.raw_line for record in alpha_records
    ]
    assert contains_stats.matching_records == len(alpha_records) == 4

    regex_records, regex_stats = indexed.query_with_stats(r"^ALP[A-Z]+$", match="regex")
    assert [record.raw_line for record in regex_records] == [
        record.raw_line for record in alpha_records
    ]
    assert regex_stats.matching_records == 4
    assert indexed.query(r"^doesnotexist$", match="regex") == []
    assert indexed.query("", match="contains") == []
    assert indexed.query("   ", match="regex") == []

    # Uppercase regex escapes retain their regex meaning instead of being
    # lowercased along with literal queries.
    assert len(indexed.query(r"\S", match="regex")) == 5

    assert (
        len(
            gai.query_index(
                fixture_paths["source"],
                fixture_paths["tbi"],
                destination,
                "ph",
                match="contains",
            )
        )
        == 4
    )
    with pytest.raises(gai.GaiInputError, match="invalid regex"):
        indexed.query("no-match[", match="regex")
    assert (
        len(
            gai.query_index(
                fixture_paths["source"],
                fixture_paths["tbi"],
                destination,
                r"^ALP[A-Z]+$",
                match="regex",
            )
        )
        == 4
    )


def test_bed_queries_support_all_match_modes(fixture_paths, tmp_path):
    destination = tmp_path / "bed.gai"
    gai.build_index(
        fixture_paths["bed_source"],
        fixture_paths["bed_tbi"],
        destination,
        ["ignored"],
    )
    indexed = gai.open_index(
        fixture_paths["bed_source"], fixture_paths["bed_tbi"], destination
    )

    exact = indexed.query("thrL")
    assert len(exact) == 1
    assert exact[0].attribute_values("name") == ["thrL"]

    prefix = indexed.query("thr", match="prefix")
    assert {
        value for record in prefix for value in record.attribute_values("name")
    } >= {
        "thrL",
        "thrA",
        "thrB",
        "thrC",
    }

    contains, contains_stats = indexed.query_with_stats(" hr ", match="contains")
    contains_names = {
        value for record in contains for value in record.attribute_values("name")
    }
    assert {"thrL", "thrA", "thrB", "thrC"} <= contains_names
    assert contains_stats.matching_records == len(contains)

    regex = indexed.query(r"^thr[ABCL]$", match="regex")
    regex_names = {
        value for record in regex for value in record.attribute_values("name")
    }
    assert regex_names == {"thrL", "thrA", "thrB", "thrC"}
    assert len(indexed.query(r"\S", match="regex")) >= len(regex)
    assert indexed.query(r"^missing$", match="regex") == []
    with pytest.raises(gai.GaiInputError, match="invalid regex"):
        indexed.query("missing[", match="regex")

    case_destination = tmp_path / "bed-case-sensitive.gai"
    gai.build_index(
        fixture_paths["bed_source"],
        fixture_paths["bed_tbi"],
        case_destination,
        case_sensitive=True,
    )
    case_indexed = gai.open_index(
        fixture_paths["bed_source"], fixture_paths["bed_tbi"], case_destination
    )
    assert len(case_indexed.query(r"^thr[ABCL]$", match="regex")) == 4
    assert case_indexed.query(r"^thr[abcl]$", match="regex") == []
    assert case_indexed.query("HR", match="contains") == []


def test_stale_corrupt_and_argument_errors(fixture_paths, tmp_path):
    destination = tmp_path / "fixture.gai"
    _build(fixture_paths, destination)

    stale_source = tmp_path / "stale.gff3.gz"
    stale_source.write_bytes(fixture_paths["source"].read_bytes() + b"\n")
    with pytest.raises(gai.GaiStaleError):
        gai.open_index(stale_source, fixture_paths["tbi"], destination)

    corrupt = tmp_path / "corrupt.gai"
    corrupt.write_bytes(
        destination.read_bytes()[:-1] + bytes([destination.read_bytes()[-1] ^ 0xFF])
    )
    with pytest.raises(gai.GaiCorruptError):
        gai.inspect_index(corrupt)

    with pytest.raises(gai.GaiInputError):
        gai.build_index(
            fixture_paths["source"], fixture_paths["tbi"], tmp_path / "empty.gai", []
        )
    with pytest.raises(gai.GaiInputError):
        gai.build_index(
            fixture_paths["source"],
            fixture_paths["tbi"],
            tmp_path / "bad-memory.gai",
            ["Name"],
            memory_budget=0,
        )
    with pytest.raises(TypeError):
        gai.build_index(
            fixture_paths["source"],
            fixture_paths["tbi"],
            tmp_path / "bad-type.gai",
            [1],
        )
