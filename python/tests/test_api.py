from __future__ import annotations

from pathlib import Path

import pytest

import gni


def test_distribution_version_is_exposed():
    assert gni.__version__ == "0.1.0"


def _build(
    paths: dict[str, Path], output: Path, *, case_sensitive: bool = False
) -> gni.BuildStats:
    return gni.build_index(
        paths["source"],
        paths["tbi"],
        output,
        ["Name", "Alias", "Name"],
        case_sensitive=case_sensitive,
        memory_budget=1,
        compression_threads=1,
        bgzf_threads=1,
    )


def test_build_query_inspect_tbi_and_source_order(fixture_paths, tmp_path):
    destination = tmp_path / "fixture.gni"
    stats = _build(fixture_paths, destination)
    assert stats.records_processed == 5
    assert stats.records_indexed == 5
    assert stats.distinct_terms >= 4
    assert stats.total_seconds >= 0

    indexed = gni.open_index(fixture_paths["source"], fixture_paths["tbi"], destination)
    metadata = indexed.metadata()
    assert metadata.attributes == ["Name", "Alias"]
    assert metadata.major_version == 1
    assert len(metadata.gff_fingerprint) == 64

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
        gni.query_index(
            fixture_paths["source"], fixture_paths["tbi"], destination, "other gene"
        )[0].start
        == 30
    )

    inspected = gni.inspect_index(destination)
    assert inspected.file_size == metadata.file_size
    assert inspected.span_data_bytes >= 0


def test_csi_case_sensitive_unknown_and_no_id(fixture_paths, tmp_path):
    destination = tmp_path / "fixture-csi.gni"
    gni.build_index(
        fixture_paths["source"],
        fixture_paths["csi"],
        destination,
        ["Name", "Alias"],
        case_sensitive=True,
        compression_threads=1,
        bgzf_threads=1,
    )
    indexed = gni.open_index(fixture_paths["source"], fixture_paths["csi"], destination)
    assert len(indexed.query("Alpha")) == 4
    assert indexed.query("alpha") == []
    assert indexed.query("missing") == []

    id_only = tmp_path / "id-only.gni"
    gni.build_index(fixture_paths["source"], fixture_paths["tbi"], id_only, ["ID"])
    assert (
        gni.open_index(fixture_paths["source"], fixture_paths["tbi"], id_only).query(
            "Alpha"
        )
        == []
    )


def test_stale_corrupt_and_argument_errors(fixture_paths, tmp_path):
    destination = tmp_path / "fixture.gni"
    _build(fixture_paths, destination)

    stale_source = tmp_path / "stale.gff3.gz"
    stale_source.write_bytes(fixture_paths["source"].read_bytes() + b"\n")
    with pytest.raises(gni.GniStaleError):
        gni.open_index(stale_source, fixture_paths["tbi"], destination)

    corrupt = tmp_path / "corrupt.gni"
    corrupt.write_bytes(
        destination.read_bytes()[:-1] + bytes([destination.read_bytes()[-1] ^ 0xFF])
    )
    with pytest.raises(gni.GniCorruptError):
        gni.inspect_index(corrupt)

    with pytest.raises(gni.GniInputError):
        gni.build_index(
            fixture_paths["source"], fixture_paths["tbi"], tmp_path / "empty.gni", []
        )
    with pytest.raises(gni.GniInputError):
        gni.build_index(
            fixture_paths["source"],
            fixture_paths["tbi"],
            tmp_path / "bad-memory.gni",
            ["Name"],
            memory_budget=0,
        )
    with pytest.raises(TypeError):
        gni.build_index(
            fixture_paths["source"],
            fixture_paths["tbi"],
            tmp_path / "bad-type.gni",
            [1],
        )
