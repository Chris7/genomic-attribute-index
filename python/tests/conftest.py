from __future__ import annotations

from pathlib import Path

import pytest


@pytest.fixture(scope="session")
def fixture_paths() -> dict[str, Path]:
    """Create a deterministic BGZF source and matching TBI/CSI indexes."""
    fixtures = Path(__file__).resolve().parents[0] / "fixtures"
    bed_fixtures = Path(__file__).resolve().parents[2] / "fixtures" / "sorting"
    return {
        "source": fixtures / "fixture.gff3.gz",
        "tbi": fixtures / "fixture.gff3.gz.tbi",
        "csi": fixtures / "fixture.gff3.gz.csi",
        "unsorted_gff": fixtures / "unsorted_fixture.gff3.gz",
        "bed_source": bed_fixtures / "ecoli_k12_mg1655.bed.gz",
        "bed_tbi": bed_fixtures / "ecoli_k12_mg1655.bed.gz.tbi",
    }
