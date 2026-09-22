from __future__ import annotations

from pathlib import Path

import pytest


@pytest.fixture(scope="session")
def fixture_paths() -> dict[str, Path]:
    """Create a deterministic BGZF source and matching TBI/CSI indexes."""
    fixtures = Path(__file__).resolve().parents[0] / "fixtures"
    return {
        "source": fixtures / "fixture.gff3.gz",
        "tbi": fixtures / "fixture.gff3.gz.tbi",
        "csi": fixtures / "fixture.gff3.gz.csi",
        "unsorted_gff": fixtures / "unsorted_fixture.gff3.gz",
    }
