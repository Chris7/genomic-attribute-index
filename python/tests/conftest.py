from __future__ import annotations

import subprocess
from pathlib import Path

import pytest


@pytest.fixture(scope="session")
def fixture_paths(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Path]:
    """Create a deterministic BGZF source and matching TBI/CSI indexes."""
    root = Path(__file__).resolve().parents[2]
    directory = tmp_path_factory.mktemp("gai-fixture")
    subprocess.run(
        [
            "cargo",
            "run",
            "--locked",
            "--quiet",
            "--example",
            "python_fixture",
            "--",
            str(directory),
        ],
        cwd=root,
        check=True,
    )
    return {
        "root": root,
        "source": directory / "fixture.gff3.gz",
        "tbi": directory / "fixture.gff3.gz.tbi",
        "csi": directory / "fixture.gff3.gz.csi",
    }
