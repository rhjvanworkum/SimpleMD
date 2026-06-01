"""Shared test fixtures."""

import json
from pathlib import Path

import pytest

DATA_DIR = Path(__file__).parent / "data"


@pytest.fixture(scope="session")
def golden_lj():
    """Golden Lennard-Jones outputs captured from the pre-refactor code."""
    with open(DATA_DIR / "golden_lj.json") as f:
        return json.load(f)


@pytest.fixture
def small_lj_kwargs():
    """A small, fast, reproducible Lennard-Jones configuration."""
    return {
        "system": "LJ",
        "integrator": "leapfrog",
        "forces": "PBP",
        "N": 3,
        "density": 0.8,
        "temperature": 1.0,
        "delta_t": 0.005,
        "steps": 50,
        "report_every": 1,
        "seed": 12345,
    }
