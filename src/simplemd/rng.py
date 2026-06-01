"""Random-number-generator seeding.

The simulation draws initial positions, velocities and orientations from the
global :mod:`numpy.random` and :mod:`random` generators. Seeding them makes a
run fully reproducible, which is what the test-suite relies on. Seeding is
*opt-in*: when no seed is configured the generators are left untouched, so the
default (non-deterministic) behavior is preserved.
"""

import random

import numpy as np


def seed_rng(seed: int) -> None:
    """Seed both the NumPy and stdlib global RNGs with ``seed``."""
    np.random.seed(seed)
    random.seed(seed)
