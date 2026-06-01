"""Regression tests that *document* known, pre-existing bugs.

These behaviors were present before the modernization refactor and are
intentionally left unfixed (out of the agreed scope). The tests pin the current
behavior so the bugs stay visible and any future fix is a deliberate, reviewed
change rather than a silent one. See the PR description for details.
"""

import pytest

from simplemd import run_simulation


def test_cell_division_force_is_broken():
    """`CD` force evaluation raises KeyError: edge particles bin to an
    out-of-range cell index (off-by-one in the cell-list construction)."""
    with pytest.raises(KeyError):
        run_simulation(system="LJ", integrator="leapfrog", forces="CD", N=3, steps=2, seed=0)


def test_predictor_corrector_is_broken_for_lennard_jones():
    """`pred-corr` assumes the higher-order history buffers (ra1, ra2, ...) that
    only the Molecule class defines, so it fails for the point-particle LJ model.
    It works for the TIP-4P model (see test_systems)."""
    with pytest.raises(AttributeError):
        run_simulation(system="LJ", integrator="pred-corr", forces="PBP", N=3, steps=2, seed=0)
