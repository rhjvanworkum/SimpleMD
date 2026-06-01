"""Tests for the high-level Python API."""

import numpy as np
import pytest

from simplemd import build_system, run_simulation
from simplemd.api import resolve_forces
from simplemd.interactions.cell_division import ComputeForcesCD
from simplemd.interactions.pair_by_pair import ComputeForcesPBP, ComputeForcesPBP_nonlinear


def test_run_simulation_returns_expected_shape(small_lj_kwargs):
    traj = run_simulation(**small_lj_kwargs)
    # 3**3 = 27 particles, 50 recorded frames, 3 dimensions
    assert traj.shape == (50, 27, 3)
    assert np.isfinite(traj).all()


def test_run_simulation_is_reproducible_with_seed(small_lj_kwargs):
    a = run_simulation(**small_lj_kwargs)
    b = run_simulation(**small_lj_kwargs)
    np.testing.assert_array_equal(a, b)


def test_run_simulation_differs_across_seeds(small_lj_kwargs):
    a = run_simulation(**{**small_lj_kwargs, "seed": 1})
    b = run_simulation(**{**small_lj_kwargs, "seed": 2})
    assert not np.array_equal(a, b)


def test_lj_conserves_momentum(small_lj_kwargs):
    sim = build_system(**{k: v for k, v in small_lj_kwargs.items() if k != "report_every"})
    sim.simulate()
    total_momentum = sum(m.rv for m in sim.mol)
    np.testing.assert_allclose(total_momentum, np.zeros(3), atol=1e-9)


@pytest.mark.parametrize("bad", [{"system": "nope"}, {"integrator": "nope"}, {"forces": "nope"}])
def test_build_system_rejects_unknown_keys(small_lj_kwargs, bad):
    kwargs = {k: v for k, v in small_lj_kwargs.items() if k != "report_every"}
    kwargs.update(bad)
    with pytest.raises(ValueError):
        build_system(**kwargs)


def test_resolve_forces_lj_uses_point_particle_function():
    assert resolve_forces("LJ", "PBP") is ComputeForcesPBP
    assert resolve_forces("LJ", "CD") is ComputeForcesCD


def test_resolve_forces_tip4p_uses_site_based_function():
    # The wiring fix: TIP-4P must route to the nonlinear (site-based) force fn.
    assert resolve_forces("TIP4P", "PBP") is ComputeForcesPBP_nonlinear


def test_resolve_forces_tip4p_rejects_cell_division():
    with pytest.raises(ValueError):
        resolve_forces("TIP4P", "CD")
