"""Tests for system construction and the TIP-4P model."""

import numpy as np
import pytest

from simplemd import build_system, run_simulation


def test_lj_system_geometry():
    sim = build_system(system="LJ", N=3, density=0.8, seed=0)
    assert sim.n_mol == 27  # 3**3 cells, one particle each
    assert sim.n_dim == 3
    assert len(sim.mol) == 27


def test_lj_initial_positions_inside_region():
    sim = build_system(system="LJ", N=3, density=0.8, seed=0)
    half = 0.5 * sim.region
    for m in sim.mol:
        assert np.all(m.r >= -half - 1e-9)
        assert np.all(m.r <= half + 1e-9)


def test_tip4p_runs_and_is_finite():
    # Regression for the wiring fix: TIP-4P used to receive the wrong force fn.
    traj = run_simulation(
        system="TIP4P",
        integrator="leapfrog",
        forces="PBP",
        N=2,
        density=0.8,
        temperature=1.0,
        delta_t=0.001,
        steps=10,
        report_every=1,
        seed=7,
    )
    # 2**3 = 8 molecules, 3 reported sites each
    assert traj.shape == (10, 24, 3)
    assert np.isfinite(traj).all()


@pytest.mark.parametrize("integrator", ["leapfrog", "verlet"])
def test_lj_runs_with_each_integrator(integrator):
    traj = run_simulation(
        system="LJ", integrator=integrator, forces="PBP", N=3, steps=10, seed=0
    )
    assert traj.shape[0] == 10
    assert np.isfinite(traj).all()


def test_tip4p_runs_with_predictor_corrector():
    traj = run_simulation(
        system="TIP4P",
        integrator="pred-corr",
        forces="PBP",
        N=2,
        density=0.8,
        temperature=1.0,
        delta_t=0.001,
        steps=5,
        seed=1,
    )
    assert traj.shape == (5, 24, 3)
    assert np.isfinite(traj).all()
