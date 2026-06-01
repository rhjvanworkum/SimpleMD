"""Regression tests for three previously-broken behaviors that are now fixed.

These guard against re-introducing the bugs documented (and previously pinned as
"known bugs") during the modernization refactor:

1. cell-division force evaluation (`CD`) used to raise KeyError for every config;
2. the predictor-corrector integrator crashed for the point-particle LJ model;
3. the TIP-4P site forces violated Newton's third law.
"""

import numpy as np
import pytest

from simplemd import build_system, run_simulation
from simplemd.interactions.cell_division import ComputeForcesCD
from simplemd.interactions.pair_by_pair import ComputeForcesPBP, ComputeForcesPBP_nonlinear

# --- Bug 1: cell-division -----------------------------------------------------


@pytest.mark.parametrize("N", [4, 5, 6])
def test_cell_division_forces_match_pair_by_pair(N):
    """CD must produce the same Lennard-Jones forces as the reference PBP method."""
    sa = build_system(system="LJ", forces="PBP", N=N, density=0.8, seed=0)
    sb = build_system(system="LJ", forces="CD", N=N, density=0.8, seed=0)
    ComputeForcesPBP(sa)
    ComputeForcesCD(sb)
    ra_pbp = np.array([m.ra for m in sa.mol])
    ra_cd = np.array([m.ra for m in sb.mol])
    np.testing.assert_allclose(ra_cd, ra_pbp, atol=1e-12)
    assert sb.u_sum == pytest.approx(sa.u_sum)
    assert sb.vir_sum == pytest.approx(sa.vir_sum)


def test_cell_division_trajectory_matches_pair_by_pair():
    a = run_simulation(system="LJ", integrator="leapfrog", forces="PBP", N=4, steps=50, seed=0)
    b = run_simulation(system="LJ", integrator="leapfrog", forces="CD", N=4, steps=50, seed=0)
    np.testing.assert_allclose(b, a, atol=1e-10)


def test_cell_division_raises_for_too_small_system():
    """Fewer than 3 cells per dimension is not a valid cell list -> clear error."""
    sim = build_system(system="LJ", forces="CD", N=3, density=0.8, seed=0)
    with pytest.raises(ValueError, match="at least 3 cells"):
        ComputeForcesCD(sim)


# --- Bug 2: predictor-corrector for the LJ model ------------------------------


def test_predictor_corrector_runs_for_lennard_jones():
    traj = run_simulation(system="LJ", integrator="pred-corr", forces="PBP", N=3, steps=50, seed=0)
    assert traj.shape == (50, 27, 3)
    assert np.isfinite(traj).all()


# --- Bug 3: Newton's third law in the TIP-4P site forces ----------------------


def test_tip4p_site_forces_conserve_momentum():
    """Equal-and-opposite pair forces => the net force on the system is ~zero."""
    sim = build_system(
        system="TIP4P", forces="PBP", N=2, density=0.8, temperature=1.0, delta_t=0.001, seed=7
    )
    ComputeForcesPBP_nonlinear(sim)
    net_force = sum(site.f for site in sim.sites)
    np.testing.assert_allclose(net_force, np.zeros(3), atol=1e-10)
