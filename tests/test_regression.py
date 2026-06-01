"""Characterization tests: outputs must match the pre-refactor code bit-for-bit.

The golden values in ``tests/data/golden_lj.json`` were captured from the
original code (before the src-layout refactor) with the global RNG seeded to
12345. These tests are the safety net that the refactor preserved behavior.
"""

import numpy as np
import pytest

from simplemd import run_simulation

CONFIGS = {
    "leapfrog_pbp": {"integrator": "leapfrog", "forces": "PBP"},
    "verlet_pbp": {"integrator": "verlet", "forces": "PBP"},
}


@pytest.mark.parametrize("name", sorted(CONFIGS))
def test_lj_trajectory_matches_golden(golden_lj, name):
    cfg = CONFIGS[name]
    params = golden_lj["_meta"]["params"]
    traj = run_simulation(
        system="LJ",
        integrator=cfg["integrator"],
        forces=cfg["forces"],
        N=params["N"],
        density=params["density"],
        temperature=params["temperature"],
        delta_t=params["delta_t"],
        steps=params["steps"],
        report_every=params["report_every"],
        seed=params["seed"],
    )
    expected = golden_lj[name]
    assert list(traj.shape) == expected["traj_shape"]
    # exact match (atol=0, rtol=0): the refactor must not perturb the numerics
    np.testing.assert_array_equal(traj[-1], np.array(expected["traj_last"]))
