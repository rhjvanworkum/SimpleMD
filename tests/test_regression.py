"""Characterization tests: outputs must match the pre-refactor code.

The golden values in ``tests/data/golden_lj.json`` were captured from the
original code (before the src-layout refactor) with the global RNG seeded to
12345. These tests are the safety net that the refactor preserved behavior.

The comparison uses a tight tolerance rather than exact equality: on a fixed
machine the refactored code reproduces the golden values bit-for-bit, but a
checked-in golden file is also compared against runs on other machines (CI),
where a different NumPy/BLAS build or CPU perturbs the last bit (~1e-16). The
tolerance here is still ~7 orders of magnitude tighter than any genuine
algorithmic regression would produce.
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
    # tight tolerance (portable across platforms; see module docstring)
    np.testing.assert_allclose(traj[-1], np.array(expected["traj_last"]), rtol=1e-9, atol=1e-12)
