"""Run a small TIP-4P water simulation and save the trajectory.

Self-contained and fast. Run with::

    uv run examples/run_tip4p_water.py

This demonstrates the rigid-body (quaternion) water model. The saved trajectory
JSON can be replayed with the optional viewer (requires the ``viz`` extra)::

    uv sync --extra viz
    uv run python -m simplemd.viz tip4p_trajectory.json
"""

import json

import numpy as np

from simplemd import run_simulation


def main() -> None:
    traj = run_simulation(
        system="TIP4P",
        integrator="leapfrog",
        forces="PBP",
        N=2,  # 2**3 = 8 water molecules
        density=0.8,
        temperature=1.0,
        delta_t=0.001,
        steps=200,
        report_every=2,
        seed=0,
    )

    print(f"trajectory shape : {traj.shape}  (frames, sites, dims)")
    print(f"all finite       : {bool(np.isfinite(traj).all())}")

    out = "tip4p_trajectory.json"
    with open(out, "w") as f:
        json.dump(traj.tolist(), f)
    print(f"saved trajectory : {out}")


if __name__ == "__main__":
    main()
