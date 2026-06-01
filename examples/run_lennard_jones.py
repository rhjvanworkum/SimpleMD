"""Run a small Lennard-Jones soft-fluid simulation.

Self-contained and fast — no external data needed. Run with::

    uv run examples/run_lennard_jones.py

It runs a seeded simulation (so the output is reproducible), prints the final
energies and pressure, and reports the trajectory shape.
"""

import numpy as np

from simplemd import build_system


def main() -> None:
    sim = build_system(
        system="LJ",
        integrator="leapfrog",
        forces="PBP",
        N=4,  # 4**3 = 64 particles
        density=0.8,
        temperature=1.0,
        delta_t=0.005,
        steps=500,
        seed=0,
    )
    sim.simulate()

    total_momentum = np.linalg.norm(sum(m.rv for m in sim.mol))
    print(f"particles            : {sim.n_mol}")
    print(f"total energy / part. : {sim.tot_energy[0]:.4f}")
    print(f"kinetic energy / part: {sim.kin_energy[0]:.4f}")
    print(f"pressure             : {sim.pressure[0]:.4f}")
    print(f"|total momentum|     : {total_momentum:.2e} (should be ~0)")


if __name__ == "__main__":
    main()
