"""Command-line interface: run a simulation and optionally save the trajectory."""

from __future__ import annotations

import argparse
import json

from simplemd.api import INTEGRATORS, SYSTEMS, run_simulation


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="simplemd", description="Run a molecular dynamics simulation."
    )
    parser.add_argument("--system", choices=sorted(SYSTEMS), default="LJ")
    parser.add_argument("--integrator", choices=sorted(INTEGRATORS), default="leapfrog")
    parser.add_argument("--forces", choices=["PBP", "CD"], default="PBP")
    parser.add_argument("-N", type=int, default=3, help="unit cells per dimension")
    parser.add_argument("--density", type=float, default=0.8)
    parser.add_argument("--temperature", type=float, default=1.0)
    parser.add_argument("--delta-t", type=float, default=0.005)
    parser.add_argument("--steps", type=int, default=1000)
    parser.add_argument("--report-every", type=int, default=1)
    parser.add_argument("--seed", type=int, default=None, help="RNG seed for a reproducible run")
    parser.add_argument("-o", "--output", default=None, help="write trajectory JSON to this file")
    args = parser.parse_args(argv)

    traj = run_simulation(
        system=args.system,
        integrator=args.integrator,
        forces=args.forces,
        N=args.N,
        density=args.density,
        temperature=args.temperature,
        delta_t=args.delta_t,
        steps=args.steps,
        report_every=args.report_every,
        seed=args.seed,
    )

    if args.output:
        with open(args.output, "w") as f:
            json.dump(traj.tolist(), f)
        print(f"wrote trajectory {traj.shape} to {args.output}")
    else:
        print(f"simulation complete: trajectory shape {traj.shape}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
