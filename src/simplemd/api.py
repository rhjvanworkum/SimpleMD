"""High-level Python API for running simulations.

This is the importable entry point::

    from simplemd import run_simulation
    traj = run_simulation(system="LJ", integrator="leapfrog", forces="PBP",
                          N=3, steps=500, seed=0)

It also exposes the registries (:data:`SYSTEMS`, :data:`INTEGRATORS`,
:data:`FORCES`) shared with the HTTP API so the two surfaces stay in sync.
"""

from __future__ import annotations

import numpy as np

from simplemd.config import Settings
from simplemd.integrators.leapfrog import LeapFrogIntegrator
from simplemd.integrators.predictor_corrector import PredictorCorrectorIntegrator
from simplemd.integrators.velocity_verlet import VerletVelocityIntegrator
from simplemd.interactions.cell_division import ComputeForcesCD
from simplemd.interactions.pair_by_pair import ComputeForcesPBP, ComputeForcesPBP_nonlinear
from simplemd.reporter import Reporter
from simplemd.systems.lennard_jones import gas_particles
from simplemd.systems.tip4p import TIP4P

SYSTEMS = {
    "LJ": gas_particles,
    "TIP4P": TIP4P,
}

INTEGRATORS = {
    "leapfrog": LeapFrogIntegrator,
    "verlet": VerletVelocityIntegrator,
    "pred-corr": PredictorCorrectorIntegrator,
}

# Point-particle force functions, keyed as exposed by the HTTP API.
FORCES = {
    "PBP": ComputeForcesPBP,
    "CD": ComputeForcesCD,
}


def resolve_forces(system: str, forces: str):
    """Return the force function appropriate for ``system``/``forces``.

    The TIP-4P water model needs the *site-based* force routine rather than the
    point-particle one. Historically the HTTP layer passed the point-particle
    function to TIP-4P, which never produced correct site forces; this mapping
    fixes that wiring while keeping the public ``forces`` keys unchanged.
    """
    if system == "TIP4P":
        if forces == "PBP":
            return ComputeForcesPBP_nonlinear
        raise ValueError(f"force method {forces!r} is not supported for the TIP4P model")
    if forces not in FORCES:
        raise ValueError(f"unknown force method {forces!r}; choose from {sorted(FORCES)}")
    return FORCES[forces]


def build_system(
    *,
    system: str = "LJ",
    integrator: str = "leapfrog",
    forces: str = "PBP",
    N: int = 3,
    density: float = 0.8,
    temperature: float = 1.0,
    delta_t: float = 0.005,
    steps: int = 1000,
    step_avg: int = 100,
    step_adjust_temp: int = 10,
    step_equilibrium: int = 20,
    show_progress: bool = False,
    show_summary: bool = False,
    seed: int | None = None,
):
    """Construct (but do not run) a configured system instance."""
    if system not in SYSTEMS:
        raise ValueError(f"unknown system {system!r}; choose from {sorted(SYSTEMS)}")
    if integrator not in INTEGRATORS:
        raise ValueError(f"unknown integrator {integrator!r}; choose from {sorted(INTEGRATORS)}")

    settings = Settings(
        N=N,
        density=density,
        temperature=temperature,
        delta_t=delta_t,
        step_limit=steps,
        step_avg=step_avg,
        step_adjust_temp=step_adjust_temp,
        step_equilibrium=step_equilibrium,
        show_progress=show_progress,
        show_summary=show_summary,
        seed=seed,
    )
    force_fn = resolve_forces(system, forces)
    return SYSTEMS[system](settings, INTEGRATORS[integrator], force_fn)


def run_simulation(*, report_every: int = 1, **kwargs) -> np.ndarray:
    """Run a simulation and return the recorded trajectory.

    Accepts the same keyword arguments as :func:`build_system` plus
    ``report_every`` (how often, in steps, to record a frame). Returns the
    trajectory as a NumPy array of shape ``(n_frames, n_bodies, n_dim)``.
    """
    system_key = kwargs.get("system", "LJ")
    sim = build_system(**kwargs)
    sim.add_reporter(Reporter, report_every, system_key)
    sim.simulate()
    return sim.reporter.traj
