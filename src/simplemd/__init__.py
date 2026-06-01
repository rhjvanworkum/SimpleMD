"""SimpleMD — a simple molecular dynamics simulator.

Two physical models are provided: a Lennard-Jones soft fluid (``"LJ"``) and a
TIP-4P water model (``"TIP4P"``). The recommended entry point is
:func:`run_simulation`.
"""

from simplemd.api import (
    FORCES,
    INTEGRATORS,
    SYSTEMS,
    build_system,
    run_simulation,
)
from simplemd.config import Settings
from simplemd.rng import seed_rng

__version__ = "0.1.0"

__all__ = [
    "FORCES",
    "INTEGRATORS",
    "SYSTEMS",
    "Settings",
    "__version__",
    "build_system",
    "run_simulation",
    "seed_rng",
]
