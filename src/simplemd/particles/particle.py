"""Point-particle data class for the Lennard-Jones model."""

import numpy as np


class Particle:
    """A single point particle: position ``r``, velocity ``rv``, acceleration ``ra``.

    The ``*1``/``*2``/``o`` buffers hold the previous-acceleration history and the
    saved position/velocity used by the predictor-corrector integrator. They are
    unused by the leapfrog and velocity-Verlet integrators.
    """

    def __init__(self, idx):
        self.idx = idx

        self.r = np.zeros(3)
        self.rv = np.zeros(3)
        self.ra = np.zeros(3)

        # predictor-corrector history (translational degrees of freedom)
        self.ra1 = np.zeros(3)
        self.ra2 = np.zeros(3)
        self.ro = np.zeros(3)
        self.rvo = np.zeros(3)
