"""Point-particle data class for the Lennard-Jones model."""

import numpy as np


class Particle:
    """A single point particle: position ``r``, velocity ``rv``, acceleration ``ra``."""

    def __init__(self, idx):
        self.idx = idx

        self.r = np.zeros(3)
        self.rv = np.zeros(3)
        self.ra = np.zeros(3)
