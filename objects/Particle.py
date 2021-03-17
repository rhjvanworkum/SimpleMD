import numpy as np

class Particle():

    def __init__(self, idx):
        self.idx = idx

        self.r = np.zeros(3)
        self.rv = np.zeros(3)
        self.ra = np.zeros(3)