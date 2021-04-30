"""
Unused for now
"""

import numpy as np
from objects.Molecule import Molecule
from models.system import System

class chain_model(System):

    def __init__(self, settings, integrator, forces):
        System.__init__(self, settings, integrator, forces)

        self.chain_dim = np.array([1, 1, 1])
        self.n_chain =  2 * self.chain_dim.prod()
        if (self.n_chain == 2): self.n_chain = 1  #in case we only study one chain
        self.chain_length = 8
        self.bond_limit = 2.1
        self.n_mol = self.cell_dim.prod() + self.n_chain * self.chain_length
        self.mol = [Molecule(self.n_dim, i) for i in range(self.n_mol)]

        self.init_coords_chain()


    def init_coords_chain(self):

        by = self.r_cut * np.cos(np.pi / 4)
        bz = self.r_cut * np.sin(np.pi / 4)
        n: int = 0

        if (self.n_chain == 1):
            for m in range(self.chain_length):
                self.mol[n].r = np.array([0, (m % 2) * by, m * bz])
                self.mol[n].r += -0.25 * self.region
                n += 1
        else:

            self.cell_size = self.region / self.chain_dim

            for z in range(self.chain_dim[2]):
                for y in range(self.chain_dim[1]):
                    for x in range(self.chain_dim[0]):
                        c = np.array([x + 0.25, y + 0.25, z + 0.25])
                        c *= self.cell_size
                        c += -0.5 * self.region
                        for j in range(2):
                            for m in range(self.chain_length):
                                self.mol[n].r = np.array([0, (m % 2) * by, m * bz])
                                self.mol[n].r += 0.5 * j * self.cell_size
                                self.mol[n].r += c
                                n += 1

            self.n_mol = n
            self.apply_boundary_conditions()

            self.cell_size = self.region / self.cell_dim

            for z in range(self.cell_dim[2]):
                for y in range(self.cell_dim[1]):
                    for x in range(self.cell_dim[0]):
                        c = np.array([x + 0.5, y + 0.5, z + 0.5])
                        c *= self.cell_size
                        c += -0.5 * self.cell_size
                        for i in range(self.n_chain * self.chain_length):
                            dir = self.mol[i].r - c
                            if (np.dot(dir, dir) < self.r_cut**2):
                                break
                        if (i == self.n_chain * self.chain_length):
                            self.mol[n].r = c
                            n += 1

            self.n_mol = n

        # assigning to chains
        n = 0
        for i in range(self.n_chain):
            for j in range(self.chain_length):
                self.mol[n].in_chain = i
                n += 1

        for n in range(self.n_chain * self.chain_length, self.n_mol):
            self.mol[n].in_chain = -1