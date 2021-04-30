"""
Currently not used
"""

from interaction_methods.pair_by_pair import ComputeForcesPBP_RK
from integrators.thermostats import apply_thermostat
import numpy as np


class RungeKuttaIntegrator():

    def __init__(self, parent, delta_t):
        self.system = parent
        self.delta_t = delta_t

    def integrate(self):

        X = np.zeros((2, self.system.n_mol, len(self.system.mol[0].r)))
        for i in range(self.system.n_mol):
            X[0][i] = self.system.mol[i].r
            X[1][i] = self.system.mol[i].rv

        k1 = self.RHS(X)
        k2 = self.RHS(X + self.delta_t / 2.0 * k1)
        k3 = self.RHS(X + self.delta_t / 2.0 * k2)
        k4 = self.RHS(X + self.delta_t * k3)

        X += self.delta_t / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4)

        for i in range(self.system.n_mol):
            self.system.mol[i].r = X[0][i]
            self.system.mol[i].rv = X[1][i]

        # make sure the periodic boundary conditions are applied
        self.system.apply_boundary_conditions()

        # apply the thermostat correctly
        apply_thermostat(self.delta_t, self.system.n_mol, self.system.mol)

        # compute forces
        self.system.calculate_forces()

    def RHS(self, X):
        output = np.zeros((2, len(X[0]), len(X[0][0])))

        # dx/dt = v
        # dv/dt = f(x)

        A = ComputeForcesPBP_RK(self.system, X[0])

        # looping over all particles
        for i in range(len(X[0])):
            output[0][i] = X[1][i]
            output[1][i] = A[i]

        return output