from integrators.thermostats import apply_thermostat
import numpy as np

class VerletVelocityIntegrator():

    def __init__(self, parent, delta_t):
        self.system = parent
        self.delta_t = delta_t

    def integrate(self):
        # calc new positions
        for i in range(self.system.n_mol):
            self.system.mol[i].r += self.delta_t * self.system.mol[i].rv + 0.5 * self.delta_t * self.system.mol[i].ra

        # make sure the periodic boundary conditions are applied
        self.system.apply_boundary_conditions()

        # calc intermediate vels
        rv_intermediate = np.zeros((self.system.n_mol, len(self.system.mol[0].rv)))
        for i in range(self.system.n_mol):
            rv_intermediate[i] = self.system.mol[i].rv + 0.5 * self.delta_t * self.system.mol[i].ra

        # compute forces
        self.system.calculate_forces()

        # calc the new velocities
        for i in range(self.system.n_mol):
            self.system.mol[i].rv = rv_intermediate[i] + 0.5 * self.delta_t * self.system.mol[i].ra

        # apply the thermostat correctly
        apply_thermostat(self.delta_t, self.system.n_mol, self.system.mol)

        # make sure the periodic boundary conditions are applied
        self.system.apply_boundary_conditions()