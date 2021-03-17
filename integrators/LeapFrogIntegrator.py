from integrators.thermostats import apply_thermostat, apply_thermostat_nonlinear

from objects.Molecule import compute_torq, compute_accels_q

class LeapFrogIntegrator():

    def __init__(self, parent, delta_t):
        self.system = parent
        self.delta_t = delta_t

    def integrate(self):
        # step 1 of the integration
        for i in range(self.system.n_mol):
            self.system.mol[i].rv += 0.5 * self.delta_t * self.system.mol[i].ra
            self.system.mol[i].r += self.delta_t * self.system.mol[i].rv

        # make sure the periodic boundary conditions are applied
        self.system.apply_boundary_conditions()

        # calculate the forces in the system now
        self.system.calculate_forces()

        # apply the thermostat correctly
        apply_thermostat(self.delta_t, self.system.n_mol, self.system.mol)

        # step 2 of te integration
        for i in range(self.system.n_mol):
            self.system.mol[i].rv += 0.5 * self.delta_t * self.system.mol[i].ra

    def integrate_nonlinear(self):
        # step 1 of the integration
        for i in range(self.system.n_mol):
            self.system.mol[i].rv += 0.5 * self.delta_t * self.system.mol[i].ra
            self.system.mol[i].r += self.delta_t * self.system.mol[i].rv

            self.system.mol[i].qv += 0.5 * self.delta_t * self.system.mol[i].qa
            self.system.mol[i].q += self.delta_t * self.system.mol[i].qv

        # make sure the periodic boundary conditions are applied
        self.system.adjust_quat()
        self.system.apply_boundary_conditions()

        # calculate the site coordinates
        self.system.gen_site_coords()

        # calculate the forces in the system now
        self.system.calculate_forces()

        for i in range(self.system.n_mol):
            compute_torq(self.system.mol[i], self.system.sites)
            compute_accels_q(self.system.mol[i])

        # apply the thermostat correctly
        apply_thermostat_nonlinear(self.delta_t, self.system.n_mol, self.system.mol)

        # step 2 of te integration
        for i in range(self.system.n_mol):
            self.system.mol[i].rv += 0.5 * self.delta_t * self.system.mol[i].ra

            self.system.mol[i].qv += 0.5 * self.delta_t * self.system.mol[i].qa