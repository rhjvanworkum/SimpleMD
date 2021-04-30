"""
Simple Lennard-Jones soft fluid disk model
"""

from models.system import System
from objects.Particle import Particle
import integrators.thermostats as thermostat

class gas_particles(System):

    def __init__(self, settings, integrator, forces):
        System.__init__(self, settings, integrator, forces)

        self.mol = [Particle(i) for i in range(self.n_mol)]

        self.init_coords()
        self.init_vels()

    """ function that actually runs the simulation experiment """
    def simulate(self):
        running = True
        while (running):
            self.step_count += 1

            self.integrator.integrate()

            self.eval_props()
            self.accum_props(1)

            if (self.step_count % self.step_avg == 0):
                self.accum_props(2)
                if self.show_summary: self.print_summary()
                self.accum_props(0)

            if (self.step_count % self.step_adjust_temp == 0):
                thermostat.adjust_temp(self)

            if (self.show_progress):
                if (self.step_count % (self.step_limit / 1000) == 0):
                    print('{:.1f}'.format(self.step_count / self.step_limit * 100) + '%')

            if self.report: self.reporter.step()

            if (self.step_count >= self.step_limit):
                # for the backend we just request this array in the API call
                if self.report and self.reporter.output_file: self.reporter.export()
                running = False