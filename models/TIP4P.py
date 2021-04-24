import numpy as np

from models.system import System
import utils as utils
import integrators.thermostats as thermostat
from objects.Molecule import Molecule, Site, MSite

class TIP4P(System):

    def __init__(self, settings, integrator, forces):
        System.__init__(self, settings, integrator, forces)

        self.mol = [Molecule(i) for i in range(self.n_mol)]

        self.n_sites = 4
        self.molecule_sites = [MSite() for i in range(self.n_sites)]
        self.sites = [Site() for i in range(self.n_sites * self.n_mol)]

        self.init_coords()
        self.init_vels()

        self.define_mol()
        self.init_ang_coords()
        self.init_ang_vels()
        self.gen_site_coords()


    def define_mol(self):
        self.molecule_sites[0].r[2] = -0.0206
        self.molecule_sites[1].r[2] = 0.0274
        self.molecule_sites[2].r[1] = 0.240
        self.molecule_sites[2].r[2] = 0.165
        self.molecule_sites[3].r[1] = - self.molecule_sites[2].r[1]
        self.molecule_sites[3].r[2] = self.molecule_sites[2].r[2]

        for n in range(self.n_mol):
            self.mol[n].m_inert = np.array([0.00980, 0.00340, 0.00640])

        ### this is a specific parameter for the water example
        self.b = 183.5

        self.molecule_sites[0].typeF = 1
        self.molecule_sites[1].typeF = 2
        self.molecule_sites[2].typeF = 3
        self.molecule_sites[3].typeF = 3

    """ rotating the molecular sites from body coordinates to world coordinates """
    def gen_site_coords(self):
        for n in range(self.n_mol):
            rMat = utils.build_rot_matrix(self.mol[n].q, 1)
            for j in range(self.n_sites):
                t = utils.Mproduct(rMat, self.molecule_sites[j].r)
                self.sites[self.n_sites * n + j].r = self.mol[n].r + t

    """ function that actually runs the simulation experiment """
    def simulate(self):
        running = True
        while (running):
            self.step_count += 1

            self.integrator.integrate_nonlinear()

            self.eval_props()
            self.accum_props(1)

            if (self.step_count % self.step_avg == 0):
                self.accum_props(2)
                if self.show_summary: self.print_summary()
                self.accum_props(0)

            if (self.step_count % self.step_adjust_temp == 0):
                thermostat.adjust_temp_nonlin(self)

            if (self.show_progress):
                if (self.step_count % (self.step_limit / 1000) == 0):
                    print('{:.1f}'.format(self.step_count / self.step_limit * 100) + '%')

            if self.report: self.reporter.step()

            if (self.step_count >= self.step_limit):
                # for the backend we just request this array in the API call
                # if self.report: self.reporter.export()
                running = False