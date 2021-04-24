import numpy as np
import random
import utils as utils


class System():

    """ initiating the system, 1. initiates the meassured the properties 2. scales the cellsize of the system """
    """ 3. sets the step count 4. initiates the positions, velocities and accelerations """
    def __init__(self, settings, integrator, forces):
        self.n_dim = 3

        self.step_limit = settings.step_limit
        self.step_avg = settings.step_avg
        self.step_adjust_temp = settings.step_adjust_temp
        self.step_equilibrium = settings.step_equilibrium

        self.show_progress = settings.show_progress
        self.show_summary = settings.show_summary

        self.integrator = integrator(self, settings.delta_t)
        self.forces = forces
        self.report = False

        self.v_sum = 0
        self.u_sum = 0
        self.vir_sum = 0

        self.tot_energy = np.zeros(3)
        self.kin_energy = np.zeros(3)
        self.pressure = np.zeros(3)

        self.cell_dim = np.ones(3).astype(int) * settings.N
        self.density = settings.density
        self.cell_size = 1 / np.sqrt(settings.density)
        self.region = self.cell_dim * self.cell_size

        self.n_mol = self.cell_dim.prod()
        self.mol = []

        self.temperature = settings.temperature
        self.r_cut = np.power(2, (1 / 6))
        self.vel_mag = np.sqrt(self.n_dim * (1 - 1 / self.n_mol) * self.temperature)

        self.step_count = 0

    """ adding a reporter to the system """
    def add_reporter(self, reporter, step_report_avg, model, output_file=None):
        self.report = True
        self.reporter = reporter(self, step_report_avg, model, output_file)

    """ Initiating cartesian coordinates """
    def init_coords(self):
        i = 0
        n_x, n_y, n_z = self.cell_dim[0], self.cell_dim[1], self.cell_dim[2]
        for x in range(n_x):
            for y in range(n_y):
                for z in range(n_z):
                    # put the position in the middle of the unit cell
                    position = np.array([x + 0.5, y + 0.5, z + 0.5])
                    # scale the position with the cell size dimensions
                    position *= self.cell_size
                    # center the position around the origin (0, 0, 0)
                    position -= 0.5 * self.region
                    self.mol[i].x = position
                    i += 1

    """ Initiating cartesian velocities """
    def init_vels(self):
        vSum = 0
        for i in range(self.n_mol):
            self.mol[i].rv = np.random.random_sample(self.n_dim)
            self.mol[i].rv *= self.vel_mag
            vSum += self.mol[i].rv
        for i in range(self.n_mol):
            self.mol[i].rv -= (vSum / self.n_mol)

    """ Initiating quaternion positions """
    def init_ang_coords(self):
        for n in range(self.n_mol):
            e = np.random.random_sample(3)
            angle = np.array([
                np.arctan2(e[0], e[1]),
                np.arccos(e[2]),
                2 * np.pi * random.random()
            ])
            self.mol[n].q = utils.euler_to_quat(angle)

    """ Initiating angular velocities """
    def init_ang_vels(self):
        for n in range(self.n_mol):
            qe = np.random.random_sample(4)
            qe[3] = 0
            self.mol[n].qv = utils.Qproduct(self.mol[n].q, qe)
            f = 0.5 * self.vel_mag / np.sqrt(
                self.mol[n].m_inert[0] * qe[0] ** 2 +
                self.mol[n].m_inert[1] * qe[1] ** 2 +
                self.mol[n].m_inert[2] * qe[2] ** 2
            )

            self.mol[n].qv *= f

    """ computing the forces acting in the system """
    def calculate_forces(self):
        self.forces(self)

    """ updating the properties of total energy, kinetic energy and pressure based on the current state of the system """
    def eval_props(self):
        # this one is just to check if the velocities sum is actually zero all the time
        self.v_sum = np.zeros(self.n_dim)
        vvSum = 0

        for i in range(self.n_mol):
            self.v_sum += self.mol[i].rv
            vv = utils.get_magnitude(self.mol[i].rv)
            vvSum += vv

        self.kin_energy[0] = 0.5 * vvSum / self.n_mol
        self.tot_energy[0] = self.kin_energy[0] + self.u_sum / self.n_mol
        self.pressure[0] = self.density * (vvSum + self.vir_sum) / (self.n_mol * self.n_dim)

    """ accumulating the properties of the system, to also store the average and standard deviation of them """
    def accum_props(self, mode):
        if (mode == 0):
            utils.set_prop_zero(self.tot_energy)
            utils.set_prop_zero(self.kin_energy)
            utils.set_prop_zero(self.pressure)
        elif (mode == 1):
            utils.accum_prop(self.tot_energy)
            utils.accum_prop(self.kin_energy)
            utils.accum_prop(self.pressure)
        elif (mode == 2):
            utils.avg_prop(self.tot_energy, self.step_avg)
            utils.avg_prop(self.kin_energy, self.step_avg)
            utils.avg_prop(self.pressure, self.step_avg)

    """ applying the periodic boundary conditions """
    def apply_boundary_conditions(self):
        for n in range(self.n_mol):
            utils.wrap_around(self.mol[n].r, self.region, self.n_dim)

    """ function to normalize the quaternions """
    def adjust_quat(self):
        for n in range(self.n_mol):
            # normalization of the quaternion
            self.mol[n].q *= 1 / np.sqrt(utils.lenSquared(self.mol[n].q))

    """ function to print out a summary of the systems curent state """
    def print_summary(self):
        print('steps: ' + str(self.step_count) + '\n' +
              'vel_sum ' + str(np.sum(self.v_sum / self.n_mol)) + '\n' +
              'E_tot: ' + str(self.tot_energy) + '\n' +
              'E_kin: ' + str(self.kin_energy) + '\n' +
              'pressure ' + str(self.pressure) + '\n' +
              '-------------------------- \n \n')