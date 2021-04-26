import numpy as np
from integrators.thermostats import apply_thermostat, apply_thermostat_nonlinear
from objects.Molecule import compute_torq, compute_accels_q


def printer(mol):
    # print('summary:')
    # for i in range(len(mol)):
    #     print(mol[i].r)
    #     print(mol[i].rv)
    print('\n')

class PredictorCorrectorIntegrator():

    def __init__(self, parent, delta_t):
        self.system = parent
        self.delta_t = delta_t

    def integrate(self):
        # predictor step of the integrator
        predictor_step(self.delta_t, self.system.n_mol, self.system.mol)

        # check for the periodic boundary conditions
        self.system.apply_boundary_conditions()

        # compute the forces in the system
        self.system.calculate_forces()

        # system apply the thermostat
        apply_thermostat(self.delta_t, self.system.n_mol, self.system.mol)

        # corrector step of the integrator
        corrector_step(self.delta_t, self.system.n_mol, self.system.mol)

        # check for the periodic boundary conditions
        self.system.apply_boundary_conditions()

    def integrate_nonlinear(self):
        # predictor step
        predictor_step(self.delta_t, self.system.n_mol, self.system.mol)
        # Q predictor step
        predictor_stepQ(self.delta_t, self.system.n_mol, self.system.mol)

        # calculate the site coordinates
        self.system.gen_site_coords()

        # calculate the forces in the system now
        self.system.calculate_forces()

        for i in range(self.system.n_mol):
            compute_torq(self.system.mol[i], self.system.sites)
            compute_accels_q(self.system.mol[i])

        # apply the thermostat correctly
        apply_thermostat_nonlinear(self.delta_t, self.system.n_mol, self.system.mol)

        # corrector step
        corrector_step(self.delta_t, self.system.n_mol, self.system.mol)
        # Q corrector step
        corrector_stepQ(self.delta_t, self.system.n_mol, self.system.mol)

        # make sure the periodic boundary conditions are applied as well as the quaternions nomralized
        self.system.adjust_quat()
        self.system.apply_boundary_conditions()

# wr = weight r = deltaT**2 / div, cr = coeffiecent r
# wv = weight v = deltaT**2 / div, cv = coeffiecent v

def pr4(mol, delta_t, wr, cr, t):
    mol.r[t] = mol.r[t] + delta_t * mol.rv[t] + wr * (cr[0] * mol.ra[t] + cr[1] * mol.ra1[t] + cr[2] * mol.ra2[t])

def pv4(mol, delta_t, wv, cv, t):
    mol.rv[t] = (mol.r[t] - mol.ro[t]) / delta_t + wv * (cv[0] * mol.ra[t] + cv[1] * mol.ra1[t] + cv[2] * mol.ra2[t])

def cr4(mol, delta_t, wr, cr, t):
    mol.r[t] = mol.ro[t] + delta_t * mol.rvo[t] + wr * (cr[0] * mol.ra[t] + cr[1] * mol.ra1[t] + cr[2] * mol.ra2[t])

def cv4(mol, delta_t, wv, cv, t):
    mol.rv[t] = (mol.r[t] - mol.ro[t]) / delta_t + wv * (cv[0] * mol.ra[t] + cv[1] * mol.ra1[t] + cv[2] * mol.ra2[t])

def prq4(mol, delta_t, wr, cr, t):
    mol.q[t] = mol.q[t] + delta_t * mol.qv[t] + wr * (cr[0] * mol.qa[t] + cr[1] * mol.qa1[t] + cr[2] * mol.qa2[t])

def pvq4(mol, delta_t, wv, cv, t):
    mol.qv[t] = (mol.q[t] - mol.qo[t]) / delta_t + wv * (cv[0] * mol.qa[t] + cv[1] * mol.qa1[t] + cv[2] * mol.qa2[t])

def crq4(mol, delta_t, wr, cr, t):
    mol.q[t] = mol.qo[t] + delta_t * mol.qvo[t] + wr * (cr[0] * mol.qa[t] + cr[1] * mol.qa1[t] + cr[2] * mol.qa2[t])

def cvq4(mol, delta_t, wv, cv, t):
    mol.qv[t] = (mol.q[t] - mol.qo[t]) / delta_t + wv * (cv[0] * mol.qa[t] + cv[1] * mol.qa1[t] + cv[2] * mol.qa2[t])

def predictor_step(delta_t, n_mol, mol):
    cr = np.array([19, -10, 3])
    cv = np.array([27, -22, 7])
    div = 24
    wr = delta_t**2 / div
    wv = delta_t / div

    for n in range(n_mol):
        mol[n].ro = mol[n].r.copy()
        mol[n].rvo = mol[n].rv.copy()
        for i in range(len(mol[0].r)):
            pr4(mol[n], delta_t, wr, cr, i)
            pv4(mol[n], delta_t, wv, cv, i)
        mol[n].ra2 = mol[n].ra1
        mol[n].ra1 = mol[n].ra

def corrector_step(delta_t, n_mol, mol):
    cr = np.array([3, 10, -1])
    cv = np.array([7, 6, -1])
    div = 24
    wr = delta_t**2 / div
    wv = delta_t / div

    for n in range(n_mol):
        for t in range(len(mol[0].r)):
            cr4(mol[n], delta_t, wr, cr, t)
            cv4(mol[n], delta_t, wv, cv, t)

def predictor_stepQ(delta_t, n_mol, mol):
    cr = np.array([19, -10, 3])
    cv = np.array([27, 22, 7])
    div = 24
    wr = delta_t**2 / div
    wv = delta_t / div

    for n in range(n_mol):
        mol[n].qo = mol[n].q.copy()
        mol[n].qvo = mol[n].qv.copy()
        for t in range(len(mol[0].r)):
            prq4(mol[n], delta_t, wr, cr, t)
            pvq4(mol[n], delta_t, wv, cv, t)
        mol[n].qa2 = mol[n].qa1
        mol[n].qa1 = mol[n].qa

def corrector_stepQ(delta_t, n_mol, mol):
    cr = np.array([3, 10, -1])
    cv = np.array([7, 6, -1])
    div = 24
    wr = delta_t**2 / div
    wv = delta_t / div

    for n in range(n_mol):
        for t in range(len(mol[0].r)):
            crq4(mol[n], delta_t, wr, cr, t)
            cvq4(mol[n], delta_t, wv, cv, t)