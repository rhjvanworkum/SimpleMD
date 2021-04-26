"""
Molecule Class Object - used for the TIP-4P model

Includes:
- Site class
- Molecular Site class
- functions to compute angular velocities, accelerations and torques
"""

import numpy as np
from utils import Qproduct, Mproduct, build_rot_matrix, lenSquared


class Molecule():
    def __init__(self, idx):
        self.idx = idx
        self.in_chain = -1

        self.r = np.zeros(3)
        self.rv = np.zeros(3)
        self.ra = np.zeros(3)
        self.ra1 = np.zeros(3)
        self.ra2 = np.zeros(3)
        self.ro = np.zeros(3)
        self.rvo = np.zeros(3)

        self.q = np.zeros(4)
        self.qv = np.zeros(4)
        self.qa = np.zeros(4)
        self.qa1 = np.zeros(4)
        self.qa2 = np.zeros(4)
        self.qo = np.zeros(4)
        self.qvo = np.zeros(4)

        self.n_sites = 0
        self.m_inert = np.zeros(3)
        self.torq = np.zeros(3)


class Site():
    def __init__(self):
        self.f = np.zeros(3)
        self.r = np.zeros(3)


class MSite():
    def __init__(self):
        self.r = np.zeros(3)
        self.typeF = 0


def compute_ang_vel(mol):
    qvt = mol.qv
    qvt[3] *= -1
    # qt is the quaternion containing the angular velocities
    # formula: 2 * W * Quaternion velocity
    # or w = 2 q_line * q_dot
    qt = 2.0 * Qproduct(qvt, mol.q)
    # put the angular velocities in a 3D vector
    return np.array([qt[0], qt[1], qt[2]])


def compute_accels_q(mol):
    w = compute_ang_vel(mol)
    # w_dot_x = (torq_x + (I_y - I_z) * w_y * w_z) / I_X
    qs = np.array([
    (mol.torq[0] + (mol.m_inert[1] - mol.m_inert[2]) * w[1] * w[2]) / mol.m_inert[0],
    (mol.torq[1] + (mol.m_inert[2] - mol.m_inert[0]) * w[2] * w[0]) / mol.m_inert[1],
    (mol.torq[2] + (mol.m_inert[0] - mol.m_inert[1]) * w[0] * w[1]) / mol.m_inert[2],
    -2. * lenSquared(mol.qv)
    ])
    # computing the quaternion accelerations
    # qa = 0.5 * W.T * (w_dot_x, w_dot_y, w_dot_z, -2 sum (q_dot_m^2))
    mol.qa = Qproduct(mol.q, qs) * 0.5


def compute_torq(mol, sites):
    mol.ra = np.zeros(3)
    torq_sum = np.zeros(3)
    # for each site we calculate the accelaration and torque
    for i in range(mol.n_sites):
        mol.ra += sites[mol.idx * mol.n_sites + i].f
        dir = sites[mol.idx * mol.n_sites + i].r - mol.r
        # total torque = sum (r_k cross f_k)
        torq_sum += np.cross(dir, sites[mol.idx * mol.n_sites + i].f)
    r_mat = build_rot_matrix(mol.q, 0)
    # torque: space-fixed coordinates -> body-fixed coordinates
    mol.torq = Mproduct(r_mat, torq_sum)  # np.dot serves as m-v multiplication here