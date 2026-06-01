import numpy as np

from simplemd.particles.molecule import compute_ang_vel
from simplemd.utils import get_magnitude


def apply_thermostat(delta_t, n_mol, mol):
    """Gaussian isokinetic constraint: remove the component of acceleration that
    would change the translational kinetic energy (point-particle version)."""
    s1 = 0
    s2 = 0

    for i in range(n_mol):
        vt = mol[i].rv + 0.5 * delta_t * mol[i].ra
        s1 += vt.dot(mol[i].ra)
        s2 += get_magnitude(vt)

    vFac = -s1 / s2

    for i in range(n_mol):
        vt = mol[i].rv + 0.5 * delta_t * mol[i].ra
        mol[i].ra += vFac * vt


def apply_thermostat_nonlinear(delta_t, n_mol, mol):
    """Isokinetic constraint including rotational degrees of freedom (rigid bodies)."""
    s1 = 0
    s2 = 0

    for n in range(n_mol):
        s1 += mol[n].rv.dot(mol[n].ra)
        s2 += get_magnitude(mol[n].rv)

    for n in range(n_mol):
        w = compute_ang_vel(mol[n])
        s1 += w.dot(mol[n].torq)
        s2 += (
            mol[n].m_inert[0] * w[0] ** 2
            + mol[n].m_inert[1] * w[1] ** 2
            + mol[n].m_inert[2] * w[2] ** 2
        )

    vFac = -s1 / s2
    for n in range(n_mol):
        mol[n].ra += vFac * mol[n].rv
        mol[n].qa += vFac * mol[n].qv


def adjust_temp(system):
    """Rescale velocity magnitudes back to the target-temperature velocity."""
    # caculate the ratio between the initial and current velocities
    vvSum = 0
    for i in range(system.n_mol):
        vvSum += get_magnitude(system.mol[i].rv)
    vFac = system.vel_mag / np.sqrt(vvSum / system.n_mol)
    # rescale the velocities
    for i in range(system.n_mol):
        system.mol[i].rv *= vFac


def adjust_temp_nonlin(system):
    """Rescale translational and rotational velocities back to the target temperature."""
    # caculate the ratio between the initial and current velocities
    vvSum = 0
    for i in range(system.n_mol):
        vvSum += get_magnitude(system.mol[i].rv)
    vFac = system.vel_mag / np.sqrt(vvSum / system.n_mol)
    # rescale the velocities
    for i in range(system.n_mol):
        system.mol[i].rv *= vFac

    vvqSum = 0
    for n in range(system.n_mol):
        w = compute_ang_vel(system.mol[n])
        vvqSum += (
            system.mol[n].m_inert[0] * w[0] ** 2
            + system.mol[n].m_inert[1] * w[1] ** 2
            + system.mol[n].m_inert[2] * w[2] ** 2
        )

    vFac = system.vel_mag / np.sqrt(vvqSum / system.n_mol)
    for n in range(system.n_mol):
        system.mol[n].qv *= vFac
