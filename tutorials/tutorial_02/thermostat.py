from utils import get_magnitude
import numpy as np

def apply_thermostat(delta_t, n_mol, mol):
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

def adjust_temp(particles, vel_mag):
    n_particles = len(particles)

    # caculate the ratio between the initial and current velocities
    vvSum = 0
    for particle in particles:
        vvSum += get_magnitude(particle.v)
    vFac = vel_mag / np.sqrt(vvSum / n_particles)
    # rescale the velocities
    for particle in particles:
        particle.v *= vFac