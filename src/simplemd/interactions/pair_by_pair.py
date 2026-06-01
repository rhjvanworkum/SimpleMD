"""Pair-by-pair force evaluation.

These functions compute interaction forces by iterating over every distinct
pair of particles (an O(N^2) scheme). They mutate the system in place, writing
accelerations onto the particles/sites and accumulating the potential energy
(``u_sum``) and virial (``vir_sum``) used for property evaluation.
"""

import numpy as np

from simplemd.utils import get_magnitude, shift_around, wrap_around


def ComputeForcesPBP(system):
    """Lennard-Jones (soft-sphere) forces between point particles.

    Used by the Lennard-Jones fluid model. Writes ``mol[i].ra`` and
    accumulates ``system.u_sum`` and ``system.vir_sum``.
    """
    rr_cut = system.r_cut**2

    # setting the acceleration to zero
    for i in range(system.n_mol):
        system.mol[i].ra = np.zeros(system.n_dim)

    system.u_sum = 0
    system.vir_sum = 0
    for i in range(system.n_mol - 1):
        for j in range(i + 1, system.n_mol):
            dir = system.mol[i].r - system.mol[j].r

            # changing the distance vector based on the position in the unit cell
            wrap_around(dir, system.region, system.n_dim)

            rr = get_magnitude(dir)
            if rr < rr_cut:
                rr_3 = (1 / rr) ** 3
                fcVal = 48 * rr_3 * (rr_3 - 0.5) * rr_3
                system.mol[i].ra += fcVal * dir
                system.mol[j].ra += -1 * fcVal * dir
                system.u_sum += 4 * rr_3 * (rr_3 - 1) + 1
                system.vir_sum += fcVal * rr


def ComputeForcesPBP_nonlinear(system):
    """Site-based forces for the TIP-4P water model.

    Iterates over molecule pairs, and for interacting pairs over their charged
    and Lennard-Jones sites. Writes per-site forces onto ``system.sites[*].f``
    and accumulates ``system.u_sum``.
    """
    rr_cut = system.r_cut**2

    # setting the forces to zero
    for i in range(system.n_mol * system.n_sites):
        system.sites[i].f = np.zeros(system.n_dim)

    system.u_sum = 0
    for i in range(system.n_mol - 1):
        for j in range(i + 1, system.n_mol):
            dir = system.mol[i].r - system.mol[j].r
            shift = np.zeros(system.n_dim)
            shift_around(dir, system.region, shift, system.n_dim)
            dir += shift

            rr = get_magnitude(dir)
            if rr < rr_cut:
                i_s = i * system.n_sites
                j_s = j * system.n_sites

                for k in range(system.n_sites):
                    for l in range(system.n_sites):
                        type_sum = system.molecule_sites[k].typeF + system.molecule_sites[l].typeF
                        if (
                            system.molecule_sites[k].typeF == system.molecule_sites[l].typeF
                            or type_sum == 5
                        ):
                            # distance vector between the two interacting sites.
                            dir = system.sites[i_s + k].r - system.sites[j_s + l].r
                            dir += shift

                            rr = get_magnitude(dir)
                            rri = 1 / rr

                            if type_sum == 2:
                                rr_3 = rri**3
                                u_val = 4 * rr_3 * (rr_3)
                                fc_val = 48 * rr_3 * (rr_3 - 0.5) * rri
                            elif type_sum == 4:
                                u_val = 4 * system.b * np.sqrt(rri)
                                fc_val = u_val * rri
                            elif type_sum == 5:
                                u_val = -2 * system.b * np.sqrt(rri)
                                fc_val = u_val * rri
                            elif type_sum == 6:
                                u_val = system.b * np.sqrt(rri)
                                fc_val = u_val * rri

                            system.sites[i_s + k].f += fc_val * dir
                            # NOTE: original code adds the same sign to both sites.
                            # This looks like a Newton's-third-law sign bug (the partner
                            # site should get the opposite force), but it is left as-is to
                            # preserve behavior — flagged for separate review.
                            system.sites[j_s + l].f += fc_val * dir
                            system.u_sum += u_val
