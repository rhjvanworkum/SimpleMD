import numpy as np
from utils import wrap_around, get_magnitude, shift_around

def ComputeForcesPBP_RK(system, positions):

    rr_cut = system.r_cut**2

    # setting the acceleration to zero
    ra = np.zeros((system.n_mol, system.n_dim))

    for i in range(system.n_mol - 1):
        for j in range(i + 1, system.n_mol):

            dir = positions[i] - positions[j]

            # changing the distance vector based on the position in the unit cell
            wrap_around(dir, system.region, system.n_dim)

            rr = get_magnitude(dir)
            if (rr < rr_cut):
                rr_3 = (1 / rr)**3
                # print(rr_3)
                fcVal = 48 * rr_3 * (rr_3 - 0.5) * rr_3
                ra[i] += fcVal * dir
                ra[j] += -1 * fcVal * dir

    return ra

def ComputeForcesPBP_nonlinear(system):

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
            if (rr < rr_cut):
                i_s = i * system.n_sites
                j_s = j * system.n_sites

                for k in range(system.n_sites):
                    for l in range(system.n_sites):
                        type_sum = system.molecule_sites[k].typeF + system.molecule_sites[l].typeF
                        if (system.molecule_sites[k].typeF == system.molecule_sites[l].typeF or type_sum == 5):
                            dir = system.sites[i_s + k].r = system.sites[j_s + l].r
                            dir += shift

                            rr = get_magnitude(dir)
                            rri = 1 / rr

                            if (type_sum == 2):
                                rr_3 = rri ** 3
                                u_val = 4 * rr_3 * (rr_3)
                                fc_val = 48 * rr_3 * (rr_3 - 0.5) * rri
                            elif (type_sum == 4):
                                u_val = 4 * system.b * np.sqrt(rri)
                                fc_val = u_val * rri
                            elif (type_sum == 5):
                                u_val = -2 * system.b * np.sqrt(rri)
                                fc_val = u_val * rri
                            elif (type_sum == 6):
                                u_val = system.b * np.sqrt(rri)
                                fc_val = u_val * rri

                            system.sites[i_s + k].f += fc_val * dir
                            system.sites[j_s + l].f += fc_val * dir
                            system.u_sum += u_val


def ComputeForcesPBP(system):

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
            if (rr < rr_cut):
                rr_3 = (1 / rr)**3
                fcVal = 48 * rr_3 * (rr_3 - 0.5) * rr_3
                system.mol[i].ra += fcVal * dir
                system.mol[j].ra += -1 * fcVal * dir
                system.u_sum += 4 * rr_3 * (rr_3 - 1) + 1
                system.vir_sum += fcVal * rr

def ComputeChainBondForces(system):

    rr_cut = system.r_cut**2

    system.u_sum = 0
    system.vir_sum = 0
    for n in range(system.n_chain):
        for i in range(system.chain_length - 1):
            j1 = n * system.chain_length + i;
            j2 = j1 + 1

            dir = system.mol[j1].r - system.mol[j2].r

            wrap_around(dir, system.region, system.n_dim)

            # here we check the soft sphere repulsive force
            rr = get_magnitude(dir)
            if (rr < rr_cut):
                rr_3 = (1 / rr) ** 3
                fcVal = 48 * rr_3 * (rr_3 - 0.5) * rr_3
                system.mol[j1].ra += fcVal * dir
                system.mol[j2].ra += -1 * fcVal * dir
                system.u_sum += 4 * rr_3 * (rr_3 - 1) + 1
                system.vir_sum += fcVal * rr

            # here we check the extra attractive interaction between
            # each pair of adjacent bonded atoms
            w = 1.0 - (system.bond_limit / np.sqrt(rr))
            if (w > 0):
                print('Bond snapped')
                break

            rr *= w**2
            if (rr < rr_cut):
                rr_3 = (1 / rr) ** 3
                fc_val = 48 * w * rr_3 * (rr_3 - 0.5) * rr_3
                system.mol[j1].ra += fc_val * dir
                system.mol[j2].ra += -1 * fc_val * dir
                system.u_sum += 4 * rr_3 * (rr_3 - 1) + 1
                system.vir_sum += fc_val * rr
