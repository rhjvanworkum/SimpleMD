import numpy as np
from utils import get_magnitude

def wrap_around_cell(cell, cells, region, shift):
    for i in range(len(cells)):
        if (cell[i] >= cells[i]):
            cell[i] = 0
            shift[i] = region[i]
        elif (cell[i] < 0):
            cell[i] = cells[i] - 1
            shift[i] = - region[i]

def ComputeForcesCD(system):
    rr_cut = system.r_cut ** 2
    # number of cells with length r_cut that fit inside region
    cells = np.ceil(system.region / system.r_cut).astype(int)

    # initializing a dictionary that contains all particles index inside each cell
    cell_dict = {}
    for i in range(cells[0]):
        for j in range(cells[1]):
            for k in range(cells[2]):
                cell_dict[(i, j, k)] = np.array([]).astype(int)

    # amount of cells we need to check for each cell
    offsets = np.array([[0,0,0], [1,0,0], [1,1,0], [0,1,0], [-1,1,0], [0,0,1],
                        [1,0,1], [1,1,1], [0,1,1], [-1,1,1], [-1,0,1],
                        [-1,-1,1], [0,-1,1], [1,-1,1]])

    # inverse of the cell length
    inverse_length = 1 / system.r_cut

    # assigning the particles to the different cells
    for i in range(system.n_mol):
        position = system.mol[i].r + 0.5 * system.region
        cell_position = position * inverse_length
        index = (np.floor(cell_position[0]), np.floor(cell_position[1]), np.floor(cell_position[2]))
        cell_dict[index] = np.append(cell_dict[index], i)

    # setting the acceleration to zero
    for mol in system.mol:
        mol.ra = np.zeros(3)

    # iterating over the cells first
    for cx in range(cells[0]):
        for cy in range(cells[1]):
            for cz in range(cells[2]):
                cell1 = np.array([cx, cy, cz])

                for offset in offsets:
                    cell2 = cell1 + offset

                    shift = np.zeros(3)
                    wrap_around_cell(cell2, cells, system.region, shift)

                    # looping over all the particles in cell1 and cell2
                    for i in cell_dict[(cell1[0], cell1[1], cell1[2])]:
                        for j in cell_dict[(cell2[0], cell2[1], cell2[2])]:
                            # perform calculation if cells are different or
                            # if cells are not different only when j2 is highter than j1
                            # this performs pair-wise checks in the same cell and not double
                            if (not np.array_equal(cell1, cell2) or j > i):

                                # calculate the distance now
                                Rij = system.mol[i].r - system.mol[j].r
                                Rij -= shift

                                # calculating the magnitude of the distance vector
                                rij = Rij[0] ** 2 + Rij[1] ** 2 + Rij[2] ** 2

                                # we only calculate the force if the distance vector is smaller than the cut off distance
                                if (rij < rr_cut):
                                    # equation 1.1
                                    rij_3 = (1 / rij) ** 3
                                    F_magnitude = 48 * (rij_3 ** 3 - 0.5 * rij_3 ** 2)

                                    # equation 1.2
                                    system.mol[i].ra += F_magnitude * Rij
                                    system.mol[j].ra += -F_magnitude * Rij
