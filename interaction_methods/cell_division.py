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
    cells = np.floor(system.region / system.r_cut).astype(int)

    # initializing a dictionary that contains all particles index inside each cell
    cell_dict = {}
    for i in range(cells[0]):
        for j in range(cells[1]):
            for k in range(cells[2]):
                cell_dict[(i, j, k)] = np.array([]).astype(int)

    # amount of cells we need to check for each cell
    n_offsets = 14
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
# def ComputeForcesCD(system):
#     # array containing the amount of cells in x and y direction
#     cells = np.floor(system.region / system.r_cut).astype(int)
#     # updating rCut based on how the cells are divided - change in future
#     system.r_cut = system.region[0] / cells[0]
#
#     cellList = {}
#     for i in range(cells[0]):
#         for j in range(cells[1]):
#             if (system.n_dim == 3):
#                 for k in range(cells[2]):
#                     cellList[(i, j, k)] = np.array([]).astype(int)
#             else:
#                 cellList[(i, j)] = np.array([]).astype(int)
#
#     # variable for the amount of neigbouring cells to check, 5 in dim=2 and 14 in dim=3
#     if (system.n_dim == 2):
#         n_offset = 5
#         offset = np.array([[0, 0], [0, 1], [1, 1], [1, 0], [1, -1]])
#     elif (system.n_dim == 3):
#         n_offset = 14
#         offset = np.array([[0,0,0], [1,0,0], [1,1,0], [0,1,0], [-1,1,0], [0,0,1],
#                         [1,0,1], [1,1,1], [0,1,1], [-1,1,1], [-1,0,1],
#                         [-1,-1,1], [0,-1,1], [1,-1,1]])
#
#     # variable for the scaling between the region size and the cell size
#     invWidth = cells / system.region
#
#     # assing cellList items
#     for i in range(system.n_mol):
#         rs = system.mol[i].r + 0.5 * system.region
#         cc = rs * invWidth
#         index = getCellIndex(cc, system.n_dim)
#         # print(cc)
#         # print(index)
#         cellList[index] = np.append(cellList[index], i)
#
#     rrCut = system.r_cut**2
#
#     # setting the acceleration to zero
#     for i in range(system.n_mol):
#         system.mol[i].ra = np.zeros(system.n_dim)
#
#     system.u_sum = 0
#     system.vir_sum = 0
#
#     if (system.n_dim == 2):
#         for cx in range(cells[0]):
#             for cy in range(cells[1]):
#                 c1v = np.array([cx, cy])
#                 for off in range(n_offset):
#                     c2v = c1v + offset[off]
#                     shift = np.zeros(system.n_dim)
#                     CellWrap(c2v, cells, system.region, shift)
#                     for j1 in cellList[(c1v[0], c1v[1])]:
#                         for j2 in cellList[(c2v[0], c2v[1])]:
#                             # perform calculation if cells are different or
#                             # if cells are not different only when j2 is highter than j1
#                             # this performs pair-wise checks in the same cell and not double
#                             if (not np.array_equal(c1v, c2v) or j2 > j1):
#
#                                 # calculate the distance now
#                                 dir = system.mol[j1].r - system.mol[j2].r
#                                 dir -= shift
#
#                                 rr = get_magnitude(dir)
#                                 if (rr < rrCut):
#                                     rr_3 = (1 / rr) ** 3
#                                     fcVal = 48 * rr_3 * (rr_3 - 0.5) * rr_3
#                                     system.mol[j1].ra += fcVal * dir
#                                     system.mol[j2].ra += -1 * fcVal * dir
#                                     system.u_sum += 4 * rr_3 * (rr_3 - 1) + 1
#                                     system.vir_sum += fcVal * rr
#     elif (system.n_dim == 3):
#         for cx in range(cells[0]):
#             for cy in range(cells[1]):
#                 for cz in range(cells[2]):
#                     c1v = np.array([cx, cy, cz])
#                     for off in range(n_offset):
#                         c2v = c1v + offset[off]
#                         shift = np.zeros(system.n_dim)
#                         CellWrap(c2v, cells, system.region, shift)
#                         for j1 in cellList[(c1v[0], c1v[1], c1v[2])]:
#                             for j2 in cellList[(c2v[0], c2v[1], c2v[2])]:
#                                 # perform calculation if cells are different or
#                                 # if cells are not different only when j2 is highter than j1
#                                 # this performs pair-wise checks in the same cell and not double
#                                 if (not np.array_equal(c1v, c2v) or j2 > j1):
#
#                                     # calculate the distance now
#                                     dir = system.mol[j1].r - system.mol[j2].r
#                                     dir -= shift
#
#                                     rr = get_magnitude(dir)
#                                     if (rr < rrCut):
#                                         rr_3 = (1 / rr) ** 3
#                                         fcVal = 48 * rr_3 * (rr_3 - 0.5) * rr_3
#                                         system.mol[j1].ra += fcVal * dir
#                                         system.mol[j2].ra += -1 * fcVal * dir
#                                         system.u_sum += 4 * rr_3 * (rr_3 - 1) + 1
#                                         system.vir_sum += fcVal * rr