"""Cell-division (cell-list) force evaluation.

Bins particles into a grid of cells of side >= ``r_cut`` so that only
neighbouring cells need to be checked, giving O(N) scaling. Produces the same
Lennard-Jones forces as the pair-by-pair method, but faster for large systems.

The cell-list half-shell stencil requires at least 3 cells per dimension; for
smaller systems use the ``PBP`` method instead.
"""

import numpy as np


def wrap_around_cell(cell, cells, region, shift):
    """Wrap a cell index across the periodic boundary, recording the image shift."""
    for i in range(len(cells)):
        if cell[i] >= cells[i]:
            cell[i] = 0
            shift[i] = region[i]
        elif cell[i] < 0:
            cell[i] = cells[i] - 1
            shift[i] = -region[i]


def ComputeForcesCD(system):
    """Lennard-Jones forces via a cell list. Equivalent to ``ComputeForcesPBP``.

    Raises:
        ValueError: if the system is too small for a valid cell list
            (fewer than 3 cells in any dimension).
    """
    rr_cut = system.r_cut**2
    # number of cells of side >= r_cut that tile the region in each dimension
    cells = np.floor(system.region / system.r_cut).astype(int)
    if np.any(cells < 3):
        raise ValueError(
            "cell-division requires at least 3 cells per dimension "
            f"(got {cells.tolist()}); increase N/density or use the 'PBP' method"
        )

    # initializing a dictionary that contains all particles index inside each cell
    cell_dict = {}
    for i in range(cells[0]):
        for j in range(cells[1]):
            for k in range(cells[2]):
                cell_dict[(i, j, k)] = np.array([]).astype(int)

    # amount of cells we need to check for each cell
    offsets = np.array(
        [
            [0, 0, 0],
            [1, 0, 0],
            [1, 1, 0],
            [0, 1, 0],
            [-1, 1, 0],
            [0, 0, 1],
            [1, 0, 1],
            [1, 1, 1],
            [0, 1, 1],
            [-1, 1, 1],
            [-1, 0, 1],
            [-1, -1, 1],
            [0, -1, 1],
            [1, -1, 1],
        ]
    )

    # Inverse cell length per dimension. Using region/cells (not r_cut) makes the
    # cells tile the region exactly, so floor(position * inverse_length) is always
    # in range [0, cells) -- this is the fix for the previous out-of-range KeyError.
    inverse_length = cells / system.region

    # assigning the particles to the different cells
    for i in range(system.n_mol):
        position = system.mol[i].r + 0.5 * system.region
        cell_position = position * inverse_length
        index = (
            int(np.floor(cell_position[0])),
            int(np.floor(cell_position[1])),
            int(np.floor(cell_position[2])),
        )
        cell_dict[index] = np.append(cell_dict[index], i)

    # setting the acceleration to zero and resetting the energy/virial accumulators
    for mol in system.mol:
        mol.ra = np.zeros(3)
    system.u_sum = 0
    system.vir_sum = 0

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
                            if not np.array_equal(cell1, cell2) or j > i:
                                # calculate the distance now
                                Rij = system.mol[i].r - system.mol[j].r
                                Rij -= shift

                                # calculating the magnitude of the distance vector
                                rij = Rij[0] ** 2 + Rij[1] ** 2 + Rij[2] ** 2

                                # only compute the force inside the cut-off radius
                                if rij < rr_cut:
                                    # equation 1.1
                                    rij_3 = (1 / rij) ** 3
                                    F_magnitude = 48 * (rij_3**3 - 0.5 * rij_3**2)

                                    # equation 1.2
                                    system.mol[i].ra += F_magnitude * Rij
                                    system.mol[j].ra += -F_magnitude * Rij

                                    # accumulate potential energy and virial
                                    # (same definitions as the pair-by-pair method)
                                    system.u_sum += 4 * rij_3 * (rij_3 - 1) + 1
                                    system.vir_sum += F_magnitude * rij
