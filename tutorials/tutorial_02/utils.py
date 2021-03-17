def get_magnitude(vec):
    mag  = 0
    ndim = len(vec)
    for i in range(ndim):
        mag += vec[i]**2
    return mag

def wrap_around(vec, region):
    for i in range(3):
        if (vec[i] >= 0.5 * region[i]):
            vec[i] -= region[i]
        elif (vec[i] < -0.5 * region[i]):
            vec[i] += region[i]

def wrap_around_cell(cell, cells, region, shift):
    for i in range(len(cells)):
        if (cell[i] >= cells[i]):
            cell[i] = 0
            shift[i] = region[i]
        elif (cell[i] < 0):
            cell[i] = cells[i] - 1
            shift[i] = - region[i]