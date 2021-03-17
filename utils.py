import numpy as np

def lenSquared(quat):
    return (quat[0])**2 + (quat[1])**2 + (quat[2])**2 + (quat[3])**2

def get_magnitude(vec):
    mag  = 0
    ndim = len(vec)
    for i in range(ndim):
        mag += vec[i]**2
    return mag

def get_dot_product(vec1, vec2):
    mag = 0
    if (len(vec1) != len(vec2)):
        raise ValueError("vectors you are trying to dot product do not have the same dim")
    else:
        ndim = len(vec1)
    for i in range(ndim):
        mag += vec1[i] * vec2[i]
    return mag

def Qproduct(q1, q2):
    return np.array([
        q1[3] * q2[0] - q1[2] * q2[1] + q1[1] * q2[2] + q1[0] * q2[3],
        q1[2] * q2[0] + q1[3] * q2[1] - q1[0] * q2[2] + q1[1] * q2[3],
       -q1[1] * q2[0] + q1[0] * q2[3] + q1[3] * q2[2] + q1[2] * q2[3],
       -q1[0] * q2[0] - q1[1] * q2[1] - q1[2] * q2[2] + q1[3] * q2[3]
    ])

def Mproduct(m, v):
    return np.array([
        m[0] * v[0] + m[3] * v[1] + m[6] * v[2],
        m[1] * v[0] + m[4] * v[1] + m[7] * v[2],
        m[2] * v[0] + m[5] * v[1] + m[8] * v[2],
    ])

def build_rot_matrix(q, transpose):
    rMat = np.zeros(9)
    p = np.zeros(10)
    k = 0
    for k2 in range(4):
        for k1 in range(k2, 4):
            p[k] = 2 * q[k1] * q[k2]
            k += 1

    rMat[0] = p[0] + p[9] - 1
    rMat[4] = p[4] + p[9] - 1
    rMat[8] = p[7] + p[9] - 1

    if (transpose):
        s = 1
    else:
        s = -1

    rMat[1] = p[1] + s * p[8]
    rMat[3] = p[1] - s * p[8]
    rMat[2] = p[2] - s * p[6]
    rMat[6] = p[2] + s * p[6]
    rMat[5] = p[5] + s * p[3]
    rMat[7] = p[5] - s * p[3]

    return rMat

def euler_to_quat(angle):
    a1 = 0.5 * angle[1]
    a2 = 0.5 * (angle[0] - angle[2])
    a3 = 0.5 * (angle[0] + angle[2])
    return np.array([
        np.sin(a1) * np.cos(a2),
        np.sin(a1) * np.sin(a2),
        np.cos(a1) * np.sin(a3),
        np.cos(a1) * np.cos(a3)
    ])

def set_prop_zero(v):
    v[1] = 0
    v[2] = 0

def accum_prop(v):
    v[1] += v[0]
    v[2] += v[0]**2

def avg_prop(v, n):
    v[1] /= n
    v[2] = np.sqrt(max(v[2] / n - v[1]**2, 0))

def wrap_around(vec, region, nDim):
    for i in range(nDim):
        if (vec[i] >= 0.5 * region[i]):
            vec[i] -= region[i]
        elif (vec[i] < -0.5 * region[i]):
            vec[i] += region[i]

def shift_around(vec, region, shift, nDim):
    for i in range(nDim):
        if (vec[i] >= 0.5 * region[i]):
            shift[i] -= region[i]
        elif (vec[i] < -0.5 * region[i]):
            shift[i] += region[i]