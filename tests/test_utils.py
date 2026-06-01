"""Unit tests for the math helpers in :mod:`simplemd.utils`."""

import numpy as np
import pytest

from simplemd import utils


def test_get_magnitude_matches_squared_norm():
    v = np.array([1.0, 2.0, 2.0])
    assert utils.get_magnitude(v) == pytest.approx(9.0)  # 1 + 4 + 4


def test_get_dot_product_matches_numpy():
    a = np.array([1.0, 2.0, 3.0])
    b = np.array([4.0, -5.0, 6.0])
    assert utils.get_dot_product(a, b) == pytest.approx(float(np.dot(a, b)))


def test_get_dot_product_raises_on_dim_mismatch():
    with pytest.raises(ValueError):
        utils.get_dot_product(np.array([1.0, 2.0]), np.array([1.0, 2.0, 3.0]))


def test_euler_to_quat_is_unit_quaternion():
    q = utils.euler_to_quat(np.array([0.3, 0.5, 0.7]))
    assert utils.lenSquared(q) == pytest.approx(1.0)


def test_wrap_around_brings_coordinate_into_box():
    region = np.array([10.0, 10.0, 10.0])
    vec = np.array([6.0, -6.0, 0.0])  # 6 >= 5 -> -=10 ; -6 < -5 -> +=10
    utils.wrap_around(vec, region, 3)
    np.testing.assert_allclose(vec, [-4.0, 4.0, 0.0])


def test_accum_and_avg_prop_compute_mean_and_std():
    # prop = [current, sum, sum_of_squares]; feed two samples of 2.0 and 4.0
    prop = np.zeros(3)
    utils.set_prop_zero(prop)
    prop[0] = 2.0
    utils.accum_prop(prop)
    prop[0] = 4.0
    utils.accum_prop(prop)
    utils.avg_prop(prop, 2)
    assert prop[1] == pytest.approx(3.0)  # mean of 2 and 4
    assert prop[2] == pytest.approx(1.0)  # std of 2 and 4
