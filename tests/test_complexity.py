import math

import numpy as np
import pytest

from kasimir import complexity, spectral


def test_gaussian_bits_of_identity():
    n = 5
    bits = complexity.gaussian_typical_bits(np.eye(n), delta=1.0)
    # (n/2) log2(2π e)
    assert bits == pytest.approx(0.5 * n * math.log2(2 * math.pi * math.e))


def test_history_rate_matches_vacuum_energy():
    omega = np.array([1.0, 2.0, 3.0])
    assert complexity.history_rate_from_frequencies(omega) == pytest.approx(3.0)


def test_1d_split_recovers_one_over_24():
    a = 2.5
    e = spectral.energy_1d_dirichlet(a)
    theta, i = complexity.split_1d(a, e)
    assert i == pytest.approx(-1.0 / 24.0)
    assert theta * i == pytest.approx(e)
    e_dn = spectral.energy_1d_dirichlet_neumann(a)
    theta_dn, i_dn = complexity.split_1d(a, e_dn)
    assert i_dn == pytest.approx(1.0 / 48.0)
    assert theta_dn * i_dn == pytest.approx(e_dn)
    assert complexity.analytic_split_1d(a, bc="dn")[1] == pytest.approx(1.0 / 48.0)


def test_3d_split_recovers_gamma():
    a, area = 0.4, 3.0
    e = spectral.energy_3d_em(a, area=area)
    theta, i = complexity.split_3d(a, e, area=area)
    assert i == pytest.approx(spectral.information_3d(a, area=area))
    assert theta * i == pytest.approx(e)


def test_non_spd_covariance_rejected():
    with pytest.raises(ValueError):
        complexity.gaussian_typical_bits(np.array([[1.0, 0.0], [0.0, -1.0]]))
