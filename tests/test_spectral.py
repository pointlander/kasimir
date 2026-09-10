import math

import numpy as np
import pytest

from kasimir import spectral


def test_1d_energy_and_force_consistent():
    a = 0.37
    e = spectral.energy_1d_dirichlet(a)
    assert e == pytest.approx(-math.pi / (24.0 * a))
    da = 1e-6
    f_num = -(spectral.energy_1d_dirichlet(a + da) - spectral.energy_1d_dirichlet(a - da)) / (
        2 * da
    )
    assert spectral.force_1d_dirichlet(a) == pytest.approx(f_num, rel=1e-8)


def test_3d_prefactors():
    a, area = 0.5, 2.0
    e = spectral.energy_3d_em(a, area=area)
    assert e == pytest.approx(-(math.pi**2) * area / (720.0 * a**3))
    f = spectral.force_3d_em(a, area=area)
    assert f == pytest.approx(-(math.pi**2) * area / (240.0 * a**4))
    # Euler: F = -dE/da = 3 E / a  (both negative; 3|E|/a = |F|)
    assert f == pytest.approx(3.0 * e / a, rel=1e-12)


def test_1d_dn_energy_and_force_repulsive():
    a = 0.37
    e = spectral.energy_1d_dirichlet_neumann(a)
    assert e == pytest.approx(math.pi / (48.0 * a))
    assert e > 0
    da = 1e-6
    f_num = -(
        spectral.energy_1d_dirichlet_neumann(a + da)
        - spectral.energy_1d_dirichlet_neumann(a - da)
    ) / (2 * da)
    assert spectral.force_1d_dirichlet_neumann(a) == pytest.approx(f_num, rel=1e-8)
    assert spectral.force_1d_dirichlet_neumann(a) > 0


def test_1d_reconstruction():
    for a in (0.2, 1.0, 7.5):
        assert spectral.reconstruct_energy(a, dim=1) == pytest.approx(
            spectral.energy_1d_dirichlet(a), rel=1e-12
        )
        assert spectral.information_1d() == pytest.approx(-1.0 / 24.0)
        assert spectral.reconstruct_energy(a, dim=1, bc="dn") == pytest.approx(
            spectral.energy_1d_dirichlet_neumann(a), rel=1e-12
        )
        assert spectral.information_1d(bc="dn") == pytest.approx(1.0 / 48.0)
        f_rec = spectral.reconstruct_force(a, dim=1, bc="dn")
        assert f_rec == pytest.approx(spectral.force_1d_dirichlet_neumann(a), rel=1e-5)


def test_3d_reconstruction_product_rule():
    area = 1.7
    for a in (0.15, 0.8, 3.0):
        rec = spectral.reconstruct_energy(a, dim=3, area=area)
        assert rec == pytest.approx(spectral.energy_3d_em(a, area=area), rel=1e-12)
        f_rec = spectral.reconstruct_force(a, dim=3, area=area)
        assert f_rec == pytest.approx(spectral.force_3d_em(a, area=area), rel=1e-5)


def test_gamma_coefficient():
    # I_K = E / Θ = [π²/720 a³] / [1/(2π a)] = π³/360 a⁻²   (magnitudes)
    a, area = 1.3, 2.0
    theta = spectral.modular_temperature(a)
    i = spectral.information_3d(a, area=area)
    assert i == pytest.approx(spectral.GAMMA_3D * area / a**2)
    assert theta * i == pytest.approx(spectral.energy_3d_em(a, area=area), rel=1e-12)


def test_extracted_I_K_is_independent_of_a():
    for a in (0.3, 1.1, 8.0):
        e_dd = spectral.energy_1d_dirichlet(a)
        e_dn = spectral.energy_1d_dirichlet_neumann(a)
        assert spectral.extracted_I_K(e_dd, a) == pytest.approx(-1.0 / 24.0)
        assert spectral.extracted_I_K(e_dn, a) == pytest.approx(1.0 / 48.0)


def test_cft_energy_linear_in_central_charge():
    a = 1.7
    for c_cft in (1.0, 2.0, 4.5):
        e = spectral.cft_energy_1d(a, c_cft, bc="dd")
        assert spectral.extracted_I_K(e, a) == pytest.approx(-c_cft / 24.0)
        e_dn = spectral.cft_energy_1d(a, c_cft, bc="dn")
        assert spectral.extracted_I_K(e_dn, a) == pytest.approx(c_cft / 48.0)


def test_rejects_nonpositive():
    with pytest.raises(ValueError):
        spectral.energy_1d_dirichlet(0.0)
    with pytest.raises(ValueError):
        spectral.energy_3d_em(-1.0)
    with pytest.raises(ValueError):
        spectral.reconstruct_energy(1.0, dim=2)  # type: ignore[arg-type]
    with pytest.raises(ValueError):
        spectral.energy_1d_dirichlet_neumann(0.0)
    with pytest.raises(ValueError):
        spectral.information_1d(bc="nn")  # type: ignore[arg-type]
    with pytest.raises(ValueError):
        spectral.cft_energy_1d(1.0, -1.0)
