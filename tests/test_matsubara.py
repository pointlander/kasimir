import math

import numpy as np
import pytest

from kasimir import lattice, spectral


def test_log_two_sinh_identity():
    x = np.array([1e-6, 0.2, 1.0, 8.0, 40.0])
    np.testing.assert_allclose(
        spectral.log_two_sinh(x), np.log(2.0 * np.sinh(x)), rtol=1e-12, atol=1e-12
    )


def test_free_energy_zero_t_is_casimir():
    a = 0.8
    assert spectral.free_energy_1d_dirichlet(a, 0.0) == pytest.approx(
        spectral.energy_1d_dirichlet(a)
    )
    assert spectral.free_energy_1d_dirichlet_neumann(a, 0.0) == pytest.approx(
        spectral.energy_1d_dirichlet_neumann(a)
    )
    assert spectral.free_energy_1d_dirichlet(a, 1e-10) == pytest.approx(
        spectral.energy_1d_dirichlet(a), rel=1e-12, abs=1e-15
    )
    assert spectral.free_energy_1d_dirichlet_neumann(a, 1e-10) == pytest.approx(
        spectral.energy_1d_dirichlet_neumann(a), rel=1e-12, abs=1e-15
    )


def test_thermal_correction_is_negative_and_vanishes_at_low_t():
    a = 1.2
    f_dd = spectral.free_energy_1d_dirichlet(a, 0.15)
    e_dd = spectral.energy_1d_dirichlet(a)
    assert f_dd < e_dd  # thermal piece T log(1-e^{-βω}) < 0
    f_dn = spectral.free_energy_1d_dirichlet_neumann(a, 0.15)
    e_dn = spectral.energy_1d_dirichlet_neumann(a)
    assert f_dn < e_dn


def test_high_t_modular_dual_dd():
    # F + π a T²/6 − (T/2) log(2 a T) → 0 for a T ≫ 1.
    a, T = 1.0, 8.0
    f = spectral.free_energy_1d_dirichlet(a, T)
    dual = f + math.pi * a * T**2 / 6.0 - 0.5 * T * math.log(2.0 * a * T)
    assert dual == pytest.approx(0.0, abs=1e-10)


def test_oscillator_helmholtz_zero_t():
    omega = np.array([0.5, 1.5, 2.5])
    assert spectral.oscillator_helmholtz(omega, 0.0) == pytest.approx(0.5 * omega.sum())
    low = spectral.oscillator_helmholtz(omega, 1e-8)
    assert low == pytest.approx(0.5 * omega.sum(), rel=1e-8)


def test_lattice_helmholtz_matches_vacuum_at_zero_t():
    n = 40
    assert lattice.helmholtz(n, 0.0, bc="dd") == pytest.approx(lattice.vacuum_energy(n, "dd"))
    assert lattice.helmholtz(n, 0.0, bc="dn") == pytest.approx(lattice.vacuum_energy(n, "dn"))


def test_lattice_low_t_equals_vacuum_energy():
    n = 50
    for bc in ("dd", "dn"):
        assert lattice.helmholtz(n, 1e-8, bc=bc) == pytest.approx(
            lattice.vacuum_energy(n, bc=bc), rel=1e-12
        )


def test_lattice_thermal_tracks_continuum_modes():
    """IR thermal piece on the chain matches continuum DD/DN frequencies.

    T is well below the lattice cutoff, so only long-wavelength modes are
    thermally populated and the discrete dispersion is close to ω = k.
    """
    n = 80
    T = 0.08
    for bc in ("dd", "dn"):
        a = lattice.interval_length(n, bc=bc)
        thermal_lat = lattice.helmholtz(n, T, bc=bc) - lattice.vacuum_energy(n, bc=bc)
        if bc == "dd":
            omega_c = np.arange(1, n + 1, dtype=float) * math.pi / a
        else:
            omega_c = (np.arange(n, dtype=float) + 0.5) * math.pi / a
        thermal_c = spectral.oscillator_helmholtz(omega_c, T) - 0.5 * float(omega_c.sum())
        assert thermal_lat == pytest.approx(thermal_c, rel=0.05)
        assert thermal_lat < 0.0


def test_bulk_free_energy_zero_t():
    assert lattice.bulk_free_energy_density(0.0) == pytest.approx(
        lattice.bulk_energy_density()
    )


def test_rejects_negative_temperature():
    with pytest.raises(ValueError):
        spectral.free_energy_1d_dirichlet(1.0, -0.1)
    with pytest.raises(ValueError):
        lattice.helmholtz(8, -1.0)
