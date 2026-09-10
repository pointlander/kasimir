import math

import pytest

from kasimir import tonks


def test_dd_matches_sum_of_squares():
    N, a = 7, 3.5
    e = tonks.energy(N, a, bc="dd")
    s = sum(n**2 for n in range(1, N + 1))
    assert e == pytest.approx((math.pi**2 / (2 * a**2)) * s)


def test_dn_matches_half_integer_sum():
    N, a = 7, 3.5
    e = tonks.energy(N, a, bc="dn")
    s = sum((n + 0.5) ** 2 for n in range(N))
    assert e == pytest.approx((math.pi**2 / (2 * a**2)) * s)


def test_expansion_reproduces_energy_at_fixed_density():
    rho = 0.8
    for bc in ("dd", "dn"):
        e_bulk, e_surf, A = tonks.expansion(rho, bc=bc)
        for N in (5, 11, 40):
            a = N / rho
            E = tonks.energy(N, a, bc=bc)
            assert E == pytest.approx(e_bulk * a + e_surf + A / a, rel=1e-12)


def test_I_K_is_one_twelfth_and_minus_one_twenty_fourth():
    rho = 1.3
    for bc, target in (("dd", 1.0 / 12.0), ("dn", -1.0 / 24.0)):
        _, _, A = tonks.expansion(rho, bc=bc)
        assert tonks.extracted_I_K(A, rho) == pytest.approx(target)
        assert tonks.I_K(bc) == pytest.approx(target)


def test_I_K_independent_of_density_mass_and_hbar():
    for rho, mass, hbar in ((0.4, 1.0, 1.0), (2.5, 3.0, 0.7)):
        for bc, target in (("dd", 1.0 / 12.0), ("dn", -1.0 / 24.0)):
            _, _, A = tonks.expansion(rho, bc=bc, hbar=hbar, mass=mass)
            assert tonks.extracted_I_K(A, rho, hbar=hbar, mass=mass) == pytest.approx(target)


def test_empirical_protocol_recovers_exact_I_K():
    rho = 0.9
    Ns = list(range(8, 41, 4))
    assert tonks.empirical_I_K(Ns, rho, bc="dd") == pytest.approx(1.0 / 12.0, rel=1e-10)
    assert tonks.empirical_I_K(Ns, rho, bc="dn") == pytest.approx(-1.0 / 24.0, rel=1e-10)


def test_not_the_relativistic_cft_values():
    # The point of the check: TG ≠ relativistic scalar.
    assert tonks.I_K("dd") != pytest.approx(-1.0 / 24.0)
    assert tonks.I_K("dn") != pytest.approx(1.0 / 48.0)


def test_rejects_bad_inputs():
    with pytest.raises(ValueError):
        tonks.energy(0, 1.0)
    with pytest.raises(ValueError):
        tonks.energy(4, -1.0)
    with pytest.raises(ValueError):
        tonks.expansion(0.0)
    with pytest.raises(ValueError):
        tonks.empirical_I_K([2, 3], 1.0)
    with pytest.raises(ValueError):
        tonks.energy(4, 1.0, bc="nn")  # type: ignore[arg-type]
