import math

import numpy as np
import pytest

from kasimir import lattice, spectral


def test_frequencies_positive_and_ordered():
    w = lattice.dirichlet_frequencies(12)
    assert w.shape == (12,)
    assert np.all(w > 0)
    assert np.all(np.diff(w) > 0)
    # lowest mode → π/a in continuum, a = 13, ω_1 = 2 sin(π/26) ≈ π/13
    assert w[0] == pytest.approx(math.pi / 13.0, rel=0.02)


def test_bulk_density_matches_large_n():
    # Dirichlet chains have a surface term O(1), so E/a = 2/π + O(1/a).
    # A two-length difference cancels the surface piece.
    n1, n2 = 200, 400
    a1, a2 = n1 + 1, n2 + 1
    bulk = (lattice.vacuum_energy(n2) - lattice.vacuum_energy(n1)) / (a2 - a1)
    assert bulk == pytest.approx(lattice.bulk_energy_density(), rel=1e-4)
    leftover = lattice.vacuum_energy(n2) - lattice.bulk_energy_density() * a2
    assert abs(leftover) < 2.0  # O(1) surface, not O(n) bulk miscount


def test_gamma_extrapolates_to_pi_over_24():
    lim = lattice.extrapolate_gamma(list(range(60, 281, 10)))
    assert lim.gamma == pytest.approx(-math.pi / 24.0, rel=5e-3)


def test_snapshot_is_logarithmic_not_casimir():
    c0, c1 = lattice.snapshot_scaling(list(range(40, 201, 8)))
    # Must have a clear log a slope, and must not be a 1/a Casimir fit.
    assert c1 != pytest.approx(0.0, abs=0.05)
    ns = list(range(40, 201, 8))
    a = np.array([n + 1.0 for n in ns])
    y = np.array([lattice.snapshot_logdet(n) for n in ns])
    # Residual of log fit should beat residual of a 1/a fit.
    log_pred = c0 + c1 * np.log(a)
    x_inv = np.column_stack([np.ones_like(a), 1.0 / a])
    inv_coeff, *_ = np.linalg.lstsq(x_inv, y, rcond=None)
    inv_pred = x_inv @ inv_coeff
    assert float(np.mean((y - log_pred) ** 2)) < float(np.mean((y - inv_pred) ** 2))


def test_lattice_casimir_tracks_continuum_formula():
    n = 80
    a = n + 1
    lim = lattice.extrapolate_gamma(list(range(40, 161, 8)))
    finite = lattice.vacuum_energy(n) - lattice.bulk_energy_density() * a - lim.intercept
    assert finite == pytest.approx(spectral.energy_1d_dirichlet(a), rel=0.05)
