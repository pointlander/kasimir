"""Computable proxies for the Kolmogorov quantities in the paper.

True K is uncomputable. What we can compute:

- Gaussian typical-sample complexity (3):  (1/2) log2 det(2π e Σ / δ²)
  This is snapshot complexity, spectral moment ζ'(0).
- History-complexity rate (4):  (1/2) Σ ω_k
  This is κ, spectral moment ζ(-1), the Casimir energy.
- Dimensionless I_K = κ / Θ, the split of §6–§7.
"""

from __future__ import annotations

import math

import numpy as np

from .spectral import information_1d, information_3d, modular_temperature, theta_1d_mode_spacing


LN2 = math.log(2.0)


def gaussian_typical_bits(cov: np.ndarray, delta: float = 1.0) -> float:
    """Typical-set Kolmogorov complexity of a Gaussian sample, in bits.

    K ≈ (1/2) log2 det(2π e Σ / δ²). `cov` must be SPD.
    """
    cov = np.asarray(cov, dtype=float)
    sign, logdet = np.linalg.slogdet(cov)
    if sign <= 0:
        raise ValueError("covariance must be positive definite")
    n = cov.shape[0]
    # (1/2) Σ log2(2π e λ_i / δ²) = (1/2) (logdet / ln 2) + (n/2) log2(2π e / δ²)
    return 0.5 * (logdet / LN2) + 0.5 * n * math.log2(2.0 * math.pi * math.e / delta**2)


def history_rate_from_frequencies(omega: np.ndarray, hbar: float = 1.0) -> float:
    """κ = (ħ/2) Σ ω. History-complexity production rate."""
    return 0.5 * hbar * float(np.asarray(omega, dtype=float).sum())


def dimensionless_from_energy(energy: float, theta: float) -> float:
    """I_K = E / Θ. Raises if Θ is not positive."""
    if theta <= 0:
        raise ValueError("conversion scale Θ must be positive")
    return energy / theta


def split_1d(a: float, energy: float, hbar: float = 1.0, c: float = 1.0) -> tuple[float, float]:
    """Return (Θ_1D, I_K) for a 1D interval, using mode-spacing conversion."""
    theta = theta_1d_mode_spacing(a, hbar=hbar, c=c)
    return theta, dimensionless_from_energy(energy, theta)


def split_3d(
    a: float, energy: float, area: float = 1.0, hbar: float = 1.0, c: float = 1.0
) -> tuple[float, float]:
    """Return (Θ_mod, I_K) for a 3D slab, using modular conversion."""
    theta = modular_temperature(a, hbar=hbar, c=c)
    return theta, dimensionless_from_energy(energy, theta)


def analytic_split_1d(a: float, hbar: float = 1.0, c: float = 1.0) -> tuple[float, float]:
    return theta_1d_mode_spacing(a, hbar=hbar, c=c), information_1d(a)


def analytic_split_3d(
    a: float, area: float = 1.0, hbar: float = 1.0, c: float = 1.0
) -> tuple[float, float]:
    return modular_temperature(a, hbar=hbar, c=c), information_3d(a, area=area)
