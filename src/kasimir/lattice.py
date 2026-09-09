"""1D harmonic chain: vacuum energy vs snapshot complexity.

A unit-mass, unit-spring Dirichlet chain of n interior sites has
    ω_j = 2 sin(j π / (2 (n+1))),  j = 1..n
and continuum limit c = 1, length a = n+1.

History-complexity rate  κ = (1/2) Σ ω_j   (vacuum energy)
Snapshot complexity      σ = (1/2) Σ log(1/ω_j)   (Gaussian entropy, up to a constant)

Only κ contains the −π/24 a^{-1} Casimir piece.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

PI = math.pi


def dirichlet_frequencies(n: int) -> np.ndarray:
    """Frequencies of an n-site Dirichlet chain, unit mass and spring."""
    if n < 1:
        raise ValueError("n must be a positive integer")
    j = np.arange(1, n + 1, dtype=float)
    return 2.0 * np.sin(j * PI / (2.0 * (n + 1)))


def vacuum_energy(n: int) -> float:
    """κ(n) = (ħ/2) Σ ω_j with ħ = 1."""
    return 0.5 * float(dirichlet_frequencies(n).sum())


def snapshot_logdet(n: int) -> float:
    """(1/2) Σ log(1/ω_j). Tracks ζ'(0), not Casimir energy.

    Ground-state position variance of oscillator j is ∝ 1/ω_j, so this is
    the N-dependent part of the typical-sample Kolmogorov complexity.
    """
    omega = dirichlet_frequencies(n)
    return 0.5 * float(np.log(1.0 / omega).sum())


def bulk_energy_density() -> float:
    """Continuum-of-modes bulk energy per unit length.

    ω(k) = 2 |sin(k/2)|. For a Dirichlet chain of n sites, length a = n+1,
    E(n)/a → (2/π). Same density as a periodic chain:
        (1/N) Σ_m (1/2) · 2 |sin(π m / N)| → (1/π) ∫_0^π sin u du = 2/π.
    """
    return 2.0 / PI


@dataclass(frozen=True)
class ContinuumLimit:
    gamma: float
    intercept: float
    slope_1_over_a2: float
    ns: tuple[int, ...]
    residuals: tuple[float, ...]


def extrapolate_gamma(ns: list[int] | None = None) -> ContinuumLimit:
    """Fit E(n) − (2/π) a = b0 + γ/a + b2 / a²  and return γ.

    Target: γ → −π/24 ≈ −0.1309.
    """
    if ns is None:
        ns = list(range(40, 241, 10))
    a = np.array([n + 1 for n in ns], dtype=float)
    y = np.array([vacuum_energy(n) - bulk_energy_density() * (n + 1) for n in ns])
    # y = b0 + γ/a + b2/a²
    x = np.column_stack([np.ones_like(a), 1.0 / a, 1.0 / a**2])
    coeff, *_ = np.linalg.lstsq(x, y, rcond=None)
    b0, gamma, b2 = (float(c) for c in coeff)
    pred = x @ coeff
    resid = tuple(float(r) for r in (y - pred))
    return ContinuumLimit(
        gamma=gamma,
        intercept=b0,
        slope_1_over_a2=b2,
        ns=tuple(ns),
        residuals=resid,
    )


def snapshot_scaling(ns: list[int] | None = None) -> tuple[float, float]:
    """Fit snapshot_logdet(n) = c0 + c1 log a. Returns (c0, c1).

    Paper §3: this tracks log a, not 1/a.
    """
    if ns is None:
        ns = list(range(40, 241, 10))
    a = np.array([n + 1 for n in ns], dtype=float)
    y = np.array([snapshot_logdet(n) for n in ns])
    x = np.column_stack([np.ones_like(a), np.log(a)])
    coeff, *_ = np.linalg.lstsq(x, y, rcond=None)
    return float(coeff[0]), float(coeff[1])
