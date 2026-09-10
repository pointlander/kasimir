"""1D harmonic chain: vacuum energy vs snapshot complexity.

A unit-mass, unit-spring chain of n interior sites.

Dirichlet–Dirichlet (DD), length a = n+1:
    ω_j = 2 sin(j π / (2 (n+1))),  j = 1..n
    continuum: ω → n π / a,  κ finite part → −π/(24 a)

Dirichlet–Neumann (DN), length a = n+1/2:
    ω_j = 2 sin((j+1/2) π / (2n+1)),  j = 0..n-1
    continuum: ω → (n+1/2) π / a,  κ finite part → +π/(48 a)

History-complexity rate  κ = (1/2) Σ ω_j   (vacuum energy)
Snapshot complexity      σ = (1/2) Σ log(1/ω_j)   (Gaussian entropy, up to a constant)

Only κ contains the Casimir 1/a piece; its sign tracks the boundary conditions.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Literal

import numpy as np

from .spectral import log_two_sinh, oscillator_helmholtz

PI = math.pi
BC = Literal["dd", "dn"]


def interval_length(n: int, bc: BC = "dd") -> float:
    """Continuum length conjugate to the discrete spectrum."""
    if n < 1:
        raise ValueError("n must be a positive integer")
    if bc == "dd":
        return float(n + 1)
    if bc == "dn":
        return float(n) + 0.5
    raise ValueError(f"unknown boundary conditions {bc!r}")


def laplacian(n: int, bc: BC = "dd") -> np.ndarray:
    """Discrete Neumann/Dirichlet Laplacian on n interior sites, unit springs.

    Dirichlet walls contribute a spring to a site held at zero (diagonal 2 at
    that end). A Neumann end has no spring beyond the last site (diagonal 1).
    """
    if n < 1:
        raise ValueError("n must be a positive integer")
    if bc not in ("dd", "dn"):
        raise ValueError(f"unknown boundary conditions {bc!r}")
    lap = np.zeros((n, n))
    for i in range(n):
        lap[i, i] = 2.0
        if i > 0:
            lap[i, i - 1] = -1.0
        if i < n - 1:
            lap[i, i + 1] = -1.0
    if bc == "dn":
        lap[n - 1, n - 1] = 1.0
    return lap


def frequencies(n: int, bc: BC = "dd") -> np.ndarray:
    """Frequencies of an n-site chain, unit mass and spring."""
    if n < 1:
        raise ValueError("n must be a positive integer")
    if bc == "dd":
        j = np.arange(1, n + 1, dtype=float)
        return 2.0 * np.sin(j * PI / (2.0 * (n + 1)))
    if bc == "dn":
        j = np.arange(n, dtype=float)
        return 2.0 * np.sin((j + 0.5) * PI / (2.0 * n + 1.0))
    raise ValueError(f"unknown boundary conditions {bc!r}")


def dirichlet_frequencies(n: int) -> np.ndarray:
    """Frequencies of an n-site Dirichlet chain (alias for frequencies(n, 'dd'))."""
    return frequencies(n, bc="dd")


def vacuum_energy(n: int, bc: BC = "dd") -> float:
    """κ(n) = (ħ/2) Σ ω_j with ħ = 1. T=0 Helmholtz."""
    return 0.5 * float(frequencies(n, bc=bc).sum())


def helmholtz(n: int, T: float, bc: BC = "dd") -> float:
    """Lattice Helmholtz free energy T Σ log(2 sinh(ω_j / 2T)).

    Equals vacuum_energy at T=0. This is (1/β) × ½ log det(−Δ_E) for a
    static chain, i.e. the Matsubara-converted complexity rate of the
    constrained ensemble.
    """
    return oscillator_helmholtz(frequencies(n, bc=bc), T)


def snapshot_logdet(n: int, bc: BC = "dd") -> float:
    """(1/2) Σ log(1/ω_j). Tracks ζ'(0), not Casimir energy.

    Ground-state position variance of oscillator j is ∝ 1/ω_j, so this is
    the N-dependent part of the typical-sample Kolmogorov complexity.
    """
    omega = frequencies(n, bc=bc)
    return 0.5 * float(np.log(1.0 / omega).sum())


def bulk_energy_density() -> float:
    """Continuum-of-modes bulk energy per unit length at T=0.

    ω(k) = 2 |sin(k/2)|. Independent of the far-end boundary condition:
        (1/N) Σ_m (1/2) · 2 |sin(π m / N)| → (1/π) ∫_0^π sin u du = 2/π.
    """
    return 2.0 / PI


def bulk_free_energy_density(T: float, n_k: int = 8000) -> float:
    """Infinite-chain Helmholtz density at temperature T.

    Average of T log(2 sinh(ω(k)/2T)) over the Brillouin zone, skipping the
    measure-zero k=0 point. At T=0 this is 2/π.
    """
    if T < 0:
        raise ValueError("temperature T must be nonnegative")
    if T == 0:
        return bulk_energy_density()
    k = np.linspace(0.0, 2.0 * PI, n_k, endpoint=False)
    omega = 2.0 * np.abs(np.sin(k / 2.0))
    omega = omega[omega > 1e-12]
    return float(T * np.mean(log_two_sinh(omega / (2.0 * T))))


@dataclass(frozen=True)
class ContinuumLimit:
    gamma: float
    intercept: float
    slope_1_over_a2: float
    ns: tuple[int, ...]
    residuals: tuple[float, ...]
    bc: str


def extrapolate_gamma(ns: list[int] | None = None, bc: BC = "dd") -> ContinuumLimit:
    """Fit E(n) − (2/π) a = b0 + γ/a + b2 / a²  and return γ.

    DD target: γ → −π/24 ≈ −0.1309.
    DN target: γ → +π/48 ≈ +0.06545.
    """
    if ns is None:
        ns = list(range(40, 241, 10))
    a = np.array([interval_length(n, bc=bc) for n in ns], dtype=float)
    y = np.array(
        [vacuum_energy(n, bc=bc) - bulk_energy_density() * interval_length(n, bc=bc) for n in ns]
    )
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
        bc=bc,
    )


def snapshot_scaling(ns: list[int] | None = None, bc: BC = "dd") -> tuple[float, float]:
    """Fit snapshot_logdet(n) = c0 + c1 log a. Returns (c0, c1).

    Paper §3: this tracks log a, not 1/a.
    """
    if ns is None:
        ns = list(range(40, 241, 10))
    a = np.array([interval_length(n, bc=bc) for n in ns], dtype=float)
    y = np.array([snapshot_logdet(n, bc=bc) for n in ns])
    x = np.column_stack([np.ones_like(a), np.log(a)])
    coeff, *_ = np.linalg.lstsq(x, y, rcond=None)
    return float(coeff[0]), float(coeff[1])
