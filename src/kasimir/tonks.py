"""Tonks–Girardeau / free-fermion box: microscopic analog of Section 12.

Hard-core bosons on an interval map to free fermions (Girardeau). The
single-particle spectrum is quadratic, E_n = ħ² k_n² / 2m, not linear.
At fixed density ρ = N/a the exact ground-state energy is a cubic polynomial
in N and decomposes exactly as

    E = e_bulk a + e_surf + A/a

with no higher orders. The Section 12 extraction I_K = A / (π ħ v),
v = π ħ ρ / m, is therefore exact at every N:

    DD (nodes at both walls):  I_K = +1/12
    DN (node / antinode):      I_K = −1/24

These are not the relativistic CFT values (−1/24, +1/48). The protocol
still isolates a dimensionless, a-independent invariant; the Hamiltonian
is different.
"""

from __future__ import annotations

import math
import numpy as np

from .spectral import BC1D

PI = math.pi
I_K_TG_DD = 1.0 / 12.0
I_K_TG_DN = -1.0 / 24.0


def fermi_velocity(rho: float, hbar: float = 1.0, mass: float = 1.0) -> float:
    """v_F = π ħ ρ / m for spinless 1D fermions / TG bosons."""
    if rho <= 0 or hbar <= 0 or mass <= 0:
        raise ValueError("rho, hbar, and mass must be positive")
    return PI * hbar * rho / mass


def energy(N: int, a: float, bc: BC1D = "dd", hbar: float = 1.0, mass: float = 1.0) -> float:
    """Ground-state energy of N TG bosons (free fermions) on an interval of length a."""
    if N < 1:
        raise ValueError("N must be a positive integer")
    if a <= 0:
        raise ValueError("length a must be positive")
    pre = (hbar**2 * PI**2) / (2.0 * mass * a**2)
    if bc == "dd":
        # k_n = n π / a, n = 1..N; Σ n² = N(N+1)(2N+1)/6
        return pre * (N * (N + 1) * (2 * N + 1) / 6.0)
    if bc == "dn":
        # k_n = (n+1/2) π / a, n = 0..N-1; Σ (n+1/2)² = N(4N²-1)/12
        return pre * (N * (4 * N**2 - 1) / 12.0)
    raise ValueError(f"unknown boundary conditions {bc!r}")


def expansion(rho: float, bc: BC1D = "dd", hbar: float = 1.0, mass: float = 1.0) -> tuple[float, float, float]:
    """Exact (e_bulk, e_surf, A) in E = e_bulk a + e_surf + A/a at density ρ."""
    if rho <= 0:
        raise ValueError("density rho must be positive")
    k = (hbar**2 * PI**2) / mass
    e_bulk = k * rho**3 / 6.0
    if bc == "dd":
        e_surf = k * rho**2 / 4.0
        A = k * rho / 12.0
    elif bc == "dn":
        e_surf = 0.0
        A = -k * rho / 24.0
    else:
        raise ValueError(f"unknown boundary conditions {bc!r}")
    return e_bulk, e_surf, A


def extracted_I_K(A: float, rho: float, hbar: float = 1.0, mass: float = 1.0) -> float:
    """I_K = A / (π ħ v_F) with v_F = π ħ ρ / m."""
    v = fermi_velocity(rho, hbar=hbar, mass=mass)
    return A / (PI * hbar * v)


def I_K(bc: BC1D = "dd") -> float:
    if bc == "dd":
        return I_K_TG_DD
    if bc == "dn":
        return I_K_TG_DN
    raise ValueError(f"unknown boundary conditions {bc!r}")


def empirical_I_K(
    Ns: list[int] | np.ndarray,
    rho: float,
    bc: BC1D = "dd",
    hbar: float = 1.0,
    mass: float = 1.0,
) -> float:
    """Section 12 protocol: fit E(a) = e1 a + e0 + A/a at fixed ρ, return A/(π ħ v)."""
    Ns_arr = np.asarray(Ns, dtype=int)
    if Ns_arr.size < 3:
        raise ValueError("need at least three particle numbers to fit bulk, surface, A")
    a = Ns_arr / rho
    E = np.array([energy(int(n), float(ai), bc=bc, hbar=hbar, mass=mass) for n, ai in zip(Ns_arr, a)])
    x = np.column_stack([a, np.ones_like(a), 1.0 / a])
    coeff, *_ = np.linalg.lstsq(x, E, rcond=None)
    A = float(coeff[2])
    return extracted_I_K(A, rho, hbar=hbar, mass=mass)
