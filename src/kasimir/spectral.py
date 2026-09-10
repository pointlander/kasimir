"""Exact spectral formulae and the Θ × I_K reconstruction.

Natural units ħ = c = 1 unless `hbar` and `c` are passed in.

1D uses the mode-spacing conversion Θ_1D = ħ π c / a. Dirichlet–Dirichlet
gives I_K = −1/24; Dirichlet–Neumann gives I_K = +1/48 (Hurwitz ζ(-1, 1/2)
= +1/24, times the oscillator ħ/2). 3D uses the modular scale
Θ = ħ c / (2 π a), which makes I_3D = −(π³/360) A / a².

The invariant in both cases is the complexity production rate κ = E_Casimir.
"""

from __future__ import annotations

import math
from typing import Literal

import numpy as np

PI = math.pi
I_1D_NATS = -1.0 / 24.0
I_1D_DN_NATS = 1.0 / 48.0
GAMMA_3D = -(PI**3) / 360.0  # I_K = GAMMA_3D * A / a^2
THETA_MODULAR = "modular"  # Θ = ħc / (2π a)
BC1D = Literal["dd", "dn"]


def energy_1d_dirichlet(a: float, hbar: float = 1.0, c: float = 1.0) -> float:
    """Casimir energy of a 1D Dirichlet scalar on an interval of length a.

    E = − ħ π c / (24 a)
    """
    if a <= 0:
        raise ValueError("separation a must be positive")
    return -(hbar * PI * c) / (24.0 * a)


def force_1d_dirichlet(a: float, hbar: float = 1.0, c: float = 1.0) -> float:
    """F = −dE/da = − ħ π c / (24 a²)  (negative = attractive)."""
    if a <= 0:
        raise ValueError("separation a must be positive")
    return -(hbar * PI * c) / (24.0 * a**2)


def energy_1d_dirichlet_neumann(a: float, hbar: float = 1.0, c: float = 1.0) -> float:
    """Casimir energy of a 1D Dirichlet–Neumann scalar on an interval of length a.

    Frequencies ω_n = (n+1/2) π c / a, n = 0,1,…  The spectral sum is the
    Hurwitz value ζ(-1, 1/2) = +1/24, and E = (ħ π c / (2a)) × that sum
    = + ħ π c / (48 a). Positive energy that falls with a is repulsive.
    """
    if a <= 0:
        raise ValueError("separation a must be positive")
    return (hbar * PI * c) / (48.0 * a)


def force_1d_dirichlet_neumann(a: float, hbar: float = 1.0, c: float = 1.0) -> float:
    """F = −dE/da = + ħ π c / (48 a²)  (positive = repulsive)."""
    if a <= 0:
        raise ValueError("separation a must be positive")
    return (hbar * PI * c) / (48.0 * a**2)


def energy_3d_em(a: float, area: float = 1.0, hbar: float = 1.0, c: float = 1.0) -> float:
    """EM Casimir energy between parallel perfect conductors.

    E = − π² ħ c A / (720 a³)
    """
    if a <= 0 or area <= 0:
        raise ValueError("separation a and area must be positive")
    return -(PI**2 * hbar * c * area) / (720.0 * a**3)


def force_3d_em(a: float, area: float = 1.0, hbar: float = 1.0, c: float = 1.0) -> float:
    """F = − π² ħ c A / (240 a⁴)  (negative = attractive)."""
    if a <= 0 or area <= 0:
        raise ValueError("separation a and area must be positive")
    return -(PI**2 * hbar * c * area) / (240.0 * a**4)


def energy_3d_scalar_dirichlet(
    a: float, area: float = 1.0, hbar: float = 1.0, c: float = 1.0
) -> float:
    """Scalar Dirichlet slab: half of the EM result (one polarisation)."""
    return 0.5 * energy_3d_em(a, area=area, hbar=hbar, c=c)


def modular_temperature(a: float, hbar: float = 1.0, c: float = 1.0) -> float:
    """Θ(a) = ħ c / (2 π a). Energy per nat of modular conversion."""
    if a <= 0:
        raise ValueError("separation a must be positive")
    return (hbar * c) / (2.0 * PI * a)


def theta_1d_mode_spacing(a: float, hbar: float = 1.0, c: float = 1.0) -> float:
    """1D conversion: lowest-mode / mode-spacing scale ħ π c / a."""
    if a <= 0:
        raise ValueError("separation a must be positive")
    return (hbar * PI * c) / a


def information_1d(a: float | None = None, bc: BC1D = "dd") -> float:
    """Dimensionless 1D complexity I_K = E/Θ in nats. Independent of a."""
    if bc == "dd":
        return I_1D_NATS
    if bc == "dn":
        return I_1D_DN_NATS
    raise ValueError(f"unknown boundary conditions {bc!r}")


def information_3d(a: float, area: float = 1.0) -> float:
    """I_K(a) = −(π³/360) A / a² nats, using modular conversion."""
    if a <= 0 or area <= 0:
        raise ValueError("separation a and area must be positive")
    return GAMMA_3D * area / a**2


def reconstruct_energy(
    a: float,
    *,
    dim: Literal[1, 3] = 3,
    area: float = 1.0,
    hbar: float = 1.0,
    c: float = 1.0,
    bc: BC1D = "dd",
) -> float:
    """E = Θ(a) I_K(a). Must match the closed 1D/3D formulae."""
    if dim == 1:
        return theta_1d_mode_spacing(a, hbar=hbar, c=c) * information_1d(a, bc=bc)
    if dim == 3:
        return modular_temperature(a, hbar=hbar, c=c) * information_3d(a, area=area)
    raise ValueError("dim must be 1 or 3")


def reconstruct_force(
    a: float,
    *,
    dim: Literal[1, 3] = 3,
    area: float = 1.0,
    hbar: float = 1.0,
    c: float = 1.0,
    da: float | None = None,
    bc: BC1D = "dd",
) -> float:
    """F = −d(Θ I_K)/da by symmetric difference, as a check of the product rule."""
    if da is None:
        da = a * 1e-6
    e_plus = reconstruct_energy(a + da, dim=dim, area=area, hbar=hbar, c=c, bc=bc)
    e_minus = reconstruct_energy(a - da, dim=dim, area=area, hbar=hbar, c=c, bc=bc)
    return -(e_plus - e_minus) / (2.0 * da)


def log_two_sinh(x: np.ndarray | float) -> np.ndarray:
    """log(2 sinh x) for x>0, stable at both ends: x + log(1 - e^{-2x})."""
    x_arr = np.asarray(x, dtype=float)
    return x_arr + np.log1p(-np.exp(-2.0 * x_arr))


def oscillator_helmholtz(
    omega: np.ndarray, T: float, hbar: float = 1.0
) -> float:
    """F = T Σ log(2 sinh(ħω / 2T)), reducing to (ħ/2) Σ ω at T=0."""
    omega = np.asarray(omega, dtype=float)
    if T < 0:
        raise ValueError("temperature T must be nonnegative")
    if T == 0:
        return 0.5 * hbar * float(omega.sum())
    return float(T * np.sum(log_two_sinh(hbar * omega / (2.0 * T))))


def _thermal_mode_sum(gaps: np.ndarray, T: float) -> float:
    """T Σ log(1 - exp(-gap_n / T)). gaps are ħ ω_n."""
    return T * float(np.log1p(-np.exp(-gaps / T)).sum())


def _n_thermal_terms(a: float, T: float, hbar: float, c: float) -> int:
    # Need n ħ π c / (a T) ≳ 40 for the tail to die.
    return max(64, int(50.0 * a * T / (PI * hbar * c)) + 32)


def free_energy_1d_dirichlet(
    a: float, T: float, hbar: float = 1.0, c: float = 1.0
) -> float:
    """Zeta-regularised Helmholtz free energy of a 1D DD scalar.

    F(a,T) = −π ħ c / (24 a) + T Σ_{n=1}^∞ log(1 − exp(−n π ħ c / (a T)))
    """
    if a <= 0:
        raise ValueError("separation a must be positive")
    if T < 0:
        raise ValueError("temperature T must be nonnegative")
    e0 = energy_1d_dirichlet(a, hbar=hbar, c=c)
    if T == 0:
        return e0
    n = np.arange(1, _n_thermal_terms(a, T, hbar, c) + 1, dtype=float)
    gaps = n * PI * hbar * c / a
    return e0 + _thermal_mode_sum(gaps, T)


def free_energy_1d_dirichlet_neumann(
    a: float, T: float, hbar: float = 1.0, c: float = 1.0
) -> float:
    """Zeta-regularised Helmholtz free energy of a 1D DN scalar.

    F(a,T) = +π ħ c / (48 a) + T Σ_{n=0}^∞ log(1 − exp(−(n+1/2) π ħ c / (a T)))
    """
    if a <= 0:
        raise ValueError("separation a must be positive")
    if T < 0:
        raise ValueError("temperature T must be nonnegative")
    e0 = energy_1d_dirichlet_neumann(a, hbar=hbar, c=c)
    if T == 0:
        return e0
    n = np.arange(_n_thermal_terms(a, T, hbar, c), dtype=float)
    gaps = (n + 0.5) * PI * hbar * c / a
    return e0 + _thermal_mode_sum(gaps, T)


def free_energy_1d(
    a: float, T: float, bc: BC1D = "dd", hbar: float = 1.0, c: float = 1.0
) -> float:
    if bc == "dd":
        return free_energy_1d_dirichlet(a, T, hbar=hbar, c=c)
    if bc == "dn":
        return free_energy_1d_dirichlet_neumann(a, T, hbar=hbar, c=c)
    raise ValueError(f"unknown boundary conditions {bc!r}")
