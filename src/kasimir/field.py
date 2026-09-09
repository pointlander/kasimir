"""Typical vacuum snapshots of a 2D Gaussian free field with Dirichlet lines.

Used only as a figure: the cavity is quieter (⟨φ²⟩ suppressed). This is the
snapshot diagnostic of paper §3 / §11.4, *not* the complexity rate that equals
Casimir energy.
"""

from __future__ import annotations

import numpy as np


def dirichlet_gff_strip(
    width: int,
    height: int,
    mass: float = 0.15,
    seed: int = 0,
    normalize: bool = True,
) -> np.ndarray:
    """Sample a massive Gaussian free field on a grid with φ=0 on the top and bottom.

    Modes: φ(x, z) = Σ_{n,k} a_{nk} sin(n π z / H) exp(2π i k x / W) / sqrt(ω),
    with ω² = k_x² + (nπ/H)² + m², a_{nk} ~ ComplexNormal, then take the real field.

    `height` is the number of interior z-points; Dirichlet walls sit on the
    boundaries so the displayed array has shape (height, width) and vanishes
    if extended by one row of zeros on each side.
    """
    if width < 4 or height < 4:
        raise ValueError("grid too small")
    rng = np.random.default_rng(seed)
    H = height + 1  # Dirichlet length in sites
    zs = np.arange(1, height + 1)
    n_modes = np.arange(1, height + 1)
    kx = 2.0 * np.pi * np.fft.fftfreq(width)

    # Draw independent real coefficients in sine-z, Fourier-x basis.
    # ω² = kx² + (nπ/H)² + m²
    n_pi_h = n_modes * np.pi / H
    omega2 = kx[None, :] ** 2 + n_pi_h[:, None] ** 2 + mass**2
    sigma = 1.0 / np.sqrt(np.sqrt(omega2))  # amplitude ~ ω^{-1/2} for a snapshot

    # Real Fourier: k and -k are conjugate. Draw real and imag in frequency.
    amp = rng.normal(size=(height, width)) * sigma
    # sine transform in z, inverse FFT in x
    # φ(z, x) = Σ_n sin(n π z / H) * field_n(x)
    sine = np.sin(np.outer(zs, n_pi_h))  # (height, n_modes) but n_modes == height
    # For each n, spatial field is ifft of amp[n, :]
    spatial_n = np.fft.ifft(amp, axis=1).real * width**0.5
    sampled = sine @ spatial_n
    if normalize:
        sampled = sampled / (np.std(sampled) + 1e-12)
    return sampled
