"""Generate the paper figures into ./figures/."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from . import field, lattice, spectral


def _save(fig, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=160, bbox_inches="tight")


def _finite_casimir(ns: list[int], bc: lattice.BC) -> tuple[np.ndarray, np.ndarray]:
    a = np.array([lattice.interval_length(n, bc=bc) for n in ns], dtype=float)
    e = np.array([lattice.vacuum_energy(n, bc=bc) for n in ns])
    finite = e - lattice.bulk_energy_density() * a
    lim = lattice.extrapolate_gamma(list(range(40, 161, 5)), bc=bc)
    return a, finite - lim.intercept


def plot_1d_lattice(out: Path) -> None:
    import matplotlib.pyplot as plt

    ns = list(range(8, 161))
    a_dd, casimir_dd = _finite_casimir(ns, "dd")
    a_dn, casimir_dn = _finite_casimir(ns, "dn")
    analytic_dd = np.array([spectral.energy_1d_dirichlet(ai) for ai in a_dd])
    analytic_dn = np.array([spectral.energy_1d_dirichlet_neumann(ai) for ai in a_dn])

    fig, axes = plt.subplots(1, 3, figsize=(11.4, 3.6))

    ax = axes[0]
    ax.plot(a_dd, casimir_dd, color="#1f4e79", lw=1.8, label="lattice, subtracted")
    ax.plot(a_dd, analytic_dd, color="#c44e52", ls="--", lw=1.6, label=r"$-\pi/(24a)$")
    ax.set_xlabel("interval length $a$")
    ax.set_ylabel("complexity production rate $\\kappa(a)$")
    ax.set_title("DD: $I_K=-1/24$, attractive")
    ax.legend(frameon=False, fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax = axes[1]
    ax.plot(a_dn, casimir_dn, color="#1f4e79", lw=1.8, label="lattice, subtracted")
    ax.plot(a_dn, analytic_dn, color="#c44e52", ls="--", lw=1.6, label=r"$+\pi/(48a)$")
    ax.axhline(0.0, color="#bbbbbb", lw=0.8)
    ax.set_xlabel("interval length $a$")
    ax.set_title(r"DN: $\zeta(-1,\frac{1}{2})=+1/24$, repulsive")
    ax.legend(frameon=False, fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax = axes[2]
    snap = np.array([lattice.snapshot_logdet(n) for n in ns])
    c0, c1 = lattice.snapshot_scaling(list(range(40, 161, 5)))
    ax.plot(a_dd, snap, color="#1f4e79", lw=1.8, label="snapshot $K$ proxy")
    ax.plot(
        a_dd, c0 + c1 * np.log(a_dd), color="#c44e52", ls="--", lw=1.6, label=r"$c_0+c_1\log a$"
    )
    ax.set_xlabel("interval length $a$")
    ax.set_ylabel(r"$\frac{1}{2}\sum\log(1/\omega_j)$")
    ax.set_title("Snapshot complexity is logarithmic")
    ax.legend(frameon=False, fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()
    _save(fig, out)
    plt.close(fig)


def plot_reconstruction(out: Path) -> None:
    import matplotlib.pyplot as plt

    a = np.linspace(0.2, 2.0, 200)
    area = 1.0
    e = np.array([spectral.energy_3d_em(ai, area=area) for ai in a])
    rec = np.array([spectral.reconstruct_energy(ai, dim=3, area=area) for ai in a])
    theta = np.array([spectral.modular_temperature(ai) for ai in a])
    info = np.array([spectral.information_3d(ai, area=area) for ai in a])
    f = np.array([spectral.force_3d_em(ai, area=area) for ai in a])
    f_rec = np.array([spectral.reconstruct_force(ai, dim=3, area=area) for ai in a])

    fig, axes = plt.subplots(1, 3, figsize=(11.2, 3.5))

    ax = axes[0]
    ax.plot(a, e, color="#1f4e79", lw=2.0, label=r"$E = -\pi^2 A/(720 a^3)$")
    ax.plot(a, rec, color="#c44e52", ls="--", lw=1.4, label=r"$\Theta(a)\,I_K(a)$")
    ax.set_xlabel("gap $a$")
    ax.set_ylabel("energy per unit area")
    ax.set_title("Energy from information × modular scale")
    ax.legend(frameon=False, fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax = axes[1]
    ax.plot(a, info, color="#1f4e79", lw=2.0, label=r"$I_K = -(\pi^3/360) A/a^2$")
    ax.set_xlabel("gap $a$")
    ax.set_ylabel("nats")
    ax.set_title("Dimensionless complexity deficit")
    ax.legend(frameon=False, fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax2 = ax.twinx()
    ax2.plot(a, theta, color="#6a9955", lw=1.4, ls=":", label=r"$\Theta=\hbar c/(2\pi a)$")
    ax2.set_ylabel("energy per nat", color="#6a9955")
    ax2.spines["top"].set_visible(False)
    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, frameon=False, fontsize=8)

    ax = axes[2]
    ax.plot(a, f, color="#1f4e79", lw=2.0, label=r"$F = -\pi^2 A/(240 a^4)$")
    ax.plot(a, f_rec, color="#c44e52", ls="--", lw=1.4, label=r"$-d(\Theta I_K)/da$")
    ax.set_xlabel("gap $a$")
    ax.set_ylabel("force per unit area")
    ax.set_title("Force from the product rule")
    ax.legend(frameon=False, fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()
    _save(fig, out)
    plt.close(fig)


def plot_quiet_cavity(out: Path) -> None:
    import matplotlib.pyplot as plt

    wide = field.dirichlet_gff_strip(
        width=256, height=96, mass=0.08, seed=1, normalize=False
    )
    narrow = field.dirichlet_gff_strip(
        width=256, height=32, mass=0.08, seed=1, normalize=False
    )
    # Shared colour scale so the narrow gap is visibly quieter.
    vmax = float(np.percentile(np.abs(wide), 99))

    fig, axes = plt.subplots(
        2, 1, figsize=(8.4, 5.2), gridspec_kw={"height_ratios": [3, 1]}
    )
    panels = (
        (axes[0], wide, "wide gap — IR modes present"),
        (axes[1], narrow, r"narrow gap — IR cut; closer to $\varphi=0$"),
    )
    for ax, data, title in panels:
        ax.imshow(
            data,
            cmap="RdBu_r",
            vmin=-vmax,
            vmax=vmax,
            aspect="auto",
            interpolation="nearest",
        )
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(title, fontsize=10, loc="left", pad=4)
        for spine in ax.spines.values():
            spine.set_color("#333333")
            spine.set_linewidth(1.4)

    fig.suptitle(
        r"Typical vacuum snapshots (Gaussian free field, Dirichlet walls). "
        r"Quietness is $\langle\varphi^2\rangle$, not $E_{\mathrm{Casimir}}$.",
        fontsize=10,
    )
    fig.tight_layout()
    _save(fig, out)
    plt.close(fig)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Write Kasimir theory figures.")
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("figures"),
        help="directory for PNG output",
    )
    args = parser.parse_args(argv)
    plot_1d_lattice(args.out / "1d-lattice.png")
    plot_reconstruction(args.out / "reconstruction.png")
    plot_quiet_cavity(args.out / "quiet-cavity.png")
    print(f"wrote figures in {args.out.resolve()}")


if __name__ == "__main__":
    main()
