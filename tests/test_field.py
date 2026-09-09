import numpy as np
import pytest

from kasimir import field


def test_gff_shape_and_finite():
    phi = field.dirichlet_gff_strip(width=64, height=24, mass=0.2, seed=3)
    assert phi.shape == (24, 64)
    assert np.isfinite(phi).all()
    assert phi.std() == pytest.approx(1.0, rel=1e-6)


def test_gff_reproducible():
    a = field.dirichlet_gff_strip(width=32, height=16, seed=9)
    b = field.dirichlet_gff_strip(width=32, height=16, seed=9)
    np.testing.assert_allclose(a, b)


def test_gff_rejects_tiny_grid():
    with pytest.raises(ValueError):
        field.dirichlet_gff_strip(width=2, height=2)


def test_narrow_gap_is_quieter():
    wide = field.dirichlet_gff_strip(
        width=128, height=64, mass=0.08, seed=2, normalize=False
    )
    narrow = field.dirichlet_gff_strip(
        width=128, height=16, mass=0.08, seed=2, normalize=False
    )
    assert float(narrow.var()) < float(wide.var())
