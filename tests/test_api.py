import re

import pytest

from molecule_grids import TooManyMolecules, draw_grid
from molecule_grids.data.examples import DEFAULT, FAMILIES, SINGLES

SMILES = [s for s, _ in DEFAULT]
NAMES = [n for _, n in DEFAULT]


def width_in(svg):
    return float(re.search(r"width='([\d.]+)pt'", svg).group(1)) / 72


@pytest.mark.parametrize("fmt, full", [("print", 7.09), ("slide", 13.0)])
def test_width_follows_format(fmt, full):
    for w in (1.0, 0.5):
        grid = draw_grid(SMILES, NAMES, format=fmt, width=w, columns=2)
        assert width_in(grid.pages[0]) == pytest.approx(w * full, abs=0.01)


@pytest.mark.parametrize("fmt, cap", [("print", 14.4), ("slide", 20.0)])
def test_bond_capped(fmt, cap):
    assert draw_grid(["CCO"], format=fmt, columns=4).bond_pt == pytest.approx(cap)


def test_save_types(tmp_path):
    grid = draw_grid(SMILES, NAMES, number=True)
    for ext in ("svg", "png", "pdf"):
        assert grid.save(tmp_path / f"grid.{ext}").stat().st_size > 0


def test_capacity_and_gif(tmp_path):
    pool = [s for fam in FAMILIES for s, _ in fam] + [s for s, _ in SINGLES]
    assert len(draw_grid(pool, format="print", width=0.5).pages) == 1  # rows grow as needed
    grid = draw_grid(pool * 3, format="slide", width=0.5)  # > 100 molecules
    assert len(grid.pages) > 1
    with pytest.raises(TooManyMolecules):
        grid.save(tmp_path / "grid.svg")
    from PIL import Image

    gif = Image.open(grid.save(tmp_path / "grid.gif"))
    assert gif.n_frames == len(grid.pages)


def test_frame():
    assert "stroke='#000'" not in draw_grid(["CCO"]).svg
    assert "stroke='#000'" in draw_grid(["CCO"], frame=True).svg


def test_errors():
    with pytest.raises(ValueError, match="invalid SMILES"):
        draw_grid(["CCO", "not-a-smiles"])
    with pytest.raises(ValueError, match="columns do not fit"):
        draw_grid(SMILES, format="print", width=0.25, columns=10)
