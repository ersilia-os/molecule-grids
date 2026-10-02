import itertools
import math
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
        grid = draw_grid(SMILES, NAMES, format=fmt, width=w)
        assert width_in(grid.pages[0]) == pytest.approx(w * full, abs=0.01)


@pytest.mark.parametrize("fmt, bond", [("print", 14.4), ("slide", 20.0)])
def test_size_presets(fmt, bond):
    for size, factor in (("small", 0.7), ("medium", 1.0), ("large", 1.4)):
        assert draw_grid(SMILES, format=fmt, size=size).bond_pt == pytest.approx(bond * factor)
    small, large = (draw_grid(SMILES, format=fmt, size=s, squeeze=False) for s in ("small", "large"))
    assert small.columns > large.columns
    tiny = draw_grid(SMILES, format=fmt, width=0.1, size="large")  # not even one column fits: shrunk
    assert tiny.bond_pt < bond * 1.4 and tiny.columns == 1
    with pytest.raises(ValueError, match="size"):
        draw_grid(SMILES, size="huge")


def test_save_types(tmp_path):
    grid = draw_grid(SMILES, NAMES, format="print", number=True)
    for ext in ("svg", "png", "pdf"):
        (path,) = grid.save(tmp_path / f"grid.{ext}")
        assert path.stat().st_size > 0


def test_pages_save_as_numbered_files(tmp_path):
    pool = [s for fam in FAMILIES for s, _ in fam] + [s for s, _ in SINGLES]
    grid = draw_grid(pool * 3, format="slide", width=0.5)
    assert len(grid.pages) > 1
    with pytest.raises(TooManyMolecules):
        _ = grid.svg
    paths = grid.save(tmp_path / "grid.png")
    assert [p.name for p in paths] == [f"grid_{k}.png" for k in range(1, len(grid.pages) + 1)]
    assert all(p.stat().st_size > 0 for p in paths)


def test_style_is_independent_of_format():
    colours = {}
    for style in ("medicinal", "computational"):
        for fmt in ("print", "slide"):
            svg = draw_grid(["OCCN"], format=fmt, style=style, frame=False).svg
            colours[style, fmt] = set(re.findall(r"#[0-9A-F]{6}", svg)) - {"#FFFFFF", "#000000"}
    assert not colours["medicinal", "print"] and not colours["medicinal", "slide"]  # black and white
    assert colours["computational", "print"] and colours["computational", "slide"]  # coloured atoms
    with pytest.raises(ValueError, match="style"):
        draw_grid(["CCO"], style="pretty")


def test_frame():
    assert "stroke='#50285A'" in draw_grid(["CCO"]).svg  # on by default; plum on slides
    assert "stroke='#000000'" in draw_grid(["CCO"], format="print").svg  # black in print
    assert "<rect x=" not in draw_grid(["CCO"], frame=False).svg


def ink_points(svg):
    """Absolute ink points of each molecule in a grid SVG."""
    groups = re.findall(r"<g transform='translate\(([\d.-]+),([\d.-]+)\)'>(.*?)</g>", svg, re.DOTALL)
    out = []
    for tx, ty, inner in groups:
        v = [float(t) for d in re.findall(r" d='([^']+)'", inner) for t in re.findall(r"-?\d+\.?\d*", d)]
        out.append([(x + float(tx), y + float(ty)) for x, y in zip(v[0::2], v[1::2])])
    return out


@pytest.mark.parametrize("mode", [True, "rows", "free"])
def test_squeeze(mode):
    square = draw_grid(SMILES, NAMES, format="print", squeeze=False)
    tight = draw_grid(SMILES, NAMES, format="print", squeeze=mode, number=True)
    assert tight.bond_pt == square.bond_pt
    assert tight.width_mm == pytest.approx(square.width_mm)  # the width never changes
    assert tight.height_mm < square.height_mm
    assert sorted(tight.order) == list(range(len(SMILES)))
    assert [tight.names[i] for i in tight.order] == [str(k + 1) for k in range(len(SMILES))]
    inks = ink_points(tight.svg)
    gap = min(math.dist(p, q) for a, b in itertools.combinations(inks, 2) for p in a for q in b)
    assert gap > 0.5 * 26  # inks never closer than half a bond (26 drawing units)
    again = draw_grid(SMILES, NAMES, format="print", squeeze=mode, number=True)
    assert again.svg == tight.svg  # deterministic


@pytest.mark.parametrize("fmt, max_in", [("print", 247 / 25.4), ("slide", 13 * 9 / 16)])
@pytest.mark.parametrize("mode", [False, "grid", "rows", "free"])
def test_page_height_cap(tmp_path, fmt, max_in, mode):
    pool = [s for fam in FAMILIES for s, _ in fam] + [s for s, _ in SINGLES]
    grid = draw_grid(pool, format=fmt, squeeze=mode)
    assert len(grid.pages) > 1
    heights = [float(re.search(r"height='([\d.]+)pt'", p).group(1)) / 72 for p in grid.pages]
    assert max(heights) <= max_in + 0.01
    assert sorted(grid.order) == list(range(len(pool)))
    assert len(grid.save(tmp_path / "grid.svg")) == len(grid.pages)


def test_one_row_shrinks_and_columns_are_used_ones():
    for mode in ("rows", False, "grid", "free"):
        one = draw_grid(["CC(=O)Oc1ccccc1C(=O)O"], format="slide", squeeze=mode)
        assert one.width_mm < 100 and one.columns == 1  # not the full 330 mm
    full = draw_grid(SMILES * 3, format="slide")
    assert full.width_mm == pytest.approx(330.2, abs=0.5)  # a full row keeps the full width


def test_max_height_paginates():
    pool = [s for fam in FAMILIES for s, _ in fam] + [s for s, _ in SINGLES]
    tall, short = draw_grid(pool, format="print"), draw_grid(pool, format="print", max_height=0.4)
    assert len(short.pages) > len(tall.pages)
    heights_mm = [float(re.search(r"height='([\d.]+)pt'", p).group(1)) * 25.4 / 72 for p in short.pages]
    assert max(heights_mm) <= 0.4 * 247 + 0.5
    with pytest.raises(ValueError, match="max_height"):
        draw_grid(["CCO"], max_height=0)


def test_errors():
    with pytest.raises(ValueError, match="invalid SMILES"):
        draw_grid(["CCO", "not-a-smiles"])
    with pytest.raises(ValueError, match="width"):
        draw_grid(SMILES, width=1.5)
