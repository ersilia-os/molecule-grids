import itertools
import math
import re

import pytest

from molecule_grids import TooManyMolecules, draw_grid
from molecule_grids.data.examples import DEFAULT, FAMILIES, SINGLES

SMILES = [s for s, _ in DEFAULT]
NAMES = [n for _, n in DEFAULT]
POOL = [m for fam in FAMILIES for m in fam] + SINGLES  # 35 drugs, more than one row


def width_in(svg):
    return float(re.search(r"width='([\d.]+)pt'", svg).group(1)) / 72


@pytest.mark.parametrize("fmt, full", [("print", 7.09), ("slide", 13.0)])
def test_width_follows_format(fmt, full):
    for w in (1.0, 0.5):
        grid = draw_grid([s for s, _ in POOL], format=fmt, width=w, fit=False)
        assert width_in(grid.pages[0]) == pytest.approx(w * full, abs=0.01)  # exactly the width
        fitted = draw_grid([s for s, _ in POOL], format=fmt, width=w)
        assert width_in(fitted.pages[0]) <= w * full + 0.01  # never wider


@pytest.mark.parametrize("fmt, labels", [("print", (5, 6, 8)), ("slide", (8, 10, 13))])
def test_size_presets(fmt, labels):
    for size, label in zip(("small", "medium", "large"), labels):  # stylia font sizes
        grid = draw_grid(SMILES, NAMES, format=fmt, size=size)
        assert grid.bond_pt == pytest.approx(label * 14.4 / 10)  # ACS proportions
        size_pt = float(re.search(r"font-size='([\d.]+)'", grid.pages[0]).group(1)) * grid.bond_pt / 26
        assert size_pt == pytest.approx(label, abs=0.05)  # captions match the labels
    bond = labels[2] * 1.44
    small, large = (draw_grid(SMILES, format=fmt, width=1, size=s, squeeze=False) for s in ("small", "large"))
    assert small.columns > large.columns
    tiny = draw_grid(SMILES, format=fmt, width=0.1, size="large")  # not even one column fits: shrunk
    assert tiny.bond_pt < bond  # shrunk so the widest molecule fits
    with pytest.raises(ValueError, match="size"):
        draw_grid(SMILES, size="huge")


def test_save_types(tmp_path):
    grid = draw_grid(SMILES, NAMES, format="print", number=True)
    for ext in ("svg", "png", "pdf"):
        (path,) = grid.save(tmp_path / f"grid.{ext}")
        assert path.stat().st_size > 0
    from PIL import Image

    png = Image.open(tmp_path / "grid.png")  # same pixel size as the SVG at 300 dpi
    assert abs(png.width - grid.width_mm / 25.4 * 300) <= 2
    assert (tmp_path / "grid.pdf").read_bytes().startswith(b"%PDF")


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
    smiles = [s for s, _ in POOL[:20]]
    square = draw_grid(smiles, format="print", width=1, squeeze=False)
    tight = draw_grid(smiles, format="print", width=1, squeeze=mode, number=True)
    assert tight.bond_pt == square.bond_pt
    assert tight.width_mm <= 180.1 + 0.01  # fitted, never wider than the 180 mm asked for
    assert tight.width_mm * tight.height_mm < square.width_mm * square.height_mm
    exact = draw_grid(smiles, format="print", width=1, squeeze=mode, fit=False)
    assert exact.width_mm == pytest.approx(180.1, abs=0.2)  # fit=False keeps the width
    assert sorted(tight.order) == list(range(len(smiles)))
    assert [tight.names[i] for i in tight.order] == [str(k + 1) for k in range(len(smiles))]
    inks = ink_points(tight.svg)
    gap = min(math.dist(p, q) for a, b in itertools.combinations(inks, 2) for p in a for q in b)
    assert gap > 0.5 * 26  # inks never closer than half a bond (26 drawing units)
    again = draw_grid(smiles, format="print", width=1, squeeze=mode, number=True)
    assert again.svg == tight.svg  # deterministic


@pytest.mark.parametrize("fmt, max_in", [("print", 247 / 25.4), ("slide", 13 * 9 / 16)])
@pytest.mark.parametrize("mode", [False, "grid", "rows", "free"])
def test_page_height_cap(tmp_path, fmt, max_in, mode):
    pool = [s for s, _ in POOL] * 2  # 70 molecules: more than one page in every mode
    grid = draw_grid(pool, format=fmt, squeeze=mode, size="large")
    assert len(grid.pages) > 1
    heights = [float(re.search(r"height='([\d.]+)pt'", p).group(1)) / 72 for p in grid.pages]
    assert max(heights) <= max_in + 0.01
    assert sorted(grid.order) == list(range(len(pool)))
    assert len(grid.save(tmp_path / "grid.svg")) == len(grid.pages)


def test_one_row_shrinks_and_columns_are_used_ones():
    for mode in ("rows", False, "grid", "free"):
        one = draw_grid(["CC(=O)Oc1ccccc1C(=O)O"], format="slide", squeeze=mode)
        assert one.width_mm < 100 and one.columns == 1  # not the full 330 mm
    full = draw_grid(SMILES * 3, format="slide", width=1)
    assert full.width_mm > 300  # full rows keep (nearly) the full width


def test_max_height_paginates():
    pool = [s for fam in FAMILIES for s, _ in fam] + [s for s, _ in SINGLES]
    tall, short = draw_grid(pool, format="print"), draw_grid(pool, format="print", max_height=0.4)
    assert len(short.pages) > len(tall.pages)
    heights_mm = [float(re.search(r"height='([\d.]+)pt'", p).group(1)) * 25.4 / 72 for p in short.pages]
    assert max(heights_mm) <= 0.4 * 247 + 0.5
    with pytest.raises(ValueError, match="max_height"):
        draw_grid(["CCO"], max_height=0)


@pytest.mark.parametrize("fmt, aspect", [("print", 3 / 2), ("slide", 16 / 9)])
def test_auto_width_aims_at_aspect(fmt, aspect):
    smiles = [s for s, _ in POOL[:20]]
    auto = draw_grid(smiles, format=fmt)
    assert 0.25 <= auto.width <= 1
    shapes = [
        abs(math.log(g.width_mm / g.height_mm / aspect))
        for g in (draw_grid(smiles, format=fmt, width=w) for w in (0.25, 1.0))
    ]
    assert abs(math.log(auto.width_mm / auto.height_mm / aspect)) <= min(shapes) + 0.15


@pytest.mark.parametrize("mode", ["rows", "grid", False])
def test_group_keeps_families_together(mode):
    smiles = [s for s, _ in POOL]
    grid = draw_grid(smiles, format="print", squeeze=mode, group=True)
    where = {i: k for k, i in enumerate(grid.order)}
    for fam in FAMILIES:
        pos = sorted(where[smiles.index(s)] for s, _ in fam)
        assert pos == list(range(pos[0], pos[0] + len(pos)))  # consecutive in reading order


def test_theme_copies_match():
    from pathlib import Path

    root = Path(__file__).parents[1]
    repo = (root / ".streamlit" / "config.toml").read_text()  # used by Streamlit Community Cloud
    package = (root / "src" / "molecule_grids" / "app" / ".streamlit" / "config.toml").read_text()
    assert repo == package  # used by `molecule-grids app`


def test_errors():
    with pytest.raises(ValueError, match="invalid SMILES"):
        draw_grid(["CCO", "not-a-smiles"])
    with pytest.raises(ValueError, match="width"):
        draw_grid(SMILES, width=1.5)
