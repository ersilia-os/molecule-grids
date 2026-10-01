"""Molecule grids: the main entry point of the package."""

import html
import math
import re
from collections import Counter
from pathlib import Path

from rdkit import Chem, RDLogger

from molecule_grids import sizing
from molecule_grids.depict import acs, slide
from molecule_grids.depict.draw import draw_one, text_width
from molecule_grids.layout.scaffold import arrange, prepare
from molecule_grids.utils.logging import logger

# A caption that is only a compound number (1, 2, 12, 3a, 10b ...) is set in bold (ACS style).
COMPOUND_NUMBER = re.compile(r"^\d+[a-z]{0,2}$")

MAX_COLUMNS = 30
MAX_MOLECULES = 100  # per figure; more are split into GIF frames
SQUEEZE_TOLERANCE = 0.10  # squeeze keeps more columns unless fewer save over 10 % of the area


class TooManyMolecules(ValueError):
    """More molecules than ``MAX_MOLECULES`` for a single figure."""


class Grid:
    """A drawn grid of molecules, possibly split over several pages.

    Build it with :func:`draw_grid`. In Jupyter, a grid displays itself inline.

    Attributes
    ----------
    format : str
        ``"print"`` or ``"slide"``.
    names : list of str
        Caption of each molecule, in input order.
    order : list of int
        Input indices in reading order (left to right, top to bottom, page by page).
    columns : int
        Number of columns (with ``squeeze``, possibly fewer than requested).
    max_columns : int
        Most columns that keep atom labels legible at this width, for these molecules.
    capacity : int
        Most molecules in one figure (``MAX_MOLECULES``, rounded down to whole rows).
    bond_pt : float
        Bond length in pt when the figure is placed at 100 %.
    width_in, height_in : float
        Figure size in inches (height of the first page).
    pages : list of str
        One SVG document per figure. More than one only when the molecules exceed ``capacity``.
    """

    def __init__(self, fmt, names, order, geometry, columns, max_columns, capacity, pages):
        self.format = fmt.name
        self._fmt = fmt
        self.names = names
        self.order = order
        self.columns = columns
        self.max_columns = max_columns
        self.capacity = capacity
        self.bond_pt = geometry.bond_pt
        self.width_in, self.height_in = (
            float(re.search(rf"{dim}='([\d.]+)pt'", pages[0]).group(1)) / sizing.PT_PER_IN
            for dim in ("width", "height")
        )
        self.pages = pages

    def __repr__(self):
        return (
            f"Grid(format={self.format!r}, molecules={len(self.names)}, columns={self.columns}, "
            f"pages={len(self.pages)}, size={self.width_in:.2f}x{self.height_in:.2f} in, "
            f"bond={self.bond_pt:.1f} pt)"
        )

    def _repr_svg_(self):
        return self.pages[0]

    def _single_page(self, what):
        if len(self.pages) > 1:
            raise TooManyMolecules(
                f"{len(self.names)} molecules exceed the {self.capacity} that fit in one figure. "
                f"Save as .gif to get {len(self.pages)} frames, or draw fewer molecules; "
                f"`.pages` holds each frame as {what}."
            )

    @property
    def svg(self):
        """The figure as an SVG string.

        Raises
        ------
        TooManyMolecules
            If the molecules exceed ``capacity``.
        """
        self._single_page("SVG")
        return self.pages[0]

    def save(self, path):
        """Save the figure; the file type follows the extension.

        Parameters
        ----------
        path : str or Path
            ``.svg``, ``.png``, ``.pdf`` or ``.gif``. A GIF has one frame per ``capacity`` molecules and is the
            only type that accepts more molecules than ``capacity``.

        Returns
        -------
        Path
            The written file.
        """
        from molecule_grids.io import export

        path = Path(path)
        ext = path.suffix.lower()
        if ext == ".gif":
            path.write_bytes(export.to_gif(self.pages, self._fmt.png_dpi))
        elif ext == ".svg":
            path.write_text(self.svg)
        elif ext == ".png":
            path.write_bytes(export.to_png(self.svg, self._fmt.png_dpi))
        elif ext == ".pdf":
            path.write_bytes(export.to_pdf(self.svg))
        else:
            raise ValueError(f"unsupported file type {ext!r}: use .svg, .png, .pdf or .gif")
        logger.success(f"Saved {path}")
        return path


def _to_mols(molecules):
    RDLogger.DisableLog("rdApp.*")
    try:
        mols = [Chem.Mol(m) if isinstance(m, Chem.Mol) else Chem.MolFromSmiles(str(m)) for m in molecules]
    finally:
        RDLogger.EnableLog("rdApp.*")
    bad = [i for i, m in enumerate(mols) if m is None]
    if bad:
        shown = ", ".join(f"{i} ({molecules[i]!r})" for i in bad[:5])
        raise ValueError(f"invalid SMILES at position {shown}" + (" ..." if len(bad) > 5 else ""))
    return mols


def draw_grid(
    molecules,
    names=None,
    format="slide",
    width=1.0,
    columns=None,
    number=False,
    group=False,
    frame=False,
    squeeze=False,
):
    """Draw molecules as a grid of square cells (or a tight grid with ``squeeze``).

    The figure width is ``width`` times the format's full width (7.09 in for print,
    13 in for slides), as in stylia. All molecules share one bond length: the largest
    at which the biggest molecule fits its cell, never above 14.4 pt in print (ACS)
    or 20 pt on slides.

    Parameters
    ----------
    molecules : sequence of str or rdkit.Chem.Mol
        SMILES (CXSMILES extensions such as ``|&1:1|`` are kept) or RDKit molecules.
    names : sequence of str, optional
        One caption per molecule.
    format : {"slide", "print"}
        ``"slide"``: RDKit colours. ``"print"``: ChemDraw ACS Document 1996 style.
    width : float
        Fraction of the format's full width, in (0, 1].
    columns : int, optional
        Number of columns. Defaults to a sensible number for the width. Rows are added as
        needed.
    number : bool
        Replace captions with bold compound numbers 1, 2, 3... in reading order.
    group : bool
        Place analogues that share a scaffold side by side. They share one orientation
        regardless. Ignored when ``squeeze`` is true.
    frame : bool
        Draw a thin black outline around the whole figure, in the format's bond line width.
    squeeze : bool
        Drop the square cells: each column is as wide as its widest molecule and each row
        as tall as its tallest, with molecules reordered to fill the figure. Same bond
        length; fewer columns are used if that saves over 10 % of the area. The figure comes
        out narrower than ``width``.

    Returns
    -------
    Grid

    Raises
    ------
    ValueError
        On invalid SMILES, a width outside (0, 1], or more columns than ``max_columns``.
    """
    fmt = sizing.get_format(format)
    if not 0 < width <= 1:
        raise ValueError(f"width is a fraction of the full {fmt.name} width and must be in (0, 1], not {width}")
    molecules = list(molecules)
    if not molecules:
        raise ValueError("no molecules to draw")
    if names is None:
        names = [m.GetProp("_Name") if isinstance(m, Chem.Mol) and m.HasProp("_Name") else "" for m in molecules]
    names = [str(k + 1) for k in range(len(molecules))] if number else [str(n) for n in names]
    if len(names) != len(molecules):
        raise ValueError(f"got {len(names)} names for {len(molecules)} molecules")

    mols, buckets = prepare(_to_mols(molecules))
    drawn = [draw_one(m, fmt.name) for m in mols]
    bold = [bool(COMPOUND_NUMBER.match(n)) for n in names]
    ink = max(max(d[3], d[4]) for d in drawn)
    ratio = acs.LABEL_RATIO if fmt.name == "print" else slide.CAPTION_RATIO

    def caption_widths(size):
        return max(text_width(n, size, b) for n, b in zip(names, bold))

    def solve(c):
        return sizing.solve(fmt, width, c, ink, caption_widths, ratio)

    max_columns = next((c for c in range(MAX_COLUMNS, 0, -1) if solve(c).bond_pt >= fmt.min_bond_pt), 1)
    if columns is None:
        target = round(width * fmt.size_in / fmt.target_cell_in)
        columns = max(1, min(target, max_columns, len(molecules)))
    elif not 1 <= columns <= max_columns:
        raise ValueError(
            f"{columns} columns do not fit: at width={width} ({width * fmt.size_in:.2f} in) at most "
            f"{max_columns} keep atom labels at {fmt.min_label_pt:g} pt or larger. "
            f"Use fewer columns or a larger width."
        )
    geo = solve(columns)
    max_rows = max(1, MAX_MOLECULES // columns)
    if squeeze:
        chunks = [
            list(range(k, min(k + columns * max_rows, len(mols)))) for k in range(0, len(mols), columns * max_rows)
        ]
        layouts = [_tight_layout(geo, drawn, names, bold, idx, columns) for idx in chunks]
    else:
        rows = arrange(buckets, len(mols), columns, group)
        layouts = [_square_layout(geo, drawn, rows[p : p + max_rows]) for p in range(0, len(rows), max_rows)]
        if len(layouts) > 1:  # GIF frames share the canvas of a full page
            layouts = [(W, geo.page_height(max_rows), place) for W, _, place in layouts]
    order = [i for _, _, place in layouts for i, *_ in place]
    if number:  # numbers follow reading order
        names = [None] * len(order)
        for k, i in enumerate(order):
            names[i] = str(k + 1)
    line = acs.BOND * (acs.LINE_RATIO if fmt.name == "print" else slide.LINE_RATIO)
    pages = _render(geo, layouts, drawn, names, bold, line if frame else None)
    used = max(Counter(round(b, 2) for *_, b in place).most_common(1)[0][1] for _, _, place in layouts)  # per row
    return Grid(fmt, names, order, geo, used, max_columns, columns * max_rows, pages)


def _square_layout(geo, drawn, rows):
    """Square cells of equal size; returns ``(W, H, [(i, cx, cy, baseline)])``."""
    W, gap, cell = geo.width, geo.gap, geo.cell
    place, y = [], gap
    for r in rows:
        rh = max(drawn[i][4] for i in r)  # tallest ink in the row
        x = (W - (len(r) * cell + (len(r) - 1) * gap)) / 2  # centre partial rows
        for i in r:
            cy = y + cell / 2
            place.append((i, x + cell / 2, cy, cy + rh / 2 + geo.caption_base))  # shared caption baseline
            x += cell + gap
        y += geo.row_pitch
    return W, geo.page_height(len(rows)), place


def _tight_layout(geo, drawn, names, bold, idx, columns):
    """Aligned grid with each column as wide as its widest molecule and each row as tall as its
    tallest. Tries every column count up to ``columns`` and two orderings, and keeps the
    smallest (fullest) figure; dropping columns must save more than ``SQUEEZE_TOLERANCE`` of the area. Returns ``(W, H, [(i, cx, cy, baseline)])``."""
    gap, strip = geo.gap, geo.caption_base + 0.3 * geo.caption
    cw = {i: max(drawn[i][3], text_width(names[i], geo.caption, bold[i]) if names[i] else 0) for i in idx}
    ch = {i: drawn[i][4] for i in idx}
    by_width = sorted(idx, key=lambda i: -cw[i])  # similar widths share a column
    by_height = sorted(idx, key=lambda i: -ch[i])  # similar heights share a row

    def candidates(C):
        R = math.ceil(len(idx) / C)
        cols = [sorted(by_width[c * R : (c + 1) * R], key=lambda i: -ch[i]) for c in range(C)]
        yield [[col[r] for col in cols if r < len(col)] for r in range(R)]
        yield [sorted(by_height[r * C : (r + 1) * C], key=lambda i: -cw[i]) for r in range(R)]

    def size(grid):
        colw = [max(cw[row[c]] for row in grid if c < len(row)) for c in range(max(map(len, grid)))]
        rowh = [max(ch[i] for i in row) for row in grid if row]
        return colw, rowh, (2 * gap + sum(colw) + (len(colw) - 1) * gap) * (gap + sum(h + strip + gap for h in rowh))

    grid = None
    for C in range(min(columns, len(idx)), 0, -1):
        best = min(candidates(C), key=lambda g: size(g)[2])
        if grid is None or size(best)[2] < (1 - SQUEEZE_TOLERANCE) * size(grid)[2]:
            grid = best  # fewer columns only when clearly fuller
    colw, rowh, _ = size(grid)
    place, y = [], gap
    for row, rh in zip((r for r in grid if r), rowh):
        x = gap
        for c, i in enumerate(row):
            place.append((i, x + colw[c] / 2, y + rh / 2, y + rh + geo.caption_base))
            x += colw[c] + gap
        y += rh + strip + gap
    return 2 * gap + sum(colw) + (len(colw) - 1) * gap, y, place


def _render(geo, layouts, drawn, names, bold, frame=None):
    """Write one SVG per layout on a shared canvas; ``frame`` is the outline width, or None."""
    W, H = max(lay[0] for lay in layouts), max(lay[1] for lay in layouts)
    outline = ""
    if frame:  # stroke drawn inside the canvas edge
        outline = (
            f"\n<rect x='{frame / 2:.2f}' y='{frame / 2:.2f}' width='{W - frame:.2f}' height='{H - frame:.2f}' "
            f"fill='none' stroke='#000' stroke-width='{frame:.2f}'/>"
        )
    pages = []
    for _, _, place in layouts:
        body = []
        for i, cx, cy, baseline in place:
            inner, x0, y0, w, h = drawn[i]
            tx, ty = cx - w / 2 - x0, cy - h / 2 - y0  # ink centred on (cx, cy)
            body.append(f"<g transform='translate({tx:.2f},{ty:.2f})'>{inner}</g>")
            if names[i]:
                weight = " font-weight='bold'" if bold[i] else ""
                body.append(
                    f"<text x='{cx:.2f}' y='{baseline:.2f}' "
                    f"font-family=\"{acs.FONT_STACK}\" font-size='{geo.caption:.2f}' "
                    f"text-anchor='middle' fill='#000'{weight}>{html.escape(names[i])}</text>"
                )
        pages.append(
            "<?xml version='1.0' encoding='utf-8'?>\n"
            f"<svg xmlns='http://www.w3.org/2000/svg' width='{W * geo.scale:.2f}pt' "
            f"height='{H * geo.scale:.2f}pt' viewBox='0 0 {W:.2f} {H:.2f}'>\n"
            "<rect width='100%' height='100%' fill='#FFFFFF'/>\n" + "\n".join(body) + outline + "\n</svg>\n"
        )
    return pages
