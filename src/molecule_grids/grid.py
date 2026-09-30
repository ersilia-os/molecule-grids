"""Molecule grids: the main entry point of the package."""

import html
import re
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


class TooManyMolecules(ValueError):
    """More molecules than fit on one page at the requested width and columns."""


class Grid:
    """A drawn grid of molecules, possibly split over several pages.

    Build it with :func:`draw_grid`. In Jupyter, a grid displays itself inline.

    Attributes
    ----------
    format : str
        ``"print"`` or ``"slide"``.
    names : list of str
        Caption of each molecule, in input order.
    columns : int
        Number of columns.
    max_columns : int
        Most columns that keep atom labels legible at this width, for these molecules.
    capacity : int
        Most molecules that fit on one page.
    bond_pt : float
        Bond length in pt when the figure is placed at 100 %.
    width_in, height_in : float
        Figure size in inches (height of the first page).
    pages : list of str
        One SVG document per page. More than one only when the molecules exceed ``capacity``.
    """

    def __init__(self, fmt, names, geometry, max_columns, capacity, pages):
        self.format = fmt.name
        self._fmt = fmt
        self.names = names
        self.columns = geometry.columns
        self.max_columns = max_columns
        self.capacity = capacity
        self.bond_pt = geometry.bond_pt
        self.width_in = geometry.width * geometry.scale / sizing.PT_PER_IN
        self.height_in = float(re.search(r"height='([\d.]+)pt'", pages[0]).group(1)) / sizing.PT_PER_IN
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
            more = f", or use more columns (up to {self.max_columns})" if self.columns < self.max_columns else ""
            raise TooManyMolecules(
                f"{len(self.names)} molecules need {len(self.pages)} pages at {self.columns} columns "
                f"(one page holds {self.capacity}). Save as .gif to get one frame per page{more}; "
                f"`.pages` holds each page as {what}."
            )

    @property
    def svg(self):
        """The figure as an SVG string.

        Raises
        ------
        TooManyMolecules
            If the grid spans several pages.
        """
        self._single_page("SVG")
        return self.pages[0]

    def save(self, path):
        """Save the figure; the file type follows the extension.

        Parameters
        ----------
        path : str or Path
            ``.svg``, ``.png``, ``.pdf`` or ``.gif``. A GIF has one frame per page and is the
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


def draw_grid(molecules, names=None, format="slide", width=1.0, columns=None, number=False, group=False):
    """Draw molecules as a grid of square cells.

    The figure width is ``width`` times the format's full width (7.09 in for print,
    13 in for slides), as in stylia. All molecules share one bond length: the largest
    at which the biggest molecule fits its cell, never above ACS 14.4 pt in print.

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
        Number of columns. By default, a sensible number for the width, widened (up to
        ``max_columns``) if the molecules would otherwise not fit on one page.
    number : bool
        Replace captions with bold compound numbers 1, 2, 3... in input order.
    group : bool
        Place analogues that share a scaffold side by side. They share one orientation
        regardless.

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
        while columns < max_columns and len(arrange(buckets, len(mols), columns, group)) > solve(columns).max_rows(fmt):
            columns += 1  # widen before paginating
    elif not 1 <= columns <= max_columns:
        raise ValueError(
            f"{columns} columns do not fit: at width={width} ({width * fmt.size_in:.2f} in) at most "
            f"{max_columns} keep atom labels at {fmt.min_label_pt:g} pt or larger. "
            f"Use fewer columns or a larger width."
        )
    geo = solve(columns)
    rows = arrange(buckets, len(mols), columns, group)
    max_rows = geo.max_rows(fmt)
    pages = _pages(geo, drawn, names, bold, rows, max_rows)
    return Grid(fmt, names, geo, max_columns, columns * max_rows, pages)


def _pages(geo, drawn, names, bold, rows, max_rows):
    """Lay out rows of drawn molecules into SVG pages of at most ``max_rows`` rows."""
    chunks = [rows[p : p + max_rows] for p in range(0, len(rows), max_rows)]
    H = geo.page_height(max_rows if len(chunks) > 1 else len(rows))  # frames share one canvas
    W, gap, cell = geo.width, geo.gap, geo.cell
    pages = []
    for page_rows in chunks:
        body, y = [], gap
        for r in page_rows:
            rh = max(drawn[i][4] for i in r)  # tallest ink in the row
            x = (W - (len(r) * cell + (len(r) - 1) * gap)) / 2  # centre partial rows
            for i in r:
                inner, x0, y0, w, h = drawn[i]
                cx, cy = x + cell / 2, y + cell / 2  # cell centre
                tx, ty = cx - w / 2 - x0, cy - h / 2 - y0  # ink centred in its square
                body.append(f"<g transform='translate({tx:.2f},{ty:.2f})'>{inner}</g>")
                if names[i]:  # captions share a baseline per row, just under the ink
                    weight = " font-weight='bold'" if bold[i] else ""
                    body.append(
                        f"<text x='{cx:.2f}' y='{cy + rh / 2 + geo.caption_base:.2f}' "
                        f"font-family=\"{acs.FONT_STACK}\" font-size='{geo.caption:.2f}' "
                        f"text-anchor='middle' fill='#000'{weight}>{html.escape(names[i])}</text>"
                    )
                x += cell + gap
            y += geo.row_pitch
        pages.append(
            "<?xml version='1.0' encoding='utf-8'?>\n"
            f"<svg xmlns='http://www.w3.org/2000/svg' width='{W * geo.scale:.2f}pt' "
            f"height='{H * geo.scale:.2f}pt' viewBox='0 0 {W:.2f} {H:.2f}'>\n"
            "<rect width='100%' height='100%' fill='#FFFFFF'/>\n" + "\n".join(body) + "\n</svg>\n"
        )
    return pages
