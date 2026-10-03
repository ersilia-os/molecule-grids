"""Figure sizing, harmonised with stylia.

A figure's width is a fraction of the format's base width (stylia's ``SIZE``): 7.09 in
for print (Nature two-column) and 13 in for slides. Molecules are drawn at a fixed size:
atom labels and captions use stylia's font sizes (print 5, 6 or 8 pt; slides 8, 10 or 13 pt
for small, medium or large), and the bond follows with ChemDraw ACS proportions
(bond = 1.44 x label). The width holds as
many square cells as fit the largest molecule at that size, one molecule per cell, and the
height follows from the number of rows, up to a page: 247 mm (9.72 in, Nature's maximum
figure height) in print and a 16:9 slide (7.31 in) on slides. stylia defines no maximum
height, so these come from the same sources as its widths.
"""

import math
from dataclasses import dataclass

from molecule_grids.depict import acs

PT_PER_IN = 72.0
MM_PER_PT = 25.4 / 72
SIZES = ("small", "medium", "large")


@dataclass(frozen=True)
class Format:
    """Output format.

    Attributes
    ----------
    name : str
        ``"print"`` or ``"slide"``.
    size_in : float
        Full figure width in inches (stylia ``SIZE``).
    label_pt : dict
        Atom label and caption size for each of ``SIZES``, in pt (stylia ``FONTSIZE_SMALL``,
        ``FONTSIZE`` and ``FONTSIZE_BIG``).
    max_height_in : float
        Tallest figure, in inches; taller grids are split into pages.
    aspect : float
        Width / height that ``width="auto"`` aims for.
    png_dpi : int
        Resolution of PNG exports.
    """

    name: str
    size_in: float
    label_pt: dict
    max_height_in: float
    aspect: float
    png_dpi: int


FORMATS = {
    # Nature two-column width and page height; stylia print font sizes.
    "print": Format("print", 7.09, dict(zip(SIZES, (5, 6, 8))), 247 / 25.4, 3 / 2, 300),
    # Wide 16:9 slide; stylia slide font sizes.
    "slide": Format("slide", 13.0, dict(zip(SIZES, (8, 10, 13))), 13.0 * 9 / 16, 16 / 9, 150),
}

GAP = 1.6  # gap between cells and page margin, in bond lengths


def get_format(name):
    """Return the :class:`Format` called ``name``."""
    try:
        return FORMATS[name]
    except KeyError:
        raise ValueError(f"format must be one of {sorted(FORMATS)}, not {name!r}") from None


@dataclass(frozen=True)
class Geometry:
    """Solved grid geometry, in internal drawing units unless stated otherwise."""

    columns: int
    cell: float  # square cell side
    pitch: float  # column pitch (cell plus its share of the spare width and gaps)
    caption: float  # caption font size
    scale: float  # pt per drawing unit
    width: float  # figure width (the most it may take)
    fit: bool = True  # narrow the figure to its content

    @property
    def gap(self):
        return GAP * acs.BOND

    @property
    def bond_pt(self):
        return acs.BOND * self.scale

    @property
    def caption_base(self):
        """Distance from the row's ink bottom to the caption baseline."""
        return 2.0 * self.caption

    @property
    def strip(self):
        """Caption band under each row."""
        return self.caption_base + 0.3 * self.caption

    def span(self, size):
        """Columns that a molecule needing ``size`` drawing units of width spans."""
        return min(self.columns, max(1, math.ceil((size + self.gap) / self.pitch - 1e-9)))


def solve(fmt, width, bond_pt, largest, typical, caption):
    """Solve the grid geometry at a fixed molecule size.

    Parameters
    ----------
    fmt : Format
        Output format.
    width : float
        Figure width as a fraction of ``fmt.size_in``.
    bond_pt : float
        Bond length for the chosen size, in pt.
    largest : float
        Width of the widest molecule (or caption), in drawing units; it must fit the figure.
    typical : float
        Cell side for a typical molecule, in drawing units; it sets the columns.
    caption : float
        Caption font size, in drawing units.

    Returns
    -------
    Geometry
        As many ``typical`` cells as fit; spare width is spread between columns. If the
        largest molecule does not fit the width, the scale is reduced until it does
        (check ``bond_pt``).
    """
    gap = GAP * acs.BOND
    width_pt = width * fmt.size_in * PT_PER_IN
    scale = bond_pt / acs.BOND
    scale = min(scale, width_pt / (2 * gap + largest))  # the largest molecule must fit
    total = width_pt / scale
    columns = max(1, int((total - gap) // (typical + gap)))
    pitch = (total - gap) / columns  # column centres evenly spaced; a cell is pitch - gap wide
    return Geometry(columns, typical, pitch, caption, scale, total)
