"""Figure sizing, harmonised with stylia.

A figure's width is a fraction of the format's base width (stylia's ``SIZE``): 7.09 in
for print (Nature two-column) and 13 in for slides. Molecules are drawn at a fixed size:
at ``size="medium"`` the bond is 14.4 pt in print (ACS) and 20 pt on slides. The width holds as
many square cells as fit the largest molecule at that size, one molecule per cell, and the
height follows from the number of rows, up to a page: 247 mm (9.72 in, Nature's maximum
figure height) in print and a 16:9 slide (7.31 in) on slides. stylia defines no maximum
height, so these come from the same sources as its widths.
"""

from dataclasses import dataclass

from molecule_grids.depict import acs

PT_PER_IN = 72.0
SIZES = {"small": 0.7, "medium": 1.0, "large": 1.4}  # molecule size, relative to medium


@dataclass(frozen=True)
class Format:
    """Output format.

    Attributes
    ----------
    name : str
        ``"print"`` or ``"slide"``.
    size_in : float
        Full figure width in inches (stylia ``SIZE``).
    bond_pt : float
        Bond length at ``size="medium"``, in pt.
    max_height_in : float
        Tallest figure, in inches; taller grids are split into pages.
    png_dpi : int
        Resolution of PNG and GIF exports.
    """

    name: str
    size_in: float
    bond_pt: float
    max_height_in: float
    png_dpi: int


FORMATS = {
    "print": Format("print", 7.09, acs.ACS["bond"], 247 / 25.4, 300),  # Nature two-column width, page height
    "slide": Format("slide", 13.0, 20.0, 13.0 * 9 / 16, 150),  # 16:9 slide; labels about 14 pt
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
    width: float  # figure width

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
    def row_pitch(self):
        return self.cell + self.caption_base + 0.3 * self.caption + self.gap

    def page_height(self, rows):
        return 2 * self.gap + rows * self.row_pitch - self.gap


def solve(fmt, width, factor, needed, caption):
    """Solve the grid geometry at a fixed molecule size.

    Parameters
    ----------
    fmt : Format
        Output format.
    width : float
        Figure width as a fraction of ``fmt.size_in``.
    factor : float
        Molecule size relative to ``fmt.bond_pt`` (a value of ``SIZES``).
    needed : float
        Smallest cell side that fits every molecule and caption, in drawing units.
    caption : float
        Caption font size, in drawing units.

    Returns
    -------
    Geometry
        As many columns as fit; spare width is spread between columns, so cells stay just
        large enough. If not even one cell fits, the scale is reduced until it does
        (check ``bond_pt``).
    """
    gap = GAP * acs.BOND
    width_pt = width * fmt.size_in * PT_PER_IN
    scale = fmt.bond_pt * factor / acs.BOND
    scale = min(scale, width_pt / (2 * gap + needed))  # one cell must fit
    total = width_pt / scale
    columns = max(1, int((total - gap) // (needed + gap)))
    pitch = (total - 2 * gap + gap) / columns  # column centres are evenly spaced across the width
    return Geometry(columns, needed, pitch, caption, scale, total)
