"""Figure sizing, harmonised with stylia.

A figure's width is a fraction of the format's base width (stylia's ``SIZE``): 7.09 in
for print (Nature two-column) and 13 in for slides. The width is split into square cells,
one per molecule, so the height follows from the number of rows. All molecules share one
bond length: the largest at which the biggest molecule fits its cell, capped at ACS
14.4 pt in print. Legibility sets a floor on the bond length (atom labels no smaller than
stylia's ``FONTSIZE_SMALL``), which in turn bounds the columns; the page height bounds
the rows.
"""

from dataclasses import dataclass

from molecule_grids.depict import acs

PT_PER_IN = 72.0


@dataclass(frozen=True)
class Format:
    """Output format.

    Attributes
    ----------
    name : str
        ``"print"`` or ``"slide"``.
    size_in : float
        Full figure width in inches (stylia ``SIZE``).
    max_height_in : float
        Tallest figure that fits a page or slide, in inches.
    min_label_pt : float
        Smallest atom label and caption, in pt (stylia ``FONTSIZE_SMALL``).
    max_bond_pt : float or None
        Largest bond length, in pt (ACS 14.4 pt for print; none for slides).
    target_cell_in : float
        Cell size used to pick a default number of columns, in inches.
    png_dpi : int
        Resolution of PNG and GIF exports.
    """

    name: str
    size_in: float
    max_height_in: float
    min_label_pt: float
    max_bond_pt: float | None
    target_cell_in: float
    png_dpi: int

    @property
    def min_bond_pt(self):
        """Shortest legible bond: atom labels are ``LABEL_RATIO`` of the bond length."""
        return self.min_label_pt / acs.LABEL_RATIO


FORMATS = {
    # Nature two-column width; 247 mm maximum page height.
    "print": Format("print", 7.09, 9.7, 5.0, acs.ACS["bond"], 1.75, 300),
    # Wide slide (13 in); height left under a title on a 16:9 slide (7.3 in tall).
    "slide": Format("slide", 13.0, 6.5, 8.0, None, 1.9, 150),
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

    def max_rows(self, fmt):
        return max(1, int((fmt.max_height_in * PT_PER_IN / self.scale - self.gap) // self.row_pitch))


def solve(fmt, width, columns, ink, caption_widths, caption_ratio):
    """Solve the grid geometry for a given number of columns.

    Parameters
    ----------
    fmt : Format
        Output format.
    width : float
        Figure width as a fraction of ``fmt.size_in``.
    columns : int
        Number of columns.
    ink : float
        Largest ink extent (width or height) over all molecules, in drawing units.
    caption_widths : callable
        Maps a caption font size to the widest caption, in drawing units.
    caption_ratio : float
        Caption font size relative to the bond length.

    Returns
    -------
    Geometry
    """
    gap = GAP * acs.BOND
    width_pt = width * fmt.size_in * PT_PER_IN
    caption = caption_ratio * acs.BOND
    for _ in range(4):  # caption size and cell size depend on each other; converges fast
        cell = max(ink, caption_widths(caption))
        scale = width_pt / (2 * gap + columns * cell + (columns - 1) * gap)
        if fmt.max_bond_pt is not None:
            scale = min(scale, fmt.max_bond_pt / acs.BOND)
        caption = max(caption_ratio * acs.BOND, fmt.min_label_pt / scale)
    total = width_pt / scale
    cell = (total - 2 * gap - (columns - 1) * gap) / columns  # fill the width exactly
    return Geometry(columns, cell, caption, scale, total)
