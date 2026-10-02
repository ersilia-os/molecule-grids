"""Computational drawing style: RDKit's colour convention with the ACS label proportions."""

from molecule_grids.depict import acs

# Line width and caption size relative to the bond, as in the original screen style
# (1.6 px lines and 14 px captions on a 25 px bond).
LINE_RATIO = 1.6 / 25
CAPTION_RATIO = 14 / 25

# RDKit's default hues, with the two that fail on white (cyan F, lime Cl) darkened.
PALETTE = {
    7: (0.13, 0.30, 0.85),
    8: (0.85, 0.10, 0.10),
    9: (0.05, 0.55, 0.55),
    15: (0.85, 0.45, 0.00),
    16: (0.75, 0.55, 0.00),
    17: (0.10, 0.55, 0.10),
    35: (0.60, 0.15, 0.10),
    53: (0.45, 0.10, 0.60),
}


def computational_options(opts):
    """Apply computational drawing options in place: RDKit colours, Arial, ACS label proportions."""
    opts.fontFile = acs.ARIAL
    opts.baseFontSize = acs.BASE_FONT_SIZE
    opts.bondLineWidth = LINE_RATIO * acs.BOND
    opts.scaleBondWidth = False
    opts.multipleBondOffset = acs.DOUBLE_OFFSET
    opts.additionalAtomLabelPadding = acs.LABEL_PADDING
    opts.annotationFontScale = 0.65
    opts.updateAtomPalette(PALETTE)
