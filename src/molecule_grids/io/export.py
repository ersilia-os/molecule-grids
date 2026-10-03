"""Render SVG pages to PNG (resvg) and vector PDF (typst).

Both renderers ship as self-contained wheels, so no system graphics library (cairo) is
needed. They use the bundled Liberation Sans only, so output is identical on every machine.
"""

import re
from pathlib import Path

from molecule_grids.depict import acs

FONTS = [acs.ARIAL, acs.ARIAL_BOLD]
MAX_PIXELS = 16000  # longest side of a raster, to keep memory in check


def _size_pt(svg):
    """``(width, height)`` of an SVG sized in pt."""
    head = svg.split(">", 2)[1]
    return tuple(float(re.search(rf"{dim}='([\d.]+)pt'", head).group(1)) for dim in ("width", "height"))


def to_png(svg, dpi):
    """Render an SVG string (sized in pt) to PNG bytes at ``dpi`` (lowered for very large figures)."""
    import resvg_py

    dpi = min(dpi, MAX_PIXELS * 72 / max(_size_pt(svg)))
    png = resvg_py.svg_to_bytes(
        svg_string=svg,
        dpi=dpi,
        background="#ffffff",
        font_files=FONTS,
        sans_serif_family="Liberation Sans",
        skip_system_fonts=True,
    )
    return bytes(png)


def to_pdf(svg):
    """Render an SVG string (sized in pt) to vector PDF bytes, one page of the same size."""
    import typst

    w, h = _size_pt(svg)
    page = f'#set page(width: {w}pt, height: {h}pt, margin: 0pt)\n#image("figure.svg", width: {w}pt, height: {h}pt)\n'
    return typst.compile(
        {"main.typ": page.encode(), "figure.svg": svg.encode()},
        format="pdf",
        font_paths=[str(Path(acs.ARIAL).parent)],
        ignore_system_fonts=True,
    )
