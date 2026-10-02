"""Rasterise SVG pages to PNG and PDF."""

import re


def _cairosvg():
    try:
        import cairosvg
    except OSError as e:  # the Python package is there, the cairo C library is not
        raise OSError(
            "PNG and PDF export need the cairo library: `brew install cairo` (macOS) "
            "or `apt install libcairo2` (Debian/Ubuntu)."
        ) from e
    return cairosvg


MAX_PIXELS = 16000  # longest side of a raster; cairo fails beyond 32767


def _dpi(svg, dpi):
    """``dpi``, lowered if needed so that the longest side stays within ``MAX_PIXELS``."""
    side_pt = max(float(v) for v in re.findall(r"(?:width|height)='([\d.]+)pt'", svg.split(">", 2)[1]))
    return min(dpi, MAX_PIXELS * 72 / side_pt)


def to_png(svg, dpi):
    """Render an SVG string (sized in pt) to PNG bytes at ``dpi`` (lowered for very large figures)."""
    return _cairosvg().svg2png(bytestring=svg.encode(), dpi=_dpi(svg, dpi), background_color="white")


def to_pdf(svg):
    """Render an SVG string to vector PDF bytes."""
    return _cairosvg().svg2pdf(bytestring=svg.encode())
