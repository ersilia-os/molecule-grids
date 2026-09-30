"""Rasterise SVG pages to PNG, PDF and animated GIF."""

import io


def _cairosvg():
    try:
        import cairosvg
    except OSError as e:  # the Python package is there, the cairo C library is not
        raise OSError(
            "PNG, PDF and GIF export need the cairo library: `brew install cairo` (macOS) "
            "or `apt install libcairo2` (Debian/Ubuntu)."
        ) from e
    return cairosvg


def to_png(svg, dpi):
    """Render an SVG string (sized in pt) to PNG bytes at ``dpi``."""
    return _cairosvg().svg2png(bytestring=svg.encode(), dpi=dpi, background_color="white")


def to_pdf(svg):
    """Render an SVG string to vector PDF bytes."""
    return _cairosvg().svg2pdf(bytestring=svg.encode())


def to_gif(pages, dpi, seconds=2.0):
    """Render SVG pages to an animated GIF, one frame per page, looping forever."""
    from PIL import Image

    frames = [Image.open(io.BytesIO(to_png(p, dpi))).convert("RGB") for p in pages]
    frames = [f.convert("P", palette=Image.Palette.ADAPTIVE, colors=64) for f in frames]
    out = io.BytesIO()
    frames[0].save(out, format="GIF", save_all=True, append_images=frames[1:], duration=int(seconds * 1000), loop=0)
    return out.getvalue()
