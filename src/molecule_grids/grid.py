"""Molecule grids: the main entry point of the package."""

import html
import re
from pathlib import Path

from rdkit import Chem, RDLogger

from molecule_grids import sizing
from molecule_grids.depict import acs, slide
from molecule_grids.depict.draw import draw_one, ink_mask, rotate, text_width
from molecule_grids.layout import pack
from molecule_grids.layout.scaffold import arrange, prepare
from molecule_grids.utils.logging import logger

# A caption that is only a compound number (1, 2, 12, 3a, 10b ...) is set in bold (ACS style).
COMPOUND_NUMBER = re.compile(r"^\d+[a-z]{0,2}$")

MAX_MOLECULES = 100  # per figure; more are split into GIF frames
SQUEEZE = {False: None, None: None, True: "grid", "grid": "grid", "rows": "rows", "free": "free"}


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
    columns : int or None
        Number of columns, set by the width (when squeezing, possibly fewer; None for
        ``squeeze="free"``).
    capacity : int
        Most molecules in one figure: up to ``MAX_MOLECULES`` (100), fewer when a slide's
        height (7.31 in, 16:9) is reached.
    bond_pt : float
        Bond length in pt when the figure is placed at 100 %.
    width_in, height_in : float
        Figure size in inches (height of the first page).
    pages : list of str
        One SVG document per figure. More than one only when the molecules exceed ``capacity``.
    """

    def __init__(self, fmt, names, order, geometry, columns, capacity, pages):
        self.format = fmt.name
        self._fmt = fmt
        self.names = names
        self.order = order
        self.columns = columns
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
            ``.svg``, ``.png``, ``.pdf`` or ``.gif``. A GIF has one frame per page (``capacity`` molecules) and is the
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
    size="medium",
    number=False,
    group=False,
    frame=False,
    squeeze=False,
):
    """Draw molecules as a grid of square cells, or squeezed to fill the figure.

    The figure width is ``width`` times the format's full width (7.09 in for print,
    13 in for slides), as in stylia. Molecules are drawn at a fixed size, 14.4 pt bonds
    in print (ACS) and 20 pt on slides at ``size="medium"``; the width holds as many columns
    as fit at that size.

    Parameters
    ----------
    molecules : sequence of str or rdkit.Chem.Mol
        SMILES (CXSMILES extensions such as ``|&1:1|`` are kept) or RDKit molecules.
    names : sequence of str, optional
        One caption per molecule.
    format : {"slide", "print"}
        ``"slide"``: RDKit colours. ``"print"``: ChemDraw ACS Document 1996 style.
    width : float
        Fraction of the format's full width, in (0, 1]. Together with ``size`` it sets the
        number of columns. Rows are added as needed.
    size : {"small", "medium", "large"}
        Molecule size: 70, 100 or 140 % of the format's bond length (print 10, 14.4 or 20 pt;
        slide 14, 20 or 28 pt). If the largest molecule would not fit the width even in one
        column, it is shrunk to fit (see ``Grid.bond_pt``).
    number : bool
        Replace captions with bold compound numbers 1, 2, 3... in reading order.
    group : bool
        Place analogues that share a scaffold side by side. They share one orientation
        regardless. Ignored when squeezing.
    frame : bool
        Draw a thin black outline around the whole figure, in the format's bond line width.
    squeeze : bool or {"grid", "rows", "free"}
        Drop the square cells and reorder (and turn by quarter turns) molecules to fill the
        figure, at the same bond length; the figure comes out smaller than ``width``.
        Analogues sharing a scaffold turn together. ``"grid"`` (or True): aligned columns
        and rows fitted to their largest molecule; fewer columns only if that saves over
        10 % of the area. ``"rows"``: rows flow like justified text. ``"free"``: molecules
        interlock by their outlines.

    Returns
    -------
    Grid

    Raises
    ------
    ValueError
        On invalid SMILES, a width out of range, or an unknown size.
    """
    fmt = sizing.get_format(format)
    if not 0 < width <= 1:
        raise ValueError(f"width is a fraction of the full {fmt.name} width and must be in (0, 1], not {width}")
    if size not in sizing.SIZES:
        raise ValueError(f"size must be one of {list(sizing.SIZES)}, not {size!r}")
    factor = sizing.SIZES[size]
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
    caption = (acs.LABEL_RATIO if fmt.name == "print" else slide.CAPTION_RATIO) * acs.BOND
    needed = max(ink, max(text_width(n, caption, b) for n, b in zip(names, bold)))
    geo = sizing.solve(fmt, width, factor, needed, caption)
    if geo.bond_pt < fmt.bond_pt * factor - 1e-6:
        logger.warning(
            f"The largest molecule does not fit width={width} at size={size!r}; "
            f"drawn at {geo.bond_pt:.1f} pt bonds instead of {fmt.bond_pt * factor:.1f}"
        )
    columns = geo.columns
    max_rows = max(1, MAX_MOLECULES // columns)
    max_h = fmt.max_height_in * sizing.PT_PER_IN / geo.scale if fmt.max_height_in else float("inf")
    if fmt.max_height_in:  # rows that fit the page height
        max_rows = max(1, min(max_rows, int((max_h - geo.gap) // geo.row_pitch)))
    variants = {(i, 0): d for i, d in enumerate(drawn)}

    def get(i, r):
        """Drawing of molecule ``i`` turned ``r`` quarter turns (drawn on first use)."""
        if (i, r) not in variants:
            variants[i, r] = draw_one(rotate(mols[i], r), fmt.name)
        return variants[i, r]

    mode = SQUEEZE.get(squeeze, "bad") if isinstance(squeeze, (bool, str, type(None))) else "bad"
    if mode == "bad":
        raise ValueError(f"squeeze must be False, True, 'grid', 'rows' or 'free', not {squeeze!r}")
    if mode:
        group = {i: idxs[0] for key, idxs in buckets.items() if len(idxs) > 1 and key != "(acyclic)" for i in idxs}
        group = [group.get(i, i) for i in range(len(mols))]
        capw = [text_width(n, geo.caption, b) if n else 0.0 for n, b in zip(names, bold)]
        squeezer = {"grid": _grid_layout, "rows": _rows_layout, "free": _free_layout}[mode]
        layouts = _paginate(len(mols), max_h, lambda idx: squeezer(geo, get, idx, columns, capw, group))
    else:
        rows = arrange(buckets, len(mols), columns, group)
        layouts = [_square_layout(geo, drawn, rows[p : p + max_rows]) for p in range(0, len(rows), max_rows)]
        if len(layouts) > 1:  # GIF frames share the canvas of a full page
            layouts = [(W, geo.page_height(max_rows), place, c) for W, _, place, c in layouts]
    order = [i for _, _, place, _ in layouts for i, *_ in place]
    if number:  # numbers follow reading order
        names = [None] * len(order)
        for k, i in enumerate(order):
            names[i] = str(k + 1)
    line = acs.BOND * (acs.LINE_RATIO if fmt.name == "print" else slide.LINE_RATIO)
    pages = _render(geo, layouts, get, names, bold, line if frame else None)
    used = None if mode == "free" else max(c for *_, c in layouts)
    if len(layouts) > 1:  # what the first page holds
        capacity = len(layouts[0][2])
    else:  # squeezed pages may hold more than the square grid
        capacity = max(columns * max_rows, len(mols)) if mode else columns * max_rows
    return Grid(fmt, names, order, geo, used, capacity, pages)


def _paginate(n, max_h, build):
    """Split molecules 0..n-1, in input order, into squeezed pages no taller than ``max_h``.

    Each page takes as many molecules as fit (at most ``MAX_MOLECULES``), found by bisection.
    """
    layouts, start = [], 0
    while start < n:
        k = min(n - start, MAX_MOLECULES)
        lay = build(list(range(start, start + k)))
        if lay[1] > max_h:
            lo, hi, lay_lo = 1, k - 1, build([start])  # largest k that fits (lay_lo is its layout)
            while lo < hi:
                mid = (lo + hi + 1) // 2
                trial = build(list(range(start, start + mid)))
                lo, hi, lay_lo = (mid, hi, trial) if trial[1] <= max_h else (lo, mid - 1, lay_lo)
            k, lay = lo, lay_lo
        layouts.append(lay)
        start += k
    return layouts


# Layouts return (W, H, [(i, r, cx, cy, baseline)], columns): figure size, and per molecule its
# quarter turns, ink centre and caption baseline, in reading order.


def _square_layout(geo, drawn, rows):
    """Square cells of equal size."""
    W, gap, cell, pitch = geo.width, geo.gap, geo.cell, geo.pitch
    place, y = [], gap
    for r in rows:
        rh = max(drawn[i][4] for i in r)  # tallest ink in the row
        x = (W - len(r) * pitch) / 2  # centre partial rows
        for i in r:
            cy = y + cell / 2
            place.append((i, 0, x + pitch / 2, cy, cy + rh / 2 + geo.caption_base))  # shared caption baseline
            x += pitch
        y += geo.row_pitch
    return W, geo.page_height(len(rows)), place, max(map(len, rows))


def _strip(geo):
    return geo.caption_base + 0.3 * geo.caption  # caption band under each row


def _dims(get, idx, turns):
    return {(i, r): get(i, r)[3:5] for i in idx for r in turns}


def _grid_layout(geo, get, idx, columns, capw, group):
    """Tight aligned grid (see :func:`molecule_grids.layout.pack.grid_pack`)."""
    gap, strip = geo.gap, _strip(geo)
    dims = _dims(get, idx, (0, 1))
    rows = pack.grid_pack(idx, columns, dims, capw, group, gap, strip)

    def box(i, r):
        return max(dims[i, r][0], capw[i]), dims[i, r][1]

    ncol = len(rows[0])
    colw = [max(box(i, r)[0] for row in rows for i, r in [row[c]] if i is not None) for c in range(ncol)]
    place, y = [], gap
    for row in rows:
        rh = max(box(i, r)[1] for i, r in row if i is not None)
        x = gap
        for c, (i, r) in enumerate(row):
            if i is not None:
                place.append((i, r, x + colw[c] / 2, y + rh / 2, y + rh + geo.caption_base))
            x += colw[c] + gap
        y += rh + strip + gap
    return gap + sum(w + gap for w in colw), y, place, ncol


def _rows_layout(geo, get, idx, columns, capw, group):
    """Rows like justified text (see :func:`molecule_grids.layout.pack.rows_pack`)."""
    gap, strip = geo.gap, _strip(geo)
    dims = _dims(get, idx, (0, 1))
    rows = pack.rows_pack(idx, geo.width, dims, capw, group, gap, strip)

    def box(i, r):
        return max(dims[i, r][0], capw[i]), dims[i, r][1]

    natural = [gap + sum(box(i, r)[0] + gap for i, r in row) for row in rows]
    W = max(natural)
    place, y = [], gap
    for k, row in enumerate(rows):
        rh = max(box(i, r)[1] for i, r in row)
        last = k == len(rows) - 1
        if len(row) > 1 and not last:  # justify: spread the slack over the inner gaps
            x, step = gap, gap + (W - natural[k]) / (len(row) - 1)
        else:  # last or lone molecule: centred, natural gaps
            x, step = (W - natural[k]) / 2 + gap, gap
        for i, r in row:
            bw = box(i, r)[0]
            place.append((i, r, x + bw / 2, y + rh / 2, y + rh + geo.caption_base))
            x += bw + step
        y += rh + strip + gap
    return W, y, place, max(map(len, rows))


def _free_layout(geo, get, idx, columns, capw, group):
    """Interlocking outlines (see :func:`molecule_grids.layout.pack.free_pack`)."""
    gap, bond = geo.gap, acs.BOND
    res, clearance = bond / 2, bond / 2  # mask pixel and half the minimum distance between inks
    caption = {i: (capw[i], geo.caption_base - 0.75 * geo.caption, geo.caption_base + 0.3 * geo.caption) for i in idx}
    masks, offset = {}, {}
    for i in idx:
        for r in range(4):
            masks[i, r], *offset[i, r] = ink_mask(get(i, r), res, clearance, caption[i] if capw[i] else None)
    placed = pack.free_pack(idx, geo.width - 2 * gap, masks, res, group)
    items = []
    for i, r, x, y in placed:
        _, _, _, w, h = get(i, r)
        left, top = x - offset[i, r][0], y - offset[i, r][1]  # ink top-left
        items.append((i, r, left + w / 2, top + h / 2, top + h + geo.caption_base, left, top, w))
    # Shift the content to the margins; canvas = content extent.
    x0 = min(min(left, cx - capw[i] / 2) for i, _, cx, _, _, left, _, _ in items)
    y0 = min(top for *_, top, _ in items)
    x1 = max(max(left + w, cx + capw[i] / 2) for i, _, cx, _, _, left, _, w in items)
    y1 = max(b + 0.3 * geo.caption if capw[i] else cy + get(i, r)[4] / 2 for i, r, _, cy, b, *_ in items)
    dx, dy = gap - x0, gap - y0
    place = [(i, r, cx + dx, cy + dy, b + dy) for i, r, cx, cy, b, *_ in items]
    place.sort(key=lambda p: (round(p[3] / (3 * bond)), p[2]))  # reading order: bands, then left to right
    return x1 - x0 + 2 * gap, y1 - y0 + 2 * gap, place, None


def _render(geo, layouts, get, names, bold, frame=None):
    """Write one SVG per layout on a shared canvas; ``frame`` is the outline width, or None."""
    W, H = max(lay[0] for lay in layouts), max(lay[1] for lay in layouts)
    outline = ""
    if frame:  # stroke drawn inside the canvas edge
        outline = (
            f"\n<rect x='{frame / 2:.2f}' y='{frame / 2:.2f}' width='{W - frame:.2f}' height='{H - frame:.2f}' "
            f"fill='none' stroke='#000' stroke-width='{frame:.2f}'/>"
        )
    pages = []
    for _, _, place, _ in layouts:
        body = []
        for i, r, cx, cy, baseline in place:
            inner, x0, y0, w, h = get(i, r)
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
