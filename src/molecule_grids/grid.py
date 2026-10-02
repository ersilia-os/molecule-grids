"""Molecule grids: the main entry point of the package."""

import dataclasses
import html
import math
import re
from pathlib import Path

from rdkit import Chem, RDLogger

from molecule_grids import sizing
from molecule_grids.depict import acs, computational
from molecule_grids.depict.draw import STYLES, draw_one, ink_mask, rotate, text_width
from molecule_grids.layout import pack
from molecule_grids.layout.scaffold import arrange, prepare
from molecule_grids.utils.logging import logger

# A caption that is only a compound number (1, 2, 12, 3a, 10b ...) is set in bold (ACS style).
COMPOUND_NUMBER = re.compile(r"^\d+[a-z]{0,2}$")

MAX_MOLECULES = 100  # per page
TYPICAL = 0.9  # columns are sized so this share of molecules fits one cell; wider ones span
FRAME_COLOR = {"print": "#000000", "slide": "#50285A"}  # black for journals, Ersilia plum on slides
SQUEEZE = {False: None, None: None, True: "rows", "grid": "grid", "rows": "rows", "free": "free"}


class TooManyMolecules(ValueError):
    """The molecules need more than one page, so there is no single SVG to return."""


class Grid:
    """A drawn grid of molecules, possibly split over several pages.

    Build it with :func:`draw_grid`. In Jupyter, a grid displays itself inline.

    Attributes
    ----------
    format : str
        ``"print"`` or ``"slide"``.
    style : str
        ``"medicinal"`` or ``"computational"``.
    width : float
        Fraction of the full width used (the one picked, with ``width="auto"``).
    names : list of str
        Caption of each molecule, in input order.
    order : list of int
        Input indices in reading order (left to right, top to bottom, page by page).
    columns : int or None
        Number of columns, set by the width (when squeezing, possibly fewer; None for
        ``squeeze="free"``).
    capacity : int
        Most molecules on one page: up to ``MAX_MOLECULES`` (100), fewer when the page
        height is reached (247 mm in print, 186 mm on slides).
    bond_pt : float
        Bond length in pt when the figure is placed at 100 %.
    width_mm, height_mm : float
        Figure size in millimetres (height of the first page).
    pages : list of str
        One SVG document per figure. More than one only when the molecules exceed ``capacity``.
    """

    def __init__(self, fmt, style, width, names, order, geometry, columns, capacity, pages):
        self.format = fmt.name
        self.style = style
        self.width = width
        self._fmt = fmt
        self.names = names
        self.order = order
        self.columns = columns
        self.capacity = capacity
        self.bond_pt = geometry.bond_pt
        self.width_mm, self.height_mm = (
            float(re.search(rf"{dim}='([\d.]+)pt'", pages[0]).group(1)) * sizing.MM_PER_PT
            for dim in ("width", "height")
        )
        self.pages = pages

    def __repr__(self):
        return (
            f"Grid(format={self.format!r}, style={self.style!r}, molecules={len(self.names)}, columns={self.columns}, "
            f"pages={len(self.pages)}, size={self.width_mm:.0f}x{self.height_mm:.0f} mm, "
            f"bond={self.bond_pt:.1f} pt)"
        )

    def _repr_svg_(self):
        return self.pages[0]

    @property
    def svg(self):
        """The figure as an SVG string.

        Raises
        ------
        TooManyMolecules
            If the molecules need several pages; use ``pages`` or :meth:`save` instead.
        """
        if len(self.pages) > 1:
            raise TooManyMolecules(
                f"{len(self.names)} molecules need {len(self.pages)} pages ({self.capacity} fit on one). "
                "Use `.pages`, or `.save()`, which writes one file per page."
            )
        return self.pages[0]

    def save(self, path):
        """Save the figure; the file type follows the extension.

        Parameters
        ----------
        path : str or Path
            ``.svg``, ``.png`` or ``.pdf``. With several pages, one file per page is written,
            numbered ``name_1.svg``, ``name_2.svg``, ...

        Returns
        -------
        list of Path
            The written files, one per page.
        """
        from molecule_grids.io import export

        path = Path(path)
        ext = path.suffix.lower()
        write = {
            ".svg": lambda p, svg: p.write_text(svg),
            ".png": lambda p, svg: p.write_bytes(export.to_png(svg, self._fmt.png_dpi)),
            ".pdf": lambda p, svg: p.write_bytes(export.to_pdf(svg)),
        }
        if ext not in write:
            raise ValueError(f"unsupported file type {ext!r}: use .svg, .png or .pdf")
        paths = (
            [path]
            if len(self.pages) == 1
            else [path.with_stem(f"{path.stem}_{k}") for k in range(1, len(self.pages) + 1)]
        )
        for p, svg in zip(paths, self.pages):
            write[ext](p, svg)
            logger.success(f"Saved {p}")
        return paths


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
    style="medicinal",
    width="auto",
    max_height=1.0,
    size="medium",
    number=False,
    group=False,
    frame=True,
    squeeze="rows",
    fit=True,
):
    """Draw molecules as a grid, squeezed to fill the figure (or in square cells).

    The figure width is ``width`` times the format's full width (180 mm for print,
    330 mm for slides), as in stylia. Molecules are drawn at a fixed size, with atom labels
    and captions at stylia's font sizes; the width holds as many columns
    as fit at that size.

    Parameters
    ----------
    molecules : sequence of str or rdkit.Chem.Mol
        SMILES (CXSMILES extensions such as ``|&1:1|`` are kept) or RDKit molecules.
    names : sequence of str, optional
        One caption per molecule.
    format : {"slide", "print"}
        The page: its full width (330 or 180 mm), maximum height, standard bond length and
        frame colour.
    style : {"medicinal", "computational"}
        How molecules are drawn: ``"medicinal"`` (medicinal chemist), ChemDraw ACS Document
        1996 style in black and white; ``"computational"``, RDKit colours.
    width : float or "auto"
        Fraction of the format's full width, in (0, 1]. Together with ``size`` it sets the
        number of columns; rows are added as needed. ``"auto"`` (default) picks the fraction
        whose figure is closest to a pleasing shape: 16:9 on slides, 3:2 in print.
    max_height : float
        Tallest figure, as a fraction of the page height (247 mm in print, 186 mm on slides),
        in (0, 1]. Molecules that do not fit go to the next page; nothing is cut.
    size : {"small", "medium", "large"}
        Molecule size, set by the atom label and caption size, which follows stylia: print
        5, 6 or 8 pt; slides 8, 10 or 13 pt. Bonds are 1.44 x that (ChemDraw ACS proportions;
        print medium is 8.6 pt). If the largest molecule would not fit the width even in one
        column, it is shrunk to fit (see ``Grid.bond_pt``).
    number : bool
        Replace captions with bold compound numbers 1, 2, 3... in reading order.
    group : bool
        Keep analogues that share a scaffold side by side, in every layout. They share one
        orientation regardless.
    frame : bool
        Draw a thin outline around the whole figure, in the format's bond line width: black
        in print, Ersilia plum (#50285A) on slides. On by default.
    squeeze : {"rows", "grid", "free"} or bool
        Reorder (and turn by quarter turns) molecules to fill the figure at the same bond
        length and width; only the height shrinks. Analogues sharing a scaffold turn
        together. ``"rows"`` (default, or True): rows flow like justified text. ``"grid"``:
        aligned columns fitted to their largest molecule. ``"free"``: molecules interlock by
        their outlines. False: square cells in input order.
    fit : bool
        Make the figure only as wide as its content needs, up to ``width`` (default). False:
        exactly ``width``, with the content centred.

    Returns
    -------
    Grid

    Raises
    ------
    ValueError
        On invalid SMILES, a width or height out of range, or an unknown size or style.
    """
    fmt = sizing.get_format(format)
    if width != "auto" and not (isinstance(width, (int, float)) and 0 < width <= 1):
        raise ValueError(f"width is a fraction of the full {fmt.name} width and must be in (0, 1], not {width}")
    if not 0 < max_height <= 1:
        raise ValueError(f"max_height is a fraction of the page height and must be in (0, 1], not {max_height}")
    if style not in STYLES:
        raise ValueError(f"style must be one of {list(STYLES)}, not {style!r}")
    if size not in sizing.SIZES:
        raise ValueError(f"size must be one of {list(sizing.SIZES)}, not {size!r}")
    bond_pt = fmt.label_pt[size] / acs.LABEL_RATIO  # ACS: label 10 pt on a 14.4 pt bond
    molecules = list(molecules)
    if not molecules:
        raise ValueError("no molecules to draw")
    if names is None:
        names = [m.GetProp("_Name") if isinstance(m, Chem.Mol) and m.HasProp("_Name") else "" for m in molecules]
    names = [str(k + 1) for k in range(len(molecules))] if number else [str(n) for n in names]
    if len(names) != len(molecules):
        raise ValueError(f"got {len(names)} names for {len(molecules)} molecules")

    mols, buckets = prepare(_to_mols(molecules))
    drawn = [draw_one(m, style) for m in mols]
    bold = [bool(COMPOUND_NUMBER.match(n)) for n in names]
    caption = acs.LABEL_RATIO * acs.BOND  # captions match the atom labels, in both styles
    capw = [text_width(n, caption, b) if n else 0.0 for n, b in zip(names, bold)]
    wide = [max(d[3], c) for d, c in zip(drawn, capw)]  # width each molecule needs
    sides = sorted(max(w, d[4]) for w, d in zip(wide, drawn))
    typical = sides[math.ceil(TYPICAL * len(sides)) - 1]  # cell side that fits most molecules
    if width == "auto":
        width = _auto_width(fmt, bond_pt, wide, [d[4] for d in drawn], typical, caption, max_height, fit)
    geo = dataclasses.replace(sizing.solve(fmt, width, bond_pt, max(wide), typical, caption), fit=fit)
    if geo.bond_pt < bond_pt - 1e-6:
        logger.warning(
            f"The largest molecule does not fit width={width} at size={size!r}; "
            f"drawn at {geo.bond_pt:.1f} pt bonds instead of {bond_pt:.1f}"
        )
    columns = geo.columns
    max_h = max_height * fmt.max_height_in * sizing.PT_PER_IN / geo.scale  # tallest page, in drawing units
    full_rows = max(1, int((max_h - geo.gap) // (geo.cell + geo.strip + geo.gap)))  # typical rows per page
    variants = {(i, 0): d for i, d in enumerate(drawn)}

    def get(i, r):
        """Drawing of molecule ``i`` turned ``r`` quarter turns (drawn on first use)."""
        if (i, r) not in variants:
            variants[i, r] = draw_one(rotate(mols[i], r), style)
        return variants[i, r]

    mode = SQUEEZE.get(squeeze, "bad") if isinstance(squeeze, (bool, str, type(None))) else "bad"
    if mode == "bad":
        raise ValueError(f"squeeze must be False, True, 'grid', 'rows' or 'free', not {squeeze!r}")
    span = [geo.span(w) for w in wide]
    cells = [s * geo.cell + (s - 1) * geo.gap for s in span]  # square cell width of each molecule
    if mode:
        fam = {i: idxs[0] for key, idxs in buckets.items() if len(idxs) > 1 and key != "(acyclic)" for i in idxs}
        fam = [fam.get(i, i) for i in range(len(mols))]  # family id: analogues turn (and group) together
        together = group
        squeezer = {"grid": _grid_layout, "rows": _rows_layout, "free": _free_layout}[mode]
        layouts = _paginate(len(mols), max_h, lambda idx: squeezer(geo, get, idx, columns, capw, fam, together))
    else:
        rows = arrange(buckets, len(mols), columns, group, span)
        layouts = _square_pages(geo, drawn, span, rows, max_h)

    def slot(i, r):
        """Width a molecule takes in its row: its cell, or (squeezed) its ink or caption."""
        return max(get(i, r)[3], capw[i]) if mode else cells[i]

    layouts = _compact(geo, layouts, slot)
    tallest = max(lay[1] for lay in layouts)
    if tallest > max_h + 1e-6:
        logger.warning(
            f"max_height={max_height} is less than the tallest molecule needs; pages are up to "
            f"{tallest * geo.scale * sizing.MM_PER_PT:.0f} mm tall (nothing is cut)"
        )
    order = [i for _, _, place, _ in layouts for i, *_ in place]
    if number:  # numbers follow reading order
        names = [None] * len(order)
        for k, i in enumerate(order):
            names[i] = str(k + 1)
    line = acs.BOND * (acs.LINE_RATIO if style == "medicinal" else computational.LINE_RATIO)  # as the bonds
    pages = _render(geo, layouts, get, names, bold, (line, FRAME_COLOR[fmt.name]) if frame else None)
    used = None if any(c is None for *_, c in layouts) else max(c for *_, c in layouts)
    if len(layouts) > 1:  # what the first page holds
        capacity = len(layouts[0][2])
    else:  # estimate for typical molecules; a squeezed page may hold more
        capacity = min(MAX_MOLECULES, columns * full_rows)
        capacity = max(capacity, len(mols)) if mode else capacity
    return Grid(fmt, style, width, names, order, geo, used, capacity, pages)


def _paginate(n, max_h, build):
    """Split molecules 0..n-1, in input order, into squeezed pages no taller than ``max_h``.

    Each page takes about as many molecules as fit (at most ``MAX_MOLECULES``): the count is
    estimated from the height of a trial packing and corrected in a few steps.
    """
    layouts, start = [], 0
    while start < n:
        k = min(n - start, MAX_MOLECULES)
        lay = build(list(range(start, start + k)))
        if lay[1] > max_h:
            fit = None  # (k, layout) of the largest packing known to fit
            for _ in range(6):
                k = max(1, min(k - 1, int(k * max_h / lay[1]))) if lay[1] > max_h else k
                lay = build(list(range(start, start + k)))
                if lay[1] <= max_h:
                    fit = (k, lay)
                    more = min(n - start, MAX_MOLECULES, int(k * max_h / lay[1]))
                    if more <= k:
                        break
                    trial = build(list(range(start, start + more)))
                    if trial[1] <= max_h:
                        fit = (more, trial)
                    break
                if k == 1:
                    break
            k, lay = fit or (k, lay)
        layouts.append(lay)
        start += k
    return layouts


# Layouts return (W, H, [(i, r, cx, cy, baseline)], columns): figure size, and per molecule its
# quarter turns, ink centre and caption baseline, in reading order.


def _auto_width(fmt, bond_pt, wide, tall, typical, caption, max_height, fit):
    """Width fraction (25 to 100 %, in 5 % steps) whose figure is closest to ``fmt.aspect``.

    Each candidate is estimated by filling rows with the molecules in order, which is quick
    and close to the real layouts; ties go to the wider figure.
    """
    best = None
    for k in range(20, 4, -1):
        width = k / 20
        geo = sizing.solve(fmt, width, bond_pt, max(wide), typical, caption)
        gap, strip = geo.gap, geo.strip
        max_h = max_height * fmt.max_height_in * sizing.PT_PER_IN / geo.scale
        rows, row, x = [], [], gap
        for w, h in zip(wide, tall):
            if row and x + w + gap > geo.width:
                rows.append((x, row))
                row, x = [], gap
            row.append(h)
            x += w + gap
        rows.append((x, row))
        W = max(x for x, _ in rows) if fit else geo.width
        H = min(max_h, gap + sum(max(r) + strip + gap for _, r in rows))
        score = abs(math.log(W / H / fmt.aspect))
        if best is None or score < best[0] - 1e-9:
            best = (score, width)
    return best[1]


def _square_pages(geo, drawn, span, rows, max_h):
    """Square cells on aligned columns, split into pages no taller than ``max_h``.

    Each cell is one column wide (``span`` columns for wide molecules) and each row is as
    tall as its tallest molecule.
    """
    W, gap, pitch = geo.width, geo.gap, geo.pitch
    if geo.fit:  # cells just large enough; the figure as wide as the widest row
        pitch = geo.cell + gap
        W = gap + max(sum(span[i] for i in r) for r in rows) * pitch
    layouts, place, y, count, used = [], [], gap, 0, 0
    for r in rows:
        rh = max(drawn[i][4] for i in r)  # each row is as tall as its tallest molecule
        if place and (y + rh + geo.strip + gap > max_h or count + len(r) > MAX_MOLECULES):
            layouts.append((W, y, place, used))
            place, y, count, used = [], gap, 0, 0
        x = (W - sum(span[i] for i in r) * pitch - gap) / 2 + gap  # centre partial rows
        for i in r:
            w = span[i] * pitch - gap  # cell width
            cy = y + rh / 2
            place.append((i, 0, x + w / 2, cy, cy + rh / 2 + geo.caption_base))  # shared caption baseline
            x += span[i] * pitch
        y += rh + geo.strip + gap
        count += len(r)
        used = max(used, sum(span[i] for i in r))  # columns actually used
    layouts.append((W, y, place, used))
    return layouts


def _compact(geo, layouts, slot):
    """Pack single-row pages tightly and give all pages one width.

    A page holding a single row is re-packed without spread-out cells. With ``geo.fit``,
    pages are as wide as the widest page's content needs; otherwise exactly ``geo.width``.
    Each page's content is centred in that width.
    """
    gap, out = geo.gap, []
    for W, H, place, cols in layouts:
        if len({round(b, 1) for *_, b in place}) == 1:  # one row: pack it tightly
            row = sorted(place, key=lambda p: p[2])
            x, packed = gap, []
            for i, r, _, cy, b in row:
                packed.append((i, r, x + slot(i, r) / 2, cy, b))
                x += slot(i, r) + gap
            W, place, cols = x, packed, len(row)
        out.append((W, H, place, cols))
    common = max(W for W, *_ in out) if geo.fit else geo.width
    return [
        (common, H, [(i, r, cx + (common - W) / 2, cy, b) for i, r, cx, cy, b in place], cols)
        for W, H, place, cols in out
    ]


def _dims(get, idx, turns):
    return {(i, r): get(i, r)[3:5] for i in idx for r in turns}


def _grid_layout(geo, get, idx, columns, capw, group, together):
    """Tight aligned grid (see :func:`molecule_grids.layout.pack.grid_pack`)."""
    gap, strip = geo.gap, geo.strip
    dims = _dims(get, idx, (0, 1))
    rows = pack.grid_pack(idx, 2 * columns, geo.width, dims, capw, group, gap, strip, together)

    def box(i, r):
        return max(dims[i, r][0], capw[i]), dims[i, r][1]

    ncol = len(rows[0])
    colw = [max(box(i, r)[0] for row in rows for i, r in [row[c]] if i is not None) for c in range(ncol)]
    spare = 0.0 if geo.fit else max(0.0, geo.width - (gap + sum(w + gap for w in colw)))
    colw = [w + spare / ncol for w in colw]  # not fitting: the spare width is shared by the columns
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


def _rows_layout(geo, get, idx, columns, capw, group, together):
    """Rows like justified text (see :func:`molecule_grids.layout.pack.rows_pack`)."""
    gap, strip = geo.gap, geo.strip
    dims = _dims(get, idx, (0, 1))
    rows = pack.rows_pack(idx, geo.width, dims, capw, group, gap, strip, together)

    def box(i, r):
        return max(dims[i, r][0], capw[i]), dims[i, r][1]

    natural = [gap + sum(box(i, r)[0] + gap for i, r in row) for row in rows]
    W = max(natural) if geo.fit else geo.width  # rows are justified to this width
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


def _free_layout(geo, get, idx, columns, capw, group, together):
    """Interlocking outlines (see :func:`molecule_grids.layout.pack.free_pack`)."""
    gap, bond = geo.gap, acs.BOND
    res, clearance = bond / 2, bond / 2  # mask pixel and half the minimum distance between inks
    caption = {i: (capw[i], geo.caption_base - 0.75 * geo.caption, geo.caption_base + 0.3 * geo.caption) for i in idx}
    masks, offset = {}, {}
    for i in idx:
        for r in range(4):
            masks[i, r], *offset[i, r] = ink_mask(get(i, r), res, clearance, caption[i] if capw[i] else None)
    placed = pack.free_pack(idx, geo.width - 2 * gap, masks, res, group, together)
    items = []
    for i, r, x, y in placed:
        _, _, _, w, h = get(i, r)
        left, top = x - offset[i, r][0], y - offset[i, r][1]  # ink top-left
        items.append((i, r, left + w / 2, top + h / 2, top + h + geo.caption_base, left, top, w))
    # Fit: the width is the content's; otherwise centre it across the full width.
    x0 = min(min(left, cx - capw[i] / 2) for i, _, cx, _, _, left, _, _ in items)
    y0 = min(top for *_, top, _ in items)
    x1 = max(max(left + w, cx + capw[i] / 2) for i, _, cx, _, _, left, _, w in items)
    y1 = max(b + 0.3 * geo.caption if capw[i] else cy + get(i, r)[4] / 2 for i, r, _, cy, b, *_ in items)
    W = x1 - x0 + 2 * gap if geo.fit else geo.width
    dx, dy = (W - (x1 - x0)) / 2 - x0, gap - y0
    place = [(i, r, cx + dx, cy + dy, b + dy) for i, r, cx, cy, b, *_ in items]
    place.sort(key=lambda p: (round(p[3] / (3 * bond)), p[2]))  # reading order: bands, then left to right
    return W, y1 - y0 + 2 * gap, place, None


def _render(geo, layouts, get, names, bold, frame=None):
    """Write one SVG per layout (page); ``frame`` is ``(line width, colour)`` or None."""
    pages = []
    for W, H, place, _ in layouts:
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
        if frame:  # stroke drawn inside the canvas edge
            lw, color = frame
            body.append(
                f"<rect x='{lw / 2:.2f}' y='{lw / 2:.2f}' width='{W - lw:.2f}' height='{H - lw:.2f}' "
                f"fill='none' stroke='{color}' stroke-width='{lw:.2f}'/>"
            )
        pages.append(
            "<?xml version='1.0' encoding='utf-8'?>\n"
            f"<svg xmlns='http://www.w3.org/2000/svg' width='{W * geo.scale:.2f}pt' "
            f"height='{H * geo.scale:.2f}pt' viewBox='0 0 {W:.2f} {H:.2f}'>\n"
            "<rect width='100%' height='100%' fill='#FFFFFF'/>\n" + "\n".join(body) + "\n</svg>\n"
        )
    return pages
