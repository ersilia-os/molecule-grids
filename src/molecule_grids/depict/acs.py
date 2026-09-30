"""ChemDraw "ACS Document 1996" drawing style, calibrated against RDKit.

Every drawing parameter is derived from the ACS style sheet (see ``ACS``). RDKit's font
size, label padding and double-bond offset do not map linearly to points, so they are
solved numerically by :func:`calibrate` on first use.
"""

import functools
import math
import re
from importlib import resources

from rdkit import Chem
from rdkit.Chem import rdCoordGen
from rdkit.Chem.Draw import rdMolDraw2D

_FONTS = resources.files("molecule_grids") / "fonts"
ARIAL = str(_FONTS / "LiberationSans-Regular.ttf")  # Arial-metric
ARIAL_BOLD = str(_FONTS / "LiberationSans-Bold.ttf")
FONT_STACK = "Arial,'Liberation Sans',sans-serif"

# ChemDraw "ACS Document 1996" style sheet, in points (Revvity/PerkinElmer support docs).
ACS = {
    "bond": 14.4,
    "line": 0.6,
    "bold": 2.0,
    "margin": 1.6,
    "hash": 2.5,
    "label": 10.0,
    "caption": 10.0,
    "double_spacing": 0.18,
}

BOND = 26  # internal drawing units per bond
LINE_RATIO = ACS["line"] / ACS["bond"]
BOLD_RATIO = ACS["bold"] / ACS["bond"]
HASH_RATIO = ACS["hash"] / ACS["bond"]
MARGIN_RATIO = ACS["margin"] / ACS["bond"]
LABEL_RATIO = ACS["label"] / ACS["bond"]  # font em size / bond

# Solved by calibrate(); the values here are only starting points.
BASE_FONT_SIZE = 0.65
LABEL_PADDING = 0.08
DOUBLE_OFFSET = 0.18


def acs_options(opts):
    """Apply ACS Document 1996 drawing options in place."""
    opts.fontFile = ARIAL
    opts.fixedBondLength = 24
    opts.baseFontSize = BASE_FONT_SIZE  # calibrated: label em = 10 pt on 14.4 pt bond
    opts.bondLineWidth = BOND * LINE_RATIO  # line 0.6 pt @ 14.4 pt bond
    opts.multipleBondOffset = DOUBLE_OFFSET  # calibrated: spacing = 18 % of bond
    opts.additionalAtomLabelPadding = LABEL_PADDING  # calibrated: gap = 1.6 pt margin
    opts.useBWAtomPalette()
    opts.padding = 0.06
    opts.clearBackground = True


def _numbers(d):
    return [float(t) for t in re.findall(r"-?\d+\.?\d*", d)]


def _median(z):
    return sorted(z)[len(z) // 2]


def probe(base, pad):
    """Draw a probe molecule; return (cap height / bond, label gap / bond)."""
    global BASE_FONT_SIZE, LABEL_PADDING
    m = Chem.MolFromSmiles("Cc1cc(NS(=O)(=O)c2ccc(N)cc2)no1")
    rdCoordGen.AddCoords(m)
    d = rdMolDraw2D.MolDraw2DSVG(300, 300)
    o = d.drawOptions()
    keep = BASE_FONT_SIZE, LABEL_PADDING
    BASE_FONT_SIZE, LABEL_PADDING = base, pad
    acs_options(o)
    o.fixedBondLength = BOND
    BASE_FONT_SIZE, LABEL_PADDING = keep
    d.DrawMolecule(m)
    d.FinishDrawing()
    svg = d.GetDrawingText()
    glyph, bonds = {}, []
    for cls, dd in re.findall(r"<path class='([^']+)' d='([^']+)'", svg):
        v = _numbers(dd)
        pts = list(zip(v[0::2], v[1::2]))
        if cls.startswith("atom-"):
            glyph.setdefault(int(cls.split("-")[1]), []).extend(pts)
        elif cls.startswith("bond-") and len(pts) == 2:
            ab = [int(x.split("-")[1]) for x in cls.split() if x.startswith("atom-")]
            bonds.append((ab, pts))
    L = max(math.dist(*p) for ab, p in bonds if not set(ab) & set(glyph))
    caps = [
        max(y for _, y in glyph[a.GetIdx()]) - min(y for _, y in glyph[a.GetIdx()])
        for a in m.GetAtoms()
        if a.GetSymbol() in "NOF" and not a.GetTotalNumHs() and a.GetIdx() in glyph
    ]
    gaps = []
    for ab, pts in bonds:
        for a in ab:
            if a in glyph:
                xs = [x for x, _ in glyph[a]]
                ys = [y for _, y in glyph[a]]
                x0, x1, y0, y1 = min(xs), max(xs), min(ys), max(ys)
                e = min(pts, key=lambda q: math.dist(q, ((x0 + x1) / 2, (y0 + y1) / 2)))
                gaps.append(math.hypot(max(x0 - e[0], 0, e[0] - x1), max(y0 - e[1], 0, e[1] - y1)))
    return _median(caps) / L, _median(gaps) / L


def _probe_double(offset):
    """Double-bond line separation / bond length for a given RDKit multipleBondOffset."""
    global DOUBLE_OFFSET
    keep, DOUBLE_OFFSET = DOUBLE_OFFSET, offset
    m = Chem.MolFromSmiles("C=C")  # acyclic: symmetric lines
    rdCoordGen.AddCoords(m)
    d = rdMolDraw2D.MolDraw2DSVG(300, 300)
    o = d.drawOptions()
    acs_options(o)
    o.fixedBondLength = BOND
    DOUBLE_OFFSET = keep
    d.DrawMolecule(m)
    d.FinishDrawing()
    segs = [_numbers(dd) for dd in re.findall(r"<path class='bond-0[^']*' d='([^']+)'", d.GetDrawingText())]
    segs = [v for v in segs if len(v) == 4]
    (a, b), (c, e) = [((v[0] + v[2]) / 2, (v[1] + v[3]) / 2) for v in segs[:2]]
    ref = rdMolDraw2D.MolDraw2DSVG(300, 300)
    o2 = ref.drawOptions()
    acs_options(o2)
    o2.fixedBondLength = BOND
    s1 = Chem.MolFromSmiles("CC")
    rdCoordGen.AddCoords(s1)
    ref.DrawMolecule(s1)
    ref.FinishDrawing()
    v = _numbers(re.search(r"<path class='bond-0[^']*' d='([^']+)'", ref.GetDrawingText()).group(1))
    return math.dist((a, b), (c, e)) / math.dist(v[:2], v[2:4])


def cap_height_em():
    """Cap height of the Arial-metric font, as a fraction of the em size."""
    from PIL import ImageFont

    b = ImageFont.truetype(ARIAL, 1000).getbbox("H")
    return (b[3] - b[1]) / 1000


def _bisect(f, target, lo, hi):
    """Solve f(x) = target for increasing f on [lo, hi]."""
    for _ in range(25):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if f(mid) < target else (lo, mid)
    return (lo + hi) / 2


@functools.cache
def calibrate():
    """Solve RDKit font size, label padding and double-bond offset against ACS 1996.

    Runs once per process (about a second) and updates the module-level settings.

    Returns
    -------
    tuple of float
        ``(BASE_FONT_SIZE, LABEL_PADDING, DOUBLE_OFFSET)``.
    """
    global BASE_FONT_SIZE, LABEL_PADDING, DOUBLE_OFFSET
    BASE_FONT_SIZE = _bisect(lambda x: probe(x, LABEL_PADDING)[0], LABEL_RATIO * cap_height_em(), 0.3, 1.2)
    LABEL_PADDING = _bisect(lambda x: probe(BASE_FONT_SIZE, x)[1], MARGIN_RATIO, 0.0, 0.5)
    DOUBLE_OFFSET = _bisect(_probe_double, ACS["double_spacing"], 0.05, 0.6)
    return BASE_FONT_SIZE, LABEL_PADDING, DOUBLE_OFFSET


def respace_hashes(inner, hashed_bonds):
    """Rebuild RDKit hashed wedges with ChemDraw spacing and line width.

    Keeps each wedge's outline (first and last hash line), regenerating the
    lines in between at HASH_RATIO x bond spacing and full line width.
    """
    pat = r"<path class='bond-{i} [^']*' d='M ([\d.-]+),([\d.-]+) L ([\d.-]+),([\d.-]+)' style='([^']*)' />"

    def mid(line):
        return (line[0] + line[2]) / 2, (line[1] + line[3]) / 2

    for i in hashed_bonds:
        found = list(re.finditer(pat.format(i=i), inner))
        if len(found) < 2:
            continue
        lines = [tuple(map(float, f.groups()[:4])) for f in found]
        lines.sort(key=lambda s: math.hypot(s[2] - s[0], s[3] - s[1]))  # narrow -> wide
        (ax, ay), (bx, by) = mid(lines[0]), mid(lines[-1])
        nx, ny = -(by - ay), bx - ax  # wedge-axis normal

        def orient(s, nx=nx, ny=ny):  # consistent sides
            mx, my = mid(s)
            return s if (s[0] - mx) * nx + (s[1] - my) * ny < 0 else (s[2], s[3], s[0], s[1])

        f0, f1 = orient(lines[0]), orient(lines[-1])
        span = math.hypot(bx - ax, by - ay)
        step = HASH_RATIO * BOND / max(span, 1e-6)  # exact ACS pitch, as a fraction
        ts, t = [], 1.0  # start at the wide end
        while t >= -1e-9:
            ts.append(max(t, 0.0))
            t -= step
        style = re.sub(r"stroke-width:[\d.]+px", f"stroke-width:{BOND * LINE_RATIO:.2f}px", found[0].group(5))
        cls = re.search(r"class='([^']+)'", found[0].group(0)).group(1)
        new = []
        for t in ts:
            p = [f0[j] + (f1[j] - f0[j]) * t for j in range(4)]
            new.append(f"<path class='{cls}' d='M {p[0]:.1f},{p[1]:.1f} L {p[2]:.1f},{p[3]:.1f}' style='{style}' />")
        first = found[0].group(0)
        for f in found[1:]:
            inner = inner.replace(f.group(0), "", 1)
        inner = inner.replace(first, "\n".join(new), 1)
    return inner


def set_wedge_width(inner):
    """Rescale RDKit solid wedges so the wide end is ACS bold width (2 pt @ 14.4 pt)."""
    target = BOLD_RATIO * BOND

    def fix(mt):
        v = _numbers(mt.group(2))
        P = list(zip(v[0::2], v[1::2]))
        if len(P) < 3:
            return mt.group(0)
        P = P[:3]
        k = min(range(3), key=lambda i: math.dist(P[i], P[(i + 1) % 3]))  # base = shortest side
        a, b, tip = P[k], P[(k + 1) % 3], P[(k + 2) % 3]
        mx, my = (a[0] + b[0]) / 2, (a[1] + b[1]) / 2
        f = target / max(math.dist(a, b), 1e-6)
        a = (mx + (a[0] - mx) * f, my + (a[1] - my) * f)
        b = (mx + (b[0] - mx) * f, my + (b[1] - my) * f)
        d = f"M {tip[0]:.1f},{tip[1]:.1f} L {a[0]:.1f},{a[1]:.1f} L {b[0]:.1f},{b[1]:.1f} Z"
        return f"{mt.group(1)}d='{d}'{mt.group(3)}"

    return re.sub(r"(<path class='bond-\d+[^']*' )d='([^']+)'( style='fill:#000000[^']*'[^>]*>)", fix, inner)
