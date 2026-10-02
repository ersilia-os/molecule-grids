"""Draw single molecules as SVG fragments with a true ink bounding box."""

import math
import re

from rdkit import Chem
from rdkit.Chem import rdAbbreviations
from rdkit.Chem.Draw import rdMolDraw2D

from molecule_grids.depict import acs
from molecule_grids.depict.slide import slide_options

# Groups drawn as condensed labels (journal convention). Labels auto-flip (CF3 -> F3C).
ABBREVIATE = ("CF3", "NO2")


def condense(mol, labels=ABBREVIATE):
    """Replace whitelisted groups with condensed labels (e.g. CF3)."""
    abbs = [a for a in rdAbbreviations.GetDefaultAbbreviations() if a.label in labels]
    return rdAbbreviations.CondenseMolAbbreviations(mol, abbs) if abbs else mol


def stereo_hs(mol):
    """Add the explicit stereo H atoms RDKit would draw (e.g. ring-fusion H) before layout.

    If added at draw time instead, they are placed after the 2D layout and can land
    inside crowded ring systems (bridged cages). Uses RDKit's own selection rule.
    """
    m = rdMolDraw2D.PrepareMolForDrawing(mol, kekulize=False, addChiralHs=True, wedgeBonds=False)
    m.RemoveAllConformers()  # coords come from CoordGen later
    return m


def chemdraw_stereo_labels(mol):
    """Drawing copy with enhanced-stereo groups shown in ChemDraw notation (&1, or1, abs).

    RDKit hardcodes 'and1'; here each group becomes a per-atom note instead and the
    groups are dropped from the drawing copy only. Wedges come from the chirality
    tags, so the depicted configuration is unchanged; the input mol keeps its groups.
    """
    groups = mol.GetStereoGroups()
    if not groups:
        return mol
    T = Chem.StereoGroupType
    counters, m = {T.STEREO_AND: 0, T.STEREO_OR: 0}, Chem.RWMol(mol)
    for g in groups:
        kind = g.GetGroupType()
        if kind == T.STEREO_ABSOLUTE:
            label = "abs"
        else:
            n = g.GetReadId() or (counters[kind] + 1)
            counters[kind] = max(counters[kind], n)
            label = f"&{n}" if kind == T.STEREO_AND else f"or{n}"
        for a in g.GetAtoms():
            m.GetAtomWithIdx(a.GetIdx()).SetProp("atomNote", label)
    m.SetStereoGroups([])
    return m.GetMol()


def rotate(mol, quarter_turns):
    """Copy of ``mol`` with its 2D coordinates turned by ``quarter_turns`` x 90 degrees.

    Wedges are recomputed from the chirality tags at draw time, so stereo is unchanged.
    """
    m = Chem.Mol(mol)
    conf = m.GetConformer()
    for _ in range(quarter_turns % 4):
        for k in range(m.GetNumAtoms()):
            p = conf.GetAtomPosition(k)
            conf.SetAtomPosition(k, (-p.y, p.x, p.z))
    return m


def text_width(text, size, bold=False):
    """Caption width in drawing units, measured with the Arial-metric font."""
    from PIL import ImageFont

    f = ImageFont.truetype(acs.ARIAL_BOLD if bold else acs.ARIAL, 100)
    return f.getlength(text) * size / 100


def ink_mask(drawn, res, clearance, caption=None):
    """Rasterise a drawn molecule (and its caption) into a dilated occupancy mask.

    Parameters
    ----------
    drawn : tuple
        Output of :func:`draw_one`.
    res : float
        Drawing units per pixel.
    clearance : float
        Free margin kept around the ink, in drawing units.
    caption : tuple, optional
        ``(width, top, bottom)`` of the caption box, the last two measured down from the
        ink's bottom edge.

    Returns
    -------
    tuple
        ``(mask, dx, dy)``: boolean array, and the mask origin relative to the ink's
        top-left corner, in drawing units.
    """
    import numpy as np
    from PIL import Image, ImageDraw, ImageFilter

    inner, x0, y0, w, h = drawn
    cx, capw, cap_bottom = x0 + w / 2, 0.0, 0.0
    if caption:
        capw, cap_top, cap_bottom = caption
    left = min(x0, cx - capw / 2) - clearance
    top = y0 - clearance
    right = max(x0 + w, cx + capw / 2) + clearance
    bottom = y0 + h + cap_bottom + clearance
    img = Image.new("L", (math.ceil((right - left) / res) + 1, math.ceil((bottom - top) / res) + 1), 0)
    pen = ImageDraw.Draw(img)

    def px(points):
        return [((x - left) / res, (y - top) / res) for x, y in points]

    for path in re.findall(r"<path [^>]*>", inner):  # bonds (stroked) and label glyphs (filled)
        v = [float(t) for t in re.findall(r"-?\d+\.?\d*", re.search(r" d='([^']+)'", path).group(1))]
        pts = px(list(zip(v[0::2], v[1::2])))
        if "fill:none" in path or len(pts) < 3:
            pen.line(pts, fill=255, width=1)
        else:
            pen.polygon(pts, fill=255)
    if caption:
        pen.rectangle(px([(cx - capw / 2, y0 + h + cap_top), (cx + capw / 2, y0 + h + cap_bottom)]), fill=255)
    k = max(1, round(clearance / res))
    img = img.filter(ImageFilter.MaxFilter(2 * k + 1))
    return np.asarray(img) > 0, left - x0, top - y0


def draw_one(mol, fmt="print"):
    """Draw one molecule with 2D coordinates at ``acs.BOND`` units per bond.

    Returns
    -------
    tuple
        ``(inner_svg, x0, y0, w, h)``, where ``(x0, y0, w, h)`` is the ink bounding box.
        RDKit pads asymmetrically, so centring must use the ink, not the canvas.
    """
    acs.calibrate()
    # Fixed canvas sized from the coordinates, with generous room. (RDKit's flexicanvas,
    # MolDraw2DSVG(-1, -1), does not honour fixedBondLength: bonds come out ~27 % short.)
    conf = mol.GetConformer()
    xy = [(conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y) for i in range(mol.GetNumAtoms())]
    bl = sorted(math.dist(xy[b.GetBeginAtomIdx()], xy[b.GetEndAtomIdx()]) for b in mol.GetBonds()) or [1.0]
    unit = acs.BOND / bl[len(bl) // 2]
    cw = (max(x for x, _ in xy) - min(x for x, _ in xy)) * unit + 8 * acs.BOND
    ch = (max(y for _, y in xy) - min(y for _, y in xy)) * unit + 8 * acs.BOND
    d = rdMolDraw2D.MolDraw2DSVG(int(cw), int(ch))
    o = d.drawOptions()
    if fmt == "print":
        acs.acs_options(o)
    else:
        slide_options(o)
    o.fixedBondLength = acs.BOND
    o.clearBackground = False
    pm = rdMolDraw2D.PrepareMolForDrawing(chemdraw_stereo_labels(mol), addChiralHs=False)
    hashed = [b.GetIdx() for b in pm.GetBonds() if b.GetBondDir() == Chem.BondDir.BEGINDASH]
    o.prepareMolsBeforeDrawing = False  # wedging done above, identically
    d.DrawMolecule(pm)
    d.FinishDrawing()
    svg = d.GetDrawingText()
    inner = svg[svg.index("<!-- END OF HEADER -->") : svg.rindex("</svg>")]
    if fmt == "print":  # ACS hash spacing and wedge width
        inner = acs.respace_hashes(inner, hashed)
        inner = acs.set_wedge_width(inner)
    xs, ys = [], []
    for dstr in re.findall(r" d='([^']+)'", inner):  # RDKit paths use absolute coords
        v = [float(t) for t in re.findall(r"-?\d+\.?\d*", dstr)]
        xs += v[0::2]
        ys += v[1::2]
    x0, y0 = min(xs), min(ys)
    return inner, x0, y0, max(xs) - x0, max(ys) - y0
