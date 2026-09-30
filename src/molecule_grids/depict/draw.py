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


def text_width(text, size, bold=False):
    """Caption width in drawing units, measured with the Arial-metric font."""
    from PIL import ImageFont

    f = ImageFont.truetype(acs.ARIAL_BOLD if bold else acs.ARIAL, 100)
    return f.getlength(text) * size / 100


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
