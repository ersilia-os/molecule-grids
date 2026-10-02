"""Measure rendered print geometry and compare it with ChemDraw ACS Document 1996."""

import math
import re
from collections import defaultdict

from rdkit import Chem
from rdkit.Chem import rdCoordGen
from rdkit.Chem.Draw import rdMolDraw2D

from molecule_grids import FORMATS, draw_grid
from molecule_grids.depict import acs
from molecule_grids.depict.draw import condense, draw_one, stereo_hs

ACS = acs.ACS


def nums(d):
    return [float(t) for t in re.findall(r"-?\d+\.?\d*", d)]


def med(z):
    return sorted(z)[len(z) // 2]


def drawn(smi):
    m = stereo_hs(condense(Chem.MolFromSmiles(smi)))
    rdCoordGen.AddCoords(m)
    return m, draw_one(m)[0]


def drawn_bond(svg=None):
    """Longest plain single-bond stroke (n-butane by default)."""
    if svg is None:
        _, svg = drawn("CCCC")
    pat = r"<path class='bond-\d+ atom-\d+ atom-\d+' d='([^']+)'"
    return max(math.dist(v[:2], v[2:]) for d in re.findall(pat, svg) for v in [nums(d)] if len(v) == 4)


rows = []


def row(name, got, want):
    rows.append((name, got, want, 100 * (got - want) / want))


# artemisinin: lines, wedges, hashes
m, svg = drawn("C[C@@H]1CC[C@H]2[C@H](C(=O)O[C@H]3[C@@]24[C@H]1CC[C@](O3)(OO4)C)C")
to_pt = ACS["bond"] / drawn_bond()
row("line width", max(float(w) for w in re.findall(r"stroke-width:([\d.]+)px", svg)) * to_pt, ACS["line"])
bases = []
for d in re.findall(r"d='([^']+)' style='fill:#000000", svg):
    P = list(zip(nums(d)[0::2], nums(d)[1::2]))[:3]
    if len(P) == 3:
        bases.append(min(math.dist(P[i], P[(i + 1) % 3]) for i in range(3)))
row("bold (wedge) width", med(bases) * to_pt, ACS["bold"])
pm = rdMolDraw2D.PrepareMolForDrawing(m, addChiralHs=False)
sp = []
for i in [b.GetIdx() for b in pm.GetBonds() if b.GetBondDir() == Chem.BondDir.BEGINDASH]:
    pat = rf"<path class='bond-{i} [^']*' d='M ([\d.-]+),([\d.-]+) L ([\d.-]+),([\d.-]+)'"
    mids = [((a + c) / 2, (b + d) / 2) for a, b, c, d in (map(float, g) for g in re.findall(pat, svg))]
    sp += [math.dist(mids[k], mids[k + 1]) for k in range(len(mids) - 1)]
row("hash spacing", med(sp) * to_pt, ACS["hash"])

# benzene: double-bond spacing
m, svg = drawn("c1ccccc1")
lines = defaultdict(list)
for cls, d in re.findall(r"<path class='(bond-\d+)[^']*' d='([^']+)'", svg):
    v = nums(d)
    if len(v) == 4:
        lines[cls].append(((v[0] + v[2]) / 2, (v[1] + v[3]) / 2))
dbl = [math.dist(*p) for p in lines.values() if len(p) == 2]
row("double-bond spacing", med(dbl) * ACS["bond"] / drawn_bond(), ACS["double_spacing"] * ACS["bond"])

# labels: cap height -> implied font size, and margin
cap, gap = acs.probe(acs.BASE_FONT_SIZE, acs.LABEL_PADDING)
row("atom label size", cap / acs.cap_height_em() * ACS["bond"], ACS["label"])
row("label margin", gap * ACS["bond"], ACS["margin"])

# page: at each size, labels follow stylia and bonds keep ACS proportions (bond = 1.44 x label)
for size in ("small", "medium", "large"):
    label = FORMATS["print"].label_pt[size]
    page = draw_grid(["CCCC"], ["Butane"], format="print", size=size).svg
    wpt = float(re.search(r"width='([\d.]+)pt'", page).group(1))
    vbw = float(re.search(r"viewBox='0 0 ([\d.]+)", page).group(1))
    row(f"{size}: bond (placed at 100 %)", drawn_bond(page) * wpt / vbw, label * ACS["bond"] / ACS["label"])
    row(f"{size}: caption size", float(re.search(r"font-size='([\d.]+)'", page).group(1)) * wpt / vbw, label)

print(f"{'metric':32s}{'rendered':>10s}{'expected':>10s}{'diff':>8s}")
for n, g, w, e in rows:
    print(f"{n:32s}{g:9.2f}pt{w:9.2f}pt{e:+7.1f}%")
