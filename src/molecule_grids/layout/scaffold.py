"""Scaffold-aware 2D layout: analogues sharing a Bemis-Murcko scaffold get one orientation."""

from collections import defaultdict

from rdkit import Chem
from rdkit.Chem import rdCoordGen, rdMolAlign
from rdkit.Chem.Scaffolds import MurckoScaffold

from molecule_grids.depict.draw import condense, stereo_hs


def bucket_by_scaffold(mols):
    """Return ``{scaffold_smiles: [indices]}`` keyed by Bemis-Murcko scaffold."""
    buckets = defaultdict(list)
    for i, m in enumerate(mols):
        scaf = MurckoScaffold.GetScaffoldForMol(m)
        key = Chem.MolToSmiles(scaf) if scaf.GetNumAtoms() else "(acyclic)"
        buckets[key].append(i)
    return buckets


def layout(mols, buckets):
    """Set 2D coords in place, aligning multi-member buckets to their shared scaffold.

    Each molecule is laid out in full by CoordGen (keeping clean side-chain geometry),
    then rigidly rotated/translated so the shared scaffold lands in a common orientation.
    This avoids the distortion of a constrained re-layout of the non-scaffold atoms.
    """
    for key, idxs in buckets.items():
        if len(idxs) > 1 and key != "(acyclic)":
            ref = MurckoScaffold.GetScaffoldForMol(mols[idxs[0]])
            rdCoordGen.AddCoords(ref)  # one clean core layout
            for i in idxs:
                rdCoordGen.AddCoords(mols[i])  # full clean molecule layout
                match = mols[i].GetSubstructMatch(ref)
                if match:  # rigid orient only, no re-layout
                    amap = list(zip(match, range(ref.GetNumAtoms())))
                    rdMolAlign.AlignMol(mols[i], ref, atomMap=amap)
        else:
            for i in idxs:
                rdCoordGen.AddCoords(mols[i])  # singleton: free layout


def prepare(mols):
    """Condense groups, add stereo H and lay out; returns new mols with 2D coords."""
    buckets = bucket_by_scaffold(mols)
    mols = [stereo_hs(condense(m)) for m in mols]
    layout(mols, buckets)
    return mols, buckets


def arrange(buckets, n, per_row, group):
    """Assign molecule indices to grid rows.

    ``group=False`` keeps input order. ``group=True`` uses first-fit packing: families
    (largest first) each go into the first row with room for the whole family, so a family
    is never split unless it is wider than ``per_row``; singletons then fill the gaps.
    """
    if not group:
        order = list(range(n))
        return [order[k : k + per_row] for k in range(0, n, per_row)]
    families = sorted((b for b in buckets.values() if len(b) > 1), key=len, reverse=True)
    singles = [[i] for b in buckets.values() if len(b) == 1 for i in b]
    rows = []
    for fam in families + singles:
        for chunk in (fam[k : k + per_row] for k in range(0, len(fam), per_row)):
            row = next((r for r in rows if len(r) + len(chunk) <= per_row), None)
            if row is None:
                rows.append(list(chunk))
            else:
                row += chunk
    return rows
