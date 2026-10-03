"""Read molecules from text and files."""

import csv
import re
from pathlib import Path

from rdkit import Chem, RDLogger

from molecule_grids.utils.logging import logger

# SMILES, optional CXSMILES block "|...|", then an optional free-text name.
LINE = re.compile(r"^(\S+(?:\s+\|[^|]*\|)?)(?:\s+(.*))?$")


def parse_text(text):
    """Parse one molecule per line: SMILES (or CXSMILES), then an optional name.

    Blank lines and lines starting with ``#`` are ignored.

    Returns
    -------
    tuple
        ``(smiles, names, errors)``, where ``errors`` lists ``(line_no, line, reason)``.
    """
    smiles, names, errors = [], [], []
    RDLogger.DisableLog("rdApp.*")
    try:
        for no, raw in enumerate(text.splitlines(), 1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            m = LINE.match(line)
            smi, name = m.group(1), (m.group(2) or "").strip()
            if Chem.MolFromSmiles(smi) is None:
                errors.append((no, line, "not a valid SMILES"))
                continue
            smiles.append(smi)
            names.append(name)
    finally:
        RDLogger.EnableLog("rdApp.*")
    return smiles, names, errors


def _pick(header, candidates):
    lower = {h.lower(): h for h in header}
    return next((lower[c] for c in candidates if c in lower), None)


def read_molecules(path):
    """Read SMILES and names from a file.

    ``.csv``/``.tsv`` files need a ``smiles`` column; a ``name`` (or ``id``) column is
    optional. Any other file is read as one ``SMILES [name]`` per line. Invalid SMILES
    are skipped with a warning.

    Parameters
    ----------
    path : str or Path
        Input file.

    Returns
    -------
    tuple of list
        ``(smiles, names)``.
    """
    path = Path(path)
    if path.suffix.lower() in (".csv", ".tsv"):
        with path.open(newline="") as f:
            rows = list(csv.DictReader(f, delimiter="\t" if path.suffix.lower() == ".tsv" else ","))
        header = list(rows[0]) if rows else []
        scol = _pick(header, ("smiles", "canonical_smiles", "input", "smi"))
        if scol is None:
            raise ValueError(f"{path} has no 'smiles' column (columns: {header})")
        ncol = _pick(header, ("name", "names", "id", "identifier", "key"))
        text = "\n".join(f"{r[scol]} {r[ncol] if ncol else ''}".strip() for r in rows)
    else:
        text = path.read_text()
    smiles, names, errors = parse_text(text)
    for no, line, reason in errors:
        logger.warning(f"Skipping line {no} of {path.name}: {line[:60]!r} is {reason}")
    logger.info(f"Read {len(smiles)} molecules from {path}")
    return smiles, names
