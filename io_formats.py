"""
bbsplit.io_formats
==================
Readers that turn an input file into a list of (mol_id, smiles) tuples.

Supported formats:
    .csv         - delimited; caller picks the SMILES (and optional ID) column.
    .smi / .txt  - one SMILES per line; an optional second whitespace- or
                   tab-separated token is treated as the molecule name/ID.
    .sdf / .mol  - MDL molfile / SD file; molecules are read with RDKit and
                   converted to canonical SMILES. The '_Name' property (the
                   molfile title line) is used as the ID when present.
"""

from __future__ import annotations

import csv
import os

from rdkit import Chem
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")


def read_smiles_lines(path: str) -> list[tuple[str, str]]:
    """Read a .smi/.txt file: one SMILES per line, optional name as 2nd token."""
    items: list[tuple[str, str]] = []
    with open(path, "r", encoding="utf-8") as fh:
        for i, raw in enumerate(fh):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split(None, 1)
            smi = parts[0]
            name = parts[1].strip() if len(parts) > 1 else f"mol_{i+1}"
            items.append((name or f"mol_{i+1}", smi))
    return items


def read_sdf(path: str) -> list[tuple[str, str]]:
    """Read an .sdf/.mol file into (id, canonical_smiles) tuples."""
    items: list[tuple[str, str]] = []
    supplier = Chem.SDMolSupplier(path, sanitize=True, removeHs=False)
    for i, mol in enumerate(supplier):
        if mol is None:
            continue
        name = ""
        if mol.HasProp("_Name"):
            name = mol.GetProp("_Name").strip()
        if not name:
            name = f"mol_{i+1}"
        try:
            smi = Chem.MolToSmiles(mol)
        except Exception:
            continue
        items.append((name, smi))
    return items


def read_csv(path: str, smiles_col: str,
             id_col: str | None = None) -> list[tuple[str, str]]:
    """Read a CSV given the SMILES column (and optional ID column)."""
    items: list[tuple[str, str]] = []
    with open(path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for i, r in enumerate(reader):
            smi = (r.get(smiles_col) or "").strip()
            if not smi:
                continue
            mid = (r.get(id_col) or "").strip() if id_col else f"mol_{i+1}"
            items.append((mid or f"mol_{i+1}", smi))
    return items


def csv_columns(path: str) -> list[str]:
    """Return the header columns of a CSV (for column pickers)."""
    with open(path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        return list(reader.fieldnames or [])


def detect_and_read(path: str) -> list[tuple[str, str]]:
    """Auto-read by extension for non-CSV formats (.smi/.txt/.sdf/.mol).
    CSV is handled separately because it needs column selection."""
    ext = os.path.splitext(path)[1].lower()
    if ext in (".smi", ".txt"):
        return read_smiles_lines(path)
    if ext in (".sdf", ".mol"):
        return read_sdf(path)
    raise ValueError(f"Use read_csv() for CSV; unsupported extension: {ext}")
