# utils/scaffold.py
# Scaffold-related utilities shared across all benchmark methods.
# Using a single module ensures that scaffold retention is judged identically
# for Scaffold-tuner, Random local editing, CReM, and REINVENT4.

from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold


def get_murcko_scaffold(mol):
    """Return the Murcko scaffold molecule for *mol*.

    Returns None if the scaffold cannot be computed.
    """
    try:
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
        if scaffold is None or scaffold.GetNumAtoms() == 0:
            return None
        Chem.SanitizeMol(scaffold)
        return scaffold
    except Exception:
        return None


def get_murcko_scaffold_smiles(mol):
    """Return the canonical SMILES of the Murcko scaffold, or None."""
    scaffold = get_murcko_scaffold(mol)
    if scaffold is None:
        return None
    return Chem.MolToSmiles(scaffold)


def get_scaffold_atom_ids(mol):
    """Return the atom indices in *mol* that belong to its Murcko scaffold.

    Useful for:
    - Random local editing: restrict edits to non-scaffold atoms.
    - CReM: pass these indices as protected_ids to prevent scaffold edits.

    Returns an empty set if the scaffold match fails.
    """
    scaffold = get_murcko_scaffold(mol)
    if scaffold is None:
        return set()

    matches = mol.GetSubstructMatches(scaffold)
    if not matches:
        return set()

    # Use the first (largest) match.
    return set(matches[0])


def scaffold_is_retained(parent_mol, candidate_mol):
    """Return True if the parent Murcko scaffold is a substructure of *candidate_mol*.

    This is the single shared scaffold-retention check used for all methods.
    Criterion: candidate_mol.HasSubstructMatch(parent_scaffold).
    """
    parent_scaffold = get_murcko_scaffold(parent_mol)
    if parent_scaffold is None:
        # If the parent has no meaningful scaffold, we cannot check retention.
        return False

    return candidate_mol.HasSubstructMatch(parent_scaffold)
