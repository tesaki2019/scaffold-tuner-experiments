# utils/parents.py
# Parent molecule preparation shared across benchmark notebooks.

from rdkit import Chem

from utils.descriptors import calc_descriptors
from utils.scaffold import get_murcko_scaffold_smiles


def prepare_parent(row):
    """Parse one row of the parent CSV and return structured data.

    Parameters
    ----------
    row : pandas Series  must have 'smiles' and either 'chembl_id' or 'parent_id'

    Returns
    -------
    tuple (parent_id, parent_smiles, parent_desc, scaffold_smiles)
    or None if the SMILES is invalid.
    """
    smiles = row.get("smiles") or row.get("parent_smiles") or row.get("canonical_smiles")
    parent_id = row.get("parent_id") or row.get("chembl_id")

    mol = Chem.MolFromSmiles(smiles) if smiles else None
    if mol is None:
        return None

    canonical_smiles = Chem.MolToSmiles(mol)
    parent_desc      = calc_descriptors(mol)
    scaffold_smiles  = get_murcko_scaffold_smiles(mol)

    return parent_id, canonical_smiles, parent_desc, scaffold_smiles
