from rdkit import Chem
from utils.descriptors import calc_descriptors


def prepare_parent(row):
    parent_smiles = row["smiles"]
    chembl_id = row["chembl_id"]

    parent = Chem.MolFromSmiles(parent_smiles)
    if parent is None:
        return None

    parent_desc = calc_descriptors(parent)

    return chembl_id, parent_smiles, parent_desc