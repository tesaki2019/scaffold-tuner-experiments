from rdkit import Chem
from utils.descriptors import calc_descriptors

def make_candidate_record(
    chembl_id,
    parent_smiles,
    generated_smiles,
    parent_desc,
    mode=None,
):
    new_mol = Chem.MolFromSmiles(generated_smiles)

    if new_mol is None:
        return None

    new_desc = calc_descriptors(new_mol)

    record = {
        "chembl_id": chembl_id,
        "parent_smiles": parent_smiles,
        "generated_smiles": generated_smiles,

        "hbd_parent": parent_desc["hbd"],
        "hbd_generated": new_desc["hbd"],
        "delta_hbd": new_desc["hbd"] - parent_desc["hbd"],

        "hba_parent": parent_desc["hba"],
        "hba_generated": new_desc["hba"],
        "delta_hba": new_desc["hba"] - parent_desc["hba"],

        "ar_parent": parent_desc["ar"],
        "ar_generated": new_desc["ar"],
        "delta_ar": new_desc["ar"] - parent_desc["ar"],

        "rotb_parent": parent_desc["rotb"],
        "rotb_generated": new_desc["rotb"],
        "delta_rotb": new_desc["rotb"] - parent_desc["rotb"],
    }

    if mode is not None:
        record["mode"] = mode

    return record

def append_candidate_record(
    records,
    chembl_id,
    parent_smiles,
    generated_smiles,
    parent_desc,
    mode=None,
):
    record = make_candidate_record(
        chembl_id=chembl_id,
        parent_smiles=parent_smiles,
        generated_smiles=generated_smiles,
        parent_desc=parent_desc,
        mode=mode,
    )

    if record is not None:
        records.append(record)

