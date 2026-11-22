# utils/descriptors.py
# Whole-molecule descriptor calculation using RDKit.
# Always recalculate on the full molecule — never estimate from fragments.

from rdkit.Chem import rdMolDescriptors


def calc_descriptors(mol):
    """Return the four benchmark descriptors for a molecule.

    Parameters
    ----------
    mol : rdkit.Chem.Mol

    Returns
    -------
    dict with keys: hbd, hba, ar, rotb
    """
    return {
        "hbd":  rdMolDescriptors.CalcNumHBD(mol),
        "hba":  rdMolDescriptors.CalcNumHBA(mol),
        "ar":   rdMolDescriptors.CalcNumAromaticRings(mol),
        "rotb": rdMolDescriptors.CalcNumRotatableBonds(mol),
    }
