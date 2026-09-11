from rdkit.Chem import rdMolDescriptors

def calc_descriptors(mol):
    return {
        "hbd": rdMolDescriptors.CalcNumHBD(mol),
        "hba": rdMolDescriptors.CalcNumHBA(mol),
        "ar": rdMolDescriptors.CalcNumAromaticRings(mol),
        "rotb": rdMolDescriptors.CalcNumRotatableBonds(mol),
    }