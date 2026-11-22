# utils/libinvent.py
# LibInvent scaffold SMILES preparation.
# Converts a parent molecule into a scaffold SMILES with attachment points
# ([*:0], [*:1], ...) for use as LibInvent input.

import re
from rdkit import Chem

from utils.scaffold import get_murcko_scaffold, get_scaffold_atom_ids


def make_libinvent_scaffold(parent_smiles):
    """Build a LibInvent-compatible scaffold SMILES from a parent molecule.

    The Murcko scaffold is extracted from the parent.  Each bond that
    crosses the scaffold / R-group boundary is cut, and an attachment point
    ([*:0], [*:1], ...) is placed on the scaffold atom at that cut.

    Returns
    -------
    (scaffold_smiles, None)  on success
    (None, reason_str)       when scaffold cannot be built

    Failure reasons
    ---------------
    invalid_parent      : parent SMILES does not parse
    no_scaffold         : Murcko scaffold is empty (e.g. pure aliphatic chain)
    no_attachment_bonds : scaffold has no bonds to any R-group
    scaffold_match_failed: substructure match between parent and scaffold failed
    sanitize_failed     : RDKit sanitisation of the annotated scaffold failed
    """
    mol = Chem.MolFromSmiles(parent_smiles)
    if mol is None:
        return None, "invalid_parent"

    scaffold = get_murcko_scaffold(mol)
    if scaffold is None:
        return None, "no_scaffold"

    scaffold_ids = get_scaffold_atom_ids(mol)
    if not scaffold_ids:
        return None, "no_scaffold"

    # Find bonds that cross the scaffold / R-group boundary.
    attachment_bonds = []  # (bond_idx, scaffold_atom_in_parent)
    for bond in mol.GetBonds():
        a1, a2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if (a1 in scaffold_ids) != (a2 in scaffold_ids):
            scaffold_atom = a1 if a1 in scaffold_ids else a2
            attachment_bonds.append((bond.GetIdx(), scaffold_atom))

    if not attachment_bonds:
        return None, "no_attachment_bonds"

    # Map: parent atom index → scaffold atom index
    matches = mol.GetSubstructMatches(scaffold)
    if not matches:
        return None, "scaffold_match_failed"

    parent_to_scaffold = {
        parent_idx: sc_idx for sc_idx, parent_idx in enumerate(matches[0])
    }

    # Build an editable scaffold with attachment-point dummies.
    rw = Chem.RWMol(scaffold)
    for ap_num, (_, scaffold_atom_in_parent) in enumerate(attachment_bonds):
        sc_idx = parent_to_scaffold.get(scaffold_atom_in_parent)
        if sc_idx is None:
            continue
        dummy = rw.AddAtom(Chem.Atom(0))
        rw.GetAtomWithIdx(dummy).SetAtomMapNum(ap_num)   # 0 → [*:0], 1 → [*:1]
        rw.AddBond(sc_idx, dummy, Chem.BondType.SINGLE)

    try:
        Chem.SanitizeMol(rw)
        raw_smi = Chem.MolToSmiles(rw.GetMol())
    except Exception as e:
        return None, f"sanitize_failed: {e}"

    # RDKit writes map_num=0 as bare '*'; restore it to '[*:0]'.
    scaffold_smi = _restore_zero_attachment(raw_smi)
    return scaffold_smi, None


def _restore_zero_attachment(smi):
    """Replace bare '*' with '[*:0]' in a SMILES string.

    RDKit outputs atoms with atom-map number 0 as '*' rather than '[*:0]'.
    LibInvent expects the explicit '[*:0]' notation.
    """
    # Temporarily protect already-annotated attachment points [*:N] (N>=1).
    protected = re.sub(r'\[\*:(\d+)\]', r'__AP\1__', smi)
    # Replace bare '*' with [*:0].
    replaced = protected.replace('*', '[*:0]')
    # Restore protected entries.
    return re.sub(r'__AP(\d+)__', r'[*:\1]', replaced)
