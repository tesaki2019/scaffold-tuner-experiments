# utils/records.py
# Candidate record creation shared across all benchmark methods.

from rdkit import Chem
from rdkit.Chem import DataStructs
from rdkit.Chem import AllChem

from utils.descriptors import calc_descriptors
from utils.scaffold import scaffold_is_retained
from utils.evaluation import calculate_descriptor_delta, is_success, is_strict_success
from configs.benchmark_config import TASKS, EVAL_DESCRIPTORS


def _tanimoto(mol_a, mol_b, radius=2, n_bits=2048):
    """Morgan fingerprint Tanimoto similarity between two molecules."""
    try:
        fp_a = AllChem.GetMorganFingerprintAsBitVect(mol_a, radius, nBits=n_bits)
        fp_b = AllChem.GetMorganFingerprintAsBitVect(mol_b, radius, nBits=n_bits)
        return DataStructs.TanimotoSimilarity(fp_a, fp_b)
    except Exception:
        return None


def _collateral_change(parent_desc, candidate_desc, task_name):
    """Sum of |delta| for all non-target descriptors."""
    target_desc = TASKS[task_name]["descriptor"]
    return sum(
        abs(candidate_desc[d] - parent_desc[d])
        for d in EVAL_DESCRIPTORS
        if d != target_desc
        and candidate_desc[d] is not None
        and parent_desc[d] is not None
    )


def make_candidate_record(
    method,
    parent_id,
    parent_smiles,
    task_name,
    candidate_rank,
    candidate_smiles,
    parent_desc=None,
    generation_strategy=None,
):
    """Build a single candidate record dict.

    Parameters
    ----------
    method              : str   e.g. "scaffold_tuner", "random", "crem"
    parent_id           : str   e.g. ChEMBL ID
    parent_smiles       : str
    task_name           : str   key in TASKS, e.g. "hbd_plus1"
    candidate_rank      : int   1-based rank in the order returned by the method
    candidate_smiles    : str
    parent_desc         : dict  optional — recalculated from parent_smiles if None
    generation_strategy : str   optional — "tier1" / "tier2" for Scaffold-tuner

    Returns
    -------
    dict or None if candidate_smiles is not a valid molecule
    """
    parent_mol = Chem.MolFromSmiles(parent_smiles)
    if parent_mol is None:
        return None

    candidate_mol = Chem.MolFromSmiles(candidate_smiles)
    valid = candidate_mol is not None

    if parent_desc is None:
        parent_desc = calc_descriptors(parent_mol)

    if valid:
        cand_desc      = calc_descriptors(candidate_mol)
        retained       = scaffold_is_retained(parent_mol, candidate_mol)
        deltas         = calculate_descriptor_delta(parent_desc, cand_desc)
        success        = is_success(parent_desc, cand_desc, task_name)
        strict_success = is_strict_success(parent_desc, cand_desc, task_name)
        tanimoto       = _tanimoto(parent_mol, candidate_mol)
        collateral     = _collateral_change(parent_desc, cand_desc, task_name)
    else:
        cand_desc      = {d: None for d in EVAL_DESCRIPTORS}
        retained       = False
        deltas         = {d: None for d in EVAL_DESCRIPTORS}
        success        = False
        strict_success = False
        tanimoto       = None
        collateral     = None

    record = {
        # --- identity ---
        "method":          method,
        "parent_id":       parent_id,
        "chembl_id":       parent_id,          # backward-compat alias
        "parent_smiles":   parent_smiles,
        "task":            task_name,
        "candidate_rank":  candidate_rank,
        "candidate_smiles": candidate_smiles,
        # --- QC ---
        "valid":             valid,
        "scaffold_retained": retained,
        # --- descriptors: parent ---
        "hbd_parent":  parent_desc["hbd"],
        "hba_parent":  parent_desc["hba"],
        "ar_parent":   parent_desc["ar"],
        "rotb_parent": parent_desc["rotb"],
        # --- descriptors: candidate ---
        "hbd_generated":  cand_desc["hbd"],
        "hba_generated":  cand_desc["hba"],
        "ar_generated":   cand_desc["ar"],
        "rotb_generated": cand_desc["rotb"],
        # --- deltas ---
        "delta_hbd":  deltas["hbd"],
        "delta_hba":  deltas["hba"],
        "delta_ar":   deltas["ar"],
        "delta_rotb": deltas["rotb"],
        # --- outcome ---
        "success":           success,
        "strict_success":    strict_success,
        # --- secondary metrics (pre-computed for efficiency) ---
        "tanimoto_to_parent":   tanimoto,
        "collateral_change":    collateral,
    }

    if generation_strategy is not None:
        record["generation_strategy"] = generation_strategy

    return record


def append_candidate_record(records, **kwargs):
    """Call make_candidate_record and append the result to *records* if valid."""
    record = make_candidate_record(**kwargs)
    if record is not None:
        records.append(record)


def compute_first_success_rank(df_candidate, all_parent_ids):
    """Compute the first candidate_rank at which success occurs per (parent_id, task).

    Parameters
    ----------
    df_candidate  : DataFrame with columns [parent_id, task, candidate_rank, success]
    all_parent_ids: iterable of all parent_ids in the benchmark

    Returns
    -------
    DataFrame with columns [parent_id, task, first_success_rank]
    NaN means no success candidate was found (candidacy 0 or all failed).
    """
    import pandas as pd
    from configs.benchmark_config import TASK_NAMES

    full_index = pd.MultiIndex.from_product(
        [sorted(all_parent_ids), TASK_NAMES],
        names=["parent_id", "task"],
    )

    success_only = df_candidate[df_candidate["success"] == True]
    first_rank = (
        success_only
        .groupby(["parent_id", "task"])["candidate_rank"]
        .min()
        .reindex(full_index)   # NaN for no-success pairs
    )
    return first_rank.reset_index().rename(columns={"candidate_rank": "first_success_rank"})
