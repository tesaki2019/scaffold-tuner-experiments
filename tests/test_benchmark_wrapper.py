"""
Verification script for scaffold_tuner_benchmark.py on Scaffold-tuner v1.0.3.
Run on the calculation server:
  cd scaffold-tuner-experiments
  python tests/test_benchmark_wrapper.py
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")
import pandas as pd

print("Testing imports...")
from scaffold_tuner.scaffold_intervention import propose_structures
print("scaffold-tuner OK")
from utils.scaffold_tuner_benchmark import generate_all_candidates
print("benchmark wrapper OK")

parent = "NCC(=O)Nc1ccc(Cl)cc1C(=O)c1ccc[nH]1"
print(f"\nParent: {parent}")

for mode in ["add_hbd", "add_hba", "remove_hbd", "remove_hba", "add_ar", "add_rb", "remove_rb"]:
    try:
        cands, stats = generate_all_candidates(parent, mode=mode, max_candidates=20)
        ok = stats["n_examined"] == sum(
            stats[k] for k in ["n_invalid","n_scaffold_failed","n_duplicate","n_valid_unique"]
        )
        print(f"  {mode:15s}: {len(cands):3d} candidates  identity={ok}")
        assert ok
    except Exception as e:
        print(f"  {mode:15s}: ERROR -- {e}")
        raise

cands, stats = generate_all_candidates(parent, mode="add_hbd", max_candidates=100)
failure_count = sum(1 for c in cands if c["actual_delta"] != 1)
print(f"\nadd_hbd: {len(cands)} total, {failure_count} non-success also ranked (expected >0)")

print("\nRunning 10-parent smoke test...")
from configs.benchmark_config import TASKS, ACTIVE_PARENT_DATASET, N_CANDIDATES_PER_PARENT
from utils.parents import prepare_parent
from utils.records import make_candidate_record

df = pd.read_csv(ACTIVE_PARENT_DATASET)
if "canonical_smiles" not in df.columns and "smiles" in df.columns:
    df = df.rename(columns={"smiles":"canonical_smiles"})
if "parent_id" not in df.columns and "chembl_id" in df.columns:
    df = df.rename(columns={"chembl_id":"parent_id"})

total, errors = 0, 0
for _, row in df.head(10).iterrows():
    prepared = prepare_parent(row)
    if not prepared: continue
    parent_id, parent_smiles, parent_desc, _ = prepared
    for task_name, task_cfg in TASKS.items():
        try:
            cands, stats = generate_all_candidates(
                parent_smiles,
                target=task_cfg["descriptor"].upper(),
                delta=task_cfg["delta"],
                max_candidates=N_CANDIDATES_PER_PARENT,
            )
        except Exception as e:
            print(f"  ERROR {parent_id}/{task_name}: {e}")
            errors += 1; continue
        for cand in cands:
            rec = make_candidate_record(
                method="scaffold_tuner", parent_id=parent_id,
                parent_smiles=parent_smiles, task_name=task_name,
                candidate_rank=cand["candidate_rank"],
                candidate_smiles=cand["generated_smiles"],
                parent_desc=parent_desc,
                generation_strategy=cand.get("generation_strategy"),
            )
            if rec: total += 1

print(f"10-parent smoke test: {total} records, {errors} errors")
assert errors == 0
print("ALL TESTS PASSED")
