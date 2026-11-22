# =============================================================================
# Benchmark configuration — shared by all notebooks and utils
# =============================================================================

# ---------------------------------------------------------------------------
# Parent datasets
# ---------------------------------------------------------------------------
N_PARENTS   = 1000
RANDOM_SEED = 42

# Development / initial benchmark set (1000 parents)
PARENT_DATASET = "../data/chembl_1000_parents.csv"

# Independent test set — created by 00_prepare_independent_test.ipynb
# Used for the primary results reported in the manuscript.
INDEPENDENT_TEST_DATASET = "../data/chembl_independent_test.csv"

# To run the benchmark on the independent test set, set:
#   ACTIVE_PARENT_DATASET = INDEPENDENT_TEST_DATASET
# For development runs, use PARENT_DATASET.
ACTIVE_PARENT_DATASET = INDEPENDENT_TEST_DATASET   # ← primary results use this

# ---------------------------------------------------------------------------
# Results output
# ---------------------------------------------------------------------------
# New benchmark results (04 methods, independent test set)
RESULTS_DIR = "../results/new_benchmark"

# Legacy results from previous benchmark are in results_legacy/
# Do not write to results_legacy/ — read-only reference only.

# ---------------------------------------------------------------------------
# ChEMBL collection filters  (used in 00_ and 01_prepare_parents)
# ---------------------------------------------------------------------------
MW_MIN   = 100
MW_MAX   = 500
HBD_MAX  = 10
HBA_MAX  = 20
ROTB_MAX = 10
ALOGP_MAX = 5
N_PARENT_CANDIDATES_MULTIPLIER = 5
SAVE_EVERY = 500

# ---------------------------------------------------------------------------
# Candidate generation
# ---------------------------------------------------------------------------
N_CANDIDATES_PER_PARENT = 100
MAX_GENERATION_ATTEMPTS = 1000

# ---------------------------------------------------------------------------
# Evaluation: candidate budgets
# ---------------------------------------------------------------------------
K_VALUES = [1, 5, 10, 20, 50, 100]

# ---------------------------------------------------------------------------
# Task definitions — 7 tasks; AR −1 is excluded from all analyses.
# ---------------------------------------------------------------------------
TASKS = {
    "hbd_plus1":  {"descriptor": "hbd",  "delta":  1, "scaffold_tuner_mode": "add_hbd"},
    "hbd_minus1": {"descriptor": "hbd",  "delta": -1, "scaffold_tuner_mode": "remove_hbd"},
    "hba_plus1":  {"descriptor": "hba",  "delta":  1, "scaffold_tuner_mode": "add_hba"},
    "hba_minus1": {"descriptor": "hba",  "delta": -1, "scaffold_tuner_mode": "remove_hba"},
    "rotb_plus1": {"descriptor": "rotb", "delta":  1, "scaffold_tuner_mode": "add_rb"},
    "rotb_minus1":{"descriptor": "rotb", "delta": -1, "scaffold_tuner_mode": "remove_rb"},
    "ar_plus1":   {"descriptor": "ar",   "delta":  1, "scaffold_tuner_mode": "add_ar"},
    # Aromatic ring decrease is not included in this benchmark.
}

TASK_NAMES = list(TASKS.keys())

# Descriptors evaluated (whole-molecule RDKit recalculation for every candidate).
EVAL_DESCRIPTORS = ["hbd", "hba", "ar", "rotb"]
