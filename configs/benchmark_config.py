# =========================
# Benchmark settings
# =========================
# Number of parent compounds
N_PARENTS = 1000

# Maximum number of generated compounds per parent
N_CANDIDATES_PER_PARENT = 100

# Random seed
RANDOM_SEED = 42

# Maximum number of attempts to generate a valid compound
MAX_GENERATION_ATTEMPTS = 1000

# ChEMBL filtering
MW_MIN = 100
MW_MAX = 500

HBD_MAX = 10
HBA_MAX = 20
ROTB_MAX = 10
ALOGP_MAX = 5

# Data collection
N_PARENT_CANDIDATES_MULTIPLIER = 5
SAVE_EVERY = 500