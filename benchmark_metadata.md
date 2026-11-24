# Benchmark Metadata

## Software versions

| Software | Version |
|---|---|
| REINVENT4 | 4.9.34 |
| PyTorch | 2.13.0 |
| Scaffold-tuner | 1.0.3 (commit 5bdfc8c) |
| CReM | (see requirements.txt) |
| RDKit | 2026.03.1 |

## REINVENT4 LibInvent prior

| Item | Value |
|---|---|
| **Filename** | `libinvent.prior` |
| **Source** | Zenodo record 15641296 |
| **DOI** | https://doi.org/10.5281/zenodo.15641296 |
| **Download URL** | https://zenodo.org/records/15641296/files/libinvent.prior?download=1 |
| **Compatible with** | REINVENT4 4.9.x |
| **SHA-256** | *(to be recorded after download — run `sha256sum external/reinvent4_priors/libinvent.prior`)* |

### Download instructions

```bash
bash scripts/download_priors.sh
```

Or manually:

```bash
curl -L -o external/reinvent4_priors/libinvent.prior \
  "https://zenodo.org/records/15641296/files/libinvent.prior?download=1"
```

### Notes on prior compatibility

- Zenodo record **15641296** is the official prior repository for REINVENT4 4.x.
- The older record 15641297 contains priors trained with an earlier version and
  fails the hash validation introduced in REINVENT4 ≥ 4.9.x.
- After download, record the SHA-256 in this file to enable full reproducibility.

## Benchmark datasets

| Dataset | Description | File |
|---|---|---|
| Development set | 1000 drug-like parents from ChEMBL | `data/chembl_1000_parents.csv` |
| Independent test set | 1000 drug-like parents from ChEMBL (no overlap with development set) | `data/chembl_independent_test.csv` |

Primary results are reported on the **independent test set**.
Both sets use the same physicochemical filters (MW 100–500, HBD ≤ 10, HBA ≤ 20, RotB ≤ 10, ALogP ≤ 5).

## Benchmark tasks

7 descriptor-modification tasks (AR −1 excluded):
HBD +1, HBD −1, HBA +1, HBA −1, RotB +1, RotB −1, AR +1

## Candidate budget

k = 1, 5, 10, 20, 50, 100

## CReM database

| Item | Value |
|---|---|
| **Database** | `chembl33_sa2_f5.db` |
| **Source** | CReM public release (ChEMBL 33) |
| **Note** | The independent test set is derived from the same ChEMBL release. Test set compounds may overlap with CReM fragment statistics. This limitation is noted in the manuscript Methods. |
