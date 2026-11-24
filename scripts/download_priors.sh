#!/usr/bin/env bash
# Download REINVENT4 LibInvent prior from Zenodo
# Official version record: https://zenodo.org/records/20701824
# (Concept DOI: https://doi.org/10.5281/zenodo.15641296)
# Compatible with REINVENT4 4.9.x
# Run from repository root: bash scripts/download_priors.sh

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PRIOR_DIR="${REPO_ROOT}/external/reinvent4_priors"
PRIOR_FILE="${PRIOR_DIR}/libinvent.prior"
ZENODO_URL="https://zenodo.org/records/20701824/files/libinvent.prior?download=1"

mkdir -p "${PRIOR_DIR}"

if [ -f "${PRIOR_FILE}" ]; then
    echo "Prior already exists: ${PRIOR_FILE}  (delete to re-download)"
    exit 0
fi

echo "Downloading LibInvent prior from Zenodo version record 20701824 ..."
echo "  URL: ${ZENODO_URL}"
curl -L --fail --progress-bar -o "${PRIOR_FILE}" "${ZENODO_URL}"

SIZE=$(du -h "${PRIOR_FILE}" | cut -f1)
echo "Downloaded: ${PRIOR_FILE}  (${SIZE})"

echo "Verifying with REINVENT4 ..."
python3 -c "
import torch, sys
prior = torch.load(sys.argv[1], map_location='cpu', weights_only=False)
print('  model_type :', prior.get('model_type','unknown'))
print('  version    :', prior.get('version','unknown'))
print('  OK')
" "${PRIOR_FILE}"
