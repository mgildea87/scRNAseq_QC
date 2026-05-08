#!/bin/bash
set -euo pipefail

CONDA_INIT="/gpfs/data/cvrcbioinfolab/gildem01/conda_envs/anaconda3/condaload_r.sh"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Some environment init scripts assume variables may be unset.
set +u
source "${CONDA_INIT}"
set -u

exec Rscript "${SCRIPT_DIR}/qc_batch_runner.R" "$@"
