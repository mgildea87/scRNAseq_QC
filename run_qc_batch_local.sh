#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONDA_INIT_DEFAULT="/gpfs/data/cvrcbioinfolab/gildem01/conda_envs/anaconda3/condaload_r.sh"
CONDA_INIT="${CONDA_INIT:-${CONDA_INIT_DEFAULT}}"

if [[ -n "${CONDA_INIT}" && -f "${CONDA_INIT}" ]]; then
	# Some environment init scripts assume variables may be unset.
	set +u
	source "${CONDA_INIT}"
	set -u
else
	echo "WARN: CONDA_INIT not found (${CONDA_INIT}). Running with Rscript from PATH." >&2
fi

exec Rscript "${SCRIPT_DIR}/qc_batch_runner.R" "$@"
