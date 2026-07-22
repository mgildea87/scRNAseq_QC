#!/bin/bash
set -euo pipefail

# Prevent RStudio-inherited environment variables from overriding batch R libs.
unset R_HOME R_LIBS R_LIBS_USER R_PROFILE R_PROFILE_USER R_ENVIRON R_ENVIRON_USER || true

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONDA_INIT_DEFAULT="/gpfs/data/cvrcbioinfolab/gildem01/conda_envs/anaconda3/condaload_r.sh"
CONDA_INIT="${CONDA_INIT:-${CONDA_INIT_DEFAULT}}"

if [[ -n "${CONDA_INIT}" && -f "${CONDA_INIT}" ]]; then
	# Some environment init scripts assume variables may be unset.
	set +u
	source /gpfs/data/cvrcbioinfolab/gildem01/conda_envs/anaconda3/condaload_r.sh
	set -u
else
	echo "ERROR: CONDA_INIT not found (${CONDA_INIT}). Refusing to run without the conda R environment." >&2
	exit 1
fi

exec Rscript --vanilla "${SCRIPT_DIR}/qc_batch_runner.R" "$@"
