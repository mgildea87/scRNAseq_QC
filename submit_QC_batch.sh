#!/bin/bash
# =============================================================================
# submit_QC_batch.sh
#
# Coordinator script — reads a sample sheet and submits ONE independent
# SLURM job per sample on the cpu_short partition.
#
# Run with plain bash (not sbatch):
#
#   bash /path/to/submit_QC_batch.sh \
#       --sample_sheet /abs/path/to/samples.csv \
#       --output_dir   /abs/path/to/results/QC
#
# OPTIONAL FLAGS (forwarded to qc_batch_runner.R):
#   --mem         <GB>    memory per job in GB              (default: 32)
#   --merge_mem   <GB>    memory for merge job in GB         (default: 64)
#   --integration_mem <GB> memory for integration job in GB   (default: 64)
#   --time        <HH:MM> wall time per job                 (default: 2:00:00)
#   --merge_only  <TRUE/FALSE> skip per-sample QC and run merge_analysis from existing sample RDS files (default: FALSE)
#   --skip_merge  <TRUE/FALSE> skip submitting merge job     (default: FALSE)
#   --run_integration <TRUE/FALSE> render integrate_RNA.Rmd after merge (default: FALSE)
#   --integration_level <Batch|Sample> required when integration runs
#   --integration_only <TRUE/FALSE> run only integration from merged_QC.rds (default: FALSE)
#   --use_cellbender <TRUE/FALSE> use cellbender_filtered.h5 as filtered input (default: FALSE)
# =============================================================================

set -euo pipefail
ORIGINAL_ARGS=("$@")

# ── Fixed paths ───────────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEMPLATE_DIR="${TEMPLATE_DIR:-${SCRIPT_DIR}}"
CONDA_INIT_DEFAULT="/gpfs/data/cvrcbioinfolab/gildem01/conda_envs/anaconda3/condaload_r.sh"
CONDA_INIT="${CONDA_INIT:-${CONDA_INIT_DEFAULT}}"

# ── Parse arguments ───────────────────────────────────────────────────────────
SAMPLE_SHEET=""
OUTPUT_DIR="QC"
MEM_GB=32
MERGE_MEM_GB=64
INTEGRATION_MEM_GB=64
WALL_TIME="2:00:00"
MERGE_ONLY="FALSE"
SKIP_MERGE="FALSE"
RUN_INTEGRATION="FALSE"
INTEGRATION_LEVEL=""
INTEGRATION_ONLY="FALSE"
USE_CELLBENDER="FALSE"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --help|-h)
      sed -n '1,26p' "$0" # prints header block lines
      exit 0
      ;;
    --sample_sheet)
      SAMPLE_SHEET="${2-}"
      if [[ -z "${SAMPLE_SHEET}" || "${SAMPLE_SHEET}" == --* ]]; then
        echo "ERROR: --sample_sheet requires a value." >&2
        exit 1
      fi
      shift 2
      ;;
    --output_dir)
      OUTPUT_DIR="${2-}"
      if [[ -z "${OUTPUT_DIR}" || "${OUTPUT_DIR}" == --* ]]; then
        echo "ERROR: --output_dir requires a value." >&2
        exit 1
      fi
      shift 2
      ;;
    --mem)
      MEM_GB="${2-}"
      if [[ -z "${MEM_GB}" || "${MEM_GB}" == --* ]]; then
        echo "ERROR: --mem requires a value." >&2
        exit 1
      fi
      shift 2
      ;;
    --merge_mem)
      MERGE_MEM_GB="${2-}"
      if [[ -z "${MERGE_MEM_GB}" || "${MERGE_MEM_GB}" == --* ]]; then
        echo "ERROR: --merge_mem requires a value." >&2
        exit 1
      fi
      shift 2
      ;;
    --integration_mem)
      INTEGRATION_MEM_GB="${2-}"
      if [[ -z "${INTEGRATION_MEM_GB}" || "${INTEGRATION_MEM_GB}" == --* ]]; then
        echo "ERROR: --integration_mem requires a value." >&2
        exit 1
      fi
      shift 2
      ;;
    --time)
      WALL_TIME="${2-}"
      if [[ -z "${WALL_TIME}" || "${WALL_TIME}" == --* ]]; then
        echo "ERROR: --time requires a value." >&2
        exit 1
      fi
      shift 2
      ;;
    --merge_only)
      MERGE_ONLY="${2-}"
      if [[ -z "${MERGE_ONLY}" || "${MERGE_ONLY}" == --* ]]; then
        echo "ERROR: --merge_only must be TRUE or FALSE." >&2
        exit 1
      fi
      shift 2
      ;;
    --skip_merge)
      SKIP_MERGE="${2-}"
      if [[ -z "${SKIP_MERGE}" || "${SKIP_MERGE}" == --* ]]; then
        echo "ERROR: --skip_merge must be TRUE or FALSE." >&2
        exit 1
      fi
      shift 2
      ;;
    --run_integration)
      RUN_INTEGRATION="${2-}"
      if [[ -z "${RUN_INTEGRATION}" || "${RUN_INTEGRATION}" == --* ]]; then
        echo "ERROR: --run_integration must be TRUE or FALSE." >&2
        exit 1
      fi
      shift 2
      ;;
    --integration_level)
      INTEGRATION_LEVEL="${2-}"
      if [[ -z "${INTEGRATION_LEVEL}" || "${INTEGRATION_LEVEL}" == --* ]]; then
        echo "ERROR: --integration_level requires a value." >&2
        exit 1
      fi
      shift 2
      ;;
    --integration_only)
      INTEGRATION_ONLY="${2-}"
      if [[ -z "${INTEGRATION_ONLY}" || "${INTEGRATION_ONLY}" == --* ]]; then
        echo "ERROR: --integration_only must be TRUE or FALSE." >&2
        exit 1
      fi
      shift 2
      ;;
    --use_cellbender)
      USE_CELLBENDER="${2-}"
      if [[ -z "${USE_CELLBENDER}" || "${USE_CELLBENDER}" == --* ]]; then
        echo "ERROR: --use_cellbender must be TRUE or FALSE." >&2
        exit 1
      fi
      shift 2
      ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
done

SKIP_MERGE_LOWER="$(echo "${SKIP_MERGE}" | tr '[:upper:]' '[:lower:]')"
if [[ "${SKIP_MERGE_LOWER}" =~ ^(true|t|1|yes|y)$ ]]; then
  SKIP_MERGE="TRUE"
elif [[ "${SKIP_MERGE_LOWER}" =~ ^(false|f|0|no|n)$ ]]; then
  SKIP_MERGE="FALSE"
else
  echo "ERROR: --skip_merge must be TRUE or FALSE." >&2
  exit 1
fi

RUN_INTEGRATION_LOWER="$(echo "${RUN_INTEGRATION}" | tr '[:upper:]' '[:lower:]')"
if [[ "${RUN_INTEGRATION_LOWER}" =~ ^(true|t|1|yes|y)$ ]]; then
  RUN_INTEGRATION="TRUE"
elif [[ "${RUN_INTEGRATION_LOWER}" =~ ^(false|f|0|no|n)$ ]]; then
  RUN_INTEGRATION="FALSE"
else
  echo "ERROR: --run_integration must be TRUE or FALSE." >&2
  exit 1
fi

INTEGRATION_ONLY_LOWER="$(echo "${INTEGRATION_ONLY}" | tr '[:upper:]' '[:lower:]')"
if [[ "${INTEGRATION_ONLY_LOWER}" =~ ^(true|t|1|yes|y)$ ]]; then
  INTEGRATION_ONLY="TRUE"
elif [[ "${INTEGRATION_ONLY_LOWER}" =~ ^(false|f|0|no|n)$ ]]; then
  INTEGRATION_ONLY="FALSE"
else
  echo "ERROR: --integration_only must be TRUE or FALSE." >&2
  exit 1
fi

if [[ "${INTEGRATION_ONLY}" == "TRUE" ]]; then
  RUN_INTEGRATION="TRUE"
fi

MERGE_ONLY_LOWER="$(echo "${MERGE_ONLY}" | tr '[:upper:]' '[:lower:]')"
if [[ "${MERGE_ONLY_LOWER}" =~ ^(true|t|1|yes|y)$ ]]; then
  MERGE_ONLY="TRUE"
elif [[ "${MERGE_ONLY_LOWER}" =~ ^(false|f|0|no|n)$ ]]; then
  MERGE_ONLY="FALSE"
else
  echo "ERROR: --merge_only must be TRUE or FALSE." >&2
  exit 1
fi

USE_CELLBENDER_LOWER="$(echo "${USE_CELLBENDER}" | tr '[:upper:]' '[:lower:]')"
if [[ "${USE_CELLBENDER_LOWER}" =~ ^(true|t|1|yes|y)$ ]]; then
  USE_CELLBENDER="TRUE"
elif [[ "${USE_CELLBENDER_LOWER}" =~ ^(false|f|0|no|n)$ ]]; then
  USE_CELLBENDER="FALSE"
else
  echo "ERROR: --use_cellbender must be TRUE or FALSE." >&2
  exit 1
fi

INTEGRATION_LEVEL_LOWER="$(echo "${INTEGRATION_LEVEL}" | tr '[:upper:]' '[:lower:]')"
if [[ "${RUN_INTEGRATION}" == "TRUE" ]]; then
  if [[ "${INTEGRATION_LEVEL_LOWER}" == "batch" ]]; then
    INTEGRATION_LEVEL="Batch"
  elif [[ "${INTEGRATION_LEVEL_LOWER}" == "sample" ]]; then
    INTEGRATION_LEVEL="Sample"
  else
    echo "ERROR: --integration_level is required when integration runs and must be Batch or Sample." >&2
    exit 1
  fi
fi

if [[ "${SKIP_MERGE}" == "TRUE" && "${RUN_INTEGRATION}" == "TRUE" && "${INTEGRATION_ONLY}" != "TRUE" ]]; then
  echo "ERROR: --run_integration TRUE requires merge; do not combine with --skip_merge TRUE." >&2
  exit 1
fi

if [[ "${MERGE_ONLY}" == "TRUE" && "${INTEGRATION_ONLY}" == "TRUE" ]]; then
  echo "ERROR: --merge_only and --integration_only are mutually exclusive." >&2
  exit 1
fi

if [[ "${MERGE_ONLY}" == "TRUE" && "${SKIP_MERGE}" == "TRUE" ]]; then
  echo "ERROR: --merge_only cannot be combined with --skip_merge TRUE." >&2
  exit 1
fi

if [[ -z "${SAMPLE_SHEET}" ]]; then
  echo "ERROR: --sample_sheet is required." >&2
  exit 1
fi

if [[ ! -f "${SAMPLE_SHEET}" ]]; then
  echo "ERROR: sample sheet not found: ${SAMPLE_SHEET}" >&2
  exit 1
fi

if [[ ! -f "${TEMPLATE_DIR}/qc_batch_runner.R" ]]; then
  echo "ERROR: TEMPLATE_DIR does not contain qc_batch_runner.R: ${TEMPLATE_DIR}" >&2
  exit 1
fi

# Resolve to absolute paths so jobs running in any working directory find them
SAMPLE_SHEET="$(realpath "${SAMPLE_SHEET}")"
OUTPUT_DIR="$(realpath -m "${OUTPUT_DIR}")"

# Ensure the main QC output directory exists even when using defaults
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${OUTPUT_DIR}/logs"
mkdir -p "${OUTPUT_DIR}/run_metadata"

# ── Record coordinator invocation metadata ───────────────────────────────────
INVOCATION_ID="$(date -u +%Y%m%dT%H%M%SZ)_$$"
RUN_CMD="bash $(realpath "$0")"
for arg in "${ORIGINAL_ARGS[@]}"; do
  RUN_CMD+=" $(printf '%q' "$arg")"
done

GIT_ROOT="$(git -C "${TEMPLATE_DIR}" rev-parse --show-toplevel 2>/dev/null || echo "${TEMPLATE_DIR}")"
REPO_URL="$(git -C "${GIT_ROOT}" remote get-url origin 2>/dev/null || echo "<unknown>")"
COMMIT_SHA="$(git -C "${GIT_ROOT}" rev-parse HEAD 2>/dev/null || echo "<unknown>")"
BRANCH_NAME="$(git -C "${GIT_ROOT}" rev-parse --abbrev-ref HEAD 2>/dev/null || echo "<unknown>")"
EXACT_TAG="$(git -C "${GIT_ROOT}" describe --tags --exact-match 2>/dev/null || echo "<none>")"
DESCRIBE_TAG="$(git -C "${GIT_ROOT}" describe --tags --always 2>/dev/null || echo "<none>")"
if [[ -n "$(git -C "${GIT_ROOT}" status --porcelain 2>/dev/null || true)" ]]; then
  WORKTREE_STATUS="dirty"
else
  WORKTREE_STATUS="clean"
fi

{
  echo "recorded_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "git_root=${GIT_ROOT}"
  echo "repo_url=${REPO_URL}"
  echo "commit_sha=${COMMIT_SHA}"
  echo "branch=${BRANCH_NAME}"
  echo "exact_tag=${EXACT_TAG}"
  echo "describe=${DESCRIBE_TAG}"
  echo "worktree_status=${WORKTREE_STATUS}"
} > "${OUTPUT_DIR}/run_metadata/release_info_submit.txt"

{
  echo "#!/bin/bash"
  echo "set -euo pipefail"
  echo "${RUN_CMD}"
} > "${OUTPUT_DIR}/run_metadata/run_command_submit_${INVOCATION_ID}.sh"

{
  echo "recorded_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "CONDA_DEFAULT_ENV=${CONDA_DEFAULT_ENV:-<unset>}"
  echo
  echo "conda info --envs:"
  conda info --envs 2>/dev/null || echo "<conda unavailable>"
} > "${OUTPUT_DIR}/run_metadata/conda_info_submit.txt"

if command -v conda >/dev/null 2>&1 && [[ -n "${CONDA_DEFAULT_ENV:-}" ]]; then
  conda env export -n "${CONDA_DEFAULT_ENV}" > "${OUTPUT_DIR}/run_metadata/conda_env_export_submit.yml" 2>/dev/null || true
fi

# ── Read sample names from CSV (skip header, column 1) ────────────────────────
mapfile -t SAMPLES < <(tail -n +2 "${SAMPLE_SHEET}" | cut -d',' -f1 | tr -d '\r')

if [[ ${#SAMPLES[@]} -eq 0 ]]; then
  echo "ERROR: no samples found in ${SAMPLE_SHEET}" >&2
  exit 1
fi

echo "Submitting ${#SAMPLES[@]} job(s)..."
echo "  Sample sheet : ${SAMPLE_SHEET}"
echo "  Output dir   : ${OUTPUT_DIR}"
echo "  Partition    : cpu_short  |  Per-sample Mem: ${MEM_GB}G  |  Merge Mem: ${MERGE_MEM_GB}G  |  Integration Mem: ${INTEGRATION_MEM_GB}G  |  Time: ${WALL_TIME}  |  Merge only: ${MERGE_ONLY}  |  Skip merge: ${SKIP_MERGE}  |  Run integration: ${RUN_INTEGRATION}  |  Integration level: ${INTEGRATION_LEVEL:-NA}  |  Integration only: ${INTEGRATION_ONLY}  |  Use CellBender: ${USE_CELLBENDER}"
echo "----------------------------------------------"

if [[ "${INTEGRATION_ONLY}" == "TRUE" ]]; then
  INTEGRATE_JOB_ID=$(sbatch \
    --job-name="QC_integrate" \
    --partition=cpu_short \
    --ntasks=1 \
    --cpus-per-task=1 \
    --mem="${INTEGRATION_MEM_GB}G" \
    --time="${WALL_TIME}" \
    --output="${OUTPUT_DIR}/logs/QC_integrate_%j.log" \
    --error="${OUTPUT_DIR}/logs/QC_integrate_%j.log" \
    --parsable \
    --wrap="
      set -euo pipefail
      if [[ -n \"${CONDA_INIT}\" && -f \"${CONDA_INIT}\" ]]; then
        set +u
        source /gpfs/data/cvrcbioinfolab/gildem01/conda_envs/anaconda3/condaload_r.sh
        set -u
      else
        echo \"ERROR: CONDA_INIT not found (${CONDA_INIT}). Refusing to run without the conda R environment.\" >&2
        exit 1
      fi
      echo \"Job ID   : \${SLURM_JOB_ID}\"
      echo \"Node     : \${SLURMD_NODENAME}\"
      echo \"Task     : integration only\"
      echo \"use_cellbender: ${USE_CELLBENDER}\"
      echo \"Start    : \$(date)\"
      Rscript \"${TEMPLATE_DIR}/qc_batch_runner.R\" \
        --sample_sheet \"${SAMPLE_SHEET}\" \
        --output_dir   \"${OUTPUT_DIR}\" \
        --integration_only TRUE \
        --integration_level \"${INTEGRATION_LEVEL}\" \
        --use_cellbender \"${USE_CELLBENDER}\"
    ")
  echo "  Submitted: integration-only job  (job ${INTEGRATE_JOB_ID})"
  echo "----------------------------------------------"
  echo "All jobs submitted. Monitor with: squeue -u \$USER"
  exit 0
fi

if [[ "${MERGE_ONLY}" == "TRUE" ]]; then
  MERGE_ONLY_JOB_ID=$(sbatch \
    --job-name="QC_merge_only" \
    --partition=cpu_short \
    --ntasks=1 \
    --cpus-per-task=1 \
    --mem="${MERGE_MEM_GB}G" \
    --time="${WALL_TIME}" \
    --output="${OUTPUT_DIR}/logs/QC_merge_only_%j.log" \
    --error="${OUTPUT_DIR}/logs/QC_merge_only_%j.log" \
    --parsable \
    --wrap="
      set -euo pipefail
      if [[ -n \"${CONDA_INIT}\" && -f \"${CONDA_INIT}\" ]]; then
        set +u
        source /gpfs/data/cvrcbioinfolab/gildem01/conda_envs/anaconda3/condaload_r.sh
        set -u
      else
        echo \"ERROR: CONDA_INIT not found (${CONDA_INIT}). Refusing to run without the conda R environment.\" >&2
        exit 1
      fi
      echo \"Job ID   : \${SLURM_JOB_ID}\"
      echo \"Node     : \${SLURMD_NODENAME}\"
      echo \"Task     : merge only\"
      echo \"Start    : \$(date)\"
      Rscript \"${TEMPLATE_DIR}/qc_batch_runner.R\" \
        --sample_sheet \"${SAMPLE_SHEET}\" \
        --output_dir   \"${OUTPUT_DIR}\" \
        --merge_only   TRUE \
        --run_integration \"${RUN_INTEGRATION}\" \
        --integration_level \"${INTEGRATION_LEVEL}\"
      echo \"Finished : \$(date)\"
    ")
  echo "  Submitted: merge-only job  (job ${MERGE_ONLY_JOB_ID})"
  echo "----------------------------------------------"
  echo "All jobs submitted. Monitor with: squeue -u \$USER"
  exit 0
fi

# ── Submit one job per sample ─────────────────────────────────────────────────
SAMPLE_JOB_IDS=()

for SAMPLE in "${SAMPLES[@]}"; do
  JOB_ID=$(sbatch \
    --job-name="QC_${SAMPLE}" \
    --partition=cpu_short \
    --ntasks=1 \
    --cpus-per-task=1 \
    --mem="${MEM_GB}G" \
    --time="${WALL_TIME}" \
    --output="${OUTPUT_DIR}/logs/QC_${SAMPLE}_%j.log" \
    --error="${OUTPUT_DIR}/logs/QC_${SAMPLE}_%j.log" \
    --parsable \
    --wrap="
      set -euo pipefail
      if [[ -n \"${CONDA_INIT}\" && -f \"${CONDA_INIT}\" ]]; then
        # Some environment init scripts assume vars like PYTHONPATH may be unset.
        # Temporarily relax nounset while sourcing, then restore strict mode.
        set +u
        source /gpfs/data/cvrcbioinfolab/gildem01/conda_envs/anaconda3/condaload_r.sh
        set -u
      else
        echo \"ERROR: CONDA_INIT not found (${CONDA_INIT}). Refusing to run without the conda R environment.\" >&2
        exit 1
      fi
      echo \"Job ID   : \${SLURM_JOB_ID}\"
      echo \"Node     : \${SLURMD_NODENAME}\"
      echo \"Sample   : ${SAMPLE}\"
      echo \"use_cellbender: ${USE_CELLBENDER}\"
      echo \"Start    : \$(date)\"
      Rscript \"${TEMPLATE_DIR}/qc_batch_runner.R\" \
        --template     \"${TEMPLATE_DIR}/sample_QC.Rmd\" \
        --sample_sheet \"${SAMPLE_SHEET}\" \
        --sample_name  \"${SAMPLE}\" \
        --output_dir   \"${OUTPUT_DIR}\" \
        --use_cellbender \"${USE_CELLBENDER}\"
      echo \"Finished : \$(date)\"
    ")
  echo "  Submitted: ${SAMPLE}  (job ${JOB_ID})"
  SAMPLE_JOB_IDS+=("${JOB_ID}")
done

# ── Submit merge job only after sample QC outputs are present ────────────────
if [[ ${#SAMPLES[@]} -gt 1 && "${SKIP_MERGE}" != "TRUE" ]]; then
  DEPENDENCY="afterok:$(IFS=:; echo "${SAMPLE_JOB_IDS[*]}")"

  MERGE_SUBMITTER_SCRIPT="${OUTPUT_DIR}/run_metadata/merge_submitter_${INVOCATION_ID}.sh"
  MERGE_JOB_SCRIPT="${OUTPUT_DIR}/run_metadata/merge_job_${INVOCATION_ID}.sh"
  INTEGRATION_JOB_SCRIPT="${OUTPUT_DIR}/run_metadata/integration_job_${INVOCATION_ID}.sh"

  SAMPLE_CHECK_LINES=""
  for SAMPLE in "${SAMPLES[@]}"; do
    SAMPLE_LABEL=$(printf '%q' "${SAMPLE}")
    SAMPLE_RDS=$(printf '%q' "${OUTPUT_DIR}/${SAMPLE}_QC.rds")
    SAMPLE_HTML=$(printf '%q' "${OUTPUT_DIR}/${SAMPLE}_QC.html")
    SAMPLE_CHECK_LINES+=$'      if [[ ! -f '"${SAMPLE_RDS}"' || ! -f '"${SAMPLE_HTML}"' ]]; then\n'
    SAMPLE_CHECK_LINES+=$'        echo "ERROR: expected QC outputs missing for sample '"${SAMPLE_LABEL}"'" >&2\n'
    SAMPLE_CHECK_LINES+=$'        exit 1\n'
    SAMPLE_CHECK_LINES+=$'      fi\n'
  done

  cat > "${INTEGRATION_JOB_SCRIPT}" <<EOF
#!/bin/bash
set -euo pipefail
if [[ -n "${CONDA_INIT}" && -f "${CONDA_INIT}" ]]; then
  set +u
  source /gpfs/data/cvrcbioinfolab/gildem01/conda_envs/anaconda3/condaload_r.sh
  set -u
else
  echo "ERROR: CONDA_INIT not found (${CONDA_INIT}). Refusing to run without the conda R environment." >&2
  exit 1
fi
echo "Job ID   : \${SLURM_JOB_ID}"
echo "Node     : \${SLURMD_NODENAME}"
echo "Task     : run integration"
echo "use_cellbender: ${USE_CELLBENDER}"
echo "Start    : \$(date)"
Rscript "${TEMPLATE_DIR}/qc_batch_runner.R" \
  --sample_sheet "${SAMPLE_SHEET}" \
  --output_dir   "${OUTPUT_DIR}" \
  --integration_only TRUE \
  --integration_level "${INTEGRATION_LEVEL}"
echo "Finished : \$(date)"
EOF
  chmod +x "${INTEGRATION_JOB_SCRIPT}"

  cat > "${MERGE_JOB_SCRIPT}" <<EOF
#!/bin/bash
set -euo pipefail
if [[ -n "${CONDA_INIT}" && -f "${CONDA_INIT}" ]]; then
  set +u
  source /gpfs/data/cvrcbioinfolab/gildem01/conda_envs/anaconda3/condaload_r.sh
  set -u
else
  echo "ERROR: CONDA_INIT not found (${CONDA_INIT}). Refusing to run without the conda R environment." >&2
  exit 1
fi
echo "Job ID   : \${SLURM_JOB_ID}"
echo "Node     : \${SLURMD_NODENAME}"
echo "Task     : merge all samples"
echo "Start    : \$(date)"
Rscript "${TEMPLATE_DIR}/qc_batch_runner.R" \
  --sample_sheet "${SAMPLE_SHEET}" \
  --output_dir   "${OUTPUT_DIR}" \
  --merge_only   TRUE
echo "Finished : \$(date)"

MERGED_QC_RDS="${OUTPUT_DIR}/merged_QC.rds"
MERGE_HTML="${OUTPUT_DIR}/merge_analysis.html"
if [[ ! -f "${MERGED_QC_RDS}" ]]; then
  echo "ERROR: merged output missing after merge job: ${MERGED_QC_RDS}" >&2
  exit 1
fi
if [[ -f "${TEMPLATE_DIR}/merge_analysis.Rmd" && ! -f "${MERGE_HTML}" ]]; then
  echo "ERROR: merge report missing after merge job: ${MERGE_HTML}" >&2
  exit 1
fi

if [[ "${RUN_INTEGRATION}" == "TRUE" ]]; then
  INTEGRATE_JOB_ID=\$(sbatch \
    --job-name="QC_integrate" \
    --partition=cpu_short \
    --ntasks=1 \
    --cpus-per-task=1 \
    --mem="${INTEGRATION_MEM_GB}G" \
    --time="${WALL_TIME}" \
    --output="${OUTPUT_DIR}/logs/QC_integrate_%j.log" \
    --error="${OUTPUT_DIR}/logs/QC_integrate_%j.log" \
    --parsable \
    "${INTEGRATION_JOB_SCRIPT}")
  echo "  Submitted: integration job  (job ${INTEGRATE_JOB_ID})"
fi
EOF
  chmod +x "${MERGE_JOB_SCRIPT}"

  cat > "${MERGE_SUBMITTER_SCRIPT}" <<EOF
#!/bin/bash
set -euo pipefail
if [[ -n "${CONDA_INIT}" && -f "${CONDA_INIT}" ]]; then
  set +u
  source /gpfs/data/cvrcbioinfolab/gildem01/conda_envs/anaconda3/condaload_r.sh
  set -u
else
  echo "ERROR: CONDA_INIT not found (${CONDA_INIT}). Refusing to run without the conda R environment." >&2
  exit 1
fi
${SAMPLE_CHECK_LINES}
echo "Submitting merge job..."
MERGE_JOB_ID=\$(sbatch \
  --job-name="QC_merge" \
  --partition=cpu_short \
  --ntasks=1 \
  --cpus-per-task=1 \
  --mem="${MERGE_MEM_GB}G" \
  --time="${WALL_TIME}" \
  --output="${OUTPUT_DIR}/logs/QC_merge_%j.log" \
  --error="${OUTPUT_DIR}/logs/QC_merge_%j.log" \
  --parsable \
  "${MERGE_JOB_SCRIPT}")
echo "  Submitted: merge job  (job ${MERGE_JOB_ID})"
EOF
  chmod +x "${MERGE_SUBMITTER_SCRIPT}"

  MERGE_SUBMITTER_JOB_ID=$(sbatch \
    --job-name="QC_merge_submit" \
    --partition=cpu_short \
    --ntasks=1 \
    --cpus-per-task=1 \
    --mem="1G" \
    --time="${WALL_TIME}" \
    --dependency="${DEPENDENCY}" \
    --output="${OUTPUT_DIR}/logs/QC_merge_submit_%j.log" \
    --error="${OUTPUT_DIR}/logs/QC_merge_submit_%j.log" \
    --parsable \
    "${MERGE_SUBMITTER_SCRIPT}")
  echo "  Submitted: merge submitter job  (job ${MERGE_SUBMITTER_JOB_ID}, depends on: ${SAMPLE_JOB_IDS[*]})"
  if [[ "${RUN_INTEGRATION}" == "TRUE" ]]; then
    echo "  Integration job will be submitted after the merge job confirms merged outputs."
  fi
elif [[ ${#SAMPLES[@]} -gt 1 ]]; then
  echo "  Skipping merge job submission because --skip_merge is TRUE"
fi

echo "----------------------------------------------"
echo "All jobs submitted. Monitor with: squeue -u \$USER"
