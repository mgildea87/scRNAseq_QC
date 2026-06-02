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
#   --outs_subdir <dir>   subdirectory inside cellranger_dir (default: outs)
#   --mem         <GB>    memory per job in GB              (default: 32)
#   --merge_mem   <GB>    memory for merge job in GB         (default: 64)
#   --integration_mem <GB> memory for integration job in GB   (default: 64)
#   --time        <HH:MM> wall time per job                 (default: 2:00:00)
#   --skip_merge  <TRUE/FALSE> skip submitting merge job     (default: FALSE)
#   --run_integration <TRUE/FALSE> render integrate_RNA.Rmd after merge (default: FALSE)
#   --integration_level <Batch|Sample> required when integration runs
#   --integration_only <TRUE/FALSE> run only integration from merged_QC.rds (default: FALSE)
#   --use_cellbender <TRUE/FALSE> use cellbender_filtered.h5 as filtered input (default: FALSE)
# =============================================================================

set -euo pipefail

# ── Fixed paths ───────────────────────────────────────────────────────────────
TEMPLATE_DIR="/gpfs/data/cvrcbioinfolab/gildem01/analysis_Rmd_templates_and_scripts/scRNAseq/Seurat_V5/QC"
CONDA_INIT="/gpfs/data/cvrcbioinfolab/gildem01/conda_envs/anaconda3/condaload_r.sh"

# ── Parse arguments ───────────────────────────────────────────────────────────
SAMPLE_SHEET=""
OUTPUT_DIR="QC"
OUTS_SUBDIR="outs"
MEM_GB=32
MERGE_MEM_GB=64
INTEGRATION_MEM_GB=64
WALL_TIME="2:00:00"
SKIP_MERGE="FALSE"
RUN_INTEGRATION="FALSE"
INTEGRATION_LEVEL=""
INTEGRATION_ONLY="FALSE"
USE_CELLBENDER="FALSE"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample_sheet) SAMPLE_SHEET="$2"; shift 2 ;;
    --output_dir)   OUTPUT_DIR="$2";   shift 2 ;;
    --outs_subdir)  OUTS_SUBDIR="$2";  shift 2 ;;
    --mem)          MEM_GB="$2";       shift 2 ;;
    --merge_mem)    MERGE_MEM_GB="$2"; shift 2 ;;
    --integration_mem) INTEGRATION_MEM_GB="$2"; shift 2 ;;
    --time)         WALL_TIME="$2";    shift 2 ;;
    --skip_merge)   SKIP_MERGE="$2";   shift 2 ;;
    --run_integration) RUN_INTEGRATION="$2"; shift 2 ;;
    --integration_level) INTEGRATION_LEVEL="$2"; shift 2 ;;
    --integration_only) INTEGRATION_ONLY="$2"; shift 2 ;;
    --use_cellbender) USE_CELLBENDER="$2"; shift 2 ;;
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

if [[ -z "${SAMPLE_SHEET}" ]]; then
  echo "ERROR: --sample_sheet is required." >&2
  exit 1
fi

if [[ ! -f "${SAMPLE_SHEET}" ]]; then
  echo "ERROR: sample sheet not found: ${SAMPLE_SHEET}" >&2
  exit 1
fi

# Resolve to absolute paths so jobs running in any working directory find them
SAMPLE_SHEET="$(realpath "${SAMPLE_SHEET}")"
OUTPUT_DIR="$(realpath -m "${OUTPUT_DIR}")"

# Ensure the main QC output directory exists even when using defaults
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${OUTPUT_DIR}/logs"

# ── Read sample names from CSV (skip header, column 1) ────────────────────────
mapfile -t SAMPLES < <(tail -n +2 "${SAMPLE_SHEET}" | cut -d',' -f1 | tr -d '\r')

if [[ ${#SAMPLES[@]} -eq 0 ]]; then
  echo "ERROR: no samples found in ${SAMPLE_SHEET}" >&2
  exit 1
fi

echo "Submitting ${#SAMPLES[@]} job(s)..."
echo "  Sample sheet : ${SAMPLE_SHEET}"
echo "  Output dir   : ${OUTPUT_DIR}"
echo "  Partition    : cpu_short  |  Per-sample Mem: ${MEM_GB}G  |  Merge Mem: ${MERGE_MEM_GB}G  |  Integration Mem: ${INTEGRATION_MEM_GB}G  |  Time: ${WALL_TIME}  |  Skip merge: ${SKIP_MERGE}  |  Run integration: ${RUN_INTEGRATION}  |  Integration level: ${INTEGRATION_LEVEL:-NA}  |  Integration only: ${INTEGRATION_ONLY}  |  Use CellBender: ${USE_CELLBENDER}"
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
      set +u
      source \"${CONDA_INIT}\"
      set -u
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
        --outs_subdir  \"${OUTS_SUBDIR}\" \
        --use_cellbender \"${USE_CELLBENDER}\"
    ")
  echo "  Submitted: integration-only job  (job ${INTEGRATE_JOB_ID})"
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
      # Some environment init scripts assume vars like PYTHONPATH may be unset.
      # Temporarily relax nounset while sourcing, then restore strict mode.
      set +u
      source \"${CONDA_INIT}\"
      set -u
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
        --outs_subdir  \"${OUTS_SUBDIR}\" \
        --use_cellbender \"${USE_CELLBENDER}\"
      echo \"Finished : \$(date)\"
    ")
  echo "  Submitted: ${SAMPLE}  (job ${JOB_ID})"
  SAMPLE_JOB_IDS+=("${JOB_ID}")
done

# ── Submit merge job (runs only if ALL sample jobs succeed) ───────────────────
if [[ ${#SAMPLES[@]} -gt 1 && "${SKIP_MERGE}" != "TRUE" ]]; then
  # Build colon-separated dependency string: afterok:id1:id2:...
  DEPENDENCY="afterok:$(IFS=:; echo "${SAMPLE_JOB_IDS[*]}")"

  MERGE_JOB_ID=$(sbatch \
    --job-name="QC_merge" \
    --partition=cpu_short \
    --ntasks=1 \
    --cpus-per-task=1 \
    --mem="${MERGE_MEM_GB}G" \
    --time="${WALL_TIME}" \
    --dependency="${DEPENDENCY}" \
    --output="${OUTPUT_DIR}/logs/QC_merge_%j.log" \
    --error="${OUTPUT_DIR}/logs/QC_merge_%j.log" \
    --parsable \
    --wrap="
      set -euo pipefail
      set +u
      source \"${CONDA_INIT}\"
      set -u
      echo \"Job ID   : \${SLURM_JOB_ID}\"
      echo \"Node     : \${SLURMD_NODENAME}\"
      echo \"Task     : merge all samples\"
      echo \"Start    : \$(date)\"
      Rscript \"${TEMPLATE_DIR}/qc_batch_runner.R\" \
        --sample_sheet \"${SAMPLE_SHEET}\" \
        --output_dir   \"${OUTPUT_DIR}\" \
        --outs_subdir  \"${OUTS_SUBDIR}\" \
        --merge_only   TRUE
      echo \"Finished : \$(date)\"
    ")
  echo "  Submitted: merge job  (job ${MERGE_JOB_ID}, depends on: ${SAMPLE_JOB_IDS[*]})"

  if [[ "${RUN_INTEGRATION}" == "TRUE" ]]; then
    INTEGRATE_JOB_ID=$(sbatch \
      --job-name="QC_integrate" \
      --partition=cpu_short \
      --ntasks=1 \
      --cpus-per-task=1 \
      --mem="${INTEGRATION_MEM_GB}G" \
      --time="${WALL_TIME}" \
      --dependency="afterok:${MERGE_JOB_ID}" \
      --output="${OUTPUT_DIR}/logs/QC_integrate_%j.log" \
      --error="${OUTPUT_DIR}/logs/QC_integrate_%j.log" \
      --parsable \
      --wrap="
        set -euo pipefail
        set +u
        source \"${CONDA_INIT}\"
        set -u
        echo \"Job ID   : \${SLURM_JOB_ID}\"
        echo \"Node     : \${SLURMD_NODENAME}\"
        echo \"Task     : run integration\"
        echo \"Start    : \$(date)\"
        Rscript \"${TEMPLATE_DIR}/qc_batch_runner.R\" \
          --sample_sheet \"${SAMPLE_SHEET}\" \
          --output_dir   \"${OUTPUT_DIR}\" \
          --integration_only TRUE \
          --integration_level \"${INTEGRATION_LEVEL}\"
        echo \"Finished : \$(date)\"
      ")
    echo "  Submitted: integration job  (job ${INTEGRATE_JOB_ID}, depends on merge job ${MERGE_JOB_ID})"
  fi
elif [[ ${#SAMPLES[@]} -gt 1 ]]; then
  echo "  Skipping merge job submission because --skip_merge is TRUE"
fi

echo "----------------------------------------------"
echo "All jobs submitted. Monitor with: squeue -u \$USER"
