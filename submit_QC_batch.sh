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
# OPTIONAL FLAGS (forwarded to run_QC_batch.R):
#   --outs_subdir <dir>   subdirectory inside cellranger_dir (default: outs)
#   --mem         <GB>    memory per job in GB              (default: 32)
#   --merge_mem   <GB>    memory for merge job in GB         (default: 64)
#   --time        <HH:MM> wall time per job                 (default: 2:00:00)
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
WALL_TIME="2:00:00"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample_sheet) SAMPLE_SHEET="$2"; shift 2 ;;
    --output_dir)   OUTPUT_DIR="$2";   shift 2 ;;
    --outs_subdir)  OUTS_SUBDIR="$2";  shift 2 ;;
    --mem)          MEM_GB="$2";       shift 2 ;;
    --merge_mem)    MERGE_MEM_GB="$2"; shift 2 ;;
    --time)         WALL_TIME="$2";    shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
done

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
echo "  Partition    : cpu_short  |  Per-sample Mem: ${MEM_GB}G  |  Merge Mem: ${MERGE_MEM_GB}G  |  Time: ${WALL_TIME}"
echo "----------------------------------------------"

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
      echo \"Start    : \$(date)\"
      Rscript \"${TEMPLATE_DIR}/run_QC_batch.R\" \
        --template     \"${TEMPLATE_DIR}/sample_QC.Rmd\" \
        --sample_sheet \"${SAMPLE_SHEET}\" \
        --sample_name  \"${SAMPLE}\" \
        --output_dir   \"${OUTPUT_DIR}\" \
        --outs_subdir  \"${OUTS_SUBDIR}\"
      echo \"Finished : \$(date)\"
    ")
  echo "  Submitted: ${SAMPLE}  (job ${JOB_ID})"
  SAMPLE_JOB_IDS+=("${JOB_ID}")
done

# ── Submit merge job (runs only if ALL sample jobs succeed) ───────────────────
if [[ ${#SAMPLES[@]} -gt 1 ]]; then
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
      Rscript \"${TEMPLATE_DIR}/run_QC_batch.R\" \
        --sample_sheet \"${SAMPLE_SHEET}\" \
        --output_dir   \"${OUTPUT_DIR}\" \
        --outs_subdir  \"${OUTS_SUBDIR}\" \
        --merge_only   TRUE
      echo \"Finished : \$(date)\"
    ")
  echo "  Submitted: merge job  (job ${MERGE_JOB_ID}, depends on: ${SAMPLE_JOB_IDS[*]})"
fi

echo "----------------------------------------------"
echo "All jobs submitted. Monitor with: squeue -u \$USER"
