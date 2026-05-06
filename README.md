# Single-sample 10x scRNA-seq QC pipeline

Generates a self-contained HTML QC report and a filtered Seurat RDS object for
each sample processed by CellRanger.  
Can be run on a single sample interactively or batched across many samples via
a sample sheet.

---

## Files

| File | Purpose |
|------|---------|
| `sample_QC.Rmd` | Parameterised R Markdown template — do not edit paths directly |
| `run_QC_batch.R` | Launcher script — renders the template for every row in a sample sheet |
| `run_QC_batch_with_r.sh` | Wrapper that sources the CVRC R conda environment and runs `run_QC_batch.R` |
| `samples.csv` | Sample sheet — edit this to point to your data |

---

## Requirements

The following R packages must be installed:

```r
install.packages(c("rmarkdown", "ggplot2", "patchwork", "viridis",
                   "pheatmap", "dplyr", "tidyr", "reshape2", "Hmisc", "knitr"))

BiocManager::install(c("Seurat", "scater", "scDblFinder"))

# From GitHub / internal
remotes::install_github("immunogenomics/presto")
# CVRCFunc — install from internal source as needed
```

---

## Sample sheet (`samples.csv`)

The sample sheet is a CSV file with the following columns:

| Column | Required | Description |
|--------|----------|-------------|
| `sample_name` | Yes | Unique label used for output filenames |
| `cellranger_dir` | Yes | Path to the CellRanger output directory (usually the folder that **contains** the `outs/` subfolder) |
| `batch` | No | Batch label stored in the per-sample Seurat object metadata (defaults to `A` if omitted or blank) |
| `min_nCount_RNA` | No | Hard lower bound on UMI count per cell |
| `max_nCount_RNA` | No | Hard upper bound on UMI count per cell |
| `min_nFeature_RNA` | No | Hard lower bound on genes detected per cell |
| `max_nFeature_RNA` | No | Hard upper bound on genes detected per cell |
| `max_percent_mt` | No | Hard upper bound on mitochondrial read fraction (%) |
| `min_malat1` | No | Minimum normalised MALAT1 expression (default 1) |

Set any threshold column to `NA` (or omit the column entirely) to use the
automatic **MAD-based default** for that metric.  
The MAD defaults are:

- `min_nCount_RNA` / `min_nFeature_RNA` → median − 4 × MAD  
- `max_nCount_RNA` / `max_nFeature_RNA` → no upper cap (Inf)  
- `max_percent_mt` → median + 5 × MAD  
- `min_malat1` → 1  

### Example

```csv
sample_name,cellranger_dir,batch,min_nCount_RNA,max_nCount_RNA,min_nFeature_RNA,max_nFeature_RNA,max_percent_mt,min_malat1
Control_1,/path/to/cellranger/count-Control-1,A,NA,NA,NA,NA,NA,1
Treatment_1,/path/to/cellranger/count-Treatment-1,A,500,25000,250,6000,20,1
```

`Control_1` uses all MAD-based defaults. `Treatment_1` uses hardcoded thresholds.

Matrix folder names are auto-detected from the chosen base directory as long as
they contain valid 10x MEX files (`matrix.mtx(.gz)`, `barcodes.tsv(.gz)`, and
`features.tsv(.gz)` or `genes.tsv(.gz)`).

> **Tip:** Run the pipeline on each sample with default thresholds first,
> inspect the QC reports, then rerun with hardcoded thresholds for any samples
> that need manual adjustment.

---

## Running the pipeline

### Multiple samples (direct shell, non-SLURM)

Use the wrapper script when you want to run directly from a shell while
ensuring the expected R environment is loaded first.

```bash
bash /path/to/templates/QC/run_QC_batch_with_r.sh \
  --sample_sheet /abs/path/to/samples.csv \
  --output_dir   /abs/path/to/results/QC
```

### Single sample (interactive / RStudio)

Open `sample_QC.Rmd` and knit with custom parameters, or run from
the R console:

```r
rmarkdown::render(
  "sample_QC.Rmd",
  params = list(
    sample_name   = "Control_1",
    filtered_path = "/path/to/outs/filtered_feature_bc_matrix",
    raw_path      = "/path/to/outs/raw_feature_bc_matrix",
    output_dir    = "results/QC"
  )
)
```

### Multiple samples — SLURM (`submit_QC_batch.sh`)

The recommended way to run across many samples. The coordinator script reads
your sample sheet and submits **one independent SLURM job per sample** on the
`cpu_short` partition. Run it with plain `bash` — do not use `sbatch`.

```bash
bash /path/to/templates/QC/submit_QC_batch.sh \
  --sample_sheet /abs/path/to/samples.csv \
  --output_dir   /abs/path/to/results/QC
```

> Use **absolute paths** for `--sample_sheet` and `--output_dir` — each sample
> runs as a separate job on a compute node where relative paths may not resolve.

Each job writes its log to `<output_dir>/logs/QC_<sample_name>_<jobid>.log`.

#### Options

| Flag | Default | Description |
|------|---------|-------------|
| `--sample_sheet` | *(required)* | Absolute path to CSV sample sheet |
| `--output_dir` | `QC` | Directory for HTML reports, RDS files, and logs |
| `--outs_subdir` | `outs` | Subdirectory inside `cellranger_dir` used as the search base for filtered/raw matrix folders. Set to `""` if `cellranger_dir` already points to the base directory containing those folders |
| `--mem` | `32` | Memory per job in GB |
| `--merge_mem` | `64` | Memory for the merge job in GB |
| `--time` | `12:00:00` | Wall time per job — max on `cpu_short` is `12:00:00` |

```bash
# Example: larger samples needing more memory
bash /path/to/submit_QC_batch.sh \
  --sample_sheet /abs/path/to/samples.csv \
  --output_dir   /abs/path/to/results/QC \
  --mem          32 \
  --merge_mem    96
```

Monitor submitted jobs with `squeue -u $USER`.

---

## Outputs

For each sample `<sample_name>`, the pipeline writes:

| File | Description |
|------|-------------|
| `<output_dir>/<sample_name>_QC.html` | Full interactive QC report |
| `<output_dir>/<sample_name>_QC.rds` | Filtered Seurat object with doublet scores, ready for integration |

When the sample sheet contains more than one sample, the merge step also writes:

| File | Description |
|------|-------------|
| `<output_dir>/merged_QC.rds` | Merged Seurat object produced from all `<sample_name>_QC.rds` files |
| `<output_dir>/merge_analysis.html` | Post-merge analysis report rendered from `merge_analysis.Rmd` |

---

## QC workflow summary

The Rmd template performs the following steps in order:

1. Load CellRanger filtered and raw count matrices
2. Barcode rank (knee) plot
3. Remove genes detected in fewer than 0.1% of barcodes
4. Compute per-cell QC metrics: UMI count, genes detected, % mitochondrial, % haemoglobin
5. Build MAD-based threshold tables for each metric
6. Pre-filter QC plots: violins, histograms, scatter pairs
7. MALAT1 scatter (proxy for cell viability)
8. Apply filters (hard thresholds from params, or MAD defaults)
9. Post-filter QC plots
10. Normalise → variable features → PCA → clustering → UMAP
11. Mean–variance plot and SCTransform model assessment
12. Per-cluster QC metric distributions
13. Top marker heatmap (Wilcoxon, top 20 per cluster)
14. Doublet detection with `scDblFinder` (annotated, not removed)
15. Save filtered + annotated Seurat object to RDS



To-do:
1. start a git repo
2. add integration step
3. use local pathas for the TF and hemoglobin gene files. These will eventually be pulled when the git repo is cloned into wherever the working directory for the project is. 