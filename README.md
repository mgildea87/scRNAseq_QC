# 10x scRNA-seq QC and processing pipeline

Generates a self-contained HTML QC report and a filtered Seurat RDS object for
each sample processed by CellRanger or similar.  
Can be run on a single sample interactively or batched across many samples via
a sample sheet.

GitHub repository: https://github.com/mgildea87/scRNAseq_QC

---

## QC workflow summary

### `sample_QC.Rmd` (per-sample QC)

- Load filtered and raw count matrices (10x MEX folders or .h5 files)
- Barcode rank (knee) plot
- Remove genes detected in fewer than 0.1% of barcodes
- Compute per-cell QC metrics: UMI count, genes detected, % mitochondrial, % haemoglobin
- Build MAD-based threshold tables for each metric
- Pre-filter QC plots: violins, histograms, scatter pairs
- MALAT1 scatter (proxy for cell viability)
- Apply filters (hard thresholds from params, or MAD defaults)
- Post-filter QC plots
- Normalise → variable features → PCA → clustering → UMAP
- Mean–variance plot and SCTransform model assessment
- Per-cluster QC metric distributions
- Top marker heatmap (Wilcoxon, top 20 per cluster)
- Doublet detection with `scDblFinder` (annotated, not removed)
- Save filtered + annotated Seurat object to RDS

### `merge_analysis.Rmd` (optional post-merge analysis)

- Load `merged_QC.rds`, detect species from feature naming, and require the merged metadata fields `batch` and `sample_name`
- Build RNA merged reductions/clusters when missing: normalize, variable features, scale, PCA, neighbors, clustering, UMAP
- Build SCT merged reductions/clusters when missing: PCA, neighbors, clustering, UMAP on the SCT assay
- Save updated merged object back to `merged_QC.rds`
- Generate RNA and SCT dimensional reduction plots by batch/sample and by merged cluster labels
- Compute and plot cluster abundance distributions across samples and batches (RNA and SCT cluster assignments)
- Compute LISI scores for batch/sample_name in RNA and SCT PCA spaces, save to `lisi.rds`, and visualize histograms/FeaturePlots/violins
- Plot QC metric distributions by cluster (for example `nCount_RNA`, `nFeature_RNA`, `percent.mt`, cell-cycle scores)
- Run marker discovery for merged RNA and SCT clusters and render top-marker heatmaps
- Render merged analysis report (`merge_analysis.html`)

### `integrate_RNA.Rmd` (optional post-merge integration)

- Load `merged_QC.rds` and validate required integration grouping metadata (`Batch` or `Sample` level)
- RNA integration workflow: split assay layers by integration group, normalize, find variable features, scale, PCA, RPCA integration, neighbors, clustering, UMAP
- SCT integration workflow: apply `SCTransform`, split SCT layers by integration group, PCA, RPCA integration, neighbors, clustering, UMAP
- Save the integrated object (`integrated.rds`, containing RNA, SCT, RPCA and Harmony results)
- Plot integrated UMAPs colored by sample, batch, and integrated cluster IDs (RNA and SCT)
- Compute and visualize cluster abundance distributions across samples and batches (RNA and SCT integrated clusters)
- Compute LISI batch-mixing scores for integrated embeddings (RNA and SCT), save RDS outputs, and visualize distributions
- Plot QC metric distributions in integrated objects (for example `nCount_RNA`, `nFeature_RNA`, `percent.mt`, `Malat1`)
- Run marker discovery on RNA/SCT integrated clusters and render top-marker heatmaps
- Render integration report (`integrate_RNA.html`)

---

## Files

| File | Purpose |
|------|---------|
| `sample_QC.Rmd` | Parameterised R Markdown template for per-sample QC |
| `merge_analysis.Rmd` | R Markdown template for post-merge analysis across samples |
| `integrate_RNA.Rmd` | Optional post-merge RNA integration report template |
| `qc_batch_runner.R` | Main batch runner — renders `sample_QC.Rmd` for each sample in the sample sheet |
| `run_qc_batch_local.sh` | Local shell launcher — loads the CVRC R conda environment, then runs `qc_batch_runner.R` |
| `submit_QC_batch.sh` | Cluster launcher — submits one SLURM job per sample (plus a dependent merge job) |
| `samples.csv` | Sample sheet — edit this to point to your data |
| `support_files/` | Bundled TF and haemoglobin reference files used by the templates |

## Requirements

The following R packages must be installed:

```r
install.packages(c("rmarkdown", "ggplot2", "patchwork", "viridis",
                   "pheatmap", "dplyr", "tidyr", "reshape2", "Hmisc", "knitr"))

BiocManager::install(c("Seurat", "scater", "scDblFinder"))

# For graph modularity diagnostics (pairwiseModularity)
BiocManager::install("bluster")

# For optional hashtag (HTO) demultiplexing / CITE-seq (ADT) support
install.packages(c("dsb", "gridExtra"))

# From GitHub / internal
remotes::install_github("immunogenomics/presto")
# CVRCFunc — install from internal source as needed
```

For local development tests, also install `testthat`:

```r
install.packages("testthat")
```

Run the focused local suite from the repository root with
`Rscript --vanilla tests/testthat.R`. Full CellRanger report renders and SLURM
submission are manual HPC checks and are not part of this fast local suite.
See [tests/README.md](tests/README.md) for per-file coverage, prerequisites,
and the manual cluster smoke checks.

---

## Sample sheet (`samples.csv`)

The sample sheet is a transposed CSV: the first column is `field`, each
remaining column is one sample, and each following row is a setting. This
keeps each sample's settings in a vertical spreadsheet column. The `sample_name`
row must exactly match the corresponding sample column header.

| Field row | Required | Description |
|-----------|----------|-------------|
| `sample_name` | Yes | Must match the sample column header; used for output filenames |
| `filtered_path` | Yes | Path to the filtered matrix directory or `.h5` file for that sample |
| `raw_path` | Yes | Path to the raw matrix directory or `.h5` file for that sample |
| `batch` | No | Batch label in Seurat metadata (defaults to `A` if omitted or blank) |
| `min_nCount_RNA` / `max_nCount_RNA` | No | Lower/upper RNA UMI bounds; `NA` uses MAD minimum or no upper cap |
| `min_nFeature_RNA` / `max_nFeature_RNA` | No | Lower/upper detected-gene bounds; `NA` uses MAD minimum or no upper cap |
| `max_percent_mt` | No | Upper mitochondrial fraction bound; `NA` uses median + 5 × MAD |
| `max_percent_rbc` | No | Upper haemoglobin-gene percentage bound; `NA` applies no cap |
| `min_malat1` | No | Minimum normalised MALAT1 expression (default 1) |
| `use_hashtag` | No | `TRUE` enables HTO loading and demultiplexing for this sample (default `FALSE`) |
| `use_adt` | No | `TRUE` enables ADT loading and DSB normalization for this sample (default `FALSE`) |
| `feature_reference_path` | Required when HTO or ADT enabled with "name" column specifying feature names | CellRanger `feature_reference.csv`; required even when separate HTO/ADT paths are provided |
| `hto_path` | No | Direct path to a separate filtered HTO matrix (10x MEX directory or `.h5`); takes priority over matrix splitting |
| `adt_path` | No | Direct path to a separate filtered ADT matrix; takes priority over matrix splitting |
| `raw_adt_path` | No | Optional direct path override for the raw ADT matrix used for DSB empty-droplet background; if blank, ADT is resolved from the `Antibody Capture` matrix in `raw_path`; barcodes must match filtered ADT |
| `hto_features` | For combined HTO/ADT matrix | Exact HTO names from the reference `name` column, separated by `;` or `,` |
| `min_nCount_ADT` / `max_nCount_ADT` | No | Lower/upper ADT UMI bounds; `NA` uses MAD minimum or no upper cap; only applies when ADT is enabled |
| `dsb_background_rna_max` | No | Override for maximum log10(RNA UMI) among empty-droplet background barcodes |
| `dsb_background_prot_min` / `dsb_background_prot_max` | No | Overrides for the log10(ADT UMI) band defining empty-droplet background barcodes |

Set a threshold to `NA` (or omit its field row) to use the automatic default:

- `min_nCount_RNA`, `min_nFeature_RNA`, and `min_nCount_ADT`: median − 4 × MAD
- `max_nCount_RNA`, `max_nFeature_RNA`, and `max_nCount_ADT`: no upper cap
- `max_percent_mt`: median + 5 × MAD
- `max_percent_rbc`: no upper cap
- `min_malat1`: 1

### Example

```csv
field,Control_1,Treatment_1
sample_name,Control_1,Treatment_1
filtered_path,/path/to/Control/filtered_feature_bc_matrix,/path/to/Treatment/filtered_feature_bc_matrix
raw_path,/path/to/Control/raw_feature_bc_matrix,/path/to/Treatment/raw_feature_bc_matrix
batch,A,A
min_nCount_RNA,NA,500
max_nCount_RNA,NA,25000
min_nFeature_RNA,NA,250
max_nFeature_RNA,NA,6000
max_percent_mt,NA,20
max_percent_rbc,NA,NA
min_malat1,1,1
use_hashtag,FALSE,FALSE
use_adt,FALSE,FALSE
min_nCount_ADT,NA,NA
max_nCount_ADT,NA,NA
```

`Control_1` uses automatic MAD thresholds. `Treatment_1` uses explicit RNA
thresholds. In spreadsheet software, add a sample by adding a new column and
filling in its field values.

For `.h5` inputs, `filtered_path` and `raw_path` must point to the file itself, not the folder that contains it. For 10x MEX inputs, the path must point to the matrix directory that contains `matrix.mtx(.gz)`, `barcodes.tsv(.gz)`, and `features.tsv(.gz)` or `genes.tsv(.gz)`.

Matrix inputs are no longer auto-detected from a parent directory. Each sample must provide direct paths to the filtered and raw matrix inputs.
Supported input formats:

- 10x MEX directories containing `matrix.mtx(.gz)`, `barcodes.tsv(.gz)`, and
  `features.tsv(.gz)` or `genes.tsv(.gz)`
- 10x-style `.h5` or `.hdf5` files

For `.h5` inputs, the QC template currently expects one scRNA matrix per file.
CellBender filtered input is opt-in via `--use_cellbender TRUE`.

---

## Hashtag demultiplexing / CITE-seq (ADT)

Both features are opt-in per sample (`use_hashtag`, `use_adt` in `samples.csv`)
and disabled by default, leaving the RNA-only analysis path unchanged. The
sample-sheet format is now transposed; convert older row-per-sample CSV sheets
before running this version.

### Input formats

HTO and ADT data usually arrive as one CellRanger `multi`/feature-barcoded
`filtered_feature_bc_matrix` / `raw_feature_bc_matrix`, where `Read10X()`
returns a named list keyed by feature type (`Gene Expression`, `Antibody
Capture`, sometimes `Multiplexing Capture`). HTO and ADT features are
frequently combined inside one `Antibody Capture` matrix. When either feature
is enabled, `feature_reference_path` is required and loaded before matrix
resolution. For a combined matrix, the reference's `name` column is
canonical and must match matrix feature rownames: `hto_features` lists exact
HTO names, and all remaining reference names are treated as ADT. Missing
reference names or HTO names are errors; `id` and `feature_type` are not used
to classify antibody features.

When HTO counts are in a separate named matrix within the filtered input,
`Multiplexing Capture` is recognized automatically (case-insensitively). For
HTO and ADT rows combined in `Antibody Capture`, set `hto_features` to the
exact HTO names from the feature reference's `name` column; remaining
reference names are treated as ADT. Other feature-type names (for example,
`hto_matrix`) are not automatically searched. If the HTO matrix is a separate
file, provide it through `hto_path`.

Separate matrices remain supported: `hto_path` and `adt_path`, when set,
override resolution for the filtered assay and may point to a standalone
single-modality MEX or `.h5` file. CellRanger's native `Multiplexing Capture`
matrix can also supply HTO counts. When ADT is enabled, `raw_adt_path` may be
left blank to resolve ADT from an `Antibody Capture` matrix in the raw input
using the same feature-reference rules, or set to a direct raw ADT matrix path
to override resolution.
Filtered and raw ADT barcodes must match for DSB. The `feature_reference_path`
is still required when HTO or ADT processing is enabled, including when
separate paths are used.

### Hashtag (HTO) demultiplexing

When `use_hashtag` is TRUE, `sample_QC.Rmd` adds a "Hashtag demultiplexing"
section that:

- CLR-normalizes the `HTO` assay and runs `HTODemux` with a fixed
  `positive.quantile` of 0.99.
- Adds `HTO_sample` (= `hash.ID`), `HTO_class` (= classification global call),
  and `replicate_id` (= `<sample_name>_<hash.ID>`) before filtering; the saved
  object contains Singlet cells only.
- Shows diagnostic QC plots (`HTOHeatmap`, `RidgePlot`, per-HTO density plots,
  and all-cell CLR density histograms). Each histogram's red line marks the
  minimum CLR expression among cells assigned as Singlets for that HTO; it is
  diagnostic only and does not change the HTODemux classifications.
- Removes HTO `Negative` cells, runs `scDblFinder::recoverDoublets` on an RNA
  PCA/UMAP of the remaining Singlet+Doublet subset (requires
  `scDblFinder`/`scater`/`gridExtra`), then retains Singlets only before RNA
  filtering/clustering. The intra-sample doublet prediction is retained as
  metadata on the saved Singlet cells.

**`sample_name` vs `library_id` semantics:** for hashtagged samples,
`qc_batch_runner.R`'s merge step keeps `library_id` = the pooled library name
from `samples.csv`, while `sample_name` in the merged object becomes the
per-cell `HTO_sample` (biological identity from the hashtag call) instead of
being overwritten with the library name. Non-hashtagged samples are
unaffected (`sample_name` == `library_id` == the sample-sheet name).

### CITE-seq (ADT) processing

When `use_adt` is TRUE, ADT count QC runs with the pre-filter RNA metrics.
The separate "CITE-seq (ADT) processing" report section appears after RNA QC
filtering, once RNA/SCT clustering and UMAPs exist, and:

- Includes `min_nCount_ADT`/`max_nCount_ADT` in the pre-filter summary, MAD
  threshold table and QC plots, then applies them with the RNA filters. `NA`
  minimum uses median − 4 × MAD; `NA` maximum means no upper cap.
- Computes a biaxial RNA-size vs ADT-size diagnostic (colored by
  mitochondrial fraction, faceted by cell/background droplet class) to help
  judge whether the default empty-droplet background selection is
  appropriate; override with `dsb_background_rna_max`,
  `dsb_background_prot_min`/`dsb_background_prot_max` if not.
- Runs `dsb::DSBNormalizeProtein()` using the raw matrix's empty droplets as
  background, with isotype controls auto-detected by name pattern
  (`Isotype`, case-insensitive) when present.
- Stores DSB-normalized values in a **separate `ADT_DSB` assay** (`data` slot
  only) — the raw `ADT` counts assay is left untouched.
- Shows per-ADT violin/ridge plots, an ADT-vs-RNA UMI scatter, `FeaturePlot`
  of ADT markers on the existing RNA UMAP, and a DSB-normalized
  average-ADT-per-cluster heatmap.

Both `HTO` and `ADT_DSB` assays, and all associated metadata, propagate
automatically through `merge()`/`SCTransform`/RPCA — `merge_analysis.Rmd` and
`integrate_RNA.Rmd` render small additive, presence-gated summary sections
for them and are otherwise unaffected for pipelines that don't use these
features.

> **Tip:** Run the pipeline on each sample with default thresholds first,
> inspect the QC reports, then rerun with hardcoded thresholds for any samples
> that need manual adjustment.

---

## Running the pipeline

### Choose a run mode

- `run_qc_batch_local.sh`: many samples on one machine from a shell (no SLURM)
- `sample_QC.Rmd` in RStudio/R console: one sample interactively
- `submit_QC_batch.sh`: many samples on the cluster via SLURM (recommended at scale)

### Many samples on one machine (terminal, no SLURM)

Use this when running directly in a shell and you want environment setup handled
for you automatically.

By default, this mode is **sequential** (`--ncores 1`), so samples are
processed one at a time. To parallelise per-sample QC on a single machine,
set `--ncores` to a value greater than 1.

When running the full sample sheet with more than one sample, merge is
performed automatically after all samples complete successfully. To skip this,
set `--skip_merge TRUE`.

To run the optional post-merge integration report template (`integrate_RNA.Rmd`)
after a successful merge, set `--run_integration TRUE`. The integration level
defaults to `Sample`; pass `--integration_level Batch` to use batch grouping.

```bash
bash /path/to/templates/QC/run_qc_batch_local.sh \
  --sample_sheet /abs/path/to/samples.csv \
  --output_dir   /abs/path/to/results/QC
```

```bash
# Example: process samples in parallel on 8 cores
bash /path/to/templates/QC/run_qc_batch_local.sh \
  --sample_sheet /abs/path/to/samples.csv \
  --output_dir   /abs/path/to/results/QC \
  --ncores       8

# Example: render all samples but skip merged_QC outputs
bash /path/to/templates/QC/run_qc_batch_local.sh \
  --sample_sheet /abs/path/to/samples.csv \
  --output_dir   /abs/path/to/results/QC \
  --skip_merge   TRUE

# Example: run optional integration report after merge
bash /path/to/templates/QC/run_qc_batch_local.sh \
  --sample_sheet     /abs/path/to/samples.csv \
  --output_dir       /abs/path/to/results/QC \
  --run_integration  TRUE \
  --integration_level Batch

# Example: run only integration from an existing merged_QC.rds
bash /path/to/templates/QC/run_qc_batch_local.sh \
  --sample_sheet      /abs/path/to/samples.csv \
  --output_dir        /abs/path/to/results/QC \
  --integration_only  TRUE
```

### One sample manually (RStudio or R console)

Open `sample_QC.Rmd` and knit with custom parameters, or run from
the R console:

When using `rmarkdown::render(..., params = list(...))`, parameter values are
resolved as follows:

- Any parameter supplied in `params = list(...)` **overrides** the value in the
  YAML `params:` block of `sample_QC.Rmd`.
- Any parameter not supplied in `params = list(...)` falls back to the default
  value defined in the YAML `params:` block.

In other words, console-supplied params take precedence, and YAML params act as
defaults.

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

# Example: .h5 inputs (including CellBender outputs). Pass the .h5 file path itself.
rmarkdown::render(
  "sample_QC.Rmd",
  params = list(
    sample_name   = "Control_1",
    filtered_path = "/path/to/cellbender_filtered.h5",
    raw_path      = "/path/to/raw_feature_bc_matrix.h5",
    output_dir    = "results/QC"
  )
)
```

### Many samples on the cluster (SLURM, `submit_QC_batch.sh`)

The recommended way to run across many samples. The coordinator script reads
your sample sheet and submits **one independent SLURM job per sample** on the
`cpu_short` partition. Run it with plain `bash` — do not use `sbatch`.

```bash
bash /path/to/templates/QC/submit_QC_batch.sh \
  --sample_sheet /abs/path/to/samples.csv \
  --output_dir   /abs/path/to/results/QC
```

For sample sheets with more than one sample, this mode submits a dependent
merge job by default. To skip that merge job, pass `--skip_merge TRUE`.

To skip all per-sample QC jobs and run only merge/report generation from
existing `<sample_name>_QC.rds` files, pass `--merge_only TRUE`.

To run the optional integration report after merge in the merge job, pass
`--run_integration TRUE`. The integration level defaults to `Sample`; pass
`--integration_level Batch` to use batch grouping.

This also works with `--merge_only TRUE`: integration will run after merge,
but only when merge succeeds and the sample sheet contains at least 2 samples.

To run only the integration stage on an existing merged object, pass
`--integration_only TRUE`; the integration level defaults to `Sample`.

> Use **absolute paths** for `--sample_sheet` and `--output_dir` — each sample
> runs as a separate job on a compute node where relative paths may not resolve.

Each job writes its log to `<output_dir>/logs/QC_<sample_name>_<jobid>.log`.

#### Options

| Flag | Default | Description |
|------|---------|-------------|
| `--sample_sheet` | *(required)* | Absolute path to CSV sample sheet |
| `--output_dir` | `QC` | Directory for HTML reports, RDS files, and logs |
| `--use_cellbender` | `FALSE` | If `TRUE`, filtered input must resolve to `cellbender_filtered.h5` (or `.hdf5`) |
| `--mem` | `32` | Memory per job in GB |
| `--merge_mem` | `64` | Memory for the merge job in GB |
| `--integration_mem` | `64` | Memory for the integration job in GB |
| `--merge_only` | `FALSE` | Skip per-sample QC jobs and submit only merge/report generation from existing per-sample RDS files |
| `--skip_merge` | `FALSE` | Skip the post-sample merge step (`merged_QC.rds` and `merge_analysis.html`) |
| `--run_integration` | `FALSE` | Run optional post-merge integration report (`integrate_RNA.Rmd`); works with `--merge_only TRUE` after successful merge (requires at least 2 samples) |
| `--integration_level` | `Sample` | Integration grouping level: `Batch` or `Sample` |
| `--integration_only` | `FALSE` | Skip sample QC and merge; run only `integrate_RNA.Rmd` using existing `merged_QC.rds` |
| `--node_type` | `cpu_short` | SLURM partition to request for all jobs |
| `--time` | `2:00:00` | Wall time per job — max on `cpu_short` is `12:00:00` |

#### Flag compatibility

Valid combinations:

- `--merge_only TRUE --run_integration TRUE [--integration_level Batch|Sample]`: runs merge from existing per-sample RDS files, then runs integration if merge succeeds (requires at least 2 samples in the sample sheet); the level defaults to `Sample`.
- `--merge_only TRUE`: runs merge/report generation only from existing per-sample RDS files.
- `--integration_only TRUE [--integration_level Batch|Sample]`: skips sample QC and merge, runs integration from an existing `merged_QC.rds`; the level defaults to `Sample`.
- `--skip_merge TRUE`: runs per-sample QC jobs only and does not submit a merge job.

Invalid combinations:

- `--merge_only TRUE --integration_only TRUE`: mutually exclusive.
- `--merge_only TRUE --skip_merge TRUE`: mutually exclusive.
- `--run_integration TRUE --skip_merge TRUE` (without `--integration_only TRUE`): integration requires merge.

```bash
# Example: larger samples needing more memory
bash /path/to/submit_QC_batch.sh \
  --sample_sheet /abs/path/to/samples.csv \
  --output_dir   /abs/path/to/results/QC \
  --use_cellbender TRUE \
  --mem          32 \
  --merge_mem    96 \
  --integration_mem 128

# Example: skip per-sample QC and run only merge_analysis/integration on existing sample RDS files
bash /path/to/submit_QC_batch.sh \
  --sample_sheet /abs/path/to/samples.csv \
  --output_dir   /abs/path/to/results/QC \
  --merge_only   TRUE \
  --run_integration TRUE \
  --integration_level Batch
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

When `--run_integration TRUE` is used and merge succeeds, the integration step can also write:

| File | Description |
|------|-------------|
| `<output_dir>/integrate_RNA.html` | Integration report rendered from `integrate_RNA.Rmd` |
| `<output_dir>/integrated.rds` | Integrated Seurat object if produced by `integrate_RNA.Rmd` |

---

## Reproducibility logging (automatic)

Reproducibility metadata is logged automatically every time you run the pipeline
using `run_qc_batch_local.sh` or `submit_QC_batch.sh`.

Use your normal run command, for example:

```bash
bash run_qc_batch_local.sh \
  --sample_sheet /abs/path/to/projects/Project_001/samples.csv \
  --output_dir /abs/path/to/projects/Project_001/results/QC
```

Metadata is written to `<output_dir>/run_metadata/` and includes:

- Git release/source info: repository URL, commit SHA, branch, exact tag (if any), describe string, worktree status
- Run command: an executable shell script with the exact invoked command
- Conda environment info: active conda environment name and `conda info --envs`
- Conda export: full environment export when conda is available and an env is active
- Run timing/status: run start and run end files with timestamp, invocation ID, and exit status

When using `submit_QC_batch.sh`, coordinator-level metadata is also logged in the
same `run_metadata` directory.