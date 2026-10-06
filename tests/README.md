# Test Suite

The fast local suite uses `testthat` and helper functions from
`support_files/qc_helpers.R`. It does not need real CellRanger data or submit
SLURM jobs. The synthetic matrix test creates temporary 10x MEX fixtures.

## Run

From the repository root:

```bash
Rscript --vanilla tests/testthat.R
```

In the CVRC environment, initialize R first:

```bash
unset R_HOME R_LIBS R_LIBS_USER R_PROFILE R_PROFILE_USER R_ENVIRON R_ENVIRON_USER
set +u
source /gpfs/data/cvrcbioinfolab/gildem01/conda_envs/anaconda3/condaload_r.sh
set -u
Rscript --vanilla tests/testthat.R
```

The suite/report requires `testthat`, `Matrix`, `rmarkdown`, and `knitr`.
Checks that use the `Seurat` reader or verify the `scDblFinder` export are
skipped if those packages are not installed. `tests/testing.Rmd` provides a
rendered report that runs the same suite and documents cluster smoke checks.

## Automated Coverage

| Test file | Coverage |
|-----------|----------|
| `test-adt-qc.R` | ADT MAD-based defaults and explicit thresholds; ADT filter boundaries and disabled no-op behavior; HTO Singlet-only mask; `recoverDoublets` PCA/UMAP name handling; validation of the optional raw ADT path override. |
| `test-hto-adt-split.R` | Splitting combined antibody features using canonical feature-reference names; clear errors for missing or invalid reference/HTO names; per-HTO thresholds derived from the minimum CLR expression among assigned Singlets. |
| `test-sample-sheet.R` | Reconstructing and selecting samples from the transposed sheet; rejecting malformed sheet structure. |
| `test-static-checks.R` | Parsing the runner and sample QC R Markdown; checking HTO demultiplexing and Singlet logic, raw ADT override/fallback routing, and the integration-level default (`Sample`) with `Batch` override; shell syntax; sample-name extraction from the CSV header; `scDblFinder::recoverDoublets` availability when installed. |
| `test-synthetic-10x.R` | Reading generated filtered/raw 10x MEX matrices through `Seurat::Read10X()`; checking combined antibody splitting, separate HTO and ADT matrices, native `Multiplexing Capture`, and raw-only barcodes. |

`helper-fixtures.R` creates the temporary MEX matrices and sample sheets used by
the tests. `helper-load.R` locates the repository and loads the shared QC
helpers; neither helper is a standalone test.

## Not Covered Locally

The suite does not render full reports against real CellRanger data, verify
end-to-end biological results, or submit/execute SLURM jobs. Those require
cluster resources and usable data. See `testing.Rmd` for the manual cluster
smoke checks, including RNA-only rendering, hashtagged+CITE-seq rendering,
merged metadata/report checks, and progress-log validation.
