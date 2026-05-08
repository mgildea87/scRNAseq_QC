# Next Steps

## Immediate
- [x] Re-run integration-only in SLURM mode and confirm the job completes.
- [x] Verify chunk progress log reaches final integration chunks.
- [x] Confirm expected outputs exist in output directory:
  - `integrated.rds`
  - `lisi_integrated_RNA.rds`
  - `FindMarkersBulk_RNA_integrated/Top_markers.csv`
- [x] Run integration-only at both levels and compare behavior:
  - `--integration_level Batch`
  - `--integration_level Sample`

## Suggested Commands
```bash
# Integration-only SLURM run (Batch)
./submit_QC_batch.sh \
  --sample_sheet samples.csv \
  --integration_only TRUE \
  --integration_level Batch \
  --run_integration TRUE

# Integration-only SLURM run (Sample)
./submit_QC_batch.sh \
  --sample_sheet samples.csv \
  --integration_only TRUE \
  --integration_level Sample \
  --run_integration TRUE
```

## After Validation
- [ ] Add SCTransform integration to the integration workflow as an addition
- [ ] Add multiple cluster resolutions?
- [ ] Update README examples with the exact final tested integration-only command(s).
- [ ] Add a troubleshooting note in README for integration chunk failures:
  - where to check SLURM logs
  - where to check `chunk_logs/`
- [ ] Decide whether to reduce `conda-package-list` verbosity in the integration report.
- [ ] Add default reference annotation? this would be at the individual sample level likely
- [ ] Containerize?

## Open Questions
- [ ] Should integration-only require merged input precheck with a clearer error message before render?
- [ ] Should `integration_mem` default be increased for large datasets?

## Change Log
- 2026-05-08: Removed unused `tf_list_path` parameter and template path resolver from `integrate_RNA.Rmd`.
- 2026-05-08: Added optional integration flow, integration-only mode, integration level requirement, and SLURM `integration_mem` controls.
