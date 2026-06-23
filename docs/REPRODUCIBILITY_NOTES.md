# Reproducibility Notes

## Recommended Use

1. Inspect `data/source_data_csv/Source_Data_CSV_Manifest.csv` to map each source-data file to its manuscript item.
2. Use `figures/main/` and `figures/supplementary/` as the submitted visual outputs.
3. Use `scripts/analysis/` as provenance code for the local final submission workflow.
4. Run `python scripts/check_repository_integrity.py` to verify that the repository contains the expected files.

## Script Caveats

The three scripts in `scripts/analysis/` were copied from the final local workflow. They document how public scRNA-seq, TCGA, IHC ROI, CMap, MR robustness, spatial boundary, and evidence-boundary artifacts were produced and synchronized during submission assembly.

Because they were originally run in a local project workspace, they may contain local path constants and package-synchronization steps. To rerun them elsewhere, edit the path constants near the top of each script and provide the required raw public inputs.

## Claim-Boundary Caveat

The source data support candidate prioritization and evidence-boundary assessment. They should not be used to claim direct functional causality, direct PARK7-positive B-cell spatial localization, B-cell-specific PARK7 protein localization, therapeutic efficacy, or clinical biomarker validation.

