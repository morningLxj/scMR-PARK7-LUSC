# Analysis Scripts

These scripts are retained as provenance workflows from the final submission package.

- `01_public_scrna_tcga_ihc_cmap_workflow.py`: public GSE200972 scRNA context, filtered CMap signature construction, TCGA-LUSC bulk-expression context, and IHC ROI/tile source-data generation.
- `02_extended_public_validation_and_mr_sensitivity.py`: GSE148071 public scRNA validation, extended MR robustness summaries, colocalization prior sensitivity, spatial neighborhood analysis, and related supplementary source data.
- `03_evidence_boundary_audit_figures.py`: evidence-boundary and reviewer-facing audit figures/tables that define what each layer supports and does not prove.

The scripts may require editing local path constants and downloading large public raw inputs before rerun. The processed outputs used in the manuscript are included in `data/source_data_csv/`.

