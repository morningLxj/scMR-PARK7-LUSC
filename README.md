# scMR-PARK7-LUSC

Code and processed source data for the manuscript:

**Cell-type-resolved genetic prioritization of PARK7 in lung squamous cell carcinoma with indirect spatial and transcriptomic context**

This repository accompanies a hypothesis-generating candidate-prioritization study of PARK7 in lung squamous cell carcinoma (LUSC). The analysis integrates cell-type-resolved Mendelian randomization, SMR/colocalization context, spatial transcriptomic redox scoring, public single-cell RNA-seq validation, TCGA-LUSC bulk-expression context, PARK7-centered network inference, exploratory CMap annotation, and supportive tissue-level PARK7 IHC source data.

## Repository Structure

```text
data/
  source_data_csv/        Processed source-data tables used for manuscript figures and supplementary analyses
  source_workbooks/       Supplementary source-data workbooks
figures/
  main/                   Submitted main figures as TIFF files
  supplementary/          Submitted supplementary figures as TIFF files
scripts/
  analysis/               Analysis/provenance scripts from the final submission workflow
  check_repository_integrity.py
docs/
  DATA_AVAILABILITY.md
  REPRODUCIBILITY_NOTES.md
requirements.txt
environment.yml
```

## Quick Check

After cloning, run:

```bash
python scripts/check_repository_integrity.py
```

The check verifies the expected source-data files, main figures, supplementary figures, and workflow scripts.

## Reproducibility Scope

This repository includes processed source data sufficient to inspect the numerical inputs behind the submitted figures and tables. The full end-to-end regeneration of every upstream result requires large public raw resources that are not redistributed here, including GWAS/eQTL summary data, public scRNA-seq matrices, TCGA-LUSC expression data, spatial transcriptomic inputs, and CMap/L1000 outputs. See `docs/DATA_AVAILABILITY.md` for data-source notes.

The scripts in `scripts/analysis/` are retained as provenance workflows from the final local submission assembly. They may require path edits and downloaded raw data before rerunning on another machine.

## Evidence Boundary

The manuscript and this repository should be interpreted as supporting a candidate-prioritization framework. They do not establish definitive causality, high-confidence colocalization, direct spatial PARK7 localization, B-cell-specific PARK7 protein localization, validated biomarker status, therapeutic vulnerability, or drug efficacy.

## Citation

If you use this repository, cite the associated manuscript when available and reference this GitHub repository:

`morningLxj/scMR-PARK7-LUSC`

