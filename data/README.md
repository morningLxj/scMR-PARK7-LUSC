# Data

## `inputs/`

Processed inputs retained from the submission workspace:

- PARK7 per-variant eQTL and leave-one-out source tables used for the LD audit.
- GSE148071 selected cell-level score table used for threshold sensitivity.
- TCGA-LUSC selected expression modules.
- PARK7 colocalization prior-sensitivity and locus-comparison tables.
- Primary Visium spot-level score and coordinate table used for revised Fig. 4.

These files are processed inputs, not raw source datasets.

## `revision_results/`

Reported major-revision outputs:

- LD audit and strict-pruned estimates.
- Complete PP0-PP4 colocalization posteriors.
- CD74-excluded and myeloid-adjusted CosMx sensitivities.
- Single-cell threshold sensitivities.
- TCGA-LUSC Cox-model summaries.
- Exploratory network audits.
- De-identified IHC H-score summaries.

The diagnostic unpruned regional estimate is retained only to document the source discrepancy. It is not an inferential result.

No direct personal identifiers or raw pathology records are included.
