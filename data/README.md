# Data

## `inputs/`

Processed inputs retained from the submission workspace:

- PARK7 per-variant eQTL and leave-one-out source tables used for the LD audit.
- Accession-labelled rs6700772 outcome associations for GCST004750 (primary
  LUSC) and GCST004748 (overall-lung-cancer sensitivity).
- GSE148071 selected cell-level score table used for threshold sensitivity.
- TCGA-LUSC selected expression modules.
- GDC masked-somatic-mutation coverage and canonical HIGH/MODERATE
  KEAP1/NFE2L2/CUL3 event status for 490 TCGA-LUSC cases, with event-level
  and query-provenance tables.
- PARK7 colocalization prior-sensitivity and locus-comparison tables.
- Primary Visium spot-level score and coordinate table used for revised Fig. 4.

These files are processed inputs, not raw source datasets.

## `revision_results/`

Reported major-revision outputs:

- LD audit and outcome-specific strict-pruned estimates.
- Complete PP0-PP4 colocalization posteriors.
- Colocalization regional-source coordinate and genome-build provenance
  boundaries.
- CD74-excluded and myeloid-adjusted CosMx sensitivities.
- Visium alternative-threshold and spatially restricted within-block
  permutation sensitivities.
- Spatial platform, sample, unit, preprocessing, and coordinate specifications.
- Single-cell threshold sensitivities.
- The complete retained 38-pair MR source table. This table is the available
  subset, not a reconstructable complete discovery universe.
- Spatial platform feature-coverage audit confirming that PARK7 was absent
  from both the Visium and CosMx panels used in the manuscript.
- TCGA-LUSC PARK7 RNA and DJ-1 RPPA Cox-model summaries, RNA-protein
  correlation, and pathway-mutation coverage summaries.
- De-identified IHC H-score summaries.

Exploratory network files retained in repository history are provenance
artifacts and are not part of the revised evidence set or submission archive.

The diagnostic unpruned regional estimate is retained only to document the source discrepancy. It is not an inferential result.

No direct personal identifiers or raw pathology records are included.
