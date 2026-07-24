# Data Availability

## Public resources

- OneK1K immune-cell eQTL resource: Yazar et al., Science 2022.
- 1000 Genomes Phase 3 European reference panel.
- 10x Genomics Visium CytAssist FFPE LUSC demonstration dataset, Space Ranger 2.0.0.
- CosMx NSCLC spatial molecular imaging dataset described by He et al., Nature Biotechnology 2022.
- NCBI GEO series GSE148071 and GSE200972.
- TCGA-LUSC expression, reverse-phase protein array, clinical, survival, and
  GDC masked-somatic-mutation resources.
- McKay et al. lung-cancer summary association resources: GCST004750
  (LUSC-specific primary outcome) and GCST004748 (overall-lung-cancer
  sensitivity outcome).

## Included here

- Processed inputs needed to inspect or rerun the included-input analyses.
- Major-revision result tables.
- Revised figure files and figure-building code.
- De-identified IHC H-score summaries.
- Harmonized rs6700772 outcome associations and accession-level provenance for
  the strict-pruned primary and sensitivity estimates.
- Visium threshold/permutation sensitivity results and a spatial
  dataset/sample specification table.
- The retained intermediate-B-cell LUSC regional colocalization source and a
  four-analysis coordinate/build provenance summary.
- The complete retained 38-pair MR source table, with its restricted
  multiplicity boundary.
- A spatial feature-coverage audit derived from the original Visium HDF5
  feature list and five CosMx matrix headers.
- GDC-derived KEAP1/NFE2L2/CUL3 mutation status, qualifying event details, and
  query provenance for the mutation-adjusted TCGA analyses.
- TCGA-LUSC PARK7 RNA/DJ-1 RPPA correlation and survival-model summaries.

## Not redistributed

- Raw GWAS or eQTL summary-statistic archives.
- The 1000 Genomes binary reference files.
- Full CosMx cell-by-gene matrices.
- Full TCGA-LUSC expression, RPPA, clinical, or survival downloads.
- Full raw single-cell matrices.
- Raw identifiable pathology records or patient images.

The repository is therefore a bounded reproducibility package. Included-input steps can be rerun directly; external-resource steps require users to obtain the source data from the original providers.
