# Reproducibility Notes

## Included-input analyses

The default analysis command reruns:

- Visium alternative-threshold and spatially restricted within-block
  permutation sensitivity.
- GSE148071 B-lineage PARK7 threshold sensitivity.
- Export of complete PARK7 PP0-PP4 posteriors at the primary prior.

```bash
python scripts/revision/run_major_revision_analyses.py
```

## External-resource analyses

Select additional steps with `--steps`:

```bash
python scripts/revision/run_major_revision_analyses.py \
  --steps ld cosmx visium scrna tcga coloc \
  --project-root /path/to/external/project \
  --cosmx-root /path/to/cosmx/samples
```

Expected external layout under `--project-root`:

```text
plink.exe
Reference/
  g1000_eur.bed
  g1000_eur.bim
  g1000_eur.fam
TCGA/
  TCGA-LUSC.survival.tsv
  TCGA-LUSC.clinical.tsv
```

`--cosmx-root` should contain the five sample directories used in the manuscript:

```text
Lung6
Lung9_Rep2
Lung12
Lung13
Lung5_Rep3
```

Environment variables may be used instead:

- `PARK7_SOURCE_ROOT`
- `PARK7_PROJECT_ROOT`
- `PARK7_COSMX_ROOT`
- `PARK7_OUT`
- `PARK7_RESULTS_ROOT`
- `PARK7_FIGURE_OUT`

## Figure regeneration

```bash
python scripts/revision/build_revision_figures.py
```

The figure script reads `data/inputs/` and `data/revision_results/` by default.

## Interpretation

The scripts reproduce reported computations from retained processed inputs. They do not convert the correlated regional PARK7 variant set into independent instruments and must not be used to claim causality, high-confidence colocalization, B-cell-specific localization, prognosis, or therapy.

The Visium sensitivity uses strict greater-than cutoffs at the global 50th,
75th, 80th, 85th, and 90th percentiles plus a z-score-above-1 definition.
Because the B-cell score contains many ties at zero, achieved high-spot
fractions are reported rather than assumed to equal the nominal quantile.
Unrestricted overlap probabilities use the corresponding hypergeometric
distribution. Spatially restricted sensitivity conditions on the antioxidant
high-spot count within 8-, 12-, and 16-array-unit blocks and simulates 100,000
within-block permutations for each threshold. These block tests preserve broad
spatial prevalence but do not fully model spot-level spatial autocorrelation.

The retained regional variant-level colocalization source covers the
intermediate-B-cell LUSC analysis only. Its eQTL and GWAS coordinate columns
are exported separately in the provenance summary because genome-build labels
were not retained. The other three analyses retain complete posteriors and
shared-SNP counts but not variant-level regional inputs; a new cross-build
regional association plot therefore cannot be reconstructed reliably.

The LD step clumps on the exposure associations and then joins the retained
variant to accession-labelled outcome associations in
`data/inputs/PARK7_StrictPruned_Outcome_Associations.csv`. GCST004750 is the
primary LUSC outcome; GCST004748 is retained only as an
overall-lung-cancer sensitivity outcome.
