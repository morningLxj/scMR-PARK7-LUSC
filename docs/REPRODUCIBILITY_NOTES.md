# Reproducibility Notes

## Included-input analyses

The default analysis command reruns:

- GSE148071 B-lineage PARK7 threshold sensitivity.
- Export of complete PARK7 PP0-PP4 posteriors at the primary prior.

```bash
python scripts/revision/run_major_revision_analyses.py
```

## External-resource analyses

Select additional steps with `--steps`:

```bash
python scripts/revision/run_major_revision_analyses.py \
  --steps ld cosmx scrna tcga coloc \
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

The LD step clumps on the exposure associations and then joins the retained
variant to accession-labelled outcome associations in
`data/inputs/PARK7_StrictPruned_Outcome_Associations.csv`. GCST004750 is the
primary LUSC outcome; GCST004748 is retained only as an
overall-lung-cancer sensitivity outcome.
