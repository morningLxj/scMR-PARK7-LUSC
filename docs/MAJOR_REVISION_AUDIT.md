# Major Revision Audit

Revision date: 2026-07-23

## Genetic audit

- 91 submitted leave-one-out rows.
- 90 unique variants because `rs225103` was duplicated.
- 89 variants matched the 1000 Genomes European reference.
- Median pairwise r-squared: 0.645.
- Pairs with r-squared at least 0.1: 3,916 of 3,916.
- Pairs with r-squared at least 0.8: 1,559.
- Strict clumping retained `rs6700772`.
- LUSC-specific outcome (GCST004750): Wald-ratio OR 1.228; 95% CI
  1.002-1.506; nominal P=0.048.
- Overall-lung-cancer sensitivity outcome (GCST004748): Wald-ratio OR 1.124;
  95% CI 0.991-1.274; P=0.069.

The retained per-variant source table had been populated from GCST004748
(overall lung cancer) despite being labelled as the primary LUSC outcome. The
strict-pruned primary estimate was therefore recomputed from the
histology-specific GCST004750 association for the same harmonized variant. The
submitted OR 1.277 was not reproduced. The unpruned diagnostic recomputation
ignores LD and is not inferential.

## Colocalization

At p1=p2=1e-4 and p12=1e-5:

- Intermediate B cells, LUSC: PP1 0.528; PP3 0.105; PP4 0.367.
- Memory B cells, LUSC: PP1 0.514; PP3 0.099; PP4 0.377.

PP1 dominated and PP4 was prior-sensitive, so a shared causal variant was not established.

## Downstream evidence

- Primary Visium: PARK7 absent; observed/expected overlap 0.73; OR 0.62; Spearman rho -0.10.
- Visium sensitivity: observed/expected overlap remained below 1 across six
  thresholds. For the primary cutoff, lower-tail within-block permutation P
  was 0.011 for 8- and 12-unit blocks and below 0.00002 for 16-unit blocks;
  significance was weaker at some alternative thresholds.
- CosMx after CD74 exclusion and myeloid adjustment: mean partial Pearson r 0.146.
- Single-cell: score associations changed by threshold.
- TCGA-LUSC: PARK7 was not associated with overall survival in univariable or adjusted models.
- IHC: tissue-level detectability only; no PARK7/B-cell-marker co-stain.

## Current claim

PARK7 remains a provisional, hypothesis-generating candidate. The current data do not establish an LD-independent causal effect, shared causal variant, B-cell-specific spatial or protein localization, PARK7-specific stress program, prognostic value, or therapeutic vulnerability.
