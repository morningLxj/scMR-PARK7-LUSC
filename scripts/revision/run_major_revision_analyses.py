from __future__ import annotations

import argparse
import math
import os
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


REPO_ROOT = Path(__file__).resolve().parents[2]
SOURCE_ROOT = Path(
    os.environ.get("PARK7_SOURCE_ROOT", REPO_ROOT / "data" / "inputs")
).resolve()
PROJECT_ROOT = Path(
    os.environ.get("PARK7_PROJECT_ROOT", REPO_ROOT / "external")
).resolve()
COSMX_ROOT = Path(
    os.environ.get("PARK7_COSMX_ROOT", PROJECT_ROOT / "cosmx")
).resolve()
OUT = Path(
    os.environ.get("PARK7_OUT", REPO_ROOT / "data" / "revision_results")
).resolve()
OUT.mkdir(parents=True, exist_ok=True)


def bh_fdr(values: pd.Series) -> pd.Series:
    result = pd.Series(np.nan, index=values.index, dtype=float)
    valid = values.notna()
    p = values.loc[valid].astype(float).to_numpy()
    if not len(p):
        return result
    order = np.argsort(p)
    ranked = p[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    restored = np.empty_like(adjusted)
    restored[order] = np.minimum(adjusted, 1.0)
    result.loc[valid] = restored
    return result


def ivw_ratio_estimate(frame: pd.DataFrame) -> dict[str, float]:
    d = frame.copy()
    x = d["exposure_beta"].astype(float).to_numpy()
    sx = d["exposure_se"].astype(float).to_numpy()
    y = d["outcome_beta"].astype(float).to_numpy()
    sy = d["outcome_se"].astype(float).to_numpy()
    ratio = y / x
    ratio_var = sy**2 / x**2 + (y**2 * sx**2) / x**4
    weights = 1.0 / ratio_var
    beta = float(np.sum(weights * ratio) / np.sum(weights))
    se = float(math.sqrt(1.0 / np.sum(weights)))
    p = float(2 * stats.norm.sf(abs(beta / se)))
    return {
        "n_instruments": int(len(d)),
        "beta": beta,
        "se": se,
        "p": p,
        "OR": math.exp(beta),
        "CI_low": math.exp(beta - 1.96 * se),
        "CI_high": math.exp(beta + 1.96 * se),
    }


def run_ld_audit() -> None:
    plink = PROJECT_ROOT / "plink.exe"
    bfile = PROJECT_ROOT / "Reference" / "g1000_eur"
    instruments = pd.read_csv(
        SOURCE_ROOT / "PARK7_PerInstrument_F_Statistics.csv"
    )
    outcomes = pd.read_csv(
        SOURCE_ROOT / "PARK7_StrictPruned_Outcome_Associations.csv"
    )
    leave_one_out = pd.read_csv(
        SOURCE_ROOT / "PARK7_LeaveOneOut_Source.csv"
    )
    primary_snps = set(leave_one_out["SNP_removed"].astype(str))
    instruments = instruments[instruments["SNP"].astype(str).isin(primary_snps)].copy()
    instruments = instruments.drop_duplicates(subset="SNP", keep="first")

    snplist = OUT / "PARK7_primary_unique_snps.txt"
    snplist.write_text("\n".join(instruments["SNP"].astype(str)) + "\n", encoding="ascii")
    ld_prefix = OUT / "PARK7_primary_unique_EUR_LD"
    subprocess.run(
        [
            str(plink),
            "--bfile",
            str(bfile),
            "--extract",
            str(snplist),
            "--r2",
            "square",
            "--out",
            str(ld_prefix),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    ld = np.loadtxt(ld_prefix.with_suffix(".ld"))
    upper = ld[np.triu_indices_from(ld, k=1)]
    ld_summary = {
        "requested_instruments": len(instruments),
        "instruments_in_1000G_EUR": int(ld.shape[0]),
        "max_pairwise_r2": float(np.nanmax(upper)),
        "median_pairwise_r2": float(np.nanmedian(upper)),
        "pairs_r2_ge_0_001": int(np.sum(upper >= 0.001)),
        "pairs_r2_ge_0_01": int(np.sum(upper >= 0.01)),
        "pairs_r2_ge_0_1": int(np.sum(upper >= 0.1)),
        "pairs_r2_ge_0_8": int(np.sum(upper >= 0.8)),
        "total_pairs": int(len(upper)),
    }

    estimate_rows = []
    retained_rows = []
    for r2 in (0.001, 0.01, 0.1):
        assoc = OUT / f"PARK7_clump_input_r2_{r2}.tsv"
        instruments[["SNP", "exposure_p"]].rename(
            columns={"exposure_p": "P"}
        ).to_csv(assoc, sep="\t", index=False)
        prefix = OUT / f"PARK7_EUR_clump_r2_{r2}"
        subprocess.run(
            [
                str(plink),
                "--bfile",
                str(bfile),
                "--clump",
                str(assoc),
                "--clump-p1",
                "1",
                "--clump-p2",
                "1",
                "--clump-kb",
                "10000",
                "--clump-r2",
                str(r2),
                "--out",
                str(prefix),
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        clumped_path = Path(f"{prefix}.clumped")
        clumped = pd.read_csv(clumped_path, sep=r"\s+")
        retained_exposure = instruments[
            instruments["SNP"].isin(clumped["SNP"])
        ].copy()
        retained_exposure = retained_exposure.drop(
            columns=[
                "outcome_effect_allele",
                "outcome_other_allele",
                "outcome_beta",
                "outcome_se",
                "outcome_p",
                "included_in_primary_lusc_mr",
            ],
            errors="ignore",
        )
        for outcome_id, outcome_rows in outcomes.groupby("outcome_id", sort=False):
            retained = retained_exposure.merge(
                outcome_rows,
                on="SNP",
                how="inner",
                validate="one_to_one",
            )
            if len(retained) != len(retained_exposure):
                raise ValueError(
                    f"{outcome_id} is missing one or more strict-pruned variants"
                )
            allele_match = (
                retained["effect_allele"].astype(str)
                == retained["outcome_effect_allele"].astype(str)
            ) & (
                retained["other_allele"].astype(str)
                == retained["outcome_other_allele"].astype(str)
            )
            if not allele_match.all():
                raise ValueError(f"Alleles are not harmonized for {outcome_id}")

            estimate = ivw_ratio_estimate(retained)
            estimate.update(
                {
                    "clump_r2": r2,
                    "outcome_id": outcome_id,
                    "outcome_label": retained["outcome_label"].iloc[0],
                    "analysis_role": retained["analysis_role"].iloc[0],
                    "cases": int(retained["cases"].iloc[0]),
                    "controls": int(retained["controls"].iloc[0]),
                    "ancestry": retained["ancestry"].iloc[0],
                }
            )
            estimate_rows.append(estimate)
            retained.insert(0, "clump_r2", r2)
            retained_rows.append(retained)

    pd.DataFrame([ld_summary]).to_csv(OUT / "MR_PARK7_LD_Audit_Summary.csv", index=False)
    pd.DataFrame(estimate_rows).to_csv(
        OUT / "MR_PARK7_LD_Pruned_Estimates.csv", index=False
    )
    pd.concat(retained_rows, ignore_index=True).to_csv(
        OUT / "MR_PARK7_LD_Pruned_Instruments.csv", index=False
    )
    outcomes.to_csv(OUT / "MR_PARK7_Outcome_Provenance.csv", index=False)


def residualize(values: pd.Series, covariate: pd.Series) -> np.ndarray:
    x = np.column_stack([np.ones(len(covariate)), covariate.to_numpy(float)])
    beta, *_ = np.linalg.lstsq(x, values.to_numpy(float), rcond=None)
    return values.to_numpy(float) - x @ beta


def partial_corr(x: pd.Series, y: pd.Series, z: pd.Series, rank: bool = False) -> float:
    frame = pd.concat([x, y, z], axis=1).dropna()
    if rank:
        frame = frame.rank(method="average")
    rx = residualize(frame.iloc[:, 0], frame.iloc[:, 2])
    ry = residualize(frame.iloc[:, 1], frame.iloc[:, 2])
    return float(np.corrcoef(rx, ry)[0, 1])


def run_cosmx_sensitivity() -> None:
    data_root = COSMX_ROOT
    samples = ["Lung6", "Lung9_Rep2", "Lung12", "Lung13", "Lung5_Rep3"]
    b_original = ["CD79A", "MS4A1", "CD74", "CD37"]
    b_no_cd74 = ["CD79A", "MS4A1", "CD37"]
    oxidative = ["SOD1", "GPX1"]
    myeloid = ["LYZ", "TYROBP", "CD68", "FCGR3A", "S100A8", "S100A9"]
    rows = []

    for sample in samples:
        sample_dir = data_root / sample / f"{sample}-Flat_files_and_images"
        expr_path = sample_dir / f"{sample}_exprMat_file.csv"
        meta_path = sample_dir / f"{sample}_metadata_file.csv"
        header = pd.read_csv(expr_path, nrows=0).columns
        marker_cols = [
            g for g in sorted(set(b_original + oxidative + myeloid)) if g in header
        ]
        expr = pd.read_csv(expr_path, usecols=["fov", "cell_ID", *marker_cols])
        meta = pd.read_csv(
            meta_path,
            usecols=["fov", "cell_ID", "CenterX_global_px", "CenterY_global_px"],
        )
        cells = pd.merge(meta, expr, on=["fov", "cell_ID"], how="inner")
        for gene in marker_cols:
            cells[gene] = pd.to_numeric(cells[gene], errors="coerce").fillna(0.0)
            cells[f"log_{gene}"] = np.log1p(cells[gene])
        cells["b_original"] = cells[[f"log_{g}" for g in b_original]].mean(axis=1)
        cells["b_no_cd74"] = cells[[f"log_{g}" for g in b_no_cd74]].mean(axis=1)
        cells["oxidative"] = cells[[f"log_{g}" for g in oxidative]].mean(axis=1)
        myeloid_used = [g for g in myeloid if g in marker_cols]
        cells["myeloid"] = cells[[f"log_{g}" for g in myeloid_used]].mean(axis=1)
        cells["grid_x"] = (cells["CenterX_global_px"] // 320).astype(int)
        cells["grid_y"] = (cells["CenterY_global_px"] // 320).astype(int)
        bins = (
            cells.groupby(["grid_x", "grid_y"], as_index=False)
            .agg(
                b_original=("b_original", "mean"),
                b_no_cd74=("b_no_cd74", "mean"),
                oxidative=("oxidative", "mean"),
                myeloid=("myeloid", "mean"),
                n_cells=("cell_ID", "size"),
            )
        )
        bins = bins[bins["n_cells"] >= 16].copy()
        rows.append(
            {
                "sample": sample,
                "n_cells": len(cells),
                "n_spatial_bins": len(bins),
                "original_B_markers": ",".join(b_original),
                "CD74_excluded_B_markers": ",".join(b_no_cd74),
                "oxidative_markers": ",".join(oxidative),
                "myeloid_markers": ",".join(myeloid_used),
                "pearson_original_B_vs_oxidative": bins["b_original"].corr(
                    bins["oxidative"], method="pearson"
                ),
                "pearson_CD74_excluded_B_vs_oxidative": bins["b_no_cd74"].corr(
                    bins["oxidative"], method="pearson"
                ),
                "partial_pearson_CD74_excluded_adjusted_myeloid": partial_corr(
                    bins["b_no_cd74"], bins["oxidative"], bins["myeloid"]
                ),
                "spearman_original_B_vs_oxidative": bins["b_original"].corr(
                    bins["oxidative"], method="spearman"
                ),
                "spearman_CD74_excluded_B_vs_oxidative": bins["b_no_cd74"].corr(
                    bins["oxidative"], method="spearman"
                ),
                "partial_spearman_CD74_excluded_adjusted_myeloid": partial_corr(
                    bins["b_no_cd74"],
                    bins["oxidative"],
                    bins["myeloid"],
                    rank=True,
                ),
            }
        )

    summary = pd.DataFrame(rows)
    summary.to_csv(OUT / "Spatial_CosMx_CD74_Myeloid_Sensitivity.csv", index=False)
    metric_cols = [c for c in summary.columns if "pearson" in c or "spearman" in c]
    aggregate = []
    for metric in metric_cols:
        vals = summary[metric].astype(float)
        aggregate.append(
            {
                "metric": metric,
                "n_samples": vals.notna().sum(),
                "mean_correlation": vals.mean(),
                "median_correlation": vals.median(),
                "positive_samples": int((vals > 0).sum()),
                "negative_samples": int((vals < 0).sum()),
            }
        )
    pd.DataFrame(aggregate).to_csv(
        OUT / "Spatial_CosMx_CD74_Myeloid_Aggregate.csv", index=False
    )


def block_permutation_summary(
    frame: pd.DataFrame,
    b_high: np.ndarray,
    a_high: np.ndarray,
    block_width: int,
    n_permutations: int,
    seed: int,
) -> dict[str, float]:
    block_row = (frame["array_row"].to_numpy(int) // block_width).astype(int)
    block_col = (frame["array_col"].to_numpy(int) // block_width).astype(int)
    block_ids = np.char.add(
        np.char.add(block_row.astype(str), "_"), block_col.astype(str)
    )
    rng = np.random.default_rng(seed)
    null_overlap = np.zeros(n_permutations, dtype=np.int32)
    expected = 0.0
    for block_id in np.unique(block_ids):
        mask = block_ids == block_id
        n_block = int(mask.sum())
        b_count = int(b_high[mask].sum())
        a_count = int(a_high[mask].sum())
        if not n_block or not b_count or not a_count:
            continue
        expected += b_count * a_count / n_block
        null_overlap += rng.hypergeometric(
            ngood=a_count,
            nbad=n_block - a_count,
            nsample=b_count,
            size=n_permutations,
        )
    observed = int(np.logical_and(b_high, a_high).sum())
    return {
        f"block_{block_width}_expected_overlap": expected,
        f"block_{block_width}_p_lower_or_equal": (
            1 + int((null_overlap <= observed).sum())
        )
        / (n_permutations + 1),
        f"block_{block_width}_p_upper_or_equal": (
            1 + int((null_overlap >= observed).sum())
        )
        / (n_permutations + 1),
    }


def run_visium_threshold_permutation() -> None:
    frame = pd.read_csv(
        SOURCE_ROOT / "Visium_Spatial_CoLocalization_Scores.csv"
    )
    required = {
        "array_row",
        "array_col",
        "B_cell_density_raw",
        "Antioxidant_raw",
        "B_cell_density_z",
        "Antioxidant_z",
    }
    missing = sorted(required.difference(frame.columns))
    if missing:
        raise ValueError(f"Visium input is missing columns: {missing}")

    threshold_specs = [
        ("global_q50", 0.50),
        ("global_q75_primary", 0.75),
        ("global_q80", 0.80),
        ("global_q85", 0.85),
        ("global_q90", 0.90),
        ("z_gt_1", None),
    ]
    n_permutations = 100_000
    rows = []
    for index, (label, quantile) in enumerate(threshold_specs):
        if quantile is None:
            b_cutoff = 1.0
            a_cutoff = 1.0
            b_high = frame["B_cell_density_z"].to_numpy(float) > b_cutoff
            a_high = frame["Antioxidant_z"].to_numpy(float) > a_cutoff
            scale = "z score"
        else:
            b_cutoff = float(frame["B_cell_density_raw"].quantile(quantile))
            a_cutoff = float(frame["Antioxidant_raw"].quantile(quantile))
            b_high = frame["B_cell_density_raw"].to_numpy(float) > b_cutoff
            a_high = frame["Antioxidant_raw"].to_numpy(float) > a_cutoff
            scale = "mean log1p score"

        n_spots = len(frame)
        b_count = int(b_high.sum())
        a_count = int(a_high.sum())
        observed = int(np.logical_and(b_high, a_high).sum())
        expected = b_count * a_count / n_spots
        b_only = b_count - observed
        a_only = a_count - observed
        neither = n_spots - observed - b_only - a_only
        odds_ratio, fisher_p = stats.fisher_exact(
            [[observed, b_only], [a_only, neither]],
            alternative="two-sided",
        )
        random_lower = float(
            stats.hypergeom.cdf(observed, n_spots, a_count, b_count)
        )
        random_upper = float(
            stats.hypergeom.sf(observed - 1, n_spots, a_count, b_count)
        )
        row = {
            "threshold": label,
            "score_scale": scale,
            "quantile": quantile,
            "strict_greater_than": True,
            "B_cutoff": b_cutoff,
            "antioxidant_cutoff": a_cutoff,
            "n_spots": n_spots,
            "B_high": b_count,
            "antioxidant_high": a_count,
            "both_high_observed": observed,
            "both_high_expected_independence": expected,
            "observed_over_expected": observed / expected if expected else np.nan,
            "fisher_odds_ratio": odds_ratio,
            "fisher_two_sided_p": fisher_p,
            "unrestricted_permutation_p_lower_or_equal": random_lower,
            "unrestricted_permutation_p_upper_or_equal": random_upper,
            "permutations_per_block_width": n_permutations,
        }
        for block_width in (8, 12, 16):
            row.update(
                block_permutation_summary(
                    frame,
                    b_high,
                    a_high,
                    block_width=block_width,
                    n_permutations=n_permutations,
                    seed=20260724 + 100 * index + block_width,
                )
            )
        rows.append(row)

    pd.DataFrame(rows).to_csv(
        OUT / "Visium_Threshold_BlockPermutation_Sensitivity.csv",
        index=False,
    )

    cosmx_path = OUT / "Spatial_CosMx_CD74_Myeloid_Sensitivity.csv"
    if cosmx_path.exists():
        cosmx = pd.read_csv(cosmx_path)
        specification_rows = [
            {
                "resource": "10x Genomics public LUSC demonstration",
                "platform_version": "Visium CytAssist; Space Ranger 2.0.0",
                "public_identifier": (
                    "CytAssist_FFPE_Human_Lung_Squamous_Cell_Carcinoma"
                ),
                "sample_id": "One FFPE LUSC section",
                "raw_public_units": "Not separately retained",
                "matched_analysis_units": 3858,
                "matched_unit_type": "under-tissue spots",
                "retained_aggregation_units": 3858,
                "retained_unit_type": "spots",
                "quality_control": (
                    "Under-tissue spots from the public Space Ranger output; "
                    "no additional spot exclusion"
                ),
                "normalization": "log1p counts; score-wise z standardization",
                "coordinate_handling": (
                    "Space Ranger array_row/array_col; Euclidean distance in "
                    "array-coordinate units"
                ),
                "analysis_level": "spot-level; no deconvolution",
            }
        ]
        for record in cosmx.to_dict(orient="records"):
            specification_rows.append(
                {
                    "resource": "Bruker CosMx NSCLC FFPE public dataset",
                    "platform_version": (
                        "CosMx SMI; release version not stated in retained files"
                    ),
                    "public_identifier": "Public NSCLC FFPE dataset",
                    "sample_id": record["sample"],
                    "raw_public_units": "Not separately retained",
                    "matched_analysis_units": int(record["n_cells"]),
                    "matched_unit_type": "segmented cells after metadata-expression merge",
                    "retained_aggregation_units": int(record["n_spatial_bins"]),
                    "retained_unit_type": "320-pixel bins with at least 16 cells",
                    "quality_control": (
                        "Metadata-expression inner merge; bins with fewer "
                        "than 16 cells excluded"
                    ),
                    "normalization": "gene-wise log1p counts before score averaging",
                    "coordinate_handling": (
                        "CenterX_global_px and CenterY_global_px divided into "
                        "320-pixel grid bins"
                    ),
                    "analysis_level": (
                        "cell-level scoring followed by bin-level correlations"
                    ),
                }
            )
        pd.DataFrame(specification_rows).to_csv(
            OUT / "Spatial_Dataset_Sample_Specification.csv",
            index=False,
        )


def run_scrna_threshold_sensitivity() -> None:
    cells = pd.read_csv(
        SOURCE_ROOT / "GSE148071_LUSC_scRNA_cell_level_selected_scores.csv"
    )
    cells = cells[cells["assigned_cell_type"].isin(["B cells", "Plasma cells"])].copy()
    score_cols = [
        "UPR_ER_stress_score",
        "NRF2_redox_score",
        "Ribosome_translation_score",
        "Apoptosis_DDR_score",
    ]
    sample_rows = []
    groups = ["B cells", "Plasma cells", "B-lineage combined"]
    methods = [
        "within_sample_positive_median",
        "within_sample_positive_quartiles",
        "detected_vs_undetected",
    ]

    for sample, sample_df in cells.groupby("sample"):
        for cell_group in groups:
            sub = (
                sample_df
                if cell_group == "B-lineage combined"
                else sample_df[sample_df["assigned_cell_type"] == cell_group]
            )
            if len(sub) < 10:
                continue
            park7 = sub["PARK7_log2CPM"].astype(float)
            positive = park7[park7 > 0]
            for method in methods:
                if method == "within_sample_positive_median":
                    if len(positive) < 6:
                        continue
                    cutoff = float(positive.median())
                    high = sub[park7 >= cutoff]
                    low = sub[park7 < cutoff]
                    lower_cutoff = np.nan
                    upper_cutoff = cutoff
                elif method == "within_sample_positive_quartiles":
                    if len(positive) < 12:
                        continue
                    q1, q3 = positive.quantile([0.25, 0.75]).astype(float)
                    high = sub[park7 >= q3]
                    low = sub[(park7 > 0) & (park7 <= q1)]
                    lower_cutoff = float(q1)
                    upper_cutoff = float(q3)
                else:
                    high = sub[park7 > 0]
                    low = sub[park7 == 0]
                    lower_cutoff = 0.0
                    upper_cutoff = 0.0
                if len(high) < 3 or len(low) < 3:
                    continue
                for score in score_cols:
                    sample_rows.append(
                        {
                            "sample": sample,
                            "cell_group": cell_group,
                            "threshold_method": method,
                            "score": score,
                            "n_high": len(high),
                            "n_low": len(low),
                            "lower_cutoff": lower_cutoff,
                            "upper_cutoff": upper_cutoff,
                            "mean_high": high[score].mean(),
                            "mean_low": low[score].mean(),
                            "delta_high_minus_low": high[score].mean()
                            - low[score].mean(),
                            "mannwhitney_p": stats.mannwhitneyu(
                                high[score], low[score], alternative="two-sided"
                            ).pvalue,
                        }
                    )

    per_sample = pd.DataFrame(sample_rows)
    per_sample["mannwhitney_fdr"] = per_sample.groupby(
        ["sample", "cell_group", "threshold_method"]
    )["mannwhitney_p"].transform(bh_fdr)
    per_sample.to_csv(
        OUT / "scRNA_PARK7_Threshold_Sensitivity_PerSample.csv", index=False
    )

    meta_rows = []
    for keys, sub in per_sample.groupby(["cell_group", "threshold_method", "score"]):
        deltas = sub["delta_high_minus_low"].dropna().astype(float)
        p = (
            stats.wilcoxon(deltas).pvalue
            if len(deltas) >= 2 and not np.allclose(deltas, 0)
            else np.nan
        )
        meta_rows.append(
            {
                "cell_group": keys[0],
                "threshold_method": keys[1],
                "score": keys[2],
                "n_samples": len(deltas),
                "median_delta_high_minus_low": deltas.median(),
                "mean_delta_high_minus_low": deltas.mean(),
                "positive_samples": int((deltas > 0).sum()),
                "negative_samples": int((deltas < 0).sum()),
                "wilcoxon_signed_rank_p": p,
            }
        )
    meta = pd.DataFrame(meta_rows)
    meta["wilcoxon_fdr"] = meta.groupby("threshold_method")[
        "wilcoxon_signed_rank_p"
    ].transform(bh_fdr)
    meta.to_csv(OUT / "scRNA_PARK7_Threshold_Sensitivity_Meta.csv", index=False)


def parse_stage(value: object) -> str | float:
    text = str(value).strip().upper()
    if not text or text in {"NAN", "NOT REPORTED", "STAGE X", "UNKNOWN"}:
        return np.nan
    if text.startswith("STAGE IV"):
        return "IV"
    if text.startswith("STAGE III"):
        return "III"
    if text.startswith("STAGE II"):
        return "II"
    if text.startswith("STAGE I"):
        return "I"
    return np.nan


def fit_cox(frame: pd.DataFrame, covariates: list[str], model: str) -> pd.DataFrame:
    from lifelines import CoxPHFitter
    from lifelines.statistics import proportional_hazard_test

    d = frame[["OS_time", "OS_event", *covariates]].dropna().copy()
    cph = CoxPHFitter()
    cph.fit(d, duration_col="OS_time", event_col="OS_event")
    out = cph.summary.reset_index().rename(columns={"covariate": "term"})
    ph = (
        proportional_hazard_test(cph, d, time_transform="rank")
        .summary["p"]
        .rename("proportional_hazards_p")
    )
    out = out.merge(ph, left_on="term", right_index=True, how="left")
    out.insert(0, "model", model)
    out.insert(1, "n", len(d))
    out.insert(2, "events", int(d["OS_event"].sum()))
    return out[
        [
            "model",
            "n",
            "events",
            "term",
            "coef",
            "exp(coef)",
            "exp(coef) lower 95%",
            "exp(coef) upper 95%",
            "p",
            "proportional_hazards_p",
        ]
    ]


def add_stage_covariates(frame: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    result = frame.copy()
    stage_dummies = pd.get_dummies(result["stage"], prefix="stage", dtype=float)
    stage_dummies.loc[result["stage"].isna(), :] = np.nan
    result = pd.concat([result, stage_dummies], axis=1)
    covariates = [
        column
        for column in ["stage_II", "stage_III", "stage_IV"]
        if column in result.columns
    ]
    return result, covariates


def run_tcga_survival() -> None:
    modules = pd.read_csv(SOURCE_ROOT / "TCGA_LUSC_selected_expression_modules.csv")
    mutation_status = pd.read_csv(
        SOURCE_ROOT / "TCGA_LUSC_NRF2_Pathway_Mutation_Status.csv"
    )
    survival = pd.read_csv(PROJECT_ROOT / "TCGA" / "TCGA-LUSC.survival.tsv", sep="\t")
    clinical = pd.read_csv(
        PROJECT_ROOT / "TCGA" / "TCGA-LUSC.clinical.tsv",
        sep="\t",
        low_memory=False,
    )

    survival = survival[
        survival["sample"].astype(str).str.contains(r"-01[A-Z]$", regex=True)
    ].copy()
    survival = survival.sort_values("OS.time", ascending=False).drop_duplicates(
        "_PATIENT"
    )
    survival = survival.rename(
        columns={"_PATIENT": "TCGA_patient", "OS.time": "OS_time", "OS": "OS_event"}
    )

    clinical = clinical[
        clinical["sample"].astype(str).str.contains(r"-01[A-Z]$", regex=True)
    ].copy()
    clinical["TCGA_patient"] = clinical["submitter_id"].astype(str)
    clinical = clinical.drop_duplicates("TCGA_patient")
    clinical["age"] = pd.to_numeric(
        clinical["age_at_index.demographic"], errors="coerce"
    )
    clinical["pack_years"] = pd.to_numeric(
        clinical["pack_years_smoked.exposures"], errors="coerce"
    )
    clinical["male"] = (
        clinical["gender.demographic"].astype(str).str.lower().eq("male").astype(float)
    )
    clinical["stage"] = clinical["ajcc_pathologic_stage.diagnoses"].map(parse_stage)
    clinical = clinical[
        ["TCGA_patient", "age", "pack_years", "male", "stage"]
    ].copy()

    merged = (
        modules.merge(
            survival[["TCGA_patient", "OS_time", "OS_event"]], on="TCGA_patient"
        )
        .merge(
            clinical,
            on="TCGA_patient",
            how="left",
        )
        .merge(
            mutation_status,
            on="TCGA_patient",
            how="left",
        )
        .copy()
    )
    merged["PARK7_z"] = (
        merged["PARK7"] - merged["PARK7"].mean()
    ) / merged["PARK7"].std(ddof=0)
    merged["NRF2_redox_z"] = (
        merged["NRF2_redox_marker_module"]
        - merged["NRF2_redox_marker_module"].mean()
    ) / merged["NRF2_redox_marker_module"].std(ddof=0)
    merged["age_per_10y"] = merged["age"] / 10.0
    merged["pack_years_per_10"] = merged["pack_years"] / 10.0
    merged, stage_covariates = add_stage_covariates(merged)

    outputs = [
        fit_cox(merged, ["PARK7_z"], "Univariable"),
        fit_cox(
            merged,
            ["PARK7_z", "age_per_10y", "male", *stage_covariates],
            "Clinical adjusted",
        ),
        fit_cox(
            merged,
            [
                "PARK7_z",
                "age_per_10y",
                "male",
                *stage_covariates,
                "pack_years_per_10",
                "NRF2_redox_z",
            ],
            "Clinical + smoking + NRF2/redox adjusted",
        ),
        fit_cox(
            merged,
            [
                "PARK7_z",
                "age_per_10y",
                "male",
                *stage_covariates,
                "pathway_mutated",
            ],
            "Clinical + KEAP1/NFE2L2/CUL3 mutation adjusted",
        ),
        fit_cox(
            merged,
            [
                "PARK7_z",
                "age_per_10y",
                "male",
                *stage_covariates,
                "pack_years_per_10",
                "NRF2_redox_z",
                "pathway_mutated",
            ],
            "Clinical + smoking + NRF2/redox + pathway mutation adjusted",
        ),
    ]
    pd.concat(outputs, ignore_index=True).to_csv(
        OUT / "TCGA_LUSC_PARK7_Cox_Models.csv", index=False
    )

    protein_wide = pd.read_csv(
        PROJECT_ROOT / "TCGA" / "TCGA-LUSC.protein.tsv",
        sep="\t",
    )
    dj1 = protein_wide.loc[
        protein_wide["peptide_target"].astype(str).eq("DJ1")
    ]
    if len(dj1) != 1:
        raise ValueError("Expected one DJ1 row in the TCGA-LUSC RPPA table.")
    protein = dj1.drop(columns="peptide_target").T.reset_index()
    protein.columns = ["sample", "DJ1_RPPA"]
    protein["TCGA_patient"] = protein["sample"].astype(str).str.slice(0, 12)
    protein["DJ1_RPPA"] = pd.to_numeric(protein["DJ1_RPPA"], errors="coerce")
    protein = (
        protein.sort_values("sample")
        .drop_duplicates("TCGA_patient")
        .merge(
            survival[["TCGA_patient", "OS_time", "OS_event"]],
            on="TCGA_patient",
        )
        .merge(clinical, on="TCGA_patient", how="left")
        .merge(mutation_status, on="TCGA_patient", how="left")
    )
    protein["DJ1_z"] = (
        protein["DJ1_RPPA"] - protein["DJ1_RPPA"].mean()
    ) / protein["DJ1_RPPA"].std(ddof=0)
    protein["age_per_10y"] = protein["age"] / 10.0
    protein["pack_years_per_10"] = protein["pack_years"] / 10.0
    protein, protein_stage_covariates = add_stage_covariates(protein)

    protein_outputs = [
        fit_cox(protein, ["DJ1_z"], "Univariable"),
        fit_cox(
            protein,
            ["DJ1_z", "age_per_10y", "male", *protein_stage_covariates],
            "Clinical adjusted",
        ),
        fit_cox(
            protein,
            [
                "DJ1_z",
                "age_per_10y",
                "male",
                *protein_stage_covariates,
                "pathway_mutated",
            ],
            "Clinical + KEAP1/NFE2L2/CUL3 mutation adjusted",
        ),
        fit_cox(
            protein,
            [
                "DJ1_z",
                "age_per_10y",
                "male",
                *protein_stage_covariates,
                "pack_years_per_10",
                "pathway_mutated",
            ],
            "Clinical + smoking + pathway mutation adjusted",
        ),
    ]
    pd.concat(protein_outputs, ignore_index=True).to_csv(
        OUT / "TCGA_LUSC_DJ1_RPPA_Cox_Models.csv",
        index=False,
    )

    matched = merged[["TCGA_patient", "PARK7"]].merge(
        protein[["TCGA_patient", "DJ1_RPPA"]],
        on="TCGA_patient",
    )
    matched = matched.dropna(subset=["PARK7", "DJ1_RPPA"])
    pearson = stats.pearsonr(matched["PARK7"], matched["DJ1_RPPA"])
    spearman = stats.spearmanr(matched["PARK7"], matched["DJ1_RPPA"])
    pd.DataFrame(
        [
            {
                "method": "Pearson",
                "n": len(matched),
                "estimate": pearson.statistic,
                "p": pearson.pvalue,
            },
            {
                "method": "Spearman",
                "n": len(matched),
                "estimate": spearman.statistic,
                "p": spearman.pvalue,
            },
        ]
    ).to_csv(
        OUT / "TCGA_LUSC_PARK7_RNA_DJ1_RPPA_Correlation.csv",
        index=False,
    )

    mutation_rows = [
        {
            "scope": "GDC masked-somatic-mutation coverage",
            "metric": "covered_cases",
            "value": len(mutation_status),
        },
        {
            "scope": "GDC masked-somatic-mutation coverage",
            "metric": "pathway_mutation_positive_cases",
            "value": int(mutation_status["pathway_mutated"].sum()),
        },
    ]
    for gene in ("KEAP1", "NFE2L2", "CUL3"):
        mutation_rows.append(
            {
                "scope": "GDC masked-somatic-mutation coverage",
                "metric": f"{gene}_mutation_positive_cases",
                "value": int(mutation_status[f"{gene}_mutated"].sum()),
            }
        )
    for scope, frame in [
        ("PARK7 RNA survival cohort", merged),
        ("DJ1 RPPA survival cohort", protein),
    ]:
        mutation_rows.extend(
            [
                {
                    "scope": scope,
                    "metric": "survival_cases",
                    "value": len(frame),
                },
                {
                    "scope": scope,
                    "metric": "mutation_covered_cases",
                    "value": int(frame["pathway_mutated"].notna().sum()),
                },
                {
                    "scope": scope,
                    "metric": "pathway_mutation_positive_cases",
                    "value": int(frame["pathway_mutated"].fillna(0).sum()),
                },
            ]
        )
    pd.DataFrame(mutation_rows).to_csv(
        OUT / "TCGA_LUSC_NRF2_Pathway_Mutation_Summary.csv",
        index=False,
    )

def export_coloc_posteriors() -> None:
    coloc = pd.read_csv(SOURCE_ROOT / "PARK7_Coloc_ABF_Prior_Sensitivity.csv")
    selected = coloc[
        coloc["cell"].isin(["bin", "bmem"])
        & coloc["histology"].isin(["LUSC", "LUAD"])
        & np.isclose(coloc["p1"], 1e-4)
        & np.isclose(coloc["p2"], 1e-4)
        & np.isclose(coloc["p12"], 1e-5)
    ].copy()
    selected["cell_type"] = selected["cell"].map(
        {"bin": "Intermediate B cells", "bmem": "Memory B cells"}
    )
    selected[
        [
            "cell_type",
            "histology",
            "n_shared_snps",
            "PP0",
            "PP1",
            "PP2",
            "PP3",
            "PP4",
            "PP4_over_PP3_plus_PP4",
        ]
    ].to_csv(OUT / "PARK7_Coloc_Complete_Posteriors.csv", index=False)

    locus_path = SOURCE_ROOT / "PARK7_LocusCompare_Source.csv"
    locus = pd.read_csv(locus_path) if locus_path.exists() else None
    provenance_rows = []
    for record in selected.to_dict(orient="records"):
        source_available = (
            locus is not None
            and record["cell"] == "bin"
            and record["histology"] == "LUSC"
        )
        provenance_rows.append(
            {
                "cell_type": (
                    "Intermediate B cells"
                    if record["cell"] == "bin"
                    else "Memory B cells"
                ),
                "histology": record["histology"],
                "n_shared_snps": int(record["n_shared_snps"]),
                "regional_variant_source_available": source_available,
                "chromosome": (
                    int(locus["CHR"].dropna().iloc[0])
                    if source_available
                    else np.nan
                ),
                "eqtl_position_min": (
                    int(locus["POS"].min()) if source_available else np.nan
                ),
                "eqtl_position_max": (
                    int(locus["POS"].max()) if source_available else np.nan
                ),
                "gwas_position_min": (
                    int(locus["pos"].min()) if source_available else np.nan
                ),
                "gwas_position_max": (
                    int(locus["pos"].max()) if source_available else np.nan
                ),
                "genome_build_labels": "Not retained in the regional source file",
                "boundary": (
                    "Coordinate columns are reported as retained; no new "
                    "cross-build regional plot was generated"
                    if source_available
                    else "Only derived posterior and shared-SNP count retained"
                ),
            }
        )
    pd.DataFrame(provenance_rows).to_csv(
        OUT / "PARK7_Coloc_Region_Provenance.csv",
        index=False,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Reproduce the PARK7 major-revision analysis tables."
    )
    parser.add_argument("--source-root", type=Path, default=SOURCE_ROOT)
    parser.add_argument("--project-root", type=Path, default=PROJECT_ROOT)
    parser.add_argument("--cosmx-root", type=Path, default=COSMX_ROOT)
    parser.add_argument("--out", type=Path, default=OUT)
    parser.add_argument(
        "--steps",
        nargs="+",
        choices=["ld", "cosmx", "visium", "scrna", "tcga", "coloc"],
        default=["visium", "scrna", "coloc"],
        help=(
            "Analyses to run. The default steps use inputs included in the "
            "repository; LD, CosMx, and TCGA require external resources."
        ),
    )
    return parser.parse_args()


def main() -> None:
    global SOURCE_ROOT, PROJECT_ROOT, COSMX_ROOT, OUT
    args = parse_args()
    SOURCE_ROOT = args.source_root.resolve()
    PROJECT_ROOT = args.project_root.resolve()
    COSMX_ROOT = args.cosmx_root.resolve()
    OUT = args.out.resolve()
    OUT.mkdir(parents=True, exist_ok=True)

    functions = {
        "ld": run_ld_audit,
        "cosmx": run_cosmx_sensitivity,
        "visium": run_visium_threshold_permutation,
        "scrna": run_scrna_threshold_sensitivity,
        "tcga": run_tcga_survival,
        "coloc": export_coloc_posteriors,
    }
    for step in args.steps:
        functions[step]()
    print(f"Major-revision analyses saved to: {OUT}")


if __name__ == "__main__":
    main()
