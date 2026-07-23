from __future__ import annotations

import argparse
import math
import os
import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle
from scipy import stats
from scipy.spatial import cKDTree
from statsmodels.nonparametric.smoothers_lowess import lowess


REPO_ROOT = Path(__file__).resolve().parents[2]
SOURCE_ROOT = Path(
    os.environ.get("PARK7_SOURCE_ROOT", REPO_ROOT / "data" / "inputs")
).resolve()
GENERATED = Path(
    os.environ.get("PARK7_RESULTS_ROOT", REPO_ROOT / "data" / "revision_results")
).resolve()
OUT = Path(
    os.environ.get("PARK7_FIGURE_OUT", REPO_ROOT / "figures" / "main")
).resolve()
OUT.mkdir(parents=True, exist_ok=True)


plt.rcParams.update(
    {
        "font.family": "Arial",
        "font.size": 8,
        "axes.titlesize": 9,
        "axes.labelsize": 8,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 7,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)


def save_figure(fig: plt.Figure, stem: str) -> None:
    fig.savefig(OUT / f"{stem}.png", dpi=400, bbox_inches="tight", facecolor="white")
    fig.savefig(
        OUT / f"{stem}.tiff",
        dpi=400,
        bbox_inches="tight",
        facecolor="white",
        pil_kwargs={"compression": "tiff_lzw"},
    )
    plt.close(fig)


def build_figure_1() -> None:
    rows = [
        (
            "Restricted candidate screen",
            "Regional PARK7 association in a supplied 38-pair subset; complete discovery universe and effect scale were unavailable.",
            "Exploratory only",
        ),
        (
            "LD audit",
            "90 unique SNPs; median pairwise r² = 0.645; all 3,916 matched EUR pairs had r² ≥ 0.1.",
            "Contradicts independence",
        ),
        (
            "Strict LD pruning",
            "rs6700772 only; OR 1.124 (95% CI 0.991–1.274); P = 0.069.",
            "Inconclusive",
        ),
        (
            "Colocalization",
            "LUSC PP4 = 0.367–0.377; PP1 had most posterior mass; results were prior-sensitive.",
            "Inconclusive",
        ),
        (
            "Primary Visium section",
            "O/E overlap 0.73; OR 0.62; Spearman ρ = −0.10; PARK7 was not measured.",
            "Against global overlap",
        ),
        (
            "Five-sample CosMx sensitivity",
            "CD74-excluded, myeloid-adjusted mean partial r = 0.146; four positive and one negative sample.",
            "Weak and heterogeneous",
        ),
        (
            "Single-cell and TCGA",
            "Broad PARK7 expression; score associations changed by threshold; no prognostic association.",
            "Contextual / mixed",
        ),
        (
            "Single-marker IHC",
            "PARK7 detected in 15 archived LUSC specimens; no B-cell co-staining was performed.",
            "Descriptive",
        ),
    ]
    status_colors = {
        "Exploratory only": "#E69F00",
        "Contradicts independence": "#D55E00",
        "Inconclusive": "#999999",
        "Against global overlap": "#D55E00",
        "Weak and heterogeneous": "#E69F00",
        "Contextual / mixed": "#56B4E9",
        "Descriptive": "#009E73",
    }

    fig, ax = plt.subplots(figsize=(7.2, 5.4), dpi=200)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    fig.suptitle(
        "What each analysis contributes to the PARK7 assessment",
        fontsize=13,
        fontweight="bold",
        y=0.98,
    )
    ax.text(
        0.5,
        0.91,
        "The analyses address distinct questions and do not form a progressive validation chain",
        ha="center",
        va="center",
        fontsize=8.5,
        color="#4B5563",
    )

    wrapped_rows = [
        [
            textwrap.fill(row[0], width=21),
            textwrap.fill(row[1], width=55),
            textwrap.fill(row[2], width=22),
        ]
        for row in rows
    ]
    table = ax.table(
        cellText=wrapped_rows,
        colLabels=["Analysis", "Observed result and limit", "Evidence direction"],
        cellLoc="left",
        colLoc="center",
        colWidths=[0.22, 0.52, 0.26],
        bbox=[0.01, 0.09, 0.98, 0.77],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(7)
    for (row_index, col_index), cell in table.get_celld().items():
        cell.set_edgecolor("#D1D5DB")
        cell.set_linewidth(0.6)
        cell.PAD = 0.06
        if row_index == 0:
            cell.set_facecolor("#264653")
            cell.get_text().set_color("white")
            cell.get_text().set_fontweight("bold")
            cell.get_text().set_ha("center")
        else:
            cell.set_facecolor("#F8FAFC" if row_index % 2 else "#FFFFFF")
            cell.get_text().set_va("center")
            if col_index == 2:
                status = rows[row_index - 1][2]
                cell.set_facecolor(status_colors[status])
                cell.get_text().set_color("white")
                cell.get_text().set_fontweight("bold")
                cell.get_text().set_ha("center")

    ax.text(
        0.01,
        0.03,
        "Overall interpretation: PARK7 remains a provisional, hypothesis-generating candidate.",
        fontsize=8,
        fontweight="bold",
        color="#24313D",
    )
    save_figure(fig, "Figure1_Cross_Modal_Assessment")


def build_figure_2() -> None:
    ld = np.loadtxt(GENERATED / "PARK7_primary_unique_EUR_LD.ld")
    audit = pd.read_csv(GENERATED / "MR_PARK7_LD_Audit_Summary.csv").iloc[0]
    strict = pd.read_csv(GENERATED / "MR_PARK7_LD_Pruned_Estimates.csv")
    strict = strict[np.isclose(strict["clump_r2"].astype(float), 0.001)].iloc[0]

    unpruned = {
        "OR": math.exp(0.10018209650404473),
        "low": math.exp(0.10018209650404473 - 1.96 * 0.007492100887441334),
        "high": math.exp(0.10018209650404473 + 1.96 * 0.007492100887441334),
    }

    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(7.2, 3.35), dpi=200, gridspec_kw={"width_ratios": [1.15, 1]}
    )
    fig.suptitle(
        "PARK7 regional instruments are highly correlated and yield an inconclusive strict-pruned estimate",
        fontsize=11.2,
        fontweight="bold",
        y=1.01,
    )

    im = ax1.imshow(ld, cmap="viridis", vmin=0, vmax=1, interpolation="nearest")
    ax1.set_title("A. Pairwise LD among 89 EUR-matched variants", loc="left", fontweight="bold")
    ax1.set_xlabel("Variant index")
    ax1.set_ylabel("Variant index")
    cbar = fig.colorbar(im, ax=ax1, fraction=0.047, pad=0.03)
    cbar.set_label("r²")
    ax1.text(
        0.02,
        0.98,
        (
            f"Median r² = {float(audit['median_pairwise_r2']):.3f}\n"
            f"Pairs with r² ≥ 0.1: {int(audit['pairs_r2_ge_0_1']):,}/"
            f"{int(audit['total_pairs']):,}\n"
            f"Pairs with r² ≥ 0.8: {int(audit['pairs_r2_ge_0_8']):,}"
        ),
        transform=ax1.transAxes,
        va="top",
        fontsize=6.8,
        color="white",
        bbox=dict(facecolor="#111827", edgecolor="none", alpha=0.78, pad=4),
    )

    labels = [
        "90-SNP diagnostic\n(LD unaccounted)",
        "Strict-pruned rs6700772\n(r² < 0.001)",
    ]
    ors = np.array([unpruned["OR"], float(strict["OR"])])
    lows = np.array([unpruned["low"], float(strict["CI_low"])])
    highs = np.array([unpruned["high"], float(strict["CI_high"])])
    y = np.array([1, 0])
    colors = ["#999999", "#0072B2"]
    for yi, or_value, lo, hi, color in zip(y, ors, lows, highs, colors):
        ax2.errorbar(
            or_value,
            yi,
            xerr=[[or_value - lo], [hi - or_value]],
            fmt="o",
            markersize=6,
            color=color,
            ecolor=color,
            elinewidth=1.8,
            capsize=3,
        )
        ax2.text(
            hi + 0.008,
            yi,
            f"{or_value:.3f} ({lo:.3f}–{hi:.3f})",
            va="center",
            fontsize=6.8,
            color=color,
        )
    ax2.axvline(1.0, color="#4B5563", linestyle="--", linewidth=1)
    ax2.set_yticks(y, labels)
    ax2.set_xlim(0.95, 1.31)
    ax2.set_ylim(-0.65, 1.65)
    ax2.set_xlabel("Odds ratio for LUSC")
    ax2.set_title("B. PARK7 effect estimates", loc="left", fontweight="bold")
    ax2.grid(axis="x", linestyle=":", alpha=0.3)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    ax2.text(
        0.02,
        0.04,
        "Only the strict-pruned estimate is used for inference.\n"
        f"Wald P = {float(strict['p']):.3f}",
        transform=ax2.transAxes,
        fontsize=6.8,
        color="#24313D",
    )
    fig.tight_layout()
    save_figure(fig, "Figure2_PARK7_LD_Audit")


def build_figure_3() -> None:
    post = pd.read_csv(GENERATED / "PARK7_Coloc_Complete_Posteriors.csv")
    post["label"] = post["cell_type"] + " / " + post["histology"]
    order = [
        "Intermediate B cells / LUSC",
        "Memory B cells / LUSC",
        "Intermediate B cells / LUAD",
        "Memory B cells / LUAD",
    ]
    post["label"] = pd.Categorical(post["label"], categories=order, ordered=True)
    post = post.sort_values("label")

    sensitivity = pd.read_csv(
        SOURCE_ROOT / "PARK7_Coloc_ABF_Prior_Sensitivity.csv"
    )
    sensitivity = sensitivity[
        sensitivity["cell"].isin(["bin", "bmem"])
        & sensitivity["histology"].eq("LUSC")
        & np.isclose(sensitivity["p1"], 1e-4)
        & np.isclose(sensitivity["p2"], 1e-4)
    ].copy()
    sensitivity["cell_type"] = sensitivity["cell"].map(
        {"bin": "Intermediate B cells", "bmem": "Memory B cells"}
    )

    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(7.2, 3.3), dpi=200, gridspec_kw={"width_ratios": [1.35, 1]}
    )
    fig.suptitle(
        "PARK7 colocalization is inconclusive at the prespecified prior and sensitive to p12",
        fontsize=11.2,
        fontweight="bold",
        y=1.01,
    )

    cols = ["PP0", "PP1", "PP2", "PP3", "PP4"]
    colors = ["#BDBDBD", "#56B4E9", "#009E73", "#E69F00", "#D55E00"]
    left = np.zeros(len(post))
    y = np.arange(len(post))
    for col, color in zip(cols, colors):
        values = post[col].to_numpy(float)
        ax1.barh(y, values, left=left, color=color, edgecolor="white", linewidth=0.4, label=col)
        left += values
    ax1.set_yticks(y, post["label"].astype(str))
    ax1.invert_yaxis()
    ax1.set_xlim(0, 1)
    ax1.set_xlabel("Posterior probability")
    ax1.set_title("A. Complete PP0–PP4 posterior distribution", loc="left", fontweight="bold")
    ax1.legend(ncol=5, frameon=False, loc="lower center", bbox_to_anchor=(0.5, -0.31))
    for yi, (_, row) in enumerate(post.iterrows()):
        ax1.text(
            0.99,
            yi,
            f"PP1 {row['PP1']:.3f} | PP3 {row['PP3']:.3f} | PP4 {row['PP4']:.3f}",
            ha="right",
            va="center",
            fontsize=5.7,
            color="#111827",
        )

    for cell_type, group in sensitivity.groupby("cell_type"):
        group = group.sort_values("p12")
        ax2.plot(
            group["p12"],
            group["PP4"],
            marker="o",
            linewidth=1.8,
            markersize=4.5,
            label=cell_type,
        )
    ax2.axhline(0.8, color="#4B5563", linestyle="--", linewidth=1, label="PP4 = 0.8")
    ax2.set_xscale("log")
    ax2.set_xticks([1e-6, 1e-5, 1e-4], [r"$10^{-6}$", r"$10^{-5}$", r"$10^{-4}$"])
    ax2.set_ylim(0, 1)
    ax2.set_xlabel("Prior probability p12")
    ax2.set_ylabel("PP4")
    ax2.set_title("B. LUSC prior sensitivity", loc="left", fontweight="bold")
    ax2.grid(axis="y", linestyle=":", alpha=0.35)
    ax2.legend(frameon=False, loc="upper left", fontsize=6.5)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    fig.tight_layout()
    save_figure(fig, "Figure3_PARK7_Colocalization_Posteriors")


def build_figure_4() -> None:
    data = pd.read_csv(SOURCE_ROOT / "Visium_Spatial_CoLocalization_Scores.csv")
    for col in ["B_high", "A_high", "Both_high"]:
        if data[col].dtype != bool:
            data[col] = data[col].astype(str).str.lower().eq("true")

    b_n = int(data["B_high"].sum())
    a_n = int(data["A_high"].sum())
    both_n = int(data["Both_high"].sum())
    expected = b_n * a_n / len(data)
    observed_expected = both_n / expected
    table = pd.crosstab(data["B_high"], data["A_high"]).reindex(
        index=[False, True], columns=[False, True], fill_value=0
    )
    odds_ratio, _ = stats.fisher_exact(table.to_numpy())
    rho_score, p_score = stats.spearmanr(
        data["B_cell_density_z"], data["Antioxidant_z"]
    )

    if "dist_to_b" not in data:
        coords = data[["array_col", "array_row"]].to_numpy(float)
        b_coords = data.loc[data["B_high"], ["array_col", "array_row"]].to_numpy(float)
        data["dist_to_b"] = cKDTree(b_coords).query(coords, k=1)[0]
    rho_dist, p_dist = stats.spearmanr(data["dist_to_b"], data["Antioxidant_z"])
    smooth = lowess(
        data["Antioxidant_z"],
        data["dist_to_b"] + np.linspace(0, 1e-7, len(data)),
        frac=0.28,
        return_sorted=True,
    )

    fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.35), dpi=200)
    fig.suptitle(
        "Spatial B-cell-marker and antioxidant scores in one Visium CytAssist LUSC section",
        fontsize=11.5,
        fontweight="bold",
        y=0.995,
    )

    ax = axes[0, 0]
    sc1 = ax.scatter(
        data["array_col"],
        data["array_row"],
        c=data["B_cell_density_z"],
        cmap="Blues",
        s=9,
        alpha=0.9,
        edgecolors="none",
    )
    ax.invert_yaxis()
    ax.set_aspect("equal")
    ax.set_title("A. B-cell-marker score (CD79A, MS4A1)", loc="left", fontweight="bold")
    ax.set_xlabel("Array column")
    ax.set_ylabel("Array row")
    cb = fig.colorbar(sc1, ax=ax, fraction=0.047, pad=0.025)
    cb.set_label("Score (z)")

    ax = axes[0, 1]
    sc2 = ax.scatter(
        data["array_col"],
        data["array_row"],
        c=data["Antioxidant_z"],
        cmap="YlOrBr",
        s=9,
        alpha=0.82,
        edgecolors="none",
    )
    ax.invert_yaxis()
    ax.set_aspect("equal")
    ax.set_title("B. Antioxidant score (SOD1, NFE2L2)", loc="left", fontweight="bold")
    ax.set_xlabel("Array column")
    ax.set_ylabel("Array row")
    cb = fig.colorbar(sc2, ax=ax, fraction=0.047, pad=0.025)
    cb.set_label("Score (z)")
    ax.text(
        0.98,
        0.98,
        "PARK7 was not measured",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=6.4,
        bbox=dict(facecolor="white", edgecolor="#D1D5DB", alpha=0.9, pad=3),
    )

    ax = axes[1, 0]
    background = data[~(data["B_high"] | data["A_high"])]
    b_only = data[data["B_high"] & ~data["A_high"]]
    a_only = data[data["A_high"] & ~data["B_high"]]
    both = data[data["Both_high"]]
    ax.scatter(background["array_col"], background["array_row"], c="#D1D5DB", s=7, alpha=0.14)
    ax.scatter(b_only["array_col"], b_only["array_row"], c="#0072B2", s=9, alpha=0.22)
    ax.scatter(a_only["array_col"], a_only["array_row"], c="#E69F00", s=9, alpha=0.20)
    ax.scatter(
        both["array_col"],
        both["array_row"],
        c="#CC79A7",
        s=19,
        alpha=0.95,
        edgecolors="white",
        linewidths=0.15,
    )
    ax.invert_yaxis()
    ax.set_aspect("equal")
    ax.set_xlabel("Array column")
    ax.set_ylabel("Array row")
    ax.set_title("C. High-score spot overlap", loc="left", fontweight="bold")
    ax.text(
        0.02,
        0.98,
        (
            f"B-marker-high: {b_n:,}; antioxidant-high: {a_n:,}\n"
            f"Observed joint-high: {both_n:,}; expected: {expected:.1f}\n"
            f"O/E = {observed_expected:.2f}; OR = {odds_ratio:.2f}\n"
            f"Global Spearman ρ = {rho_score:.2f}"
        ),
        transform=ax.transAxes,
        va="top",
        fontsize=6.4,
        bbox=dict(facecolor="white", edgecolor="#D1D5DB", alpha=0.92, pad=3),
    )

    ax = axes[1, 1]
    ax.scatter(
        data["dist_to_b"],
        data["Antioxidant_z"],
        c="#94A3B8",
        s=7,
        alpha=0.17,
        edgecolors="none",
        rasterized=True,
    )
    ax.plot(smooth[:, 0], smooth[:, 1], color="#CC79A7", linewidth=2.2)
    ax.set_xlabel("Distance to nearest B-marker-high spot (array-coordinate units)")
    ax.set_ylabel("Antioxidant score (z)")
    ax.set_title("D. Distance association", loc="left", fontweight="bold")
    ax.text(
        0.03,
        0.97,
        f"Spearman ρ = {rho_dist:.3f}\nP = {p_dist:.2e}\nNegligible positive association",
        transform=ax.transAxes,
        va="top",
        fontsize=6.4,
        bbox=dict(facecolor="white", edgecolor="#D1D5DB", alpha=0.92, pad=3),
    )
    ax.grid(axis="y", linestyle=":", alpha=0.3)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout(rect=[0, 0, 1, 0.975])
    save_figure(fig, "Figure4_Visium_Spatial_Constraint")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Rebuild the four main PARK7 major-revision figures."
    )
    parser.add_argument("--source-root", type=Path, default=SOURCE_ROOT)
    parser.add_argument("--results-root", type=Path, default=GENERATED)
    parser.add_argument("--out", type=Path, default=OUT)
    return parser.parse_args()


def main() -> None:
    global SOURCE_ROOT, GENERATED, OUT
    args = parse_args()
    SOURCE_ROOT = args.source_root.resolve()
    GENERATED = args.results_root.resolve()
    OUT = args.out.resolve()
    OUT.mkdir(parents=True, exist_ok=True)
    build_figure_1()
    build_figure_2()
    build_figure_3()
    build_figure_4()
    print(f"Revised figures saved to {OUT}")


if __name__ == "__main__":
    main()
