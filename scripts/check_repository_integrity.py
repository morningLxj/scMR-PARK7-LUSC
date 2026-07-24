from __future__ import annotations

import csv
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

REQUIRED = [
    "README.md",
    "data/README.md",
    "data/inputs/PARK7_PerInstrument_F_Statistics.csv",
    "data/inputs/PARK7_StrictPruned_Outcome_Associations.csv",
    "data/inputs/PARK7_LeaveOneOut_Source.csv",
    "data/inputs/GSE148071_LUSC_scRNA_cell_level_selected_scores.csv",
    "data/inputs/PARK7_Coloc_ABF_Prior_Sensitivity.csv",
    "data/inputs/Visium_Spatial_CoLocalization_Scores.csv",
    "data/revision_results/MR_PARK7_LD_Audit_Summary.csv",
    "data/revision_results/MR_PARK7_LD_Pruned_Estimates.csv",
    "data/revision_results/MR_PARK7_Outcome_Provenance.csv",
    "data/revision_results/PARK7_Coloc_Complete_Posteriors.csv",
    "data/revision_results/PARK7_Coloc_Region_Provenance.csv",
    "data/revision_results/Spatial_CosMx_CD74_Myeloid_Sensitivity.csv",
    "data/revision_results/Spatial_Dataset_Sample_Specification.csv",
    "data/revision_results/Visium_Threshold_BlockPermutation_Sensitivity.csv",
    "data/revision_results/scRNA_PARK7_Threshold_Sensitivity_Meta.csv",
    "data/revision_results/TCGA_LUSC_PARK7_Cox_Models.csv",
    "figures/main/Figure1_Cross_Modal_Assessment.png",
    "figures/main/Figure2_PARK7_LD_Audit.png",
    "figures/main/Figure3_PARK7_Colocalization_Posteriors.png",
    "figures/main/Figure4_Visium_Spatial_Constraint.png",
    "figures/supplementary/Supplementary_Figure_S1_PARK7_IHC.tif",
    "scripts/revision/run_major_revision_analyses.py",
    "scripts/revision/build_revision_figures.py",
]


def read_rows(relative: str) -> list[dict[str, str]]:
    with (ROOT / relative).open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def close(actual: float, expected: float, tolerance: float = 5e-4) -> bool:
    return math.isclose(actual, expected, rel_tol=tolerance, abs_tol=tolerance)


def main() -> None:
    missing = [path for path in REQUIRED if not (ROOT / path).is_file()]
    empty = [
        path
        for path in REQUIRED
        if (ROOT / path).is_file() and (ROOT / path).stat().st_size == 0
    ]
    if missing or empty:
        raise SystemExit(f"Missing files: {missing}; empty files: {empty}")

    forbidden_suffixes = {".doc", ".docx", ".pdf", ".ppt", ".pptx", ".zip"}
    forbidden = [
        path.relative_to(ROOT).as_posix()
        for path in ROOT.rglob("*")
        if path.is_file()
        and ".git" not in path.parts
        and path.suffix.lower() in forbidden_suffixes
    ]
    if forbidden:
        raise SystemExit(f"Submission/private files found in repository: {forbidden}")

    audit = read_rows(
        "data/revision_results/MR_PARK7_LD_Audit_Summary.csv"
    )[0]
    assert int(float(audit["requested_instruments"])) == 90
    assert int(float(audit["instruments_in_1000G_EUR"])) == 89
    assert close(float(audit["median_pairwise_r2"]), 0.645)
    assert int(float(audit["pairs_r2_ge_0_1"])) == 3916

    pruned = read_rows(
        "data/revision_results/MR_PARK7_LD_Pruned_Estimates.csv"
    )
    assert all(int(float(row["n_instruments"])) == 1 for row in pruned)
    primary_lusc = [
        row for row in pruned
        if row["outcome_id"] == "GCST004750"
        and close(float(row["clump_r2"]), 0.001)
    ]
    overall_lung = [
        row for row in pruned
        if row["outcome_id"] == "GCST004748"
        and close(float(row["clump_r2"]), 0.001)
    ]
    assert len(primary_lusc) == 1
    assert len(overall_lung) == 1
    assert close(float(primary_lusc[0]["OR"]), 1.228, tolerance=1e-3)
    assert close(float(overall_lung[0]["OR"]), 1.124, tolerance=1e-3)

    coloc = read_rows(
        "data/revision_results/PARK7_Coloc_Complete_Posteriors.csv"
    )
    lusc = [row for row in coloc if row["histology"] == "LUSC"]
    assert len(lusc) == 2
    assert all(float(row["PP1"]) > float(row["PP4"]) for row in lusc)

    coloc_region = read_rows(
        "data/revision_results/PARK7_Coloc_Region_Provenance.csv"
    )
    assert len(coloc_region) == 4
    assert sum(
        row["regional_variant_source_available"].lower() == "true"
        for row in coloc_region
    ) == 1

    spatial = read_rows(
        "data/revision_results/Spatial_CosMx_CD74_Myeloid_Sensitivity.csv"
    )
    assert len(spatial) == 5

    spatial_specification = read_rows(
        "data/revision_results/Spatial_Dataset_Sample_Specification.csv"
    )
    assert len(spatial_specification) == 6
    cosmx_specification = [
        row for row in spatial_specification if row["resource"].startswith("Bruker")
    ]
    assert sum(int(row["matched_analysis_units"]) for row in cosmx_specification) == 501403
    assert sum(
        int(row["retained_aggregation_units"]) for row in cosmx_specification
    ) == 16438

    visium_sensitivity = read_rows(
        "data/revision_results/Visium_Threshold_BlockPermutation_Sensitivity.csv"
    )
    assert len(visium_sensitivity) == 6
    assert all(float(row["observed_over_expected"]) < 1 for row in visium_sensitivity)
    primary_visium = [
        row for row in visium_sensitivity
        if row["threshold"] == "global_q75_primary"
    ][0]
    assert int(primary_visium["both_high_observed"]) == 126
    assert close(
        float(primary_visium["both_high_expected_independence"]),
        172.962,
        tolerance=1e-3,
    )
    assert close(
        float(primary_visium["block_8_p_lower_or_equal"]),
        0.0111,
        tolerance=1e-3,
    )

    oversized = [
        path.relative_to(ROOT).as_posix()
        for path in ROOT.rglob("*")
        if path.is_file() and path.stat().st_size >= 100 * 1024 * 1024
    ]
    if oversized:
        raise SystemExit(f"Files at or above GitHub's 100 MB limit: {oversized}")

    print(
        "Repository integrity check passed: required files, claim-critical "
        "values, privacy exclusions, and size limits are consistent."
    )


if __name__ == "__main__":
    main()
