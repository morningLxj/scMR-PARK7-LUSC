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
    "data/inputs/TCGA_LUSC_NRF2_Pathway_Mutation_Status.csv",
    "data/inputs/TCGA_LUSC_NRF2_Pathway_Mutation_Events.csv",
    "data/inputs/TCGA_LUSC_NRF2_Pathway_Mutation_Provenance.csv",
    "data/revision_results/MR_PARK7_LD_Audit_Summary.csv",
    "data/revision_results/MR_PARK7_LD_Pruned_Estimates.csv",
    "data/revision_results/MR_PARK7_Outcome_Provenance.csv",
    "data/revision_results/PARK7_Coloc_Complete_Posteriors.csv",
    "data/revision_results/PARK7_Coloc_Region_Provenance.csv",
    "data/revision_results/Spatial_CosMx_CD74_Myeloid_Sensitivity.csv",
    "data/revision_results/Spatial_Dataset_Sample_Specification.csv",
    "data/revision_results/Spatial_Feature_Coverage_Audit.csv",
    "data/revision_results/Visium_Threshold_BlockPermutation_Sensitivity.csv",
    "data/revision_results/scRNA_PARK7_Threshold_Sensitivity_Meta.csv",
    "data/revision_results/Restricted_38_Pair_MR_Source_Table.csv",
    "data/revision_results/TCGA_LUSC_PARK7_Cox_Models.csv",
    "data/revision_results/TCGA_LUSC_DJ1_RPPA_Cox_Models.csv",
    "data/revision_results/TCGA_LUSC_PARK7_RNA_DJ1_RPPA_Correlation.csv",
    "data/revision_results/TCGA_LUSC_NRF2_Pathway_Mutation_Summary.csv",
    "figures/main/Figure1_Cross_Modal_Assessment.png",
    "figures/main/Figure2_PARK7_LD_Audit.png",
    "figures/main/Figure3_PARK7_Colocalization_Posteriors.png",
    "figures/main/Figure4_Visium_Spatial_Constraint.png",
    "figures/supplementary/Supplementary_Figure_S1_PARK7_IHC.tif",
    "scripts/revision/run_major_revision_analyses.py",
    "scripts/revision/prepare_additional_nonexperimental_evidence.py",
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

    restricted_mr = read_rows(
        "data/revision_results/Restricted_38_Pair_MR_Source_Table.csv"
    )
    assert len(restricted_mr) == 38
    assert len({row["gene"] for row in restricted_mr}) == 12
    assert len({row["cell"] for row in restricted_mr}) == 14

    spatial_coverage = read_rows(
        "data/revision_results/Spatial_Feature_Coverage_Audit.csv"
    )
    assert len(spatial_coverage) == 2
    assert all(
        row["PARK7_measured"].lower() == "false"
        for row in spatial_coverage
    )
    visium_coverage = [
        row for row in spatial_coverage if row["resource"].startswith("10x Visium")
    ][0]
    cosmx_coverage = [
        row for row in spatial_coverage if row["resource"].startswith("CosMx")
    ][0]
    assert int(visium_coverage["feature_count"]) == 18085
    assert int(visium_coverage["analysis_unit_count"]) == 3858
    assert int(cosmx_coverage["feature_count"]) == 960
    assert int(cosmx_coverage["analysis_unit_count"]) == 501403

    mutation_status = read_rows(
        "data/inputs/TCGA_LUSC_NRF2_Pathway_Mutation_Status.csv"
    )
    assert len(mutation_status) == 490
    assert sum(int(row["pathway_mutated"]) for row in mutation_status) == 134
    assert sum(int(row["KEAP1_mutated"]) for row in mutation_status) == 51
    assert sum(int(row["NFE2L2_mutated"]) for row in mutation_status) == 69
    assert sum(int(row["CUL3_mutated"]) for row in mutation_status) == 22

    rna_models = [
        row
        for row in read_rows(
            "data/revision_results/TCGA_LUSC_PARK7_Cox_Models.csv"
        )
        if row["term"] == "PARK7_z"
    ]
    assert len(rna_models) == 5
    full_rna = [
        row
        for row in rna_models
        if row["model"]
        == "Clinical + smoking + NRF2/redox + pathway mutation adjusted"
    ][0]
    assert int(full_rna["n"]) == 399
    assert close(float(full_rna["exp(coef)"]), 0.976, tolerance=1e-3)

    protein_models = [
        row
        for row in read_rows(
            "data/revision_results/TCGA_LUSC_DJ1_RPPA_Cox_Models.csv"
        )
        if row["term"] == "DJ1_z"
    ]
    assert len(protein_models) == 4
    univariable_protein = [
        row for row in protein_models if row["model"] == "Univariable"
    ][0]
    assert int(univariable_protein["n"]) == 323
    assert close(
        float(univariable_protein["exp(coef)"]),
        1.026,
        tolerance=1e-3,
    )

    rna_protein = read_rows(
        "data/revision_results/TCGA_LUSC_PARK7_RNA_DJ1_RPPA_Correlation.csv"
    )
    pearson = [row for row in rna_protein if row["method"] == "Pearson"][0]
    assert int(pearson["n"]) == 320
    assert close(float(pearson["estimate"]), 0.526, tolerance=1e-3)

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
