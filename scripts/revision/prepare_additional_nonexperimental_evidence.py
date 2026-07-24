from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[2]
INPUT_OUT = REPO_ROOT / "data" / "inputs"
RESULT_OUT = REPO_ROOT / "data" / "revision_results"
PATHWAY_GENES = ("KEAP1", "NFE2L2", "CUL3")
SPATIAL_TARGETS = (
    "PARK7",
    "CD79A",
    "MS4A1",
    "CD74",
    "CD37",
    "SOD1",
    "GPX1",
    "NFE2L2",
)
COSMX_SAMPLES = ("Lung6", "Lung9_Rep2", "Lung12", "Lung13", "Lung5_Rep3")


def export_restricted_mr_source(source: Path, result_out: Path) -> None:
    frame = pd.read_csv(source)
    required = {
        "cell",
        "gene",
        "n_iv",
        "beta_ivw",
        "se_ivw",
        "p_ivw",
        "q_bh_fdr",
    }
    if set(frame.columns) != required:
        raise ValueError(
            "The retained MR source table does not have the expected seven columns."
        )
    if len(frame) != 38:
        raise ValueError(f"Expected 38 retained rows, found {len(frame)}.")
    if frame["gene"].nunique() != 12 or frame["cell"].nunique() != 14:
        raise ValueError(
            "The retained source table must contain 12 genes and 14 cell labels."
        )
    frame.to_csv(
        result_out / "Restricted_38_Pair_MR_Source_Table.csv",
        index=False,
    )


def parse_gdc_events(occurrences: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for _, record in occurrences.iterrows():
        retained = False
        for index in range(14):
            prefix = f"ssm.consequence.{index}.transcript."
            canonical = str(record.get(prefix + "is_canonical", "")).lower()
            gene = record.get(prefix + "gene.symbol")
            impact = record.get(prefix + "annotation.vep_impact")
            consequence = record.get(prefix + "consequence_type")
            if (
                canonical == "true"
                and gene in PATHWAY_GENES
                and impact in {"HIGH", "MODERATE"}
            ):
                rows.append(
                    {
                        "TCGA_patient": record["case.submitter_id"],
                        "gene": gene,
                        "vep_impact": impact,
                        "consequence_type": consequence,
                        "ssm_id": record["ssm.ssm_id"],
                    }
                )
                retained = True
        if not retained:
            raise ValueError(
                "A supplied GDC occurrence had no qualifying canonical "
                "HIGH/MODERATE consequence."
            )
    events = pd.DataFrame(rows).drop_duplicates()
    return events.sort_values(
        ["TCGA_patient", "gene", "ssm_id"],
        kind="stable",
    ).reset_index(drop=True)


def export_gdc_mutation_status(
    occurrences_path: Path,
    coverage_path: Path,
    input_out: Path,
    access_date: str,
) -> None:
    occurrences = pd.read_csv(occurrences_path, sep="\t", low_memory=False)
    coverage = pd.read_csv(coverage_path, sep="\t", low_memory=False)
    events = parse_gdc_events(occurrences)

    covered = sorted(
        coverage["cases.0.submitter_id"].dropna().astype(str).unique()
    )
    positive = set(events["TCGA_patient"].astype(str))
    if not positive.issubset(covered):
        raise ValueError("Mutation-positive cases were absent from GDC MAF coverage.")

    status = pd.DataFrame({"TCGA_patient": covered})
    status["WXS_covered"] = 1
    for gene in PATHWAY_GENES:
        cases = set(events.loc[events["gene"].eq(gene), "TCGA_patient"])
        status[f"{gene}_mutated"] = status["TCGA_patient"].isin(cases).astype(int)
    status["pathway_mutated"] = (
        status[[f"{gene}_mutated" for gene in PATHWAY_GENES]].max(axis=1)
    )

    if len(status) != 490 or int(status["pathway_mutated"].sum()) != 134:
        raise ValueError(
            "Unexpected GDC case counts; expected 490 covered and 134 "
            "pathway-mutation-positive cases."
        )

    status.to_csv(
        input_out / "TCGA_LUSC_NRF2_Pathway_Mutation_Status.csv",
        index=False,
    )
    events.to_csv(
        input_out / "TCGA_LUSC_NRF2_Pathway_Mutation_Events.csv",
        index=False,
    )

    provenance = pd.DataFrame(
        [
            {
                "project": "TCGA-LUSC",
                "mutation_endpoint": "https://api.gdc.cancer.gov/ssm_occurrences",
                "coverage_endpoint": "https://api.gdc.cancer.gov/files",
                "genes": ";".join(PATHWAY_GENES),
                "event_definition": (
                    "Canonical-transcript somatic SSM with VEP impact "
                    "HIGH or MODERATE"
                ),
                "negative_definition": (
                    "No qualifying event among cases represented in the "
                    "official GDC masked-somatic-mutation file coverage"
                ),
                "coverage_workflow": (
                    "Aliquot Ensemble Somatic Variant Merging and Masking"
                ),
                "covered_cases": len(status),
                "mutation_positive_cases": int(status["pathway_mutated"].sum()),
                "qualifying_events": len(events),
                "access_date": access_date,
            }
        ]
    )
    provenance.to_csv(
        input_out / "TCGA_LUSC_NRF2_Pathway_Mutation_Provenance.csv",
        index=False,
    )


def decode_h5_values(values: np.ndarray) -> set[str]:
    decoded = []
    for value in values:
        if isinstance(value, bytes):
            decoded.append(value.decode("utf-8"))
        else:
            decoded.append(str(value))
    return set(decoded)


def export_spatial_feature_coverage(
    visium_h5: Path,
    cosmx_root: Path,
    result_out: Path,
) -> None:
    with h5py.File(visium_h5, "r") as handle:
        matrix = handle["matrix"]
        shape = tuple(int(value) for value in matrix["shape"][:])
        features = decode_h5_values(matrix["features"]["name"][:])
    if len(shape) != 2:
        raise ValueError(f"Unexpected Visium matrix shape: {shape}")

    rows: list[dict[str, object]] = []
    visium_row: dict[str, object] = {
        "resource": "10x Visium CytAssist FFPE LUSC",
        "sample_scope": "One public LUSC section",
        "feature_count": shape[0],
        "analysis_unit_count": shape[1],
        "analysis_unit": "under-tissue spots",
    }
    for target in SPATIAL_TARGETS:
        visium_row[f"{target}_measured"] = target in features
    rows.append(visium_row)

    panel_reference: set[str] | None = None
    for sample in COSMX_SAMPLES:
        expression = next(
            (cosmx_root / sample).rglob(f"{sample}_exprMat_file.csv"),
            None,
        )
        if expression is None:
            raise FileNotFoundError(f"CosMx expression matrix not found: {sample}")
        columns = pd.read_csv(expression, nrows=0).columns.tolist()
        panel = {
            column
            for column in columns[2:]
            if not str(column).startswith("NegPrb")
        }
        if panel_reference is None:
            panel_reference = panel
        elif panel != panel_reference:
            raise ValueError("The five retained CosMx samples use different panels.")

    assert panel_reference is not None
    cosmx_row: dict[str, object] = {
        "resource": "CosMx SMI NSCLC 960-gene panel",
        "sample_scope": ";".join(COSMX_SAMPLES),
        "feature_count": len(panel_reference),
        "analysis_unit_count": 501403,
        "analysis_unit": "matched segmented cells across five LUSC samples",
    }
    for target in SPATIAL_TARGETS:
        cosmx_row[f"{target}_measured"] = target in panel_reference
    rows.append(cosmx_row)

    audit = pd.DataFrame(rows)
    if audit["PARK7_measured"].any():
        raise ValueError("PARK7 was unexpectedly present in a spatial panel.")
    audit.to_csv(
        result_out / "Spatial_Feature_Coverage_Audit.csv",
        index=False,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Prepare retained source tables and coverage audits for the "
            "nonexperimental major-revision extension."
        )
    )
    parser.add_argument("--restricted-mr-source", type=Path, required=True)
    parser.add_argument("--visium-h5", type=Path, required=True)
    parser.add_argument("--cosmx-root", type=Path, required=True)
    parser.add_argument("--gdc-occurrences", type=Path, required=True)
    parser.add_argument("--gdc-coverage", type=Path, required=True)
    parser.add_argument("--input-out", type=Path, default=INPUT_OUT)
    parser.add_argument("--result-out", type=Path, default=RESULT_OUT)
    parser.add_argument("--access-date", default="2026-07-24")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    input_out = args.input_out.resolve()
    result_out = args.result_out.resolve()
    input_out.mkdir(parents=True, exist_ok=True)
    result_out.mkdir(parents=True, exist_ok=True)

    export_restricted_mr_source(
        args.restricted_mr_source.resolve(),
        result_out,
    )
    export_gdc_mutation_status(
        args.gdc_occurrences.resolve(),
        args.gdc_coverage.resolve(),
        input_out,
        args.access_date,
    )
    export_spatial_feature_coverage(
        args.visium_h5.resolve(),
        args.cosmx_root.resolve(),
        result_out,
    )
    print("Additional nonexperimental evidence tables prepared.")


if __name__ == "__main__":
    main()
