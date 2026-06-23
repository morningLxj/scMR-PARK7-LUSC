from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


EXPECTED = {
    "main_figures": 5,
    "supplementary_figures": 19,
    "source_data_files": 61,
    "analysis_scripts": 3,
}


REQUIRED_FILES = [
    "README.md",
    "requirements.txt",
    "environment.yml",
    "docs/DATA_AVAILABILITY.md",
    "docs/REPRODUCIBILITY_NOTES.md",
    "data/source_data_csv/Source_Data_CSV_Manifest.csv",
    "data/source_data_csv/Table1_MR_Summary_Source.csv",
    "data/source_data_csv/SupplementaryTableS8_Evidence_Reproducibility_Matrix.csv",
    "data/source_workbooks/Source_Data_and_Tables.xlsx",
    "figures/main/Figure_1.tif",
    "figures/main/Figure_5.tif",
    "figures/supplementary/S1.tif",
    "figures/supplementary/S19_Network_scRNA_Control_Context.tif",
    "scripts/analysis/01_public_scrna_tcga_ihc_cmap_workflow.py",
    "scripts/analysis/02_extended_public_validation_and_mr_sensitivity.py",
    "scripts/analysis/03_evidence_boundary_audit_figures.py",
]


def count_files(path: str, pattern: str = "*") -> int:
    return sum(1 for p in (ROOT / path).glob(pattern) if p.is_file())


def main() -> int:
    missing = [rel for rel in REQUIRED_FILES if not (ROOT / rel).exists()]
    counts = {
        "main_figures": count_files("figures/main", "*.tif"),
        "supplementary_figures": count_files("figures/supplementary", "*.tif"),
        "source_data_files": count_files("data/source_data_csv"),
        "analysis_scripts": count_files("scripts/analysis", "*.py"),
    }

    print("Repository integrity check")
    print(f"Root: {ROOT}")
    for key, value in counts.items():
        expected = EXPECTED[key]
        status = "OK" if value == expected else f"expected {expected}"
        print(f"- {key}: {value} ({status})")

    if missing:
        print("\nMissing required files:")
        for rel in missing:
            print(f"- {rel}")
        return 1

    bad_counts = {k: v for k, v in counts.items() if v != EXPECTED[k]}
    if bad_counts:
        return 1

    print("\nAll required files are present.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
