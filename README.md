# Paper Code Repository

This repository contains the source code, processed data, and results for the paper.

## Repository Structure

### `Scripts/`
Contains all Python (`.py`) and R (`.R`) scripts used for data analysis, Mendelian Randomization (MR), and figure generation.
- **Key Scripts:**
  - `run_mr_*.py/R`: Main Mendelian Randomization analysis pipelines.
  - `generate_figure*.py`: Scripts to generate the figures in the manuscript.
  - `generate_table*.py`: Scripts to generate the tables in the manuscript.
  - `run_pp4_*.R`: Colocalization analysis scripts.

### `Data/`
Contains key processed data files and classification results.
- `Final_Classified_Results.csv`: Final classification of risk genes.
- `Final_Gold_Standard.csv`: Gold standard gene set used for validation.
- `Final_Subtype_Specificity.csv`: Subtype specificity analysis results.
- `candidates.txt`: List of candidate genes analyzed.

### `Tables/`
Contains the final tables used in the manuscript.
- `Table1_MR_Summary.csv`
- `Table2_PARK7_LUSC_Specificity.csv`
- `Table3_Protein_Validation.csv`
- `Table4_CMap_Translation.csv`
- `TableS1_Full_MR_Statistics.csv`: Supplementary table with full statistics.

### `Figures/`
Contains the high-resolution figures (PDF and PNG) generated for the manuscript.

## Usage

1.  **Prerequisites**: Ensure you have Python 3.x and R 4.x installed.
2.  **Data**: Place any large raw data files (e.g., GWAS summary statistics, huge eQTL tables) in a `RawData/` directory (not included here due to size limits) if you intend to re-run the full pipeline.
3.  **Running Analysis**:
    - The scripts in `Scripts/` are generally named according to their function.
    - Start with `run_mr_python.py` for the main MR analysis.
    - Use `generate_figure*.py` to reproduce specific figures.

## License

[Specify License here, e.g., MIT]
