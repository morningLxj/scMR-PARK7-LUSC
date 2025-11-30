
import pandas as pd

# Input and Output paths
input_csv = r"D:\2026YJ\My_MR_Project\B_cell_PARK7_High_vs_Low_DEGs_Real.csv"
output_csv = r"D:\2026YJ\论文写作\spring模板\TableS4_SingleCell_DEGs.csv"

try:
    # Read the file
    df = pd.read_csv(input_csv)
    
    # Add Description column based on gene symbols
    # Mapping based on our known gene functions
    gene_desc = {
        "PARK7": "Antioxidant",
        "CD79B": "B cell receptor",
        "MS4A1": "B cell marker (CD20)",
        "CD79A": "B cell receptor",
        "GPX1": "Glutathione peroxidase",
        "SOD1": "Superoxide dismutase",
        "CAT": "Hydrogen peroxide catabolism",
        "CD40": "Costimulatory molecule",
        "CD19": "B cell marker",
        "XRCC1": "DNA repair",
        "CHEK1": "Cell cycle checkpoint",
        "BRCA1": "DNA repair",
        "GSR": "Glutathione metabolism",
        "TXN": "Redox regulation",
        "ATM": "DNA damage response",
        "ATR": "DNA damage response",
        "NFE2L2": "Oxidative stress response",
        "HMOX1": "Heme catabolism",
        "NQO1": "Quinone reductase"
    }
    
    df['Description'] = df['gene'].map(gene_desc).fillna("Uncharacterized")
    
    # Reorder and Rename columns to match Table S4 style
    # gene, avg_log2FC, p_val_adj, pct.1, pct.2, Description
    df_final = df[['gene', 'avg_log2FC', 'p_val_adj', 'pct.1', 'pct.2', 'Description']]
    df_final.columns = ['Gene', 'Avg_Log2FC', 'P_Val_Adj', 'Pct.1', 'Pct.2', 'Description']
    
    # Sort by P_Val_Adj
    df_final = df_final.sort_values('P_Val_Adj')
    
    # Save
    df_final.to_csv(output_csv, index=False)
    print(f"Successfully generated Table S4 CSV at: {output_csv}")
    
except Exception as e:
    print(f"Error generating Table S4 CSV: {e}")
