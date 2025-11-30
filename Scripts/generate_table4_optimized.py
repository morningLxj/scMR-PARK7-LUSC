import pandas as pd
import numpy as np

# Paths
INPUT_FILE = "d:/2026YJ/My_MR_Project/CMap_Results/Real_CMap_Top_Reversers.csv"
OUTPUT_FILE = "d:/2026YJ/My_MR_Project/Table4_CMap_Translation_Optimized.csv"

def generate_table4():
    # 1. Load Data
    df = pd.read_csv(INPUT_FILE)
    
    # 2. Filter Top 7 Reversers
    df_sorted = df.sort_values(by='Connectivity Score', ascending=True)
    top7_df = df_sorted.head(7).copy()
    
    # 3. Rename Columns
    cols_map = {
        'src_set_id': 'Drug/Mechanism Class',
        'cell_iname': 'Cell Line',
        'Connectivity Score': 'Connectivity Score',
        'fdr_q_nlog10': 'Significance (-log10(FDR))',
        'set_size': 'N_total', 
        'num_hiq_sig': 'N_sig'
    }
    top7_df = top7_df.rename(columns=cols_map)
    
    # 4. Select Final Columns
    final_cols = ['Drug/Mechanism Class', 'Cell Line', 'Connectivity Score', 'Significance (-log10(FDR))', 'N_total', 'N_sig']
    final_df = top7_df[final_cols]
    
    # 5. Format Numbers
    # Connectivity Score: 4 decimals (-0.9996)
    # Significance: 3 decimals (0.583)
    
    final_df['Connectivity Score'] = final_df['Connectivity Score'].apply(lambda x: f"{x:.4f}")
    final_df['Significance (-log10(FDR))'] = final_df['Significance (-log10(FDR))'].apply(lambda x: f"{x:.3f}")
    
    # 6. Clean Class Name - Use underscores as requested in example, or keep as is?
    # User example: SODIUM/CALCIUM_EXCHANGE_INHIBITOR
    # Input has underscores. Previous optimization replaced them with spaces.
    # User input shows underscores: "SODIUM/CALCIUM_EXCHANGE_INHIBITOR"
    # So we KEEP the underscores this time, or ensure they are present.
    # The raw data has underscores.
    
    # 7. Save
    final_df.to_csv(OUTPUT_FILE, index=False)
    print(f"Table 4 saved to {OUTPUT_FILE}")
    print(final_df)

if __name__ == "__main__":
    generate_table4()
