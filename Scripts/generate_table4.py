import pandas as pd
import numpy as np

# Paths
INPUT_FILE = "d:/2026YJ/My_MR_Project/CMap_Results/Real_CMap_Top_Reversers.csv"
OUTPUT_FILE = "d:/2026YJ/My_MR_Project/Table4_CMap_Translation.csv"

def generate_table4():
    # 1. Load Data
    df = pd.read_csv(INPUT_FILE)
    
    # 2. Filter Top 7 Reversers (Most Negative Connectivity Score)
    # Sort by Connectivity Score (ascending)
    df_sorted = df.sort_values(by='Connectivity Score', ascending=True)
    top7_df = df_sorted.head(7).copy()
    
    # 3. Rename Columns
    # Drug/Mechanism Class: id or src_set_id? 
    # The 'id' column contains full string "CLASS:CELL:TRT". 'src_set_id' seems to be the class name.
    # Let's use src_set_id for cleaner name.
    
    cols_map = {
        'src_set_id': 'Drug/Mechanism Class',
        'cell_iname': 'Cell Line',
        'Connectivity Score': 'Connectivity Score',
        'fdr_q_nlog10': 'Significance (-log10(FDR))',
        'set_size': 'N_total', # Total sets tested? Or size of this set? Usually set_size is number of perturbations in the class.
        'num_hiq_sig': 'N_sig' # Number of high quality significant instances
    }
    
    # Check if columns exist
    for col in cols_map.keys():
        if col not in top7_df.columns:
            print(f"Warning: Column {col} not found in input.")
            
    top7_df = top7_df.rename(columns=cols_map)
    
    # 4. Select Final Columns
    final_cols = ['Drug/Mechanism Class', 'Cell Line', 'Connectivity Score', 'Significance (-log10(FDR))', 'N_total', 'N_sig']
    final_df = top7_df[final_cols]
    
    # 5. Format Numbers
    # Connectivity Score: 3 decimals
    # Significance: 2 decimals
    
    final_df['Connectivity Score'] = final_df['Connectivity Score'].apply(lambda x: f"{x:.3f}")
    final_df['Significance (-log10(FDR))'] = final_df['Significance (-log10(FDR))'].apply(lambda x: f"{x:.2f}")
    
    # Clean Class Name (remove underscores)
    final_df['Drug/Mechanism Class'] = final_df['Drug/Mechanism Class'].str.replace('_', ' ')
    
    # 6. Save
    final_df.to_csv(OUTPUT_FILE, index=False)
    print(f"Table 4 saved to {OUTPUT_FILE}")
    print(final_df)

if __name__ == "__main__":
    generate_table4()
