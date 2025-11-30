import pandas as pd

# Paths
INPUT_FILE = "d:/2026YJ/My_MR_Project/Final_Subtype_Specificity.csv"
OUTPUT_FILE = "d:/2026YJ/My_MR_Project/Table2_PARK7_LUSC_Specificity.csv"

def generate_table2():
    # 1. Load Data
    df = pd.read_csv(INPUT_FILE)
    
    # 2. Filter for PARK7 and LUSC-Specific
    # Condition: gene == 'PARK7' AND Specificity == 'LUSC-Specific'
    # Actually, user said "mechanism stratification", so we might want to show ALL PARK7 rows?
    # User said: "精确筛选出 gene = 'PARK7' 且 Specificity = 'LUSC-Specific' 的所有数据行"
    # So we only keep those.
    
    target_df = df[(df['gene'] == 'PARK7') & (df['Specificity'] == 'LUSC-Specific')].copy()
    
    # 3. Sort by Delta (descending)
    target_df = target_df.sort_values(by='Delta', ascending=False)
    
    # 4. Rename Columns
    # gene -> Gene
    # cell -> Cell Type
    # PP4_LUAD -> Coloc PP4 (LUAD)
    # PP4_LUSC -> Coloc PP4 (LUSC)
    # Delta -> Delta (LUSC - LUAD)
    # Specificity -> Subtype Specificity
    
    rename_map = {
        'gene': 'Gene',
        'cell': 'Cell Type',
        'PP4_LUAD': 'Coloc PP4 (LUAD)',
        'PP4_LUSC': 'Coloc PP4 (LUSC)',
        'Delta': 'Delta Score',
        'Specificity': 'Classification'
    }
    target_df = target_df.rename(columns=rename_map)
    
    # 5. Format Numbers (4 decimal places)
    cols_to_format = ['Coloc PP4 (LUAD)', 'Coloc PP4 (LUSC)', 'Delta Score']
    
    for col in cols_to_format:
        target_df[col] = target_df[col].apply(lambda x: f"{x:.4f}" if pd.notnull(x) else "-")
        
    # 6. Save
    target_df.to_csv(OUTPUT_FILE, index=False)
    print(f"Table 2 saved to {OUTPUT_FILE}")
    print(target_df)

if __name__ == "__main__":
    generate_table2()
