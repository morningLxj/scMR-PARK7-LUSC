import pandas as pd
import numpy as np

# Paths
INPUT_FILE = "d:/2026YJ/My_MR_Project/Proteome_MR_Results_Final.csv"
OUTPUT_FILE = "d:/2026YJ/My_MR_Project/Table3_Protein_Validation.csv"

def generate_table3():
    # 1. Load Data
    df = pd.read_csv(INPUT_FILE)
    
    # 2. Filter: Keep only PARK7 and CTSW (Exposure)
    # Based on file content, exposure column has "PARK7 (Protein)" and "CTSW (Protein)"
    target_genes = ["PARK7", "CTSW"]
    # Use string contains to be safe
    df = df[df['exposure'].str.contains('|'.join(target_genes))]
    
    # 3. Select Columns
    # We need: Exposure (Protein), Outcome (Cohort), OR, 95% CI, P-value
    # Input cols: exposure, outcome, or, or_lci95, or_uci95, pval
    
    # 4. Format OR (95% CI)
    def format_or_ci(row):
        or_val = row['or']
        lci = row['or_lci95']
        uci = row['or_uci95']
        
        # If OR is huge (like for CTSW Adeno), keep it readable but scientific if needed?
        # User said "对于极小的 P-value 使用科学计数法，OR 和 CI 保留 3 位小数（对于极大的 OR，保持其原始精度）"
        # Let's check magnitude.
        
        if or_val > 1000:
             return f"{or_val:.1e} ({lci:.1e}-{uci:.1e})"
        else:
             return f"{or_val:.3f} ({lci:.3f}-{uci:.3f})"

    df['OR (95% CI)'] = df.apply(format_or_ci, axis=1)
    
    # 5. Format P-value
    def format_pval(p):
        if p < 0.001:
            return f"{p:.2e}"
        else:
            return f"{p:.3f}"
            
    df['P_MR'] = df['pval'].apply(format_pval)
    
    # 6. Rename and Reorder
    # Gene (Protein) | Outcome Cohort | OR (95% CI) | P_MR
    
    # Clean exposure name: "PARK7 (Protein)" -> "PARK7"
    df['Gene'] = df['exposure'].str.replace(r" \(Protein\)", "", regex=True)
    
    # Clean outcome: "Total" -> "Lung Cancer (Total)", "Adeno" -> "Lung Adenocarcinoma (LUAD)", "Squamous" -> "Lung Squamous Cell Carcinoma (LUSC)"
    outcome_map = {
        "Total": "Lung Cancer (Total)",
        "Adeno": "Lung Adenocarcinoma (LUAD)",
        "Squamous": "Lung Squamous Cell Carcinoma (LUSC)"
    }
    df['Outcome Cohort'] = df['outcome'].map(outcome_map).fillna(df['outcome'])
    
    final_cols = ['Gene', 'Outcome Cohort', 'OR (95% CI)', 'P_MR']
    final_df = df[final_cols]
    
    # Sort by Gene then Outcome
    final_df = final_df.sort_values(by=['Gene', 'Outcome Cohort'])
    
    # 7. Save
    final_df.to_csv(OUTPUT_FILE, index=False)
    print(f"Table 3 saved to {OUTPUT_FILE}")
    print(final_df)

if __name__ == "__main__":
    generate_table3()
