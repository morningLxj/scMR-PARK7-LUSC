import pandas as pd
import numpy as np

# Paths
INPUT_FILE = "d:/2026YJ/My_MR_Project/Proteome_MR_Results_Final.csv"
OUTPUT_FILE = "d:/2026YJ/My_MR_Project/Table3_Protein_Validation_Optimized.csv"

def generate_table3():
    # 1. Load Data
    df = pd.read_csv(INPUT_FILE)
    
    # 2. Filter: Keep only PARK7 and CTSW
    target_genes = ["PARK7", "CTSW"]
    df = df[df['exposure'].str.contains('|'.join(target_genes))]
    
    # 3. Columns mapping
    # Protein Exposure | Study Cohort | MR Method | NSNP | Beta | SE | OR(95%CI) | P-value
    
    # Clean exposure name
    # PARK7 (Protein) -> PARK7 (Protein)
    # Keep as is or just ensure consistent format
    
    # Clean Outcome
    outcome_map = {
        "Total": "Total",
        "Adeno": "Adeno",
        "Squamous": "Squamous"
    }
    df['Study Cohort'] = df['outcome'].map(outcome_map).fillna(df['outcome'])
    
    # Format numbers
    
    def format_beta(x):
        return f"{x:.3f}"
        
    def format_se(x):
        return f"{x:.3f}"
        
    def format_or_ci(row):
        or_val = row['or']
        lci = row['or_lci95']
        uci = row['or_uci95']
        
        # Special handling for large numbers
        if or_val > 100:
            return f"{or_val:.3f} ({lci:.3f}-{uci:.3f})"
        else:
            return f"{or_val:.3f} ({lci:.3f}-{uci:.3f})"

    def format_pval(p):
        # Use LaTeX-like scientific notation for very small numbers?
        # User requested: 2.79 x 10^-1 style or similar
        # Let's use standard scientific notation first, user can format to latex if needed.
        # Python's 'e' notation is standard for CSV.
        # But user example: 1.06 x 10^-7
        
        if p < 0.001:
            # Convert 1.06e-07 to 1.06 x 10^-7
            s = f"{p:.2e}"
            base, exponent = s.split('e')
            return f"{base} x 10^{int(exponent)}"
        else:
            # For larger P, e.g. 0.279 -> 2.79 x 10^-1 or just 0.279?
            # User example: 0.279 -> 2.79 x 10^-1
            s = f"{p:.2e}"
            base, exponent = s.split('e')
            return f"{base} x 10^{int(exponent)}"

    df['Beta'] = df['b'].apply(format_beta)
    df['SE'] = df['se'].apply(format_se)
    df['OR(95%CI)'] = df.apply(format_or_ci, axis=1)
    df['P-value'] = df['pval'].apply(format_pval)
    
    # Rename Columns
    cols_map = {
        'exposure': 'Protein Exposure',
        'Study Cohort': 'Study Cohort',
        'method': 'MR Method',
        'nsnp': 'NSNP',
        'Beta': 'Beta',
        'SE': 'SE',
        'OR(95%CI)': 'OR(95%CI)',
        'P-value': 'P-value'
    }
    
    final_df = df.rename(columns=cols_map)
    
    # Select final columns
    final_cols = ['Protein Exposure', 'Study Cohort', 'MR Method', 'NSNP', 'Beta', 'SE', 'OR(95%CI)', 'P-value']
    final_df = final_df[final_cols]
    
    # Sort
    # Order: PARK7 Total, CTSW Total, PARK7 Adeno, CTSW Adeno...
    # Custom sort
    final_df['Gene_Rank'] = final_df['Protein Exposure'].apply(lambda x: 0 if 'PARK7' in x else 1)
    cohort_order = {'Total': 0, 'Adeno': 1, 'Squamous': 2}
    final_df['Cohort_Rank'] = final_df['Study Cohort'].map(cohort_order)
    
    final_df = final_df.sort_values(by=['Cohort_Rank', 'Gene_Rank'])
    final_df = final_df.drop(columns=['Gene_Rank', 'Cohort_Rank'])
    
    # Save
    final_df.to_csv(OUTPUT_FILE, index=False)
    print(f"Table 3 saved to {OUTPUT_FILE}")
    print(final_df)

if __name__ == "__main__":
    generate_table3()
