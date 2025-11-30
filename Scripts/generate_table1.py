import pandas as pd
import time
import os

# Paths
MR_FILE = "d:/2026YJ/My_MR_Project/Final_Classified_Results.csv"
TCGA_FILE = "d:/2026YJ/My_MR_Project/TCGA_Results/TCGA_Batch_Survival_Summary.csv"
OUTPUT_FILE = "d:/2026YJ/My_MR_Project/Table1_MR_Summary.csv"

def generate_table1():
    # 1. Load Data
    print("Loading MR results...")
    mr_df = pd.read_csv(MR_FILE)
    
    print("Loading TCGA results...")
    tcga_df = pd.read_csv(TCGA_FILE)
    
    # 2. Filter & Sort MR Results to get Top Candidates
    evidence_order = {"Strong": 0, "Moderate": 1, "Suggestive": 2, "Weak": 3}
    mr_df['Evidence_Rank'] = mr_df['Evidence_Level'].map(evidence_order)
    
    # Sort: Rank (asc), p_ivw (asc), pp4_ilcco (desc)
    mr_df = mr_df.sort_values(by=['Evidence_Rank', 'p_ivw', 'pp4_ilcco'], ascending=[True, True, False])
    
    top_genes_df = mr_df.drop_duplicates(subset=['gene'], keep='first').head(15)
    
    # 3. Format Table 1 Columns
    table_rows = []
    
    for _, row in top_genes_df.iterrows():
        gene = row['gene']
        cell = row['cell']
        mr_p = row['p_ivw']
        coloc_pp4 = row['pp4_ilcco']
        
        # TCGA Lookup
        tcga_hits = tcga_df[tcga_df['Gene'] == gene]
        tcga_str = "NS"
        if not tcga_hits.empty:
            res_strs = []
            for _, t_row in tcga_hits.iterrows():
                cohort = t_row['Cohort'].replace("TCGA-", "")
                hr = t_row['HR']
                p_val = t_row['PValue']
                sig_mark = "*" if p_val < 0.05 else ""
                res_strs.append(f"{cohort}: {hr:.2f} ({p_val:.3f}){sig_mark}")
            tcga_str = " | ".join(res_strs)
        
        # Format numbers
        if mr_p == 0:
            mr_p_str = "< 1.00e-10"
        elif mr_p < 0.001:
            mr_p_str = f"{mr_p:.2e}"
        else:
            mr_p_str = f"{mr_p:.3f}"
            
        table_rows.append({
            "Gene": gene,
            "Cell Type": cell,
            "MR P-value": mr_p_str,
            "Coloc PP4": f"{coloc_pp4:.3f}",
            "TCGA HR (P-value)": tcga_str,
            "Evidence Level": row['Evidence_Level']
        })
        
    result_df = pd.DataFrame(table_rows)
    
    # Save with retry logic
    max_retries = 5
    for i in range(max_retries):
        try:
            result_df.to_csv(OUTPUT_FILE, index=False)
            print(f"Table 1 saved to {OUTPUT_FILE}")
            print(result_df)
            break
        except PermissionError:
            if i < max_retries - 1:
                print(f"File locked, retrying in 1s... ({i+1}/{max_retries})")
                time.sleep(1)
            else:
                print("Could not save file. Please close it if it's open.")
                # Save to a new file as backup
                backup_file = OUTPUT_FILE.replace(".csv", "_v2.csv")
                result_df.to_csv(backup_file, index=False)
                print(f"Saved to backup: {backup_file}")

if __name__ == "__main__":
    generate_table1()
