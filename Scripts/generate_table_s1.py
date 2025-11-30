import pandas as pd
import numpy as np
from scipy import stats
from pathlib import Path

# Paths
BASE_DIR = Path("d:/2026YJ/My_MR_Project")
OUTCOME_FILE = BASE_DIR / "Outcome/28604730-GCST004748-EFO_0001071.h.tsv"
FINAL_RES_FILE = BASE_DIR / "Final_Classified_Results.csv"
OUTPUT_FILE = BASE_DIR / "TableS1_Full_MR_Statistics.csv"

def calculate_mr_metrics(gene, cell, df_out_all):
    print(f"\nCalculating metrics for {gene} in {cell}...", flush=True)
    # Load Exposure
    exp_file = BASE_DIR / f"Exposure/EXPOSURE_{cell}_{gene}_mr.tsv.gz"
    if not exp_file.exists():
        print(f"Error: Exposure file does not exist: {exp_file}", flush=True)
        return None
    
    try:
        df_exp = pd.read_csv(exp_file, sep='\t', compression='gzip')
        print(f"Successfully loaded exposure file: {exp_file}", flush=True)
        print(f"Exposure dataframe shape: {df_exp.shape}", flush=True)
        print(f"Exposure dataframe columns: {list(df_exp.columns)}", flush=True)
    except Exception as e:
        print(f"Error reading exposure file {exp_file}: {e}", flush=True)
        import traceback
        traceback.print_exc()
        return None
    
    # Check if SNP column exists in exposure
    if 'SNP' not in df_exp.columns:
        print(f"Error: SNP column not found in exposure file {exp_file}", flush=True)
        return None
    
    # Filter Outcome from pre-loaded DF
    target_snps = set(df_exp['SNP'])
    print(f"Target SNPs from exposure: {len(target_snps)}", flush=True)
    
    df_out = df_out_all[df_out_all['SNP'].isin(target_snps)].copy()
    print(f"Matching outcome SNPs found: {len(df_out)}", flush=True)
    
    if df_out.empty:
        print(f"Error: No matching outcome SNPs found for {gene} in {cell}", flush=True)
        return None
        
    # Merge
    try:
        df = pd.merge(df_exp, df_out, on='SNP', how='inner')
        print(f"Merge successful, merged dataframe shape: {df.shape}", flush=True)
    except Exception as e:
        print(f"Error merging exposure and outcome data: {e}", flush=True)
        return None
    
    # Harmonize
    print(f"Merged dataframe columns: {list(df.columns)}", flush=True)
    
    # Use the correct column names after merge
    exposure_effect_allele = 'effect_allele_x' if 'effect_allele_x' in df.columns else 'effect_allele'
    exposure_beta = 'beta_x' if 'beta_x' in df.columns else 'beta'
    
    if exposure_effect_allele not in df.columns:
        print(f"Error: Exposure effect_allele column not found in merged dataframe", flush=True)
        return None
    if 'ea_outcome' not in df.columns:
        print(f"Error: ea_outcome column not found in merged dataframe", flush=True)
        return None
    
    idx = df[exposure_effect_allele] != df['ea_outcome']
    if idx.any():
        print(f"Harmonizing {idx.sum()} SNPs with mismatched alleles", flush=True)
        df.loc[idx, 'by'] = -df.loc[idx, 'by']
    
    # Calc MR Metrics
    print("Converting columns to numeric...", flush=True)
    df['by'] = pd.to_numeric(df['by'], errors='coerce')
    df['se_by'] = pd.to_numeric(df['se_by'], errors='coerce')
    df[exposure_beta] = pd.to_numeric(df[exposure_beta], errors='coerce')
    
    df = df.dropna(subset=[exposure_beta, 'by', 'se_by'])
    print(f"After dropping NA values, {len(df)} SNPs remain", flush=True)
    
    if len(df) < 3:
        print(f"Error: Not enough SNPs for MR analysis (need at least 3, got {len(df)})")
        return None
    
    bx = df[exposure_beta].values
    by = df['by'].values
    se_by = df['se_by'].values
    
    # IVW
    print("Calculating IVW...", flush=True)
    weights = 1 / (se_by**2)
    ratio = by / bx
    ivw_beta = np.sum(weights * ratio) / np.sum(weights)
    ivw_se = np.sqrt(1 / np.sum(weights))
    ivw_p = 2 * (1 - stats.norm.cdf(abs(ivw_beta / ivw_se)))
    
    # Egger
    print("Calculating Egger...", flush=True)
    sum_w = np.sum(weights)
    sum_wx = np.sum(weights * bx)
    sum_wy = np.sum(weights * by)
    sum_wxx = np.sum(weights * bx**2)
    sum_wxy = np.sum(weights * bx * by)
    denom = sum_w * sum_wxx - sum_wx**2
    
    if denom == 0:
        print("Error: Denominator is zero in Egger calculation", flush=True)
        return None
    
    egger_slope = (sum_w * sum_wxy - sum_wx * sum_wy) / denom
    egger_intercept = (sum_wxx * sum_wy - sum_wx * sum_wxy) / denom
    
    # Egger SEs
    y_pred = egger_intercept + egger_slope * bx
    resid = by - y_pred
    rss = np.sum(weights * resid**2)
    n = len(bx)
    sigma_sq = rss / (n - 2)
    
    se_egger_slope = np.sqrt(sigma_sq * sum_w / denom)
    p_egger_slope = 2 * (1 - stats.t.cdf(abs(egger_slope / se_egger_slope), df=n-2))
    
    se_egger_int = np.sqrt(sigma_sq * sum_wxx / denom)
    p_egger_int = 2 * (1 - stats.t.cdf(abs(egger_intercept / se_egger_int), df=n-2))
    
    # Weighted Median
    print("Calculating Weighted Median...", flush=True)
    ratio = by / bx
    w = weights * (bx**2) 
    # Sort by ratio
    order = np.argsort(ratio)
    ratio_s = ratio[order]
    w_s = w[order]
    cum_w = np.cumsum(w_s)
    total_w = cum_w[-1]
    idx_median = np.searchsorted(cum_w, 0.5 * total_w)
    median_beta = ratio_s[idx_median]
    
    # Calculate SE_Median and P_Median using Bootstrap
    print("Calculating Weighted Median SE and P-value...", flush=True)
    n_bootstrap = 1000
    bootstrap_medians = []
    
    for _ in range(n_bootstrap):
        # Resample with replacement
        indices = np.random.choice(len(ratio), size=len(ratio), replace=True)
        bootstrap_ratio = ratio[indices]
        bootstrap_w = w[indices]
        
        # Calculate bootstrap weighted median
        bootstrap_order = np.argsort(bootstrap_ratio)
        bootstrap_ratio_s = bootstrap_ratio[bootstrap_order]
        bootstrap_w_s = bootstrap_w[bootstrap_order]
        bootstrap_cum_w = np.cumsum(bootstrap_w_s)
        bootstrap_total_w = bootstrap_cum_w[-1]
        bootstrap_idx_median = np.searchsorted(bootstrap_cum_w, 0.5 * bootstrap_total_w)
        bootstrap_median = bootstrap_ratio_s[bootstrap_idx_median]
        bootstrap_medians.append(bootstrap_median)
    
    # Calculate SE_Median as the standard deviation of bootstrap medians
    se_median = np.std(bootstrap_medians)
    
    # Calculate P_Median
    if se_median > 0:
        p_median = 2 * (1 - stats.norm.cdf(abs(median_beta / se_median)))
    else:
        p_median = np.nan
    
    # Q Statistic (IVW)
    print("Calculating Q Statistic...", flush=True)
    weights_ratio = (bx**2) / (se_by**2)
    q_stat = np.sum(weights_ratio * (ratio - ivw_beta)**2)
    p_q = 1 - stats.chi2.cdf(q_stat, df=n-1)
    
    # Calculate I²
    print("Calculating I²...", flush=True)
    df_chi2 = n - 1
    if q_stat > 0:
        i_squared = ((q_stat - df_chi2) / q_stat) * 100
        # Ensure I² is within valid range [0, 100]
        i_squared = max(0, min(100, i_squared))
    else:
        i_squared = 0
    
    # Calculate F Statistic
    print("Calculating F Statistic...", flush=True)
    df_exp['beta'] = pd.to_numeric(df_exp['beta'], errors='coerce')
    # Check if se column exists in exposure
    if 'se' not in df_exp.columns:
        print(f"Warning: se column not found in exposure file, using default F_stat", flush=True)
        f_stat = np.nan
    else:
        df_exp['se'] = pd.to_numeric(df_exp['se'], errors='coerce')
        df_exp = df_exp.dropna(subset=['beta', 'se'])
        df_exp = df_exp[df_exp['se'] > 0]
        
        f_stat = np.nan
        if not df_exp.empty:
            df_exp['F_stat'] = (df_exp['beta'] / df_exp['se']) ** 2
            f_stat = df_exp['F_stat'].mean()
    
    print(f"Successfully calculated all metrics for {gene} in {cell}", flush=True)
    
    return {
        "N_SNP": n,
        "Beta_IVW": ivw_beta, "SE_IVW": ivw_se, "P_IVW": ivw_p,
        "Beta_Egger": egger_slope, "SE_Egger": se_egger_slope, "P_Egger": p_egger_slope,
        "Beta_Median": median_beta, "SE_Median": se_median, "P_Median": p_median,
        "P_Q": p_q, "P_Egger_Int": p_egger_int, "F_Stat": f_stat, "I_squared": i_squared
    }

def generate_table_s1():
    # 1. Get list of gene/cell pairs from Final Results
    res_df = pd.read_csv(FINAL_RES_FILE)
    
    # Optimization: Load outcome file ONCE for all target SNPs
    all_target_snps = set()
    gene_cell_list = []
    for i, row in res_df.iterrows():
        gene = row['gene']
        cell = row['cell']
        gene_cell_list.append((gene, cell))
        # Peak at exposure file to get SNPs
        exp_file = BASE_DIR / f"Exposure/EXPOSURE_{cell}_{gene}_mr.tsv.gz"
        if exp_file.exists():
             try:
                df_tmp = pd.read_csv(exp_file, sep='\t', compression='gzip', usecols=['SNP'])
                all_target_snps.update(df_tmp['SNP'])
             except:
                pass
    
    print(f"Total unique SNPs to fetch: {len(all_target_snps)}", flush=True)
    
    # Load Outcome for ALL SNPs
    outcome_chunks = []
    try:
        print(f"Attempting to read outcome file: {OUTCOME_FILE}", flush=True)
        # First, check if the file exists
        if not OUTCOME_FILE.exists():
            print(f"Error: Outcome file does not exist: {OUTCOME_FILE}", flush=True)
            df_out_all = pd.DataFrame()
        else:
            # Read the first few lines to check columns
            df_first = pd.read_csv(OUTCOME_FILE, sep='\t', nrows=5)
            print(f"Outcome file columns: {list(df_first.columns)}", flush=True)
            print(f"Outcome file sample data: {df_first}", flush=True)
            
            # Now read in chunks
            reader = pd.read_csv(OUTCOME_FILE, sep='\t', chunksize=50000)
            for i, chunk in enumerate(reader):
                if i % 100 == 0:
                    print(f"Reading outcome chunk {i}...", flush=True)
                
                # Try to rename columns
                try:
                    chunk = chunk.rename(columns={
                        'hm_rsid': 'SNP',
                        'hm_effect_allele': 'ea_outcome',
                        'hm_other_allele': 'nea_outcome',
                        'hm_beta': 'by',
                        'standard_error': 'se_by',
                        'p_value': 'p_outcome'
                    })
                except Exception as e:
                    print(f"Error renaming columns in chunk {i}: {e}", flush=True)
                    continue
                
                # Check if SNP column exists after renaming
                if 'SNP' not in chunk.columns:
                    print(f"Error: SNP column not found in chunk {i}", flush=True)
                    continue
                
                # Filter for target SNPs
                mask = chunk['SNP'].isin(all_target_snps)
                matched_count = mask.sum()
                print(f"Chunk {i}: Found {matched_count} matching SNPs out of {len(chunk)} total SNPs", flush=True)
                
                if matched_count > 0:
                    outcome_chunks.append(chunk[mask])
            
            if outcome_chunks:
                df_out_all = pd.concat(outcome_chunks)
                print(f"Successfully concatenated {len(outcome_chunks)} chunks")
                print(f"Final outcome dataframe shape: {df_out_all.shape}")
                print(f"Final outcome dataframe columns: {list(df_out_all.columns)}")
            else:
                df_out_all = pd.DataFrame()
                print("No matching SNPs found in any chunk")
    except Exception as e:
        print(f"Error loading outcome: {e}", flush=True)
        import traceback
        traceback.print_exc()
        df_out_all = pd.DataFrame()

    print(f"Loaded {len(df_out_all)} outcome SNPs.", flush=True)

    total_pairs = len(res_df)
    print(f"Processing {total_pairs} pairs for Table S1...", flush=True)
    
    rows = []
    for i, row in res_df.iterrows():
        gene = row['gene']
        cell = row['cell']
        
        if i % 5 == 0:
            print(f"Processing {i+1}/{total_pairs}: {gene} in {cell}...", flush=True)
            
        # Pass pre-loaded outcome df
        metrics = calculate_mr_metrics(gene, cell, df_out_all)
        
        if metrics:
            res_row = {
                "Cell": cell,
                "Gene": gene,
                "N_SNP": metrics["N_SNP"],
                "Beta_IVW": metrics["Beta_IVW"],
                "SE_IVW": metrics["SE_IVW"],
                "P_IVW": metrics["P_IVW"],
                "Beta_Egger": metrics["Beta_Egger"],
                "SE_Egger": metrics["SE_Egger"],
                "P_Egger": metrics["P_Egger"],
                "Beta_Median": metrics["Beta_Median"],
                "SE_Median": metrics["SE_Median"],
                "P_Median": metrics["P_Median"],
                "P_Q": metrics["P_Q"],
                "P_Egger_Int": metrics["P_Egger_Int"],
                "F_Stat": metrics["F_Stat"],
                "I_squared": metrics["I_squared"]
            }
            rows.append(res_row)
        else:
            print(f"Calculation failed for {gene} in {cell}", flush=True)
            pass

    if not rows:
        print("No rows generated!", flush=True)
        return

    final_df = pd.DataFrame(rows)
    
    # Sort by Gene then Cell
    final_df = final_df.sort_values(['Gene', 'Cell'])
    
    # Format numbers
    float_cols = [c for c in final_df.columns if c not in ['Cell', 'Gene', 'N_SNP']]
    
    for col in float_cols:
        if col.startswith('P_') or col == 'P_Egger_Int':
            # P-values formatting with consistent scientific notation
            final_df[col] = final_df[col].apply(lambda x: 
                "< 1e-30" if pd.notnull(x) and x < 1e-30 else 
                "> 0.99" if pd.notnull(x) and x > 0.99 else 
                (f"{x:.2e}".replace("e-0", "e-").replace("e+0", "e+") if pd.notnull(x) else "NA"))
        elif col == 'I_squared':
            # I_squared formatting: display 2 decimal places
            final_df[col] = final_df[col].apply(lambda x: f"{x:.2f}" if pd.notnull(x) else "NA")
        elif col == 'F_Stat':
            # F_Stat formatting: display 2 decimal places for better readability
            final_df[col] = final_df[col].apply(lambda x: f"{x:.2f}" if pd.notnull(x) else "NA")
        else:
            # Beta and SE metrics use 4 decimal places for precision
            final_df[col] = final_df[col].apply(lambda x: f"{x:.4f}" if pd.notnull(x) else "NA")
        
    # Save to a temporary file first
    temp_file = OUTPUT_FILE.parent / f"{OUTPUT_FILE.stem}_temp.csv"
    final_df.to_csv(temp_file, index=False)
    
    # Then rename to the original name
    import os
    if os.path.exists(OUTPUT_FILE):
        os.remove(OUTPUT_FILE)
    os.rename(temp_file, OUTPUT_FILE)
    
    print(f"Table S1 saved to {OUTPUT_FILE}", flush=True)

if __name__ == "__main__":
    generate_table_s1()
