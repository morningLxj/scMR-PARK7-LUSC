import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy import stats
from statistics import NormalDist
import matplotlib.gridspec as gridspec

# Settings
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.family'] = 'Arial'

# Paths
BASE_DIR = Path("My_MR_Project")
EXP_DIR = BASE_DIR / "Exposure"
OUTCOME_FILE = BASE_DIR / "Outcome/28604730-GCST004748-EFO_0001071.h.tsv"
OUTPUT_DIR = BASE_DIR / "FigureS2"
OUTPUT_DIR.mkdir(exist_ok=True)

# Genes and Cells of Interest (Core Genes)
# PARK7, CTSW, TMEM50A
TARGETS = {
    "PARK7": "cd8et",
    "CTSW": "cd8et",
    "TMEM50A": "cd8et"
}

def load_data(gene, cell):
    # Exposure
    exp_path = EXP_DIR / f"EXPOSURE_{cell}_{gene}_mr.tsv.gz"
    if not exp_path.exists():
        print(f"Exposure file not found: {exp_path}")
        return None
    
    df_exp = pd.read_csv(exp_path, sep='\t', compression='gzip')
    df_exp['SNP'] = df_exp['SNP'].astype(str)
    
    target_snps = set(df_exp['SNP'].unique())
    
    # Outcome - Read in chunks to save memory
    if not OUTCOME_FILE.exists():
        print(f"Outcome file not found: {OUTCOME_FILE}")
        return None
        
    outcome_chunks = []
    chunk_size = 100000 # Read 100k rows at a time
    
    try:
        # Check columns first to ensure we have the right names
        # We need to rename columns as well.
        # Standard names: hm_rsid -> SNP
        
        print(f"Reading outcome file in chunks...", flush=True)
        logging.info("Reading outcome file in chunks...")
        reader = pd.read_csv(OUTCOME_FILE, sep='\t', chunksize=10000)
        
        for i, chunk in enumerate(reader):
            if i % 100 == 0:
                msg = f"Processing chunk {i}..."
                print(msg, flush=True)
                logging.info(msg)
            
            # Rename first to match logic
            chunk = chunk.rename(columns={
                'hm_rsid': 'SNP',
                'hm_effect_allele': 'ea_outcome',
                'hm_other_allele': 'nea_outcome',
                'hm_beta': 'by',
                'standard_error': 'se_by',
                'p_value': 'p_outcome'
            })
            chunk['SNP'] = chunk['SNP'].astype(str)
            
            # Filter
            filtered_chunk = chunk[chunk['SNP'].isin(target_snps)]
            if not filtered_chunk.empty:
                outcome_chunks.append(filtered_chunk)
                
            # Optional: break early if we found all SNPs? 
            # No, SNPs might be duplicated or scattered.
            
        if outcome_chunks:
            df_out = pd.concat(outcome_chunks, ignore_index=True)
        else:
            print(f"No matching SNPs found in outcome for {gene}")
            return None
            
    except Exception as e:
        print(f"Error reading outcome file in chunks: {e}")
        return None
    
    # Merge
    df = pd.merge(df_exp, df_out, on='SNP', how='inner')
    
    # Harmonize
    idx = df['effect_allele'] != df['ea_outcome']
    if idx.any():
        df.loc[idx, 'by'] = -df.loc[idx, 'by']
    
    # Ensure numeric
    for col in ['beta', 'se', 'by', 'se_by']:
        df[col] = pd.to_numeric(df[col], errors='coerce')
        
    df = df.dropna(subset=['beta', 'se', 'by', 'se_by'])
    
    return df

def mr_ivw(bx, by, se_by):
    # Simple IVW for regression line
    weights = 1 / (se_by ** 2)
    # Fixed effect IVW (slope) - simplified for plotting line
    # In regression terms: by ~ 0 + bx, weights=weights
    # Weighted least squares
    
    # However, standard IVW uses specific formula.
    # For plotting the line: slope = sum(w * ratio) / sum(w)
    ratio = by / bx
    ivw_beta = np.sum(weights * ratio) / np.sum(weights)
    return ivw_beta

def mr_egger_stats(bx, by, se_by):
    # Egger regression: by ~ intercept + bx
    weights = 1 / (se_by ** 2)
    
    sum_w = np.sum(weights)
    sum_wx = np.sum(weights * bx)
    sum_wy = np.sum(weights * by)
    sum_wxx = np.sum(weights * bx * bx)
    sum_wxy = np.sum(weights * bx * by)
    
    denom = sum_w * sum_wxx - sum_wx * sum_wx
    
    intercept = (sum_wxx * sum_wy - sum_wx * sum_wxy) / denom
    slope = (sum_w * sum_wxy - sum_wx * sum_wy) / denom
    
    # Calculate P-value for intercept
    # Residuals
    y_pred = intercept + slope * bx
    resid = by - y_pred
    # Weighted sum of squared residuals / (n-2)
    rss = np.sum(weights * resid**2)
    n = len(by)
    if n > 2:
        sigma_sq = rss / (n - 2)
        # Variance of intercept
        # Var(beta_0) = sigma^2 * sum(w * x^2) / denom
        var_intercept = sigma_sq * sum_wxx / denom
        se_intercept = np.sqrt(var_intercept)
        t_stat = intercept / se_intercept
        p_val = 2 * (1 - stats.t.cdf(np.abs(t_stat), df=n-2))
    else:
        p_val = 1.0
        
    return intercept, slope, p_val

def plot_scatter(ax, df, gene):
    bx = df['beta'].values
    by = df['by'].values
    se_by = df['se_by'].values
    
    # Points with error bars
    ax.errorbar(bx, by, yerr=se_by, fmt='o', color='black', ecolor='lightgray', elinewidth=1, capsize=0, markersize=4, alpha=0.6)
    
    # IVW Line
    slope_ivw = mr_ivw(bx, by, se_by)
    x_range = np.linspace(min(bx), max(bx), 100)
    ax.plot(x_range, slope_ivw * x_range, color='#1f77b4', linestyle='-', linewidth=2, label='IVW')
    
    # Egger Line
    intercept_egg, slope_egg, p_intercept = mr_egger_stats(bx, by, se_by)
    ax.plot(x_range, intercept_egg + slope_egg * x_range, color='#d62728', linestyle='--', linewidth=2, label='MR-Egger')
    
    ax.set_xlabel(f'SNP effect on {gene}')
    ax.set_ylabel('SNP effect on Outcome')
    ax.set_title(f'Scatter Plot\nEgger Intercept P={p_intercept:.3f}')
    ax.legend()
    ax.grid(True, linestyle=':', alpha=0.6)

def plot_funnel(ax, df):
    # Funnel plot: IV estimate vs 1/SE_IV
    bx = df['beta'].values
    by = df['by'].values
    
    ratio = by / bx
    # Approx SE of ratio using Delta method or simple approximation se_by/|bx|
    # Standard funnel uses 1/SE as Y axis, Ratio as X axis
    # se_ratio approx = se_by / abs(bx)
    se_ratio = df['se_by'].values / np.abs(bx)
    precision = 1 / se_ratio
    
    # Global IVW estimate for center line
    ivw_est = mr_ivw(bx, by, df['se_by'].values)
    
    ax.scatter(ratio, precision, color='#2ca02c', alpha=0.6, s=30)
    ax.axvline(x=ivw_est, color='black', linestyle='--', linewidth=1.5)
    
    ax.set_xlabel('MR Estimate (beta)')
    ax.set_ylabel('Precision (1/SE)')
    ax.set_title('Funnel Plot')
    ax.grid(True, linestyle=':', alpha=0.6)

def plot_leave_one_out(ax, df):
    # Leave-one-out sensitivity
    # Re-run IVW excluding one SNP at a time
    indices = df.index.tolist()
    estimates = []
    labels = []
    ses = []
    
    bx_all = df['beta'].values
    by_all = df['by'].values
    se_by_all = df['se_by'].values
    
    # Full estimate
    full_est = mr_ivw(bx_all, by_all, se_by_all)
    
    for i in range(len(df)):
        # Exclude i
        mask = np.ones(len(df), dtype=bool)
        mask[i] = False
        
        bx = bx_all[mask]
        by = by_all[mask]
        se_by = se_by_all[mask]
        
        est = mr_ivw(bx, by, se_by)
        estimates.append(est)
        labels.append(df['SNP'].iloc[i])
        
        # Crude SE for plotting error bars (not strictly calculated here, just point estimate for LOO usually sufficient or use standard error propagation)
        # For visualization, let's just plot the point estimates. 
        # If we want error bars, we need to calc SE for each. 
        # Simplified: standard IVW SE formula
        # se_ivw = sqrt(1 / sum(weights))
        weights = 1 / (se_by**2)
        se_est = np.sqrt(1 / np.sum(weights))
        ses.append(se_est)

    # Plot
    y_pos = np.arange(len(df))
    ax.errorbar(estimates, y_pos, xerr=ses, fmt='o', color='black', ecolor='gray', elinewidth=1, capsize=2, markersize=4)
    ax.axvline(x=full_est, color='red', linestyle='--', linewidth=1)
    ax.axvline(x=0, color='gray', linestyle=':', linewidth=1)
    
    ax.set_yticks(y_pos)
    # If too many SNPs, don't label all y-axis
    if len(df) > 20:
        ax.set_yticklabels([''] * len(df))
        ax.set_ylabel(f'{len(df)} SNPs (Leave-one-out)')
    else:
        ax.set_yticklabels(labels, fontsize=8)
        
    ax.set_xlabel('MR Estimate')
    ax.set_title('Leave-one-out Sensitivity')
    ax.invert_yaxis() # Top is first SNP
    ax.grid(True, axis='x', linestyle=':', alpha=0.6)

def plot_forest(ax, df):
    # Forest Plot of individual SNP effects
    
    # Calculate individual Wald ratios
    bx = df['beta'].values
    by = df['by'].values
    se_by = df['se_by'].values
    se_bx = df['se'].values
    
    ratios = by / bx
    # SE of ratio (Delta method approx)
    se_ratios = np.sqrt((se_by**2/bx**2) + (by**2 * se_bx**2 / bx**4))
    
    y_pos = np.arange(len(df))
    
    ax.errorbar(ratios, y_pos, xerr=1.96*se_ratios, fmt='o', color='#ff7f0e', ecolor='gray', elinewidth=1, capsize=2, markersize=4)
    
    # Add IVW estimate line
    ivw_est = mr_ivw(bx, by, se_by)
    ax.axvline(x=ivw_est, color='red', linestyle='--', label='IVW')
    ax.axvline(x=0, color='gray', linestyle=':')
    
    ax.set_xlabel('Wald Ratio (beta)')
    ax.set_title('Single SNP Analysis (Forest Plot)')
    ax.set_yticks(y_pos)
    if len(df) > 20:
        ax.set_yticklabels([''] * len(df))
        ax.set_ylabel(f'{len(df)} SNPs')
    else:
        ax.set_yticklabels(df['SNP'], fontsize=8)
    ax.invert_yaxis()
    ax.grid(True, axis='x', linestyle=':', alpha=0.6)


import logging

# Setup logging
logging.basicConfig(filename='d:/2026YJ/My_MR_Project/debug_s2.log', level=logging.DEBUG, format='%(asctime)s %(message)s')

def load_all_data_and_process():
    # 1. Load all exposure data first
    exposures = {}
    all_snps = set()
    
    for gene, cell in TARGETS.items():
        print(f"Loading exposure for {gene} in {cell}...", flush=True)
        exp_path = EXP_DIR / f"EXPOSURE_{cell}_{gene}_mr.tsv.gz"
        if not exp_path.exists():
            print(f"Exposure file not found: {exp_path}", flush=True)
            continue
            
        try:
            df_exp = pd.read_csv(exp_path, sep='\t', compression='gzip')
            df_exp['SNP'] = df_exp['SNP'].astype(str)
            exposures[gene] = df_exp
            all_snps.update(df_exp['SNP'].unique())
        except Exception as e:
            print(f"Error loading exposure for {gene}: {e}", flush=True)
            
    if not exposures:
        print("No exposure data loaded.", flush=True)
        return

    print(f"Total unique SNPs to fetch: {len(all_snps)}", flush=True)

    # 2. Read Outcome file ONCE
    if not OUTCOME_FILE.exists():
        print(f"Outcome file not found: {OUTCOME_FILE}", flush=True)
        return

    outcome_rows = []
    print(f"Reading outcome file...", flush=True)
    # Columns to keep and rename
    # We only want specific columns to avoid merge conflicts (like effect_allele)
    
    try:
        reader = pd.read_csv(OUTCOME_FILE, sep='\t', chunksize=50000)
        
        for i, chunk in enumerate(reader):
            if i % 50 == 0:
                print(f"Scanning chunk {i}...", flush=True)
            
            # Rename columns we need
            # Note: standard_error might be hm_standard_error or standard_error. 
            # Based on file inspection: standard_error exists.
            chunk = chunk.rename(columns={
                'hm_rsid': 'SNP',
                'hm_effect_allele': 'ea_outcome',
                'hm_other_allele': 'nea_outcome',
                'hm_beta': 'by',
                'standard_error': 'se_by',
                'p_value': 'p_outcome'
            })
            chunk['SNP'] = chunk['SNP'].astype(str)
            
            # Filter for SNPs we need
            mask = chunk['SNP'].isin(all_snps)
            if mask.any():
                # Keep only relevant columns to avoid merge suffixes
                cols_to_keep = ['SNP', 'ea_outcome', 'nea_outcome', 'by', 'se_by', 'p_outcome']
                # Check if all cols exist
                available_cols = [c for c in cols_to_keep if c in chunk.columns]
                outcome_rows.append(chunk.loc[mask, available_cols])
                
        if outcome_rows:
            df_outcome_all = pd.concat(outcome_rows, ignore_index=True)
            print(f"Loaded {len(df_outcome_all)} matching outcome rows.", flush=True)
        else:
            print("No matching outcome rows found.", flush=True)
            return

    except Exception as e:
        print(f"Error reading outcome file: {e}", flush=True)
        return

    # 3. Process each gene
    suffix_map = {
        "PARK7": "S2A",
        "CTSW": "S2B",
        "TMEM50A": "S2C"
    }

    for gene, df_exp in exposures.items():
        if gene not in suffix_map:
            continue
            
        print(f"Processing {gene}...", flush=True)
        
        # Merge
        df = pd.merge(df_exp, df_outcome_all, on='SNP', how='inner')
        
        if df.empty:
            print(f"No overlapping SNPs for {gene}", flush=True)
            continue
            
        # Harmonize
        # Ensure effect_allele exists (from exposure) and ea_outcome exists (from outcome)
        if 'effect_allele' not in df.columns:
            print(f"Error: 'effect_allele' column missing for {gene}. Columns: {df.columns}", flush=True)
            continue
            
        idx = df['effect_allele'] != df['ea_outcome']
        if idx.any():
            df.loc[idx, 'by'] = -df.loc[idx, 'by']
        
        # Ensure numeric
        for col in ['beta', 'se', 'by', 'se_by']:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors='coerce')
            
        df = df.dropna(subset=['beta', 'se', 'by', 'se_by'])
        
        if len(df) < 3:
            print(f"Not enough valid SNPs for {gene} (n={len(df)})", flush=True)
            continue
            
        # Plot
        create_panel_plot(gene, df, suffix_map[gene])

def create_panel_plot(gene, df, suffix):
    try:
        fig = plt.figure(figsize=(12, 10))
        gs = gridspec.GridSpec(2, 2, height_ratios=[1, 1], width_ratios=[1, 1])
        
        # A: Scatter Plot (includes Egger line)
        ax1 = fig.add_subplot(gs[0, 0])
        plot_scatter(ax1, df, gene)
        ax1.text(-0.1, 1.05, 'A', transform=ax1.transAxes, fontsize=16, fontweight='bold', va='top', ha='right')
        
        # B: Funnel Plot
        ax2 = fig.add_subplot(gs[0, 1])
        plot_funnel(ax2, df)
        ax2.text(-0.1, 1.05, 'B', transform=ax2.transAxes, fontsize=16, fontweight='bold', va='top', ha='right')
        
        # C: Leave-one-out
        ax3 = fig.add_subplot(gs[1, 0])
        plot_leave_one_out(ax3, df)
        ax3.text(-0.1, 1.05, 'C', transform=ax3.transAxes, fontsize=16, fontweight='bold', va='top', ha='right')
        
        # D: Forest Plot (Single SNP)
        ax4 = fig.add_subplot(gs[1, 1])
        plot_forest(ax4, df)
        ax4.text(-0.1, 1.05, 'D', transform=ax4.transAxes, fontsize=16, fontweight='bold', va='top', ha='right')
        
        plt.tight_layout()
        
        # Save
        out_png = OUTPUT_DIR / f"Figure{suffix}_{gene}_Sensitivity.png"
        out_pdf = OUTPUT_DIR / f"Figure{suffix}_{gene}_Sensitivity.pdf"
        
        plt.savefig(out_png, dpi=600, bbox_inches='tight')
        plt.savefig(out_pdf, dpi=600, bbox_inches='tight')
        print(f"Saved {gene} panel to {out_png}", flush=True)
        plt.close()
    except Exception as e:
        print(f"Error plotting {gene}: {e}", flush=True)
        import traceback
        traceback.print_exc()

def main():
    sns.set_style("whitegrid", {'axes.grid': True, 'grid.linestyle': '--'})
    load_all_data_and_process()

if __name__ == "__main__":
    main()
