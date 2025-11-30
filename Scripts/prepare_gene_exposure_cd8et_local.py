import sys
import pandas as pd
import numpy as np
from pathlib import Path
try:
    from scipy.stats import norm
    invcdf = lambda p: norm.ppf(p)
except Exception:
    from statistics import NormalDist
    invcdf = lambda p: NormalDist().inv_cdf(p)

INPUT_FILE = r"My_MR_Project/Exposure/cd8et_eqtl_table.tsv.gz"
OUT_DIR = Path("My_MR_Project/Exposure")

def calculate_se(p_value, beta):
    p = np.where(p_value == 0, np.finfo(float).tiny, p_value)
    z = invcdf(1 - p / 2)
    se = np.where(z != 0, np.abs(beta / z), np.nan)
    return se

def main():
    if len(sys.argv) < 3:
        print("usage: python prepare_gene_exposure_cd8et_local.py <GENE_SYMBOL> <P_THRESH>")
        return
    gene = sys.argv[1]
    p_thresh = float(sys.argv[2])
    df = pd.read_csv(INPUT_FILE, sep='\t', compression='gzip')
    m = (df['GENE'] == gene) & (df['P_VALUE'] < p_thresh)
    sub = df.loc[m].copy()
    if sub.empty:
        print("no instruments found")
        return
    sub['se'] = calculate_se(sub['P_VALUE'], sub['SPEARMANS_RHO'])
    out_df = pd.DataFrame({
        'SNP': sub['RSID'],
        'effect_allele': sub['A2'],
        'other_allele': sub['A1'],
        'eaf': sub['A2_FREQ_ONEK1K'],
        'beta': sub['SPEARMANS_RHO'],
        'se': sub['se'],
        'pval': sub['P_VALUE'],
        'exposure': f'cd8et eQTL ({gene})',
        'id.exposure': f'cd8et_{gene}'
    })
    out_path = OUT_DIR / f"EXPOSURE_cd8et_{gene}_mr.tsv.gz"
    out_df.dropna(subset=['se'], inplace=True)
    out_df.to_csv(out_path, sep='\t', index=False, compression='gzip')
    print(f"saved {len(out_df)} IVs to {out_path}")

if __name__ == "__main__":
    main()
