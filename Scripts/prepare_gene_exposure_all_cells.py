import sys
import pandas as pd
import numpy as np
from pathlib import Path
try:
    from scipy.stats import norm
    def invcdf(p):
        return norm.ppf(p)
except Exception:
    from statistics import NormalDist
    def invcdf(p):
        return NormalDist().inv_cdf(p)

def calculate_se(p_value, beta):
    p = np.where(p_value == 0, np.finfo(float).tiny, p_value)
    z = invcdf(1 - p / 2)
    se = np.where(z != 0, np.abs(beta / z), np.nan)
    return se

def prepare_for_gene(gene_symbol, gene_id, p_thresh):
    base = Path("My_MR_Project/Exposure")
    files = sorted(base.glob("*_eqtl_table.tsv.gz"))
    total_iv = 0
    for f in files:
        try:
            df = pd.read_csv(f, sep='\t', compression='gzip')
        except Exception:
            continue
        has_symbol = gene_symbol in df.get('GENE', pd.Series([], dtype=str)).values
        has_id = gene_id in df.get('GENE_ID', pd.Series([], dtype=str)).values
        if not has_symbol and not has_id:
            continue
        m = (df['P_VALUE'] < p_thresh) & ((df['GENE'] == gene_symbol) | (df['GENE_ID'] == gene_id))
        sub = df.loc[m].copy()
        if sub.empty:
            continue
        sub['se'] = calculate_se(sub['P_VALUE'], sub['SPEARMANS_RHO'])
        cell = sub['CELL_ID'].iloc[0] if 'CELL_ID' in sub.columns else f.stem.split('_')[0]
        out = base / f"EXPOSURE_{cell}_{gene_symbol}_mr.tsv.gz"
        out_df = pd.DataFrame({
            'SNP': sub['RSID'],
            'effect_allele': sub['A2'],
            'other_allele': sub['A1'],
            'eaf': sub['A2_FREQ_ONEK1K'],
            'beta': sub['SPEARMANS_RHO'],
            'se': sub['se'],
            'pval': sub['P_VALUE'],
            'exposure': f'{cell} eQTL ({gene_symbol})',
            'id.exposure': f'{cell}_{gene_symbol}'
        })
        out_df.dropna(subset=['se'], inplace=True)
        if len(out_df) == 0:
            continue
        out_df.to_csv(out, sep='\t', index=False, compression='gzip')
        print(f"saved {len(out_df)} IVs to {out}")
        total_iv += len(out_df)
    if total_iv == 0:
        print("no instruments found")

def main():
    if len(sys.argv) < 4:
        print("usage: python prepare_gene_exposure_all_cells.py <GENE_SYMBOL> <GENE_ID> <P_THRESH>")
        return
    gene_symbol = sys.argv[1]
    gene_id = sys.argv[2]
    p_thresh = float(sys.argv[3])
    prepare_for_gene(gene_symbol, gene_id, p_thresh)

if __name__ == "__main__":
    main()

