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

INPUT_FILE = "My_MR_Project/Exposure/cd8et_eqtl_table.tsv.gz"
P_THRESHOLD = 1e-6
TOP_N = 10

def calculate_se(p_value, beta):
    p = np.where(p_value == 0, np.finfo(float).tiny, p_value)
    z = invcdf(1 - p / 2)
    se = np.where(z != 0, np.abs(beta / z), np.nan)
    return se

def main():
    base = Path("My_MR_Project/Exposure")
    counts = {}
    for chunk in pd.read_csv(INPUT_FILE, sep='\t', compression='gzip', chunksize=500000):
        m = chunk['P_VALUE'] < P_THRESHOLD
        sub = chunk.loc[m, ['GENE']]
        for g, c in sub['GENE'].value_counts().items():
            counts[g] = counts.get(g, 0) + c
    genes = sorted(counts.items(), key=lambda x: x[1], reverse=True)[:TOP_N]
    targets = [g for g, _ in genes]
    buffers = {g: [] for g in targets}
    for chunk in pd.read_csv(INPUT_FILE, sep='\t', compression='gzip', chunksize=500000):
        m = (chunk['P_VALUE'] < P_THRESHOLD) & (chunk['GENE'].isin(targets))
        sub = chunk.loc[m, ['GENE','RSID','A2','A1','A2_FREQ_ONEK1K','SPEARMANS_RHO','P_VALUE']]
        if sub.empty:
            continue
        for g in targets:
            gdf = sub[sub['GENE'] == g]
            if len(gdf) == 0:
                continue
            buffers[g].append(gdf)
    for gene in targets:
        if len(buffers[gene]) == 0:
            continue
        df = pd.concat(buffers[gene], ignore_index=True)
        df['se'] = calculate_se(df['P_VALUE'], df['SPEARMANS_RHO'])
        df_mr = pd.DataFrame({
            'SNP': df['RSID'],
            'effect_allele': df['A2'],
            'other_allele': df['A1'],
            'eaf': df['A2_FREQ_ONEK1K'],
            'beta': df['SPEARMANS_RHO'],
            'se': df['se'],
            'pval': df['P_VALUE'],
            'exposure': f'cd8et eQTL ({gene})',
            'id.exposure': f'cd8et_{gene}'
        })
        df_mr.dropna(subset=['se'], inplace=True)
        if len(df_mr) == 0:
            continue
        out = base / f"EXPOSURE_cd8et_{gene}_mr.tsv.gz"
        df_mr.to_csv(out, sep='\t', index=False, compression='gzip')
    print("done")

if __name__ == "__main__":
    main()
