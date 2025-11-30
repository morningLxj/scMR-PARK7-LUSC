import pandas as pd
import numpy as np
from pathlib import Path
try:
    from scipy.stats import norm
    invcdf = lambda p: norm.ppf(p)
except Exception:
    from statistics import NormalDist
    invcdf = lambda p: NormalDist().inv_cdf(p)

EQTL = Path("My_MR_Project/Exposure/cd8et_eqtl_table.tsv.gz")
LIST = Path("My_MR_Project/top_genes_fast.txt")
OUTDIR = Path("My_MR_Project/Exposure")

def se_from_p_beta(p, b):
    p = np.where(p == 0, np.finfo(float).tiny, p)
    z = invcdf(1 - p / 2)
    return np.where(z != 0, np.abs(b / z), np.nan)

genes = [x.strip() for x in LIST.read_text(encoding='utf-8').splitlines() if x.strip()]
cols = ['GENE','RSID','A1','A2','A2_FREQ_ONEK1K','SPEARMANS_RHO','P_VALUE']
df = pd.read_csv(EQTL, sep='\t', compression='gzip', usecols=cols)
for gene in genes:
    sub = df[(df['GENE']==gene) & (df['P_VALUE']<1e-6)].copy()
    if sub.empty:
        continue
    sub['se'] = se_from_p_beta(sub['P_VALUE'].values, sub['SPEARMANS_RHO'].values)
    sub = sub.dropna(subset=['se'])
    out = OUTDIR / f"EXPOSURE_cd8et_{gene}_mr.tsv.gz"
    m = pd.DataFrame({
        'SNP': sub['RSID'].astype(str),
        'effect_allele': sub['A2'],
        'other_allele': sub['A1'],
        'eaf': sub['A2_FREQ_ONEK1K'],
        'beta': sub['SPEARMANS_RHO'],
        'se': sub['se'],
        'pval': sub['P_VALUE'],
        'exposure': [f'cd8et eQTL ({gene})'] * len(sub),
        'id.exposure': [f'cd8et_{gene}'] * len(sub)
    })
    if len(m) == 0:
        continue
    m.to_csv(out, sep='\t', index=False, compression='gzip')
print('done exposures for', len(genes))

