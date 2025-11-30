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
OUTCOME_ILCCO = Path("My_MR_Project/Outcome/28604730-GCST004748-EFO_0001071.h.tsv")
OUTCOME_FG = Path("My_MR_Project/Outcome/finngen_R12_C3_BRONCHUS_LUNG_EXALLC.gz")
OUTCSV = Path("My_MR_Project/CD8ET_TopN_MR_Results.csv")
P_THRESH = 1e-6
TOP_N = 10
LOG = Path("My_MR_Project/mr_py_log.txt")

def calculate_se(p_value, beta):
    p = np.where(p_value == 0, np.finfo(float).tiny, p_value)
    z = invcdf(1 - p / 2)
    se = np.where(z != 0, np.abs(beta / z), np.nan)
    return se

def main():
    df = pd.read_csv(EQTL, sep='\t', compression='gzip', usecols=['GENE','RSID','CHR','POS','A1','A2','A2_FREQ_ONEK1K','SPEARMANS_RHO','P_VALUE'], nrows=1000000)
    sub = df[df['P_VALUE'] < P_THRESH]
    if sub.empty:
        print('no instruments at p<1e-6 in first 1e6 rows')
        return
    top = sub['GENE'].value_counts().head(TOP_N).index.tolist()
    dsub = sub[sub['GENE'].isin(top)].copy()
    dsub['se'] = calculate_se(dsub['P_VALUE'].values, dsub['SPEARMANS_RHO'].values)
    dsub = dsub.dropna(subset=['se'])
    out_ilcco = pd.read_csv(OUTCOME_ILCCO, sep='\t')
    out_ilcco = out_ilcco.rename(columns={
        'hm_rsid':'SNP',
        'hm_effect_allele':'ea_outcome',
        'hm_other_allele':'nea_outcome',
        'hm_beta':'by',
        'standard_error':'se_by',
        'hm_effect_allele_frequency':'eaf_outcome',
        'p_value':'p_outcome'
    })
    out_ilcco['SNP'] = out_ilcco['SNP'].astype(str)
    results = []
    for gene in top:
        gdf = dsub[dsub['GENE'] == gene].copy()
        gdf = gdf.rename(columns={'RSID':'SNP','CHR':'CHR','POS':'POS','A2':'effect_allele','A1':'other_allele','A2_FREQ_ONEK1K':'eaf','SPEARMANS_RHO':'beta','P_VALUE':'pval'})
        gdf['SNP'] = gdf['SNP'].astype(str)
        m = pd.merge(gdf, out_ilcco, on='SNP', how='inner')
        if m.empty:
            continue
        try:
            prev = LOG.read_text(encoding='utf-8')
        except Exception:
            prev = ""
        LOG.write_text(prev + f"merged_{gene}={len(m)}\n", encoding='utf-8')
        for col in ['beta','se','by','se_by']:
            m[col] = pd.to_numeric(m[col], errors='coerce')
        m = m.replace([np.inf, -np.inf], np.nan).dropna(subset=['beta','se','by','se_by'])
        m = m[(m['beta'] != 0) & (m['se'] > 0)]
        if len(m) < 1:
            continue
        idx = m['effect_allele'] != m['ea_outcome']
        if idx.any():
            m.loc[idx, 'by'] = -m.loc[idx, 'by']
        bx = m['beta'].values
        by = m['by'].values
        se_bx = m['se'].values
        se_by = m['se_by'].values
        ratio = by / bx
        var_ratio = (se_by**2 / (bx**2)) + ((by**2) * (se_bx**2) / (bx**4))
        w = 1 / var_ratio
        ivw = np.sum(w * ratio) / np.sum(w)
        se_ivw = np.sqrt(1 / np.sum(w))
        try:
            from statistics import NormalDist
            p_ivw = 2 * (1 - NormalDist().cdf(abs(ivw / se_ivw)))
        except Exception:
            p_ivw = float('nan')
        res = {'gene': f'cd8et_{gene}', 'n_iv': int(len(m)), 'beta_ivw': float(ivw), 'se_ivw': float(se_ivw), 'p_ivw': float(p_ivw)}
        print('RES', res)
        results.append(res)
    if results:
        pd.DataFrame(results).to_csv(OUTCSV, index=False)
        print('SAVED', OUTCSV)
    else:
        print('NO_RESULTS')

if __name__ == '__main__':
    main()
