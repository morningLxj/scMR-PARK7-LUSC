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
OUTCOME = Path("My_MR_Project/Outcome/28604730-GCST004748-EFO_0001071.h.tsv")
OUTCSV = Path("My_MR_Project/CD8ET_TopN_MR_Results.csv")
LOG = Path("My_MR_Project/mr_simple_log.txt")
P_THRESH = 1e-6
TOP_N = 10

def calculate_se(p, b):
    p = np.where(p == 0, np.finfo(float).tiny, p)
    z = invcdf(1 - p / 2)
    return np.where(z != 0, np.abs(b / z), np.nan)

def main():
    df = pd.read_csv(EQTL, sep='\t', compression='gzip', usecols=['GENE','RSID','A1','A2','A2_FREQ_ONEK1K','SPEARMANS_RHO','P_VALUE'], nrows=1000000)
    sig = df[df['P_VALUE'] < P_THRESH]
    if sig.empty:
        LOG.write_text("no_sig_eqtl_rows\n", encoding='utf-8')
        print('no_sig_eqtl_rows')
        return
    top = sig['GENE'].value_counts().head(TOP_N).index.tolist()
    LOG.write_text(f"top_genes={','.join(top)}\n", encoding='utf-8')
    out = pd.read_csv(OUTCOME, sep='\t')
    out = out.rename(columns={
        'hm_rsid':'SNP',
        'hm_effect_allele':'ea_outcome',
        'hm_other_allele':'nea_outcome',
        'hm_beta':'by',
        'standard_error':'se_by',
        'hm_effect_allele_frequency':'eaf_outcome',
        'p_value':'p_outcome'
    })
    out['SNP'] = out['SNP'].astype(str)
    results = []
    for gene in top:
        g = sig[sig['GENE'] == gene].copy()
        g['se'] = calculate_se(g['P_VALUE'].values, g['SPEARMANS_RHO'].values)
        g = g.dropna(subset=['se'])
        g = g.rename(columns={'RSID':'SNP','A2':'effect_allele','A1':'other_allele','A2_FREQ_ONEK1K':'eaf','SPEARMANS_RHO':'beta','P_VALUE':'pval'})
        g['SNP'] = g['SNP'].astype(str)
        m = pd.merge(g, out, on='SNP', how='inner')
        LOG.write_text(LOG.read_text(encoding='utf-8') + f"merged_rows_{gene}={len(m)}\n", encoding='utf-8')
        if m.empty:
            continue
        idx = m['effect_allele'] != m['ea_outcome']
        if idx.any():
            m.loc[idx, 'by'] = -m.loc[idx, 'by']
        for col in ['beta','se','by','se_by']:
            m[col] = pd.to_numeric(m[col], errors='coerce')
        m = m.replace([np.inf, -np.inf], np.nan).dropna(subset=['beta','se','by','se_by'])
        m = m[(m['beta'] != 0) & (m['se'] > 0)]
        if len(m) < 1:
            continue
        bx = m['beta'].values
        by = m['by'].values
        se_bx = m['se'].values
        se_by = m['se_by'].values
        ratio = by / bx
        var_ratio = (se_by**2 / (bx**2)) + ((by**2) * (se_bx**2) / (bx**4))
        w = 1 / var_ratio
        ivw = float(np.sum(w * ratio) / np.sum(w))
        se_ivw = float(np.sqrt(1 / np.sum(w)))
        try:
            from statistics import NormalDist
            p_ivw = float(2 * (1 - NormalDist().cdf(abs(ivw / se_ivw))))
        except Exception:
            p_ivw = float('nan')
        results.append({'gene': f'cd8et_{gene}', 'n_iv': int(len(m)), 'beta_ivw': ivw, 'se_ivw': se_ivw, 'p_ivw': p_ivw})
    if results:
        pd.DataFrame(results).to_csv(OUTCSV, index=False)
        print('saved:', OUTCSV)
    else:
        print('no MR results produced')

if __name__ == '__main__':
    main()

