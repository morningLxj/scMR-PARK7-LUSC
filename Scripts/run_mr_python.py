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
P_THRESH = 1e-6
TOP_N = 10
LOG = Path("My_MR_Project/mr_py_log.txt")

def calculate_se(p_value, beta):
    p = np.where(p_value == 0, np.finfo(float).tiny, p_value)
    z = invcdf(1 - p / 2)
    se = np.where(z != 0, np.abs(beta / z), np.nan)
    return se

def main():
    counts = {}
    cols = ['GENE','RSID','CHR','POS','A1','A2','A2_FREQ_ONEK1K','SPEARMANS_RHO','P_VALUE']
    for chunk in pd.read_csv(EQTL, sep='\t', compression='gzip', usecols=['GENE','P_VALUE'], chunksize=500000):
        m = chunk['P_VALUE'] < P_THRESH
        sub = chunk.loc[m, ['GENE']]
        vc = sub['GENE'].value_counts()
        for g, c in vc.items():
            counts[g] = counts.get(g, 0) + c
    targets = [g for g, _ in sorted(counts.items(), key=lambda x: x[1], reverse=True)[:TOP_N]]
    LOG.write_text(f"targets_count={len(targets)}\n", encoding='utf-8')
    buffers = {g: [] for g in targets}
    for chunk in pd.read_csv(EQTL, sep='\t', compression='gzip', usecols=cols, chunksize=500000):
        chunk = chunk[chunk['GENE'].isin(targets)]
        m = chunk['P_VALUE'] < P_THRESH
        sub = chunk.loc[m]
        if sub.empty:
            continue
        for g in targets:
            gdf = sub[sub['GENE'] == g]
            if len(gdf) == 0:
                continue
            buffers[g].append(gdf)
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
    # prepare position columns for fallback
    if 'hm_chrom' in out.columns and 'hm_pos' in out.columns:
        out['hm_chrom'] = out['hm_chrom'].astype(str)
        out['hm_pos'] = out['hm_pos'].astype(int)
    if 'chromosome' in out.columns and 'base_pair_location' in out.columns:
        out['chromosome'] = out['chromosome'].astype(str)
        out['base_pair_location'] = out['base_pair_location'].astype(int)
    results = []
    for gene in targets:
        if len(buffers[gene]) == 0:
            continue
        df = pd.concat(buffers[gene], ignore_index=True)
        df['se'] = calculate_se(df['P_VALUE'].values, df['SPEARMANS_RHO'].values)
        df = df.dropna(subset=['se'])
        df = df.rename(columns={'RSID':'SNP','CHR':'CHR','POS':'POS','A2':'effect_allele','A1':'other_allele','A2_FREQ_ONEK1K':'eaf','SPEARMANS_RHO':'beta','P_VALUE':'pval'})
        # primary merge by rsid
        m = pd.merge(df, out, on='SNP', how='inner')
        # fallback merge by position
        if m.empty and ('hm_chrom' in out.columns and 'hm_pos' in out.columns):
            df['CHR'] = df['CHR'].astype(str)
            df['POS'] = df['POS'].astype(int)
            m = pd.merge(df, out, left_on=['CHR','POS'], right_on=['hm_chrom','hm_pos'], how='inner')
        if m.empty and ('chromosome' in out.columns and 'base_pair_location' in out.columns):
            df['CHR'] = df['CHR'].astype(str)
            df['POS'] = df['POS'].astype(int)
            m = pd.merge(df, out, left_on=['CHR','POS'], right_on=['chromosome','base_pair_location'], how='inner')
        try:
            prev = LOG.read_text(encoding='utf-8')
        except Exception:
            prev = ""
        LOG.write_text(prev + f"merged_rows_{gene}={len(m)}\n", encoding='utf-8')
        if m.empty:
            continue
        idx = m['effect_allele'] != m['ea_outcome']
        if idx.any():
            m.loc[idx, 'by'] = -m.loc[idx, 'by']
            tmp = m.loc[idx, 'ea_outcome'].copy()
            m.loc[idx, 'ea_outcome'] = m.loc[idx, 'nea_outcome']
            m.loc[idx, 'nea_outcome'] = tmp
        bx = m['beta'].values
        by = m['by'].values
        se_bx = m['se'].values
        se_by = m['se_by'].values
        ratio = by / bx
        var_ratio = (se_by**2 / (bx**2)) + ((by**2) * (se_bx**2) / (bx**4))
        w = 1 / var_ratio
        ivw = np.sum(w * ratio) / np.sum(w)
        se_ivw = np.sqrt(1 / np.sum(w))
        z = np.abs(ivw / se_ivw)
        try:
            from math import erf
            p_ivw = 2 * (1 - 0.5 * (1 + erf(z / np.sqrt(2))))
        except Exception:
            # NormalDist fallback
            try:
                from statistics import NormalDist
                p_ivw = 2 * (1 - NormalDist().cdf(z))
            except Exception:
                p_ivw = float('nan')
        results.append({'gene': f'cd8et_{gene}', 'n_iv': int(len(m)), 'beta_ivw': float(ivw), 'se_ivw': float(se_ivw), 'p_ivw': float(p_ivw)})
    if results:
        pd.DataFrame(results).to_csv(OUTCSV, index=False)
        print("saved:", OUTCSV)
    else:
        print("no MR results produced")

if __name__ == '__main__':
    main()
