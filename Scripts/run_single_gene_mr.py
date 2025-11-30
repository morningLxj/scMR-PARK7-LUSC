import pandas as pd
import numpy as np
from pathlib import Path
from statistics import NormalDist

EQTL = Path("My_MR_Project/Exposure/cd8et_eqtl_table.tsv.gz")
OUTCOME = Path("My_MR_Project/Outcome/28604730-GCST004748-EFO_0001071.h.tsv")
GENE = "KCNAB2"
OUTTXT = Path("My_MR_Project/mr_results_preview.txt")

def main():
    df = pd.read_csv(EQTL, sep='\t', compression='gzip', usecols=['GENE','RSID','A1','A2','A2_FREQ_ONEK1K','SPEARMANS_RHO','P_VALUE'], nrows=1000000)
    sub = df[(df['GENE']==GENE) & (df['P_VALUE']<1e-6)].copy()
    if sub.empty:
        OUTTXT.write_text("no_iv_for_gene\n", encoding='utf-8')
        print('no_iv_for_gene')
        return
    out = pd.read_csv(OUTCOME, sep='\t')
    out = out.rename(columns={'hm_rsid':'SNP','hm_effect_allele':'ea_outcome','hm_other_allele':'nea_outcome','hm_beta':'by','standard_error':'se_by'})
    sub = sub.rename(columns={'RSID':'SNP','A2':'effect_allele','A1':'other_allele','A2_FREQ_ONEK1K':'eaf','SPEARMANS_RHO':'beta','P_VALUE':'pval'})
    sub['SNP'] = sub['SNP'].astype(str)
    m = pd.merge(sub, out[['SNP','ea_outcome','nea_outcome','by','se_by']], on='SNP', how='inner')
    if m.empty:
        OUTTXT.write_text("no_merge_rows\n", encoding='utf-8')
        print('no_merge_rows')
        return
    idx = m['effect_allele'] != m['ea_outcome']
    if idx.any():
        m.loc[idx, 'by'] = -m.loc[idx, 'by']
    m['se'] = (m['beta'].abs() / NormalDist().inv_cdf(1 - m['pval']/2)).replace([np.inf,-np.inf], np.nan)
    m = m.dropna(subset=['se','beta','by','se_by'])
    m = m[(m['beta'] != 0) & (m['se'] > 0)]
    if len(m) < 1:
        OUTTXT.write_text("insufficient_clean_iv\n", encoding='utf-8')
        print('insufficient_clean_iv')
        return
    bx = m['beta'].values
    by = m['by'].values
    se_bx = m['se'].values
    se_by = m['se_by'].values
    ratio = by / bx
    var_ratio = (se_by**2 / (bx**2)) + ((by**2) * (se_bx**2) / (bx**4))
    w = 1 / var_ratio
    ivw = float(np.sum(w * ratio) / np.sum(w))
    se_ivw = float(np.sqrt(1 / np.sum(w)))
    p_ivw = float(2 * (1 - NormalDist().cdf(abs(ivw / se_ivw))))
    OUTTXT.write_text(f"gene={GENE}\nN={len(m)}\nbeta={ivw}\nse={se_ivw}\np={p_ivw}\n", encoding='utf-8')
    print('written:', OUTTXT)

if __name__ == '__main__':
    main()

