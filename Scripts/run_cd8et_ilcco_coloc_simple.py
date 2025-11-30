import pandas as pd
import numpy as np
from pathlib import Path
from statistics import NormalDist

MR = Path("My_MR_Project/CD8ET_TopN_MR_Results.csv")
OUTCOME = Path("My_MR_Project/Outcome/28604730-GCST004748-EFO_0001071.h.tsv")
EXP_DIR = Path("My_MR_Project/Exposure")
OUTCSV = Path("My_MR_Project/CD8ET_ILCCO_coloc_sensitivity.csv")

def main():
    mr = pd.read_csv(MR)
    mr['p_ivw'] = pd.to_numeric(mr['p_ivw'], errors='coerce')
    mr = mr[mr['p_ivw'] < 5e-8]
    targets = [x.replace('cd8et_','') for x in mr['gene'].astype(str).tolist()]
    out = pd.read_csv(OUTCOME, sep='\t', low_memory=False)
    out = out.rename(columns={'hm_rsid':'SNP','hm_chrom':'CHR','hm_pos':'POS','hm_beta':'by','standard_error':'se_by','hm_effect_allele':'ea_outcome','hm_other_allele':'nea_outcome'})
    out['SNP'] = out['SNP'].astype(str)
    out['CHR'] = out['CHR'].astype(str)
    out['POS'] = pd.to_numeric(out['POS'], errors='coerce')
    out = out.dropna(subset=['POS'])
    out['POS'] = out['POS'].astype(int)

    rows = []
    for g in targets:
        f = EXP_DIR / f"EXPOSURE_cd8et_{g}_mr.tsv.gz"
        if not f.exists():
            continue
        d = pd.read_csv(f, sep='\t', compression='gzip')
        d['SNP'] = d['SNP'].astype(str)
        m = pd.merge(d, out[['SNP','ea_outcome','nea_outcome','by','se_by','CHR','POS']], on='SNP', how='inner')
        if len(m) < 2:
            continue
        idx = m['effect_allele'] != m['ea_outcome']
        if idx.any():
            m.loc[idx, 'by'] = -m.loc[idx, 'by']
        bx = m['beta'].astype(float).values
        by = m['by'].astype(float).values
        se_bx = m['se'].astype(float).values
        se_by = m['se_by'].astype(float).values
        ratio = by / bx
        var_ratio = (se_by**2 / (bx**2)) + ((by**2) * (se_bx**2) / (bx**4))
        w = 1 / var_ratio
        ivw = float(np.sum(w * ratio) / np.sum(w))
        Q = float(np.sum(w * (ratio - ivw)**2))
        X = np.vstack([np.ones_like(bx), bx]).T
        W = np.diag(1.0/(se_by**2))
        XtW = X.T @ W
        beta_hat = np.linalg.pinv(XtW @ X) @ XtW @ by
        eg_int = float(beta_hat[0])
        resid = by - X @ beta_hat
        sigma2 = float((resid.T @ W @ resid) / (len(bx) - 2))
        cov = np.linalg.pinv(XtW @ X) * sigma2
        eg_int_se = float(np.sqrt(cov[0,0]))
        eg_int_p = float(2 * (1 - NormalDist().cdf(abs(eg_int / eg_int_se))))
        conc = float(np.mean(bx * by > 0))

        z_eqtl_ins = bx / se_bx
        z_ilcco_ins = by / se_by
        z_corr = float(pd.Series(z_eqtl_ins).corr(pd.Series(z_ilcco_ins)))
        pe_idx = int(np.argmax(np.abs(z_eqtl_ins)))
        po_idx = int(np.argmax(np.abs(z_ilcco_ins)))
        peak_shared = int(str(m.iloc[pe_idx]['SNP']) == str(m.iloc[po_idx]['SNP']))
        try:
            peak_dist = int(abs(int(m.iloc[pe_idx]['POS']) - int(m.iloc[po_idx]['POS'])))
        except Exception:
            peak_dist = -1
        chr_mode = str(m['CHR'].astype(str).mode().iloc[0])
        center = int(np.median(m['POS'].astype(int).values))

        rows.append({'gene': g, 'chr': chr_mode, 'center_pos': center, 'n_eqtl_region': int(len(m)), 'n_ilcco_region': int(len(m)), 'n_shared_rsids': int(len(m)), 'z_corr': z_corr, 'peak_shared': peak_shared, 'peak_dist_bp': peak_dist, 'Q': Q, 'Q_p': np.nan, 'egger_intercept': eg_int, 'egger_intercept_se': eg_int_se, 'egger_intercept_p': eg_int_p, 'concordance': conc})

    if rows:
        pd.DataFrame(rows).to_csv(OUTCSV, index=False)
        print('saved:', OUTCSV)
    else:
        cols = ['gene','chr','center_pos','n_eqtl_region','n_ilcco_region','n_shared_rsids','z_corr','peak_shared','peak_dist_bp','Q','Q_p','egger_intercept','egger_intercept_se','egger_intercept_p','concordance']
        pd.DataFrame(columns=cols).to_csv(OUTCSV, index=False)
        print('saved empty:', OUTCSV)

if __name__ == '__main__':
    main()
