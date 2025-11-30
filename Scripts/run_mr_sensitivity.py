import pandas as pd
import numpy as np
from pathlib import Path
from statistics import NormalDist

GENES = ["CTSW","EIF4G3","TMEM50A","PARK7"]
EXP_DIR = Path("My_MR_Project/Exposure")
OUTCOME = Path("My_MR_Project/Outcome/28604730-GCST004748-EFO_0001071.h.tsv")
FG_OUTCOME = Path("My_MR_Project/Outcome/finngen_R12_C3_BRONCHUS_LUNG_EXALLC.gz")
OUTCSV = Path("My_MR_Project/CD8ET_MR_Sensitivity.csv")

def egger_intercept(bx, by, se_by):
    X = np.vstack([np.ones_like(bx), bx]).T
    W = np.diag(1.0/(se_by**2))
    XtW = X.T @ W
    beta_hat = np.linalg.pinv(XtW @ X) @ XtW @ by
    resid = by - X @ beta_hat
    sigma2 = float((resid.T @ W @ resid) / (len(bx) - 2))
    cov = np.linalg.pinv(XtW @ X) * sigma2
    inter = float(beta_hat[0])
    se = float(np.sqrt(cov[0,0]))
    p = float(2 * (1 - NormalDist().cdf(abs(inter / se))))
    return inter, se, p

def presso_like_clean(bx, by, se_by, threshold=3.0, max_iter=5):
    X = np.vstack([np.ones_like(bx), bx]).T
    W = np.diag(1.0/(se_by**2))
    for _ in range(max_iter):
        XtW = X.T @ W
        beta_hat = np.linalg.pinv(XtW @ X) @ XtW @ by
        resid = by - X @ beta_hat
        # leverage h = diag(X*(X'WX)^-1*X'W)
        H = X @ np.linalg.pinv(XtW @ X) @ XtW
        h = np.clip(np.diag(H), 0, 0.99)
        sigma2 = float((resid.T @ W @ resid) / max(1, len(bx) - 2))
        stud = resid / np.sqrt(sigma2 * (1 - h))
        idx_out = np.where(np.abs(stud) > threshold)[0]
        if len(idx_out) == 0:
            break
        keep = np.ones(len(bx), dtype=bool)
        keep[idx_out] = False
        bx, by, se_by = bx[keep], by[keep], se_by[keep]
        X = np.vstack([np.ones_like(bx), bx]).T
        W = np.diag(1.0/(se_by**2))
    return bx, by, se_by

def main():
    out = pd.read_csv(OUTCOME, sep='\t')
    out = out.rename(columns={'hm_rsid':'SNP','hm_effect_allele':'ea_outcome','hm_other_allele':'nea_outcome','hm_beta':'by','standard_error':'se_by'})
    out['SNP'] = out['SNP'].astype(str)
    # FinnGen outcome mapping
    fg = pd.read_csv(FG_OUTCOME, sep='\t', compression='gzip')
    # rsids may contain ';' list, take first rsid
    if 'rsids' in fg.columns:
        fg['SNP'] = fg['rsids'].astype(str).str.split(';').str[0]
    else:
        fg['SNP'] = fg.get('variant', '').astype(str)
    fg = fg.rename(columns={'#chrom':'CHR','pos':'POS','beta':'by','sebeta':'se_by','alt':'ea_outcome','ref':'nea_outcome'})
    fg = fg[['SNP','CHR','POS','by','se_by','ea_outcome','nea_outcome']].dropna(subset=['SNP'])
    rows = []
    for g in GENES:
        f = EXP_DIR / f"EXPOSURE_cd8et_{g}_mr.tsv.gz"
        if not f.exists():
            continue
        d = pd.read_csv(f, sep='\t', compression='gzip')
        d['SNP'] = d['SNP'].astype(str)
        m = pd.merge(d, out[['SNP','ea_outcome','nea_outcome','by','se_by']], on='SNP', how='inner')
        if len(m) < 3:
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
        inter, se, p = egger_intercept(bx, by, se_by)
        conc = float(np.mean(bx * by > 0))
        # MR-PRESSO-like outlier removal for TMEM50A
        Q_clean = np.nan; inter_clean = np.nan; se_clean = np.nan; p_clean = np.nan; n_clean = np.nan
        if g == 'TMEM50A':
            bx_c, by_c, se_by_c = presso_like_clean(bx.copy(), by.copy(), se_by.copy())
            if len(bx_c) >= 3:
                ratio_c = by_c / bx_c
                var_ratio_c = (se_by_c**2 / (bx_c**2))
                w_c = 1 / (var_ratio_c + 1e-12)
                ivw_c = float(np.sum(w_c * ratio_c) / np.sum(w_c))
                Q_clean = float(np.sum(w_c * (ratio_c - ivw_c)**2))
                inter_clean, se_clean, p_clean = egger_intercept(bx_c, by_c, se_by_c)
                n_clean = int(len(bx_c))

        # FinnGen replicate
        mf = pd.merge(d, fg, on='SNP', how='inner')
        fgQ = np.nan; fgInter = np.nan; fgSe = np.nan; fgP = np.nan; fgConc = np.nan; fgN = int(len(mf))
        if len(mf) >= 3:
            idxf = mf['effect_allele'] != mf['ea_outcome']
            if idxf.any():
                mf.loc[idxf, 'by'] = -mf.loc[idxf, 'by']
            bx_f = mf['beta'].astype(float).values
            by_f = mf['by'].astype(float).values
            se_bx_f = mf['se'].astype(float).values
            se_by_f = mf['se_by'].astype(float).values
            ratio_f = by_f / bx_f
            var_ratio_f = (se_by_f**2 / (bx_f**2)) + ((by_f**2) * (se_bx_f**2) / (bx_f**4))
            w_f = 1 / var_ratio_f
            ivw_f = float(np.sum(w_f * ratio_f) / np.sum(w_f))
            fgQ = float(np.sum(w_f * (ratio_f - ivw_f)**2))
            fgInter, fgSe, fgP = egger_intercept(bx_f, by_f, se_by_f)
            fgConc = float(np.mean(bx_f * by_f > 0))

        rows.append({'gene': g, 'n_iv': int(len(m)), 'Q': Q, 'Egger_intercept': inter, 'Egger_se': se, 'Egger_p': p, 'concordance': conc,
                     'n_iv_clean': n_clean, 'Q_clean': Q_clean, 'Egger_intercept_clean': inter_clean, 'Egger_se_clean': se_clean, 'Egger_p_clean': p_clean,
                     'fg_n_iv': fgN, 'fg_Q': fgQ, 'fg_Egger_intercept': fgInter, 'fg_Egger_se': fgSe, 'fg_Egger_p': fgP, 'fg_concordance': fgConc})
    pd.DataFrame(rows).to_csv(OUTCSV, index=False)
    print('saved:', OUTCSV)

if __name__ == '__main__':
    main()
