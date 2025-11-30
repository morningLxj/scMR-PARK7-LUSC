import pandas as pd
import numpy as np
from pathlib import Path
from statistics import NormalDist

MR = Path("My_MR_Project/CD8ET_TopN_MR_Results.csv")
EQTL = Path("My_MR_Project/Exposure/cd8et_eqtl_table.tsv.gz")
OUTCOME = Path("My_MR_Project/Outcome/28604730-GCST004748-EFO_0001071.h.tsv")
EXP_DIR = Path("My_MR_Project/Exposure")
OUTCSV = Path("CD8ET_ILCCO_coloc_sensitivity.csv")

def mode_chr(x):
    vals, counts = np.unique(x, return_counts=True)
    return str(vals[np.argmax(counts)])

def read_eqtl_genes(genes):
    cols = ['GENE','RSID','CHR','POS','SPEARMANS_RHO','P_VALUE']
    buf = []
    print(f"Reading eQTL file for {len(genes)} genes...", flush=True)
    try:
        for i, chunk in enumerate(pd.read_csv(EQTL, sep='\t', compression='gzip', usecols=cols, chunksize=500000)):
            if i % 1 == 0:
                print(f"Processing chunk {i}, buf size: {len(buf)}", flush=True)
            # if i > 20: break # Debug limit
            g = chunk[chunk['GENE'].isin(genes)]
            if len(g):
                buf.append(g)
    except Exception as e:
        print(f"Error reading eQTL file: {e}", flush=True)
        return None
        
    if not buf:
        return pd.DataFrame(columns=cols)
    return pd.concat(buf, ignore_index=True)

def main():
    mr = pd.read_csv(MR)
    mr['p_ivw'] = pd.to_numeric(mr['p_ivw'], errors='coerce')
    mr = mr[mr['p_ivw'] < 5e-8]
    targets = [x.replace('cd8et_','') for x in mr['gene'].astype(str).tolist()]
    print('targets_count', len(targets), flush=True)
    
    # Read all eQTL data for targets at once
    e_all = read_eqtl_genes(set(targets))
    if e_all is None:
        print("No eQTL data found for targets (None)", flush=True)
        return
    print(f"e_all rows: {len(e_all)}", flush=True)
    if len(e_all) == 0:
        print("No eQTL data found for targets (0 rows)", flush=True)
        # still write empty schema
        out_path = Path("My_MR_Project/CD8ET_ILCCO_coloc_sensitivity.csv")
        cols = ['gene','chr','center_pos','n_eqtl_region','n_ilcco_region','n_shared_rsids','z_corr','peak_shared','peak_dist_bp','Q','Q_p','egger_intercept','egger_intercept_se','egger_intercept_p','concordance']
        pd.DataFrame(columns=cols).to_csv(out_path, index=False)
        print('saved empty:', out_path)
        return
    
    out = pd.read_csv(OUTCOME, sep='\t', low_memory=False)
    out = out.rename(columns={'hm_rsid':'SNP','hm_chrom':'CHR','hm_pos':'POS','hm_beta':'by','standard_error':'se_by','hm_effect_allele':'ea_outcome','hm_other_allele':'nea_outcome'})
    out['CHR'] = out['CHR'].astype(str)
    out['POS'] = pd.to_numeric(out['POS'], errors='coerce')
    out = out.dropna(subset=['POS'])
    out['POS'] = out['POS'].astype(int)
    
    rows = []
    for g in targets:
        print(f"Processing gene {g}...", flush=True)
        eg = e_all[e_all['GENE'] == g].copy()
        if len(eg)==0:
            exp_file = EXP_DIR / f"EXPOSURE_cd8et_{g}_mr.tsv.gz"
            if not exp_file.exists():
                print(f"No data for {g}", flush=True)
                continue
            d = pd.read_csv(exp_file, sep='\t', compression='gzip')
            d['SNP'] = d['SNP'].astype(str)
            m = pd.merge(d, out[['SNP','ea_outcome','nea_outcome','by','se_by','CHR','POS']], on='SNP', how='inner')
            if len(m) < 2:
                print(f"No overlap instruments for {g}", flush=True)
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
            chr_ = str(m['CHR'].astype(str).mode().iloc[0])
            center = int(np.median(m['POS'].astype(int).values))
            rows.append({'gene': g, 'chr': chr_ , 'center_pos': center, 'n_eqtl_region': int(len(m)), 'n_ilcco_region': int(len(m)), 'n_shared_rsids': int(len(m)), 'z_corr': z_corr, 'peak_shared': peak_shared, 'peak_dist_bp': peak_dist, 'Q': Q, 'Q_p': np.nan, 'egger_intercept': eg_int, 'egger_intercept_se': eg_int_se, 'egger_intercept_p': eg_int_p, 'concordance': conc})
            continue
        # sanitize POS
        eg['POS'] = pd.to_numeric(eg['POS'], errors='coerce')
        eg = eg.dropna(subset=['POS'])
        eg['POS'] = eg['POS'].astype(int)
        if len(eg)==0:
            print(f"No valid POS for {g}", flush=True)
            continue
        chr_ = mode_chr(eg['CHR'].astype(str).values)
        center = int(np.median(eg['POS'].values))
        start = center - 1000000
        end = center + 1000000
        e_reg = eg[(eg['CHR'].astype(str)==chr_) & (eg['POS'].astype(int)>=start) & (eg['POS'].astype(int)<=end)].copy()
        if len(e_reg)==0:
            # fall back to instrument-level colocalization when region extraction empty
            e_reg = pd.DataFrame(columns=['RSID','POS','z_eqtl'])
        pvals = e_reg['P_VALUE'].astype(float).values
        pvals = np.clip(pvals, 1e-300, 1.0)
        z_eqtl = []
        rhos = e_reg['SPEARMANS_RHO'].values
        for i, p in enumerate(pvals):
            if np.isfinite(p) and 0.0 < float(1 - float(p)/2) < 1.0:
                z = NormalDist().inv_cdf(1 - float(p)/2)
            else:
                z = np.nan
            z_eqtl.append(z * np.sign(rhos[i]))
        e_reg['z_eqtl'] = np.array(z_eqtl)
        o_reg = out[(out['CHR']==chr_) & (out['POS']>=start) & (out['POS']<=end)].copy()
        if len(o_reg)==0:
            # fall back to instrument-level colocalization when region extraction empty
            o_reg = pd.DataFrame(columns=['SNP','POS','z_ilcco'])
        # compute z with guard against division by zero
        se_safe = o_reg['se_by'].replace(0, np.nan)
        o_reg['z_ilcco'] = o_reg['by'].values / se_safe.values
        j = pd.merge(e_reg[['RSID','z_eqtl','POS']].rename(columns={'RSID':'SNP'}), o_reg[['SNP','z_ilcco','POS']].rename(columns={'POS':'POS_OUT'}), on='SNP', how='inner')
        z_corr = float(j[['z_eqtl','z_ilcco']].corr().iloc[0,1]) if len(j)>=3 else np.nan
        
        # Filter out rows with NaN/inf z-scores before argmax
        e_reg_valid = e_reg[np.isfinite(e_reg['z_eqtl'])]
        o_reg_valid = o_reg[np.isfinite(o_reg['z_ilcco'])]
        
        if len(e_reg_valid) > 0 and len(o_reg_valid) > 0:
            pe = e_reg_valid.iloc[int(np.argmax(np.abs(e_reg_valid['z_eqtl'].values)))]
            po = o_reg_valid.iloc[int(np.argmax(np.abs(o_reg_valid['z_ilcco'].values)))]
            peak_shared = int(str(pe['RSID']) == str(po['SNP']))
            peak_dist = int(abs(int(pe['POS']) - int(po['POS'])))
        else:
            # instrument-level fallback if possible
            exp_file = EXP_DIR / f"EXPOSURE_cd8et_{g}_mr.tsv.gz"
            if exp_file.exists():
                d = pd.read_csv(exp_file, sep='\t', compression='gzip')
                d['SNP'] = d['SNP'].astype(str)
                d['z_eqtl_ins'] = d['beta'].astype(float) / d['se'].astype(float)
                mo = pd.merge(d[['SNP','z_eqtl_ins']], out[['SNP','by','se_by','POS']], on='SNP', how='inner')
                if len(mo) >= 2:
                    z_ilcco_ins = mo['by'].astype(float) / mo['se_by'].replace(0, np.nan).astype(float)
                    z_corr = float(pd.Series(mo['z_eqtl_ins']).corr(z_ilcco_ins))
                    pe_idx = int(np.argmax(np.abs(mo['z_eqtl_ins'].values)))
                    po_idx = int(np.argmax(np.abs(z_ilcco_ins.values)))
                    peak_shared = int(str(mo.iloc[pe_idx]['SNP']) == str(mo.iloc[po_idx]['SNP']))
                    # use outcome POS for distance when available
                    try:
                        peak_dist = int(abs(int(mo.iloc[pe_idx]['POS']) - int(mo.iloc[po_idx]['POS'])))
                    except Exception:
                        peak_dist = -1
                else:
                    peak_shared = 0
                    peak_dist = -1
            else:
                peak_shared = 0
                peak_dist = -1
        exp_file = EXP_DIR / f"EXPOSURE_cd8et_{g}_mr.tsv.gz"
        Q = np.nan; Q_p = np.nan; eg_int = np.nan; eg_int_se = np.nan; eg_int_p = np.nan; conc = np.nan
        if exp_file.exists():
            d = pd.read_csv(exp_file, sep='\t', compression='gzip')
            d['SNP'] = d['SNP'].astype(str)
            m = pd.merge(d, out[['SNP','ea_outcome','nea_outcome','by','se_by']], on='SNP', how='inner')
            if len(m)>=2:
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
                Q_p = np.nan
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
        rows.append({'gene': g, 'chr': chr_, 'center_pos': center, 'n_eqtl_region': int(len(e_reg)), 'n_ilcco_region': int(len(o_reg)), 'n_shared_rsids': int(len(j)), 'z_corr': z_corr, 'peak_shared': peak_shared, 'peak_dist_bp': peak_dist, 'Q': Q, 'Q_p': Q_p, 'egger_intercept': eg_int, 'egger_intercept_se': eg_int_se, 'egger_intercept_p': eg_int_p, 'concordance': conc})
    print('rows_len', len(rows))
    if rows:
        out_path = Path("My_MR_Project/CD8ET_ILCCO_coloc_sensitivity.csv")
        pd.DataFrame(rows).to_csv(out_path, index=False)
        print('saved:', out_path)
    else:
        cols = ['gene','chr','center_pos','n_eqtl_region','n_ilcco_region','n_shared_rsids','z_corr','peak_shared','peak_dist_bp','Q','Q_p','egger_intercept','egger_intercept_se','egger_intercept_p','concordance']
        out_path = Path("My_MR_Project/CD8ET_ILCCO_coloc_sensitivity.csv")
        pd.DataFrame(columns=cols).to_csv(out_path, index=False)
        print('saved empty:', out_path)

if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        Path("My_MR_Project/coloc_err.txt").write_text(str(e))
