import glob
import pandas as pd
import numpy as np
from pathlib import Path
from statistics import NormalDist

EXP_DIR = Path("My_MR_Project/Exposure")
OUTCOME_FG = Path("My_MR_Project/Outcome/finngen_R12_C3_BRONCHUS_LUNG_EXALLC.gz")
OUTCOME_ILCCO = Path("My_MR_Project/Outcome/28604730-GCST004748-EFO_0001071.h.tsv")
EQTL_FULL = Path("My_MR_Project/Exposure/cd8et_eqtl_table.tsv.gz")
OUTCSV = Path("My_MR_Project/AllCells_MR_Results.csv")

def compute(out):
    files = sorted(glob.glob(str(EXP_DIR / 'EXPOSURE_*_mr.tsv.gz')))
    rows = []
    for f in files:
        d = pd.read_csv(f, sep='\t', compression='gzip')
        name = Path(f).name.replace('EXPOSURE_','').replace('_mr.tsv.gz','')
        parts = name.split('_')
        gene = parts[-1]
        cell = '_'.join(parts[:-1])
        d['SNP'] = d['SNP'].astype(str)
        m = pd.merge(d, out, on='SNP', how='inner')
        if m.empty:
            eqtl_path = EXP_DIR / f"{cell}_eqtl_table.tsv.gz"
            if eqtl_path.exists():
                eq = pd.read_csv(eqtl_path, sep='\t', compression='gzip', usecols=['GENE','RSID','CHR','POS'])
                eq = eq[eq['GENE'] == gene].rename(columns={'RSID':'SNP'})
                d2 = pd.merge(d, eq[['SNP','CHR','POS']], on='SNP', how='left')
                d2['CHR'] = d2['CHR'].astype(str)
                d2['POS'] = pd.to_numeric(d2['POS'], errors='coerce')
                d2 = d2.dropna(subset=['POS'])
                d2['POS'] = d2['POS'].astype(int)
                ilc = pd.read_csv(OUTCOME_ILCCO, sep='\t', low_memory=False)
                ilc = ilc.rename(columns={'hm_chrom':'CHR','hm_pos':'POS','hm_effect_allele':'ea_outcome','hm_other_allele':'nea_outcome','hm_beta':'by','standard_error':'se_by'})
                ilc['CHR'] = ilc['CHR'].astype(str)
                ilc['POS'] = pd.to_numeric(ilc['POS'], errors='coerce')
                ilc = ilc.dropna(subset=['POS'])
                ilc['POS'] = ilc['POS'].astype(int)
                m = pd.merge(d2, ilc[['CHR','POS','ea_outcome','nea_outcome','by','se_by']], on=['CHR','POS'], how='inner')
        if m.empty:
            continue
        idx = m['effect_allele'] != m['ea_outcome']
        if idx.any():
            m.loc[idx, 'by'] = -m.loc[idx, 'by']
        for col in ['beta','se','by','se_by']:
            m[col] = pd.to_numeric(m[col], errors='coerce')
        m = m.replace([np.inf, -np.inf], np.nan).dropna(subset=['beta','se','by','se_by'])
        m = m[(m['beta'] != 0) & (m['se'] > 0)]
        if len(m) == 1:
            bx = float(m['beta'].iloc[0])
            by = float(m['by'].iloc[0])
            se_bx = float(m['se'].iloc[0])
            se_by = float(m['se_by'].iloc[0])
            ivw = by / bx
            se_ivw = abs(se_by / bx)
            p_ivw = float(2 * (1 - NormalDist().cdf(abs(ivw / se_ivw))))
        else:
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
        rows.append({'cell': cell, 'gene': gene, 'n_iv': int(len(m)), 'beta_ivw': ivw, 'se_ivw': se_ivw, 'p_ivw': p_ivw})
    return rows

def main():
    fg = pd.read_csv(OUTCOME_FG, sep='\t', compression='gzip')
    fg = fg.rename(columns={'#chrom':'chrom','pos':'pos','ref':'nea_outcome','alt':'ea_outcome','beta':'by','sebeta':'se_by','af_alt':'eaf_outcome','rsids':'SNP','pval':'p_outcome'})
    fg['SNP'] = fg['SNP'].fillna('')
    fg = fg.assign(SNP=fg['SNP'].str.split(';')).explode('SNP')
    fg['SNP'] = fg['SNP'].str.strip()
    res = compute(fg)
    if not res:
        ilc = pd.read_csv(OUTCOME_ILCCO, sep='\t', low_memory=False)
        ilc = ilc.rename(columns={'hm_rsid':'SNP','hm_effect_allele':'ea_outcome','hm_other_allele':'nea_outcome','hm_beta':'by','standard_error':'se_by'})
        ilc['SNP'] = ilc['SNP'].astype(str)
        res = compute(ilc)
    if res:
        df = pd.DataFrame(res)
        df.to_csv(OUTCSV, index=False)
        base = {}
        for g in ['PARK7','CTSW','TMEM50A']:
            r = df[(df['gene']==g) & (df['cell']=='cd8et')]
            if len(r)>0:
                base[g] = float(np.sign(r['beta_ivw'].iloc[0]))
        sel = df[df['gene'].isin(['PARK7','CTSW','TMEM50A'])].copy()
        sel['significant'] = sel['p_ivw'] < 5e-8
        def cons(row):
            b = base.get(row['gene'], np.nan)
            if np.isnan(b):
                return np.nan
            return float(np.sign(row['beta_ivw']) == b)
        sel['direction_consistent_cd8et'] = sel.apply(cons, axis=1)
        sel.to_csv(Path('My_MR_Project/AllCells_MR_Summary.csv'), index=False)
        print('saved:', OUTCSV, 'rows', len(res))
    else:
        print('no MR results produced')

if __name__ == '__main__':
    main()
