import pandas as pd
import numpy as np
try:
    from scipy.stats import norm
    def invcdf(p):
        return norm.ppf(p)
except Exception:
    from statistics import NormalDist
    def invcdf(p):
        return NormalDist().inv_cdf(p)

INPUT_FILE = "My_MR_Project/Exposure/cd8et_eqtl_table.tsv.gz"
OUTPUT_FILE = "My_MR_Project/Exposure/CD8ET_PDCD1_exposure_mr.tsv.gz"
TARGET_GENE_ID = "ENSG00000188389"
P_THRESHOLD = 5e-8

def calculate_se(p_value, beta):
    p = np.where(p_value == 0, np.finfo(float).tiny, p_value)
    z = invcdf(1 - p / 2)
    se = np.where(z != 0, np.abs(beta / z), np.nan)
    return se

def prepare_exposure_data(input_path, output_path, gene_id, p_thresh):
    df = pd.read_csv(input_path, sep='\t', compression='gzip')
    df = df[(df['GENE_ID'] == gene_id) & (df['P_VALUE'] < p_thresh)].copy()
    if df.empty:
        print(f"no significant eQTL for {gene_id} with P<{p_thresh}")
        return
    df['se'] = calculate_se(df['P_VALUE'], df['SPEARMANS_RHO'])
    df_mr = pd.DataFrame({
        'SNP': df['RSID'],
        'effect_allele': df['A2'],
        'other_allele': df['A1'],
        'eaf': df['A2_FREQ_ONEK1K'],
        'beta': df['SPEARMANS_RHO'],
        'se': df['se'],
        'pval': df['P_VALUE'],
        'exposure': 'CD8ET eQTL (PDCD1)',
        'id.exposure': 'CD8ET_PDCD1'
    })
    df_mr.dropna(subset=['se'], inplace=True)
    df_mr.to_csv(output_path, sep='\t', index=False, compression='gzip')
    print(f"saved {len(df_mr)} IVs to {output_path}")

if __name__ == "__main__":
    prepare_exposure_data(INPUT_FILE, OUTPUT_FILE, TARGET_GENE_ID, P_THRESHOLD)

