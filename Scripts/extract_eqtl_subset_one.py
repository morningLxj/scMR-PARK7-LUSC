import os
import pandas as pd
from pathlib import Path

cell = os.environ.get('CELL')
assert cell, 'CELL env required'

base = Path('My_MR_Project/Exposure')
mr_csv = Path('My_MR_Project/AllCells_MR_Results.csv')
sum_csv = Path('My_MR_Project/Complete_Analysis_Summary.csv')

mr = pd.read_csv(mr_csv)
gene_pool = set(mr.loc[mr['p_ivw'] < 0.05, 'gene'].astype(str))
try:
    sm = pd.read_csv(sum_csv)
    gene_pool |= set(sm.loc[sm['p_ivw'] < 0.05, 'gene'].astype(str))
except Exception:
    pass

src = base / f"{cell}_eqtl_table.tsv.gz"
dst = base / f"eqtl_subset_{cell}.csv"
usecols = ['RSID','CHR','POS','SPEARMANS_RHO','P_VALUE','A2_FREQ_ONEK1K','GENE']
header_written = False
print('processing', cell, 'from', src)
for chunk in pd.read_csv(src, sep='\t', compression='gzip', usecols=usecols, chunksize=300000):
    sub = chunk[chunk['GENE'].astype(str).isin(gene_pool)]
    if sub.empty:
        continue
    sub.to_csv(dst, index=False, mode=('a' if header_written else 'w'), header=not header_written)
    header_written = True
print('written', dst, 'header_written=', header_written)

