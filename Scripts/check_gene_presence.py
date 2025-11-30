import pandas as pd
from pathlib import Path
EQTL = Path("My_MR_Project/Exposure/cd8et_eqtl_table.tsv.gz")
gene = "KCNAB2"
cnt = 0
for chunk in pd.read_csv(EQTL, sep='\t', compression='gzip', usecols=['GENE','P_VALUE'], chunksize=500000):
    cnt += int(((chunk['GENE']==gene) & (chunk['P_VALUE']<1e-6)).sum())
print("count_p1e-6", gene, cnt)
