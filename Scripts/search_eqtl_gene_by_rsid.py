import pandas as pd
from pathlib import Path
EQTL = Path("My_MR_Project/Exposure/cd8et_eqtl_table.tsv.gz")
target = 'rs4951859'
df = pd.read_csv(EQTL, sep='\t', compression='gzip', usecols=['RSID','GENE','P_VALUE'], nrows=1000000)
hit = df[df['RSID'] == target]
Path('My_MR_Project/eqtl_gene_rs4951859.txt').write_text(hit.head(10).to_csv(index=False), encoding='utf-8')
print('written My_MR_Project/eqtl_gene_rs4951859.txt')
