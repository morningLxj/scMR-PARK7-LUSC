import pandas as pd
from pathlib import Path
EQTL = Path("My_MR_Project/Exposure/cd8et_eqtl_table.tsv.gz")
df = pd.read_csv(EQTL, sep='\t', compression='gzip', usecols=['GENE','P_VALUE'], nrows=1000000)
sub = df[df['P_VALUE']<1e-6]
top = sub['GENE'].value_counts().head(10)
Path('My_MR_Project/top_genes_fast.txt').write_text("\n".join(top.index.tolist()), encoding='utf-8')
print('written My_MR_Project/top_genes_fast.txt')
