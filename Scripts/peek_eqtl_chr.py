import pandas as pd
from pathlib import Path
EQTL = Path("My_MR_Project/Exposure/cd8et_eqtl_table.tsv.gz")
df = pd.read_csv(EQTL, sep='\t', compression='gzip', usecols=['CHR','POS','RSID'], nrows=10)
print(df)
