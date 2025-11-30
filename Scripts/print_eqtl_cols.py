import pandas as pd
from pathlib import Path
EQTL = Path("My_MR_Project/Exposure/cd8et_eqtl_table.tsv.gz")
OUT = Path("My_MR_Project/eqtl_cols_py.txt")
df = pd.read_csv(EQTL, sep='\t', compression='gzip', nrows=5)
OUT.write_text(",".join(df.columns.tolist()), encoding='utf-8')
print("written:", OUT)
