import pandas as pd
from pathlib import Path
OUTCOME = Path("My_MR_Project/Outcome/28604730-GCST004748-EFO_0001071.h.tsv")
OUT = Path("My_MR_Project/out_cols_py.txt")
df = pd.read_csv(OUTCOME, sep='\t', nrows=5)
OUT.write_text(",".join(df.columns.tolist()), encoding='utf-8')
print("written:", OUT)
