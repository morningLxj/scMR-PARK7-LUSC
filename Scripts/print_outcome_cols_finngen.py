import pandas as pd
from pathlib import Path
OUTCOME = Path("My_MR_Project/Outcome/finngen_R12_C3_BRONCHUS_LUNG_EXALLC.gz")
OUT = Path("My_MR_Project/out_cols_finngen_py.txt")
df = pd.read_csv(OUTCOME, sep='\t', compression='gzip', nrows=5)
OUT.write_text(",".join(df.columns.tolist()), encoding='utf-8')
print("written:", OUT)
