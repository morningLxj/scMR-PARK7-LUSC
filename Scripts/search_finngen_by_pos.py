import pandas as pd
from pathlib import Path
OUTCOME_FG = Path("My_MR_Project/Outcome/finngen_R12_C3_BRONCHUS_LUNG_EXALLC.gz")
df = pd.read_csv(OUTCOME_FG, sep='\t', compression='gzip', usecols=['#chrom','pos','rsids'], nrows=2000000)
hit = df[(df['#chrom'].astype(str)=='1') & (df['pos'].astype(int)==729679)]
Path('My_MR_Project/fg_pos_1_729679.txt').write_text(hit.head(10).to_csv(index=False), encoding='utf-8')
print('written My_MR_Project/fg_pos_1_729679.txt')
