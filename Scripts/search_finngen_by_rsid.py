import pandas as pd
from pathlib import Path
OUTCOME_FG = Path("My_MR_Project/Outcome/finngen_R12_C3_BRONCHUS_LUNG_EXALLC.gz")
target = 'rs4951859'
df = pd.read_csv(OUTCOME_FG, sep='\t', compression='gzip', usecols=['rsids'], nrows=2000000)
df['rsids'] = df['rsids'].fillna('')
mask = df['rsids'].str.contains(target)
Path('My_MR_Project/fg_rsid_rs4951859.txt').write_text(df.loc[mask].head(10).to_csv(index=False), encoding='utf-8')
print('written My_MR_Project/fg_rsid_rs4951859.txt')
