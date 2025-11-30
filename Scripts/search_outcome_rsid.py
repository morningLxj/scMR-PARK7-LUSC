import pandas as pd
from pathlib import Path
ILCCO = Path("My_MR_Project/Outcome/28604730-GCST004748-EFO_0001071.h.tsv")
rs = "rs11121518"
df = pd.read_csv(ILCCO, sep='\t', usecols=['hm_rsid'], nrows=1000000)
hit = (df['hm_rsid'].astype(str) == rs).sum()
Path('My_MR_Project/outcome_hit_rs11121518.txt').write_text(str(int(hit)), encoding='utf-8')
print('written outcome_hit_rs11121518.txt')
