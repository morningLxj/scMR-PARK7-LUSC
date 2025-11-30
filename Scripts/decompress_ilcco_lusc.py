import pandas as pd
from pathlib import Path

lusc_gz = Path('My_MR_Project/Outcome/28604730-GCST004750-EFO_0000708.h.tsv.gz')
lusc_out = Path('My_MR_Project/Outcome/28604730-GCST004750-EFO_0000708.h.tsv')

header_written = False
for chunk in pd.read_csv(lusc_gz, sep='\t', compression='gzip', chunksize=500000):
    chunk.to_csv(lusc_out, sep='\t', index=False, header=not header_written, mode=('a' if header_written else 'w'))
    header_written = True
print('written:', lusc_out)

