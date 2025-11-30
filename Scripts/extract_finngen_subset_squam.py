import pandas as pd
from pathlib import Path

squam_gz = Path("My_MR_Project/Outcome/finngen_R12_C3_NSCLC_SQUAM_EXALLC.gz")
squam_out = Path("My_MR_Project/Outcome/finngen_R12_C3_NSCLC_SQUAM_subset.csv")
usecols = ['#chrom','pos','rsids','beta','sebeta','pval']

header_written = False
for chunk in pd.read_csv(squam_gz, sep='\t', compression='gzip', usecols=usecols, chunksize=300000):
    chunk = chunk.copy()
    chunk['rsid'] = chunk['rsids'].astype(str).str.split(';').str[0]
    chunk.rename(columns={'#chrom': 'chrom', 'sebeta': 'se', 'pval': 'p'}, inplace=True)
    out = chunk[['rsid','chrom','pos','beta','se','p']]
    out.to_csv(squam_out, index=False, mode=('a' if header_written else 'w'), header=not header_written)
    header_written = True
print("written:", squam_out)

