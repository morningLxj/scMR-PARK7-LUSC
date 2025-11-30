import pandas as pd
from pathlib import Path

adeno_gz = Path("My_MR_Project/Outcome/finngen_R12_C3_NSCLC_ADENO_EXALLC.gz")
squam_gz = Path("My_MR_Project/Outcome/finngen_R12_C3_NSCLC_SQUAM_EXALLC.gz")

adeno_out = Path("My_MR_Project/Outcome/finngen_R12_C3_NSCLC_ADENO_subset.csv")
squam_out = Path("My_MR_Project/Outcome/finngen_R12_C3_NSCLC_SQUAM_subset.csv")

usecols = ['#chrom','pos','rsids','beta','sebeta','pval']

def write_subset(src_gz: Path, dst_csv: Path):
    header_written = False
    for chunk in pd.read_csv(src_gz, sep='\t', compression='gzip', usecols=usecols, chunksize=500000):
        chunk = chunk.copy()
        chunk['rsid'] = chunk['rsids'].astype(str).str.split(';').str[0]
        chunk.rename(columns={'#chrom': 'chrom', 'sebeta': 'se', 'pval': 'p'}, inplace=True)
        out = chunk[['rsid','chrom','pos','beta','se','p']]
        out.to_csv(dst_csv, index=False, mode=('a' if header_written else 'w'), header=not header_written)
        header_written = True

write_subset(adeno_gz, adeno_out)
write_subset(squam_gz, squam_out)
print("written:", adeno_out, squam_out)

