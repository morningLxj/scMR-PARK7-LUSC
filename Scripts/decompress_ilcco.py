import pandas as pd
from pathlib import Path

luad_gz = Path('My_MR_Project/Outcome/28604730-GCST004744-EFO_0000571.h.tsv.gz')
lusc_gz = Path('My_MR_Project/Outcome/28604730-GCST004750-EFO_0000708.h.tsv.gz')

luad_out = Path('My_MR_Project/Outcome/28604730-GCST004744-EFO_0000571.h.tsv')
lusc_out = Path('My_MR_Project/Outcome/28604730-GCST004750-EFO_0000708.h.tsv')

cols = None  # write full columns to preserve mapping flexibility

def gunzip_tsv(src_gz: Path, dst_tsv: Path):
    header_written = False
    for chunk in pd.read_csv(src_gz, sep='\t', compression='gzip', chunksize=500000):
        mode = 'a' if header_written else 'w'
        chunk.to_csv(dst_tsv, sep='\t', index=False, header=not header_written, mode=mode)
        header_written = True
    print('written:', dst_tsv)

gunzip_tsv(luad_gz, luad_out)
gunzip_tsv(lusc_gz, lusc_out)

