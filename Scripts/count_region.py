import gzip
import sys

fin = 'My_MR_Project/Outcome/finngen_R12_C3_BRONCHUS_LUNG_EXALLC.gz'
chrom = 11
pos0 = 116671040
win = 500000
start = pos0 - win
end = pos0 + win

n_region = 0
rsids = []
with gzip.open(fin, 'rt', encoding='utf-8', errors='replace') as g:
    header = g.readline().strip().split('\t')
    for line in g:
        parts = line.strip().split('\t')
        try:
            c = int(parts[0]); p = int(parts[1])
        except Exception:
            continue
        if c == chrom and start <= p <= end:
            n_region += 1
            rsids.append(parts[4])

open('My_MR_Project/count_region_out.txt','w').write(f'n_region={n_region}\nsample={rsids[:10]}')
