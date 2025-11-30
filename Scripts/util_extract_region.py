import gzip, csv
import sys
from pathlib import Path

def main():
    if len(sys.argv) < 6:
        print('usage: util_extract_region.py <in.gz> <chrom> <start> <end> <out.tsv>')
        sys.exit(1)
    in_path = Path(sys.argv[1])
    chrom_req = sys.argv[2]
    start = int(sys.argv[3])
    end = int(sys.argv[4])
    out_path = Path(sys.argv[5])
    with gzip.open(in_path, 'rt', encoding='utf-8', errors='replace') as f, out_path.open('w', encoding='utf-8') as out:
        reader = csv.DictReader(f, delimiter='\t')
        out.write('\t'.join(['chrom','pos','rsids','beta','sebeta','pval'])+'\n')
        for row in reader:
            chrom = row.get('#chrom') or row.get('chrom')
            if chrom != chrom_req:
                continue
            pos = int(row['pos'])
            if start <= pos <= end:
                out.write('\t'.join([chrom, row['pos'], row['rsids'], row['beta'], row['sebeta'], row['pval']])+'\n')
    print('written', out_path)

if __name__ == '__main__':
    main()
