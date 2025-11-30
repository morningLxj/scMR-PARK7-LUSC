import gzip
import pandas as pd
from pathlib import Path

MR = Path("My_MR_Project/CD8ET_TopN_MR_Results.csv")
EQTL = Path("My_MR_Project/Exposure/cd8et_eqtl_table.tsv.gz")
OUT_FILTERED = Path("My_MR_Project/cd8et_coloc_filtered.tsv")

def main():
    print("Reading MR results...")
    mr = pd.read_csv(MR)
    mr['p_ivw'] = pd.to_numeric(mr['p_ivw'], errors='coerce')
    mr = mr[mr['p_ivw'] < 5e-8]
    targets = set([x.replace('cd8et_','') for x in mr['gene'].astype(str).tolist()])
    print(f"Targets: {targets}")
    
    if not targets:
        print("No targets found.")
        return

    print(f"Scanning {EQTL}...")
    
    with gzip.open(EQTL, 'rt') as f_in, open(OUT_FILTERED, 'w') as f_out:
        header = f_in.readline()
        f_out.write(header)
        
        cols = header.strip().split('\t')
        try:
            gene_idx = cols.index('GENE')
        except ValueError:
            print("GENE column not found in header")
            return
            
        count = 0
        match_count = 0
        for line in f_in:
            count += 1
            if count % 1000000 == 0:
                print(f"Processed {count} lines, found {match_count} matches...", flush=True)
            
            parts = line.split('\t')
            if len(parts) > gene_idx:
                gene = parts[gene_idx]
                if gene in targets:
                    f_out.write(line)
                    match_count += 1
                    
    print(f"Done. Scanned {count} lines. Saved {match_count} lines to {OUT_FILTERED}")

if __name__ == '__main__':
    main()
