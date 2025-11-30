import pandas as pd
from pathlib import Path
EQTL = Path("My_MR_Project/Exposure/cd8et_eqtl_table.tsv.gz")
OUT = Path("My_MR_Project/eqtl_p1e6_stats.txt")
total = 0
sig = 0
genes = {}
# read a limited number of rows for quick stats
df = pd.read_csv(EQTL, sep='\t', compression='gzip', usecols=['GENE','P_VALUE'], nrows=500000)
total = len(df)
m = df['P_VALUE'] < 1e-6
sig = int(m.sum())
vc = df.loc[m, 'GENE'].value_counts()
for g, c in vc.items():
    genes[g] = int(c)
lines = [
    f"total_rows {total}",
    f"sig_rows_p<1e-6 {sig}",
    f"sig_genes_count {len(genes)}",
    f"top_genes {sorted(genes.items(), key=lambda x: x[1], reverse=True)[:10]}"
]
OUT.write_text("\n".join(lines), encoding='utf-8')
print('written:', OUT)
