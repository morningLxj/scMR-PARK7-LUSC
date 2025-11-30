import pandas as pd
from pathlib import Path

INPUT_FILE = "My_MR_Project/Exposure/cd8et_eqtl_table.tsv.gz"
P_THRESHOLD = 1e-6
TOP_N = 20
OUT_TXT = Path("My_MR_Project/top_genes_cd8et.txt")

def main():
    df = pd.read_csv(INPUT_FILE, sep='\t', compression='gzip', usecols=['GENE','P_VALUE'])
    sub = df[df['P_VALUE'] < P_THRESHOLD]
    cnt = sub['GENE'].value_counts().head(TOP_N)
    OUT_TXT.write_text("\n".join(cnt.index.tolist()), encoding='utf-8')
    print("written:", OUT_TXT)

if __name__ == "__main__":
    main()

