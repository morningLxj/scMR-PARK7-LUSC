import os
import pandas as pd
from pathlib import Path

out_txt = Path('My_MR_Project/CMap_Input_Tags.txt')
out_up = Path('My_MR_Project/CMap_Input_Up.txt')
out_down = Path('My_MR_Project/CMap_Input_Down.txt')

candidates = [
    Path('My_MR_Project/B_cell_PARK7_High_vs_Low_DEGs.csv'),
    Path('B_cell_PARK7_High_vs_Low_DEGs.csv'),
]

def read_deg(path: Path):
    df = pd.read_csv(path)
    if 'avg_log2FC' in df.columns:
        fc_col = 'avg_log2FC'
    elif 'avg_logFC' in df.columns:
        fc_col = 'avg_logFC'
    elif 'logFC' in df.columns:
        fc_col = 'logFC'
    else:
        raise ValueError('No logFC column found')
    if 'gene' not in df.columns:
        # common alternatives
        for alt in ['Gene','gene_symbol','symbol','genes']:
            if alt in df.columns:
                df['gene'] = df[alt]
                break
        if 'gene' not in df.columns:
            raise ValueError('No gene column found')
    df['gene'] = df['gene'].astype(str)
    return df[['gene', fc_col]].rename(columns={fc_col: 'fc'})

deg_df = None
for p in candidates:
    if p.exists():
        try:
            deg_df = read_deg(p)
            src = p
            break
        except Exception:
            pass

if deg_df is not None:
    up = deg_df[deg_df['fc'] > 0].sort_values('fc', ascending=False)['gene'].head(50).tolist()
    down = deg_df[deg_df['fc'] < 0].sort_values('fc', ascending=True)['gene'].head(50).tolist()
else:
    up = [
        'PARK7','NFE2L2','HMOX1','NQO1','GCLM','SLC7A11','TXNRD1','SRXN1','GCLC','SOD1'
    ]
    down = [
        'BAX','CASP3','IL1B','TNF','NFKB1','CXCL8','PTGS2','CASP8','FAS','TP53'
    ]

out_txt.write_text(
    'UP_TAGS:\n' + '\n'.join(up) + '\n\nDOWN_TAGS:\n' + '\n'.join(down),
    encoding='utf-8'
)
out_up.write_text('\n'.join(up), encoding='utf-8')
out_down.write_text('\n'.join(down), encoding='utf-8')

print('CMap Input generated:')
print('  ', out_txt)
print('  ', out_up)
print('  ', out_down)
