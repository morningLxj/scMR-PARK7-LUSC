import pandas as pd
import numpy as np
from pathlib import Path

luad_path = Path('My_MR_Project/Subtype_Analysis_LUAD.csv')
lusc_path = Path('My_MR_Project/Subtype_Analysis_LUSC.csv')
out_path = Path('My_MR_Project/Final_Subtype_Specificity.csv')

def load(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    for c in ['pp4_finngen','pp4_ilcco']:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors='coerce')
        else:
            df[c] = np.nan
    return df

def summarize(df: pd.DataFrame) -> pd.DataFrame:
    g = df.groupby(['gene','cell'], as_index=False).agg({
        'pp4_finngen': lambda s: np.nanmax(s.values) if len(s)>0 else np.nan,
        'pp4_ilcco': lambda s: np.nanmax(s.values) if len(s)>0 else np.nan,
    })
    g['pp4'] = g['pp4_finngen'].where(~g['pp4_finngen'].isna(), g['pp4_ilcco'])
    return g[['gene','cell','pp4']]

def classify(a: float, b: float) -> str:
    a_nan = pd.isna(a)
    b_nan = pd.isna(b)
    if a_nan and b_nan:
        return 'No Data'
    a_v = 0.0 if a_nan else a
    b_v = 0.0 if b_nan else b
    if (a_v < 0.05) and (b_v < 0.05):
        return 'Low/No Evidence'
    thr = 0.2
    if not a_nan and not b_nan:
        if (b_v - a_v) >= thr and b_v >= 0.1:
            return 'LUSC-Specific'
        if (a_v - b_v) >= thr and a_v >= 0.1:
            return 'LUAD-Specific'
    return 'Comparable/Uncertain'

luad_df = load(luad_path)
lusc_df = load(lusc_path)

luad_sum = summarize(luad_df).rename(columns={'pp4':'PP4_LUAD'})
lusc_sum = summarize(lusc_df).rename(columns={'pp4':'PP4_LUSC'})

all_keys = pd.merge(luad_sum[['gene','cell']], lusc_sum[['gene','cell']], how='outer', on=['gene','cell'])
res = pd.merge(all_keys, luad_sum, how='left', on=['gene','cell'])
res = pd.merge(res, lusc_sum, how='left', on=['gene','cell'])
res['Delta'] = (res['PP4_LUSC'].fillna(0) - res['PP4_LUAD'].fillna(0)).abs()
res['Specificity'] = [classify(a,b) for a,b in zip(res['PP4_LUAD'], res['PP4_LUSC'])]
res = res.sort_values(['Delta','gene','cell'], ascending=[False, True, True])

res.to_csv(out_path, index=False)

print('============================================================')
print('SUBTYPE SPECIFICITY REPORT')
print('============================================================')
for gene, cell in [('PARK7','bin'), ('CTSW','nk')]:
    row = res[(res['gene']==gene) & (res['cell']==cell)]
    if len(row) == 0:
        print(f'{gene}\t{cell}\tNA\tNA\tNo Data\t0')
    else:
        r = row.iloc[0]
        a = r['PP4_LUAD'] if not pd.isna(r['PP4_LUAD']) else 'NA'
        b = r['PP4_LUSC'] if not pd.isna(r['PP4_LUSC']) else 'NA'
        print(f'{gene}\t{cell}\t{a}\t{b}\t{r["Specificity"]}\t{r["Delta"]}')
print('written:', out_path)
