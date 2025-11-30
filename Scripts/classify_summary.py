import csv

def classify_evidence(pp4):
    try:
        v = float(pp4)
    except Exception:
        return "Weak"
    if v >= 0.5:
        return "Strong"
    elif v >= 0.2:
        return "Moderate"
    elif v >= 0.05:
        return "Suggestive"
    else:
        return "Weak"

in_path = 'My_MR_Project/Complete_Analysis_Summary.csv'
out_path = 'My_MR_Project/Final_Classified_Results.csv'

rows = []
with open(in_path, 'r', newline='', encoding='utf-8') as f:
    reader = csv.DictReader(f)
    for r in reader:
        # normalize numeric strings
        r['pp4_finngen'] = r.get('pp4_finngen', '')
        r['pp4_ilcco'] = r.get('pp4_ilcco', '')
        r['p_ivw'] = r.get('p_ivw', '')
        # add alias column to match requested report
        r['mr_p'] = r['p_ivw']
        r['Evidence_Level'] = classify_evidence(r['pp4_finngen'])
        rows.append(r)

# sort by FinnGen PP4 desc
def to_float(x):
    try:
        return float(x)
    except Exception:
        return -1.0

rows_sorted = sorted(rows, key=lambda r: to_float(r['pp4_finngen']), reverse=True)

# write full classified table
fieldnames = list(rows_sorted[0].keys()) if rows_sorted else []
with open(out_path, 'w', newline='', encoding='utf-8') as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows_sorted)

# build highlights (non-Weak)
highlights = [r for r in rows_sorted if r['Evidence_Level'] != 'Weak']

strong_n = sum(1 for r in rows_sorted if r['Evidence_Level'] == 'Strong')
moderate_n = sum(1 for r in rows_sorted if r['Evidence_Level'] == 'Moderate')
suggestive_n = sum(1 for r in rows_sorted if r['Evidence_Level'] == 'Suggestive')

print('='*60)
print('       FINAL SC-MR PROJECT EXECUTIVE SUMMARY')
print('='*60)
print(f'Total pairs analyzed: {len(rows_sorted)}')
print(f'Strong Evidence (PP4 >= 0.5):     {strong_n}')
print(f'Moderate Evidence (PP4 0.2-0.5):  {moderate_n}')
print(f'Suggestive Evidence (PP4 0.05-0.2): {suggestive_n}')
print('-'*60)
print('TOP FINDINGS (Sorted by Evidence Strength):')

# select columns for display
cols = ['gene', 'cell', 'mr_p', 'pp4_finngen', 'pp4_ilcco', 'Evidence_Level']
if highlights:
    # compute column widths
    widths = {c: max(len(c), max(len(str(r.get(c, ''))) for r in highlights)) for c in cols}
    # header
    print(' '.join(c.ljust(widths[c]) for c in cols))
    # rows
    for r in highlights:
        print(' '.join(str(r.get(c, '')).ljust(widths[c]) for c in cols))
else:
    print('No non-Weak findings.')

print('='*60)
print(f'\nFull classified table saved to: {out_path}')

