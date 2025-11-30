import csv
in_path = 'My_MR_Project/Complete_Analysis_Summary.csv'
out_path = 'My_MR_Project/Final_Gold_Standard.csv'
rows = []
with open(in_path, 'r', newline='', encoding='utf-8') as f:
    reader = csv.DictReader(f)
    for r in reader:
        try:
            fpp4 = float(r.get('pp4_finngen', ''))
        except Exception:
            fpp4 = -1.0
        try:
            ipp4 = float(r.get('pp4_ilcco', ''))
        except Exception:
            ipp4 = -1.0
        if fpp4 >= 0.8 or ipp4 >= 0.8:
            rows.append(r)
with open(out_path, 'w', newline='', encoding='utf-8') as f:
    if rows:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    else:
        f.write('gene,cell,pp4_finngen,pp4_ilcco\n')
        f.write('# No pairs meet PP4 >= 0.8 in current results\n')
print('Saved:', out_path)
