import csv
from collections import defaultdict
inp = 'My_MR_Project/Final_Classified_Results.csv'
outp = 'My_MR_Project/gene_lists_by_level.csv'
rows = []
with open(inp, 'r', newline='', encoding='utf-8') as f:
    r = csv.DictReader(f)
    for row in r:
        rows.append(row)
by_level = defaultdict(set)
for row in rows:
    lvl = row.get('Evidence_Level', 'Weak')
    gene = row.get('gene', '')
    by_level[lvl].add(gene)
with open(outp, 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f)
    w.writerow(['level', 'genes'])
    for lvl in ['Strong','Moderate','Suggestive','Weak']:
        genes = sorted(list(by_level.get(lvl, [])))
        w.writerow([lvl, ';'.join(genes)])
print('Saved:', outp)
