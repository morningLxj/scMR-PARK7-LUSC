import requests
import pandas as pd
import time
from pathlib import Path

genes = ['PARK7', 'CTSW', 'TMEM50A', 'RPS26', 'PNRC2']
out = Path('My_MR_Project/Drug_Repurposing_DGIdb.csv')

def req_json(session, url, params, retries=3, delay=1.0):
    for i in range(retries):
        r = session.get(url, params=params, headers={'Accept': 'application/json', 'User-Agent': 'Mozilla/5.0'}, allow_redirects=True, timeout=15)
        if r.status_code == 200 and r.headers.get('Content-Type','').startswith('application/json'):
            try:
                return r.json()
            except Exception:
                pass
        time.sleep(delay)
    return None

def query_gene(session, gene):
    data = req_json(session, 'https://dgidb.org/api/v2/interactions.json', {'genes': gene})
    res = []
    if not data:
        return res
    for m in data.get('matched_terms', []):
        g = m.get('gene_name')
        for it in m.get('interactions', []):
            res.append({
                'Target_Gene': g,
                'Drug_Name': it.get('drug_name'),
                'Interaction_Type': ', '.join(it.get('interaction_types') or []) if it.get('interaction_types') is not None else '',
                'Score': it.get('score'),
                'Source': ', '.join(it.get('source_names') or []) if it.get('source_names') is not None else ''
            })
    return res

session = requests.Session()
all_rows = []
for g in genes:
    rows = query_gene(session, g)
    if rows:
        all_rows.extend(rows)
    time.sleep(1.0)

if all_rows:
    df = pd.DataFrame(all_rows)
    df = df.sort_values(['Target_Gene', 'Score'], ascending=[True, False])
    df.to_csv(out, index=False)
    print('Success:', len(df), 'interactions')
    print('Saved to:', out)
    print(df.head(10).to_string(index=False))
else:
    pd.DataFrame(columns=['Target_Gene','Drug_Name','Interaction_Type','Score','Source']).to_csv(out, index=False)
    print('No interactions found; wrote empty file to', out)
