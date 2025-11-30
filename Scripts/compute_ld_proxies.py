import json
import gzip
import sys
import os
import numpy as np
from pysnptools.snpreader import Bed

bbase = 'My_MR_Project/Reference/g1000_eur'
fin = 'My_MR_Project/Outcome/finngen_R12_C3_BRONCHUS_LUNG_EXALLC.gz'
# allow passing anchors via CLI, e.g., python compute_ld_proxies.py rs12564169 rs7757989 1500000
args = sys.argv[1:]
win = 1500000
if args and args[-1].isdigit():
    win = int(args[-1])
    args = args[:-1]
targets = args if args else ['rs4956891','rs7757989']

bed = Bed(bbase, count_A1=True)
print('bed loaded')

sids = []
chrom = []
pos = []
with open(bbase + '.bim','r') as f:
    for line in f:
        p = line.strip().split()
        chrom.append(int(p[0]))
        sids.append(p[1])
        pos.append(int(p[3]))

info = {}
for t in targets:
    try:
        idx = sids.index(t)
        info[t] = {'idx': idx, 'chrom': chrom[idx], 'pos': pos[idx]}
    except ValueError:
        info[t] = None
print('info', info)

proxies = {}
for t, dat in info.items():
    if not dat:
        proxies[t] = []
        continue
    c = dat['chrom']
    p0 = dat['pos']
    idx = dat['idx']
    cand = [i for i,(cc,pp) in enumerate(zip(chrom,pos)) if cc==c and abs(pp-p0)<=win]
    if idx not in cand:
        cand.append(idx)
    print('cand', t, len(cand))
    X = bed.read(columns=cand).val
    print('read window', t, len(cand))
    xi = np.nan_to_num(X[:, cand.index(idx)])
    r2_list = []
    for j,i in enumerate(cand):
        if i == idx:
            continue
        xj = np.nan_to_num(X[:, j])
        r2 = 0.0
        if xj.std() != 0 and xi.std() != 0:
            r = float(np.corrcoef(xi, xj)[0,1])
            r2 = r*r
        r2_list.append({'target': t, 'rsid': sids[i], 'chrom': int(chrom[i]), 'pos': int(pos[i]), 'r2': r2})
    r2_list.sort(key=lambda x: -x['r2'])
    strong = [x for x in r2_list if x['r2'] >= 0.8]
    if not strong:
        strong = [x for x in r2_list if x['r2'] >= 0.6]
    proxies[t] = strong[:100]

out_candidates = 'My_MR_Project/ld_proxy_candidates.tsv'
with open(out_candidates,'w') as f:
    f.write('target\trsids\tchrom\tpos\tr2\n')
    for t, arr in proxies.items():
        for x in arr:
            f.write(f"{t}\t{x['rsid']}\t{x['chrom']}\t{x['pos']}\t{x['r2']:.6f}\n")

hits = []
with gzip.open(fin, 'rt', encoding='utf-8', errors='replace') as g:
    header = g.readline().strip().split('\t')
    rs_idx = header.index('rsids')
    beta_idx = header.index('beta')
    se_idx = header.index('sebeta')
    alt_idx = header.index('alt')
    ref_idx = header.index('ref')
    af_idx = header.index('af_alt')
    p_idx = header.index('pval')
    rsset = set([x['rsid'] for vals in proxies.values() for x in vals])
    coordset = set([(x['chrom'], x['pos']) for vals in proxies.values() for x in vals])
    for line in g:
        parts = line.strip().split('\t')
        if len(parts) <= max(rs_idx, beta_idx, se_idx, alt_idx, ref_idx, af_idx, p_idx):
            continue
        rsid = parts[rs_idx]
        if rsid in rsset:
            hits.append({
                'rsid': rsid,
                'beta': parts[beta_idx],
                'se': parts[se_idx],
                'alt': parts[alt_idx],
                'ref': parts[ref_idx],
                'eaf': parts[af_idx],
                'pval': parts[p_idx]
            })
        else:
            try:
                cc = int(parts[0]); pp = int(parts[1])
                if (cc, pp) in coordset:
                    hits.append({
                        'rsid': rsid,
                        'beta': parts[beta_idx],
                        'se': parts[se_idx],
                        'alt': parts[alt_idx],
                        'ref': parts[ref_idx],
                        'eaf': parts[af_idx],
                        'pval': parts[p_idx]
                    })
            except:
                pass

out_hits = 'My_MR_Project/ld_proxy_hits_in_finngen.tsv'
with open(out_hits,'w') as f:
    f.write('rsid\tbeta\tse\talt\tref\teaf\tpval\n')
    for h in hits:
        f.write('\t'.join([h['rsid'], h['beta'], h['se'], h['alt'], h['ref'], h['eaf'], h['pval']]) + '\n')
print('written', out_candidates, 'and', out_hits)
