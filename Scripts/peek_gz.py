import gzip
import sys
from pathlib import Path

p = Path(sys.argv[1])
out = Path(sys.argv[2])
with gzip.open(p, 'rt', encoding='utf-8', errors='replace') as f:
    lines = [next(f) for _ in range(10)]
out.write_text(''.join(lines), encoding='utf-8')
print('written', out)
