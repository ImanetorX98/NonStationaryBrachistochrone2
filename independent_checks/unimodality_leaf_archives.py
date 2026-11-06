"""Independent check of the archived leaves of the unimodality certificates (../interval_arithmetic/leaves_unimodality):
for each of the 20 grid cases, the leaves tile [t_lo, t_hi] exactly (shared binary end points, read from their
hexadecimal form), every printed bound has the sign required by its label (P: Delta' > 0, N: Delta' < 0,
C: Delta'' < 0), and the ordered labels satisfy the criterion [PC]*[CN]* with at least one C.
It reads only the JSON archives and imports no other file of this repository.
Written through OpenAI Codex as a second, independent implementation (see the paragraph on AI tools in the paper).
Usage: python3 unimodality_leaf_archives.py"""
# Optional progress bar (../tools/progress.py, standard error only); without it this file is self-contained.
try:
    import os as _os, sys as _sys
    _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'tools'))
    from progress import Progress
except ImportError:
    class Progress:
        def __init__(self, total, label=''): pass
        def step(self, n=1): pass
        def close(self): pass
import json, re, sys
from decimal import Decimal
from pathlib import Path

ARCH = Path(__file__).resolve().parent.parent/'interval_arithmetic'/'leaves_unimodality'
paths = sorted(ARCH.glob('*.json'))
assert len(paths) == 20, f'expected 20 archives, found {len(paths)}'
pg = Progress(len(paths), 'leaf archives')
for path in paths:
    data = json.loads(path.read_text()); leaves = data['leaves']
    assert data['certified'] is True
    lo = float.fromhex(data['t_lo']); hi = float.fromhex(data['t_hi'])
    assert float.fromhex(leaves[0][0]) == lo and float.fromhex(leaves[-1][1]) == hi
    seq = ''
    for i, l in enumerate(leaves):
        a, b, label, stage, K, direction, bound = l
        assert float.fromhex(a) < float.fromhex(b)
        if i: assert a == leaves[i-1][1]
        seq += label
        assert Decimal(bound) > 0 if label == 'P' else Decimal(bound) < 0
    assert re.fullmatch('[PC]*[CN]*', seq) and 'C' in seq
    print(f'PASS leaves {path.name}: {len(leaves)} cells, complete tiling, strict printed bounds, pattern.', flush=True)
    pg.step()
pg.close()
print(f'SUMMARY: {len(paths)}/20 archives PASS')
