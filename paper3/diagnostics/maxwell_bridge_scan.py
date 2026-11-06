"""maxwell_bridge_scan.py -- grid diagnostic, not certified: scan of (q, r0) with q < 3/2 in search of interior
bridges of the convex hull of S (true Maxwell angles below delta_c). Uses maxwell_rule_scan.analyse."""
import sys
import mpmath as mp
from multiprocessing import Pool
from maxwell_rule_scan import analyse

def one(c):
    q, r0 = c
    try:
        r = analyse(q, r0, 160)
    except Exception as e:
        return f"q={q} r0={r0:g}: ERROR {e}"
    br = r['bridges']
    inner = [b for b in br if not b['from_dive']]
    return (f"q={q} r0={r0:g}: maxima {len(r['maxima'])} "
            f"[{', '.join(mp.nstr(D/mp.pi,5) for _, D in r['maxima'])}] pi; "
            f"delta_c ~ {', '.join(mp.nstr(b['delta']/mp.pi,6) for b in br if b['from_dive'])} pi; "
            f"interior bridges {len(inner)}"
            + ("".join(f" [angle {mp.nstr(b['delta']/mp.pi,6)} pi, t {mp.nstr(b['ta'],5)}-{mp.nstr(b['tb'],5)}]" for b in inner)))

if __name__ == "__main__":
    cases = [(q, r0) for q in [1.001, 1.01, 1.05, 1.1, 1.2, 1.3, 1.4] for r0 in [10.0, 1e2, 1e3, 1e4, 1e5, 1e6]]
    with Pool(10) as p:
        for line in p.imap(one, cases):
            print(line); sys.stdout.flush()
