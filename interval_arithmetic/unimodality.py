"""unimodality.py -- interval certificate of the unimodality of Delta(t) on [t_lo, t_hi] (supplementary material,
"The certified grid", part of the proof of the proposition on the certified grid of critical separations).

Bisection in t on [t_lo, t_hi] = [2 + eps (r0-2), r0 - eps (r0-2)] with labels P (Delta' > 0), C (Delta'' < 0),
N (Delta' < 0). Criterion: the leaves tile [t_lo, t_hi] (shared binary end points), there is at least one C and no
N precedes a P ([PC]*[CN]*). With Delta' > 0 at t_lo and Delta' < 0 at t_hi (edges.py), Delta' has exactly one zero,
non-degenerate, in the block of C between the last P and the first N: a C between two P's has Delta' decreasing and
positive at its right end.
Quadrature in u: first the zeroth-order one of jets_t.py, then midpoint + Lagrange remainder with nested jets
(jets_nested.py), error O(1/K^2).
q is an interval (jets_nested.qint): for '1.2' it contains both 6/5 and the binary64 value.
Every leaf (exact end points in hexadecimal, label, stage, K, bound proved with directed rounding) is archived in
leaves_unimodality/.
Usage (from this folder): python3 unimodality.py [q,r0] ; the 20 cases run in parallel (10 processes), several hours."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'tools'))
from progress import Progress
import sys, time, re, hashlib
from multiprocessing import Pool
from jets_nested import iv, make, integrate, qint, down, up, manifest
import json, os
from jets_t import make as make1, enclose as enclose1

def certify(qs, r0, eps=1e-4, Kmax=1024, minw=1e-10):
    Qi = qint(qs)
    G = make(Qi, r0, 'angle'); G1 = make1(Qi, r0)
    tlo = float(2 + eps*(r0 - 2)); thi = float(r0 - eps*(r0 - 2))
    stack = [(tlo, thi)]; leaves = []; evals = 0
    while stack:
        a, b = stack.pop(); T = iv.mpf([a, b]); label = None   # a, b binary floats: exact endpoints
        # stage 1: fast zeroth-order quadrature of jets_t.py (also rigorous, O(1/K) width)
        for K1 in (64, 256):
            try:
                d1, d2 = enclose1(G1, T, K1); evals += K1
            except Exception:
                break
            if d1 > 0: label = 'P'; proof = ('v1', K1, "Delta'>=", down(d1.a, 12))
            elif d1 < 0: label = 'N'; proof = ('v1', K1, "Delta'<=", up(d1.b, 12))
            elif d2 < 0: label = 'C'; proof = ('v1', K1, "Delta''<=", up(d2.b, 12))
            if label: break
        # stage 2: midpoint + Lagrange remainder, O(1/K^2) width
        K = 16
        while label is None and K <= Kmax:
            try:
                _, d1, d2 = integrate(G, T, K); evals += 2*K
            except Exception:          # sqrt/division not certified on this (too wide) box
                break
            if d1 > 0: label = 'P'; proof = ('v2', K, "Delta'>=", down(d1.a, 12)); break
            if d1 < 0: label = 'N'; proof = ('v2', K, "Delta'<=", up(d1.b, 12)); break
            if d2 < 0: label = 'C'; proof = ('v2', K, "Delta''<=", up(d2.b, 12)); break
            # refine K only if the box is already narrow (otherwise width in t dominates)
            if K*4 > Kmax or (d1.delta < 4*(b - a) and d2.delta < 4*(b - a)): break
            K *= 4
        if label: leaves.append((a, b, label, proof)); continue
        if b - a < minw: return None, [['?', 0]], evals, (a, b), None, []
        m = 0.5*(a + b); stack.append((m, b)); stack.append((a, m))
    leaves.sort()
    for (a1, b1, *_), (a2, b2, *_) in zip(leaves, leaves[1:]):
        assert b1 == a2, "leaves do not tile the interval"
    assert leaves[0][0] == tlo and leaves[-1][1] == thi
    seq = ''.join(l[2] for l in leaves)
    # Delta' > 0 on P, < 0 on N, strictly decreasing on C. If no N leaf precedes a P leaf, a C leaf
    # between two P's has Delta' > 0 at its right end hence on all of it (decreasing), likewise C
    # between N's; with the certified edge signs (Delta' > 0 at t_lo, < 0 at t_hi) Delta' has exactly
    # one zero, inside the C block between the last P and the first N.
    ok = re.fullmatch(r'[PC]*[CN]*', seq) is not None and 'C' in seq
    comp = []
    for l in (x[2] for x in leaves):
        if comp and comp[-1][0] == l: comp[-1][1] += 1
        else: comp.append([l, 1])
    # t-range of the C block (contains the unique maximum)
    lastP = max([i for i, x in enumerate(leaves) if x[2] == 'P'], default=-1)
    firstN = min([i for i, x in enumerate(leaves) if x[2] == 'N'], default=len(leaves))
    cs = [(x[0], x[1]) for x in leaves[lastP + 1:firstN]]
    cblock = (cs[0][0], cs[-1][1]) if cs else None
    return ok, comp, evals, (tlo, thi), cblock, leaves

def run(case):
    qs, r0 = case; t0 = time.time()
    ok, comp, evals, rng, cb, leaves = certify(qs, r0)
    tag = "CERTIFIED unimodal" if ok else ("NOT certified" if ok is False
          else f"FAILED to resolve at width limit near t in [{rng[0]!r}, {rng[1]!r}]")
    cbs = f"; max in t in [{cb[0]:.6f}, {cb[1]:.6f}]" if ok else ""
    os.makedirs('leaves_unimodality', exist_ok=True)
    with open(f'leaves_unimodality/q{qs}_r0{r0:g}.json', 'w') as fh:
        json.dump({'q': qs, 'q_interval': 'hull(iv.mpf(str), iv.mpf(float))', 'r0': r0,
                   't_lo': rng[0].hex() if ok is not None else None,
                   't_hi': rng[1].hex() if ok is not None else None, 'certified': ok,
                   'leaves': [[a.hex(), b.hex(), l, *p] for a, b, l, p in leaves]}, fh)
    return (f"q={qs:>4} r0={r0:5.1f}: {tag} on t in [{rng[0]!r}, {rng[1]!r}]; "
            f"pattern {''.join(f'{l}x{n} ' for l, n in comp)}{cbs}; {evals} integrand evals, {time.time()-t0:.0f}s")

if __name__ == "__main__":
    print(manifest()); sys.stdout.flush()
    cases = [(qs, r0) for qs in ['1.2', '1.5', '2', '4', '10'] for r0 in [4.0, 6.0, 10.0, 20.0]]
    if len(sys.argv) > 1:
        qs, r0 = sys.argv[1].split(','); cases = [(qs, float(r0))]
    with Pool(10) as p:
        pg = Progress(len(cases), 'unimodality: (q, r0) cases')
        for line in p.imap(run, cases):
            print(line); sys.stdout.flush()
            pg.step()
        pg.close()
