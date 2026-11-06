"""critical_separation.py -- rigorous enclosure of the critical separation delta_c (midpoint quadrature with
remainder, jets_nested.py; q an interval through qint; printed end points with directed rounding; t* as exact
binary numbers in hexadecimal).
Assumption certified elsewhere: Delta is unimodal (unimodality.py + edges.py).
(i) encloses ell(t) - 2R0 at t_c - h and t_c + h (t_c a 30-digit estimate) and checks opposite signs: the root t*
    lies in (t_c - h, t_c + h);
(ii) Delta' < 0 on the box [t_c - h, t_c + h] (enclosed): the box lies on the shallow branch, Delta decreases there,
    so delta_c = Delta(t*) in [Delta(t_c + h), Delta(t_c - h)].
Uniqueness of t*: ell' = L Delta' < 0 on the shallow branch (after the maximum), and the deep branch is longer than
the dive (proof of the proposition on the certified grid), so ell - 2R0 has a single zero.
Usage (from this folder): python3 critical_separation.py ; the 20 cases run in parallel (10 processes), a few minutes."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'tools'))
from progress import Progress
import sys, time, hashlib
from multiprocessing import Pool
from mpmath import mp, mpf
mp.dps = 30
from jets_nested import iv, make, integrate, R0_enclosure, qint, down, up, manifest
from critical_separation_estimate import legs, dive

def run(case, K=512):
    qs, r0 = case; t0 = time.time(); q = float(qs); Qi = qint(qs)   # float q only for the estimate t_c
    D = dive(q, r0)
    ts = [2 + (r0 - 2)*mp.mpf(10)**(-7 + 7*k/40) for k in range(41)]; ts[-1] = r0 - (r0 - 2)*mp.mpf('1e-7')
    gaps = [legs(q, r0, t)[1] - D for t in ts]
    kc = next(i for i in range(40) if gaps[i]*gaps[i+1] < 0)
    tc = mp.findroot(lambda t: legs(q, r0, t)[1] - D, (ts[kc], ts[kc+1]), solver='anderson')
    h = float(1e-5*(r0 - 2))
    a, b = float(tc) - h, float(tc) + h
    G = make(Qi, r0, 'angle'); Gl = make(Qi, r0, 'length')
    twoR0 = R0_enclosure(Qi, r0, 4*K)
    la = integrate(Gl, iv.mpf(a), K)[0] - twoR0; lb = integrate(Gl, iv.mpf(b), K)[0] - twoR0
    dprime = integrate(G, iv.mpf([a, b]), K)[1]
    Da = integrate(G, iv.mpf(a), K)[0]; Db = integrate(G, iv.mpf(b), K)[0]
    ok = (la > 0) is True and (lb < 0) is True and (dprime < 0) is True
    lo, hi = (Db/iv.pi).a, (Da/iv.pi).b
    return (f"q={qs:>4} r0={r0:5.1f}: {'CERTIFIED' if ok else 'NOT certified'} "
            f"delta_c/pi in [{down(lo, 9)}, {up(hi, 9)}]  t* in ({a!r}, {b!r}) = ({a.hex()}, {b.hex()})  "
            f"(ell-2R0 >= {down(la.a, 9)} at t_c-h, <= {up(lb.b, 9)} at t_c+h; "
            f"Delta' <= {up(dprime.b, 6)} on box; {time.time()-t0:.0f}s)")

if __name__ == "__main__":
    print(manifest()); sys.stdout.flush()
    cases = [(qs, r0) for qs in ['1.2', '1.5', '2', '4', '10'] for r0 in [4.0, 6.0, 10.0, 20.0]]
    with Pool(10) as p:
        pg = Progress(len(cases), 'critical separation: (q, r0) cases')
        for line in p.imap(run, cases):
            print(line); sys.stdout.flush()
            pg.step()
        pg.close()
