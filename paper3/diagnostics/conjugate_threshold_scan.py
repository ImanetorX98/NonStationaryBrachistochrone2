"""Floating-point diagnostic, not a certificate: threshold of in-plane conjugate points before pi, constant q. For each q
and launch radius d, the smallest radius r_conj at which a ray from d has a zero of u with phi < pi; control q >= 3/2
(must give none). Source of the markers of Figure 2 (q = 1.05, 1.1, 1.2, 1.3). Several minutes (8 processes)."""
import numpy as np, sys
from multiprocessing import Pool
sys.path.insert(0, '.')
from conjugate_points_static import run
def Rq(q):
    if q >= 1.5: return np.inf
    return (33*q-38+np.sqrt((33*q-38)**2+480*(3-2*q)*(q-1)))/(6*(3-2*q)*(q-1))
def task(args):
    q, d = args
    best = None
    for a0 in np.concatenate([np.linspace(0.01, np.pi/2, 60), np.pi - np.geomspace(1e-4, np.pi/2, 240)]):
        zs, reach, upi, _ = run(q, d, a0, rmax=60*d, smax=1e8)
        if zs:
            rc, pc = zs[0]
            if best is None or max(d, rc) < max(d, best[0]): best = (rc, pc, a0)
    return q, d, best
if __name__ == '__main__':
    qs = [1.05, 1.1, 1.2, 1.3, 1.4, 1.45, 1.5, 1.6, 2.0]
    ds = [10., 20., 40., 80., 160., 320., 640., 1280., 2560., 5120.]
    with Pool(8) as p:
        res = p.map(task, [(q, d) for q in qs for d in ds])
    for q, d, best in res:
        s = "none" if best is None else f"r_conj={best[0]:.1f} phi={best[1]:.4f} a0={best[2]:.5f}"
        print(f"q={q} R={Rq(q):.2f} d={d}: {s}", flush=True)
