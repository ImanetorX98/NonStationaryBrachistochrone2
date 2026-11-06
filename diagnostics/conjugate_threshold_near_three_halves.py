"""Floating-point diagnostic, not a certificate: q = 1.4, 1.45 at large radii (expected threshold ~ (3/2 - q)^-3); only
nearly radial incoming rays. Source of the markers of Figure 2 at q = 1.4, 1.45."""
import numpy as np, sys
from multiprocessing import Pool
sys.path.insert(0, '.')
from conjugate_points_static import run
def task(args):
    q, d = args
    for a0 in np.pi - np.geomspace(2e-5, 0.05, 200):
        zs, reach, upi, _ = run(q, d, a0, rmax=60*d, smax=1e10)
        if zs: return q, d, (zs[0][0], zs[0][1], a0)
    return q, d, None
if __name__ == '__main__':
    jobs = [(1.4, d) for d in [1e4, 2e4, 4e4, 8e4, 1.6e5]] + [(1.45, d) for d in [4e4, 1.6e5, 6.4e5]]
    with Pool(8) as p:
        for q, d, b in p.imap(task, jobs):
            print(f"q={q} d={d:.0f}: " + ("none" if b is None else f"r_conj={b[0]:.0f} phi={b[1]:.4f} a0=pi-{np.pi-b[2]:.2e}"), flush=True)
