"""vaidya_step_branches.py -- floating-point diagnostic, not a certificate: refinement and search for other branches
for the accretion step of vaidya_step_continuation.py (main text, remark on an example with three phases: "a second,
deep family of extremals stays above the escape price"). Usage: python3 vaidya_step_branches.py <m1> <eps>"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'tools'))
from progress import Progress
import sys, numpy as np
sys.argv = [sys.argv[0]] + sys.argv[1:3]
import vaidya_step_continuation as G

def multistart(d, Ts, psis):
    sols = []
    for T in Ts:
        for psi in psis:
            x = G.solve(np.array([T, psi]), d)
            if x is None: continue
            if any(abs(x[0]-y[0]) < 1e-6 and abs(x[1]-y[1]) < 1e-6 for y, _ in sols): continue
            s = G.shoot(x[0], x[1], d); sols.append((x, s))
    return sorted(sols, key=lambda t: t[1]['C'])

print(f'# m1={G.M1} eps={G.EPS}')
# (a) multi-start on a grid of d
pg = Progress(14, 'step branches (a): d')
for d in [-50.0, 0.0, 0.25, 0.5, 0.75, 1.0, 2.0, 5.0, 10.0, 13.0, 13.5, 13.75, 14.0, 20.0]:
    B = G.B_V(d)[0]
    sols = multistart(d, np.linspace(4, 40, 13), np.linspace(0.05, 3.0, 9))
    print(f'd={d:7.3f} B_V={B:.8f} solutions={len(sols)}')
    for x, s in sols:
        print(f'     T={x[0]:10.6f} psi={x[1]:9.6f} C={s["C"]:.8f} C-B={s["C"]-B:+.3e} J={s["J"]:.6f} q0={s["q0"]:+.4e} det={G.jac(x,d):+.3e}')
    sys.stdout.flush()
    pg.step()
pg.close()
# (b) fine continuation in the transition zone, from d=14.2 to 13.2
x = multistart(14.2, [14.0], [0.22])[0][0]
print('# fine continuation d in [13.2,14.2]')
pg = Progress(len(np.arange(14.2, 13.19, -0.02)), 'step branches (b): d')
for d in np.arange(14.2, 13.19, -0.02):
    xn = G.solve(x, d)
    if xn is None: print(f'{d:.3f} LOST'); break
    x = xn; s = G.shoot(x[0], x[1], d)
    print(f'{d:7.3f} T={x[0]:.6f} psi={x[1]:.6f} C={s["C"]:.8f} J={s["J"]:.6f} q0={s["q0"]:+.4e} det={G.jac(x,d):+.3e}')
    sys.stdout.flush()
    pg.step()
pg.close()
