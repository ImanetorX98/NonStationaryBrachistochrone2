"""tau_var_late_drop.py -- floating-point diagnostic (SLSQP, reduced problem, local optimum, not a certificate):
direct optimisation for the late drop eta_d = 15 (A: 1 -> 0.5), q = 2, r0 = r_f = 6, delta = pi/2, started from the
explicit competitor (down to 2.3M, angular arc until eta = 24, up; cost 4.3427). The main text reports that the
direct optimisation drifts towards r = 2M with cost approaching B = 3.938 (tau_var_asymmetric_dive.py)."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'tools'))
from progress import Progress
import numpy as np, sys
from tau_var_direct_optimisation import solve
E = np.sqrt(2.0); Afun = lambda e: 0.75 + 0.25*np.tanh(-(e - 15))
pg = Progress(3, 'late drop: starts')
for n, T, (e1, e2) in ((160, 33.12, (10.02, 24.0)), (160, 28.0, (8.0, 20.0)), (200, 38.0, (11.0, 28.0))):
    ee = (np.arange(n) + 0.5)*T/n
    psi = np.where(ee < e1, np.pi*0.98, np.where(ee < e2, np.pi/2, 0.02))
    res, rs, ang, cost = solve(Afun, E, 6.0, 6.0, np.pi/2, n=n, x0=np.concatenate([[T], psi]))
    print(f"n={n}, start T={T}: success {res.success} ({res.message}); cost {cost:.6f}, T={res.x[0]:.3f}, "
          f"min r {rs.min():.4f}, angle {ang/np.pi:.3f}pi, time at r<2.5: {np.mean(rs[:-1] < 2.5)*res.x[0]:.2f}", flush=True)
    pg.step()
pg.close()
