"""tau_var_direct_optimisation.py -- floating-point diagnostic (scipy SLSQP, LOCAL optimum, not a certificate):
direct optimisation of the reduced radial problem of the proper-time branch (a = 0, any A(eta)).

Reduction: the cost tau = int omega(eta, r) d eta depends only on r(eta); a rail has
  dr/deta = (f/sqrt N) cos psi,  |dOmega/deta| = (sqrt f/(r sqrt N)) sin psi,  psi in [0, pi],
and the reachable angle is every delta <= int (sqrt f/(r sqrt N)) sin psi d eta (the excess is spent going back and
forth). Variables: T and psi_k (n intervals, RK2 in the radius); constraints r(T) = r_f, angular budget >= delta,
r >= 2 + 1e-4. M = 1. Normalisation of critical_separation_estimate.py: Ehat = sqrt(q0), A = 1 at launch.
Control: constant A, q = 2, r0 = r_f = 6: the value must reproduce the length of the shallow dip of angle delta
(legs() of critical_separation_estimate.py), the global minimum for delta < delta_c (main text, theorem on a single
maximum).
Usage: python3 tau_var_direct_optimisation.py [eta_d list, decreasing A] [eta_d list, increasing A]
(defaults 4,8,15 and 4,8, used for the log in ../results)."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'tools'))
from progress import Progress
import numpy as np, mpmath as mp
from scipy.optimize import minimize
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'interval_arithmetic'))
from critical_separation_estimate import legs

def solve(Afun, E, r0, rf, delta, n=120, x0=None, T0=None):
    f = lambda r: 1 - 2/r
    N = lambda e, r: 1/(1 - Afun(e)**2*f(r)/E**2)
    def traj(x):
        T, psi = x[0], x[1:]; h = T/n; r = r0; rs = [r]; ang = 0.0; cost = 0.0
        for k in range(n):
            e = k*h; c, s = np.cos(psi[k]), np.sin(psi[k])
            v = lambda e, r: f(r)/np.sqrt(N(e, r))*c
            rm = r + 0.5*h*v(e, r); em = e + 0.5*h
            r = r + h*v(em, rm); rs.append(r)
            ang += h*np.sqrt(f(rm))/(rm*np.sqrt(N(em, rm)))*s
            cost += h*Afun(em)**2*f(rm)/E
        return np.array(rs), ang, cost
    cons = [{'type': 'eq', 'fun': lambda x: traj(x)[0][-1] - rf},
            {'type': 'ineq', 'fun': lambda x: traj(x)[1] - delta},
            {'type': 'ineq', 'fun': lambda x: traj(x)[0] - 2.0001}]
    if x0 is None:
        k = np.arange(n); x0 = np.concatenate([[T0], np.where(k < n/2, np.pi*0.7, np.pi*0.3)])
    bnds = [(1e-3, 500)] + [(0, np.pi)]*n
    res = minimize(lambda x: traj(x)[2], x0, method='SLSQP', bounds=bnds, constraints=cons,
                   options=dict(maxiter=600, ftol=1e-11))
    rs, ang, cost = traj(res.x)
    return res, rs, ang, cost

def dip_length(q, r0, delta):
    """shallow dip of angle delta (normalisation A = 1, Ehat = sqrt q): find t on the shallow branch"""
    ts = [2 + (r0 - 2)*k/200 for k in range(1, 200)]
    vals = [legs(q, r0, t)[0] for t in ts]
    for a, b, va, vb in zip(ts, ts[1:], vals, vals[1:]):
        if (va - delta)*(vb - delta) < 0 and vb < va:            # decreasing (shallow) branch
            t = mp.findroot(lambda t: legs(q, r0, t)[0] - delta, (a, b), solver='anderson')
            return float(legs(q, r0, t)[1]), float(t)

if __name__ == "__main__":
    q0, r0 = 2.0, 6.0; E = np.sqrt(q0)
    print("Control: A constant, q = 2, r0 = r_f = 6 (must reproduce the shallow dip, the global minimum)")
    sols = {}
    pg = Progress(2, 'direct optimisation: constant A')
    for dl in (np.pi/4, 0.5*np.pi):
        ell, t = dip_length(q0, r0, dl)
        best = None
        for T0 in (6.0, 10.0, 16.0):
            res, rs, ang, cost = solve(lambda e: 1.0, E, r0, r0, dl, T0=T0)
            if res.success and (best is None or cost < best[3]): best = (res, rs, ang, cost)
        res, rs, ang, cost = best; sols[dl] = res.x
        print(f"  delta={dl/np.pi:.3f}pi: reduced cost {cost:.6f}, dip length {ell:.6f} (periastron {t:.4f}), "
              f"rel. diff {(cost-ell)/ell:+.2e}; min r {rs.min():.4f} (dip periastron {t:.4f}); T={res.x[0]:.3f}")
        pg.step()
    pg.close()
    dl = 0.5*np.pi
    print("Decreasing A: A = 1 until eta_d, then 0.5 (tanh, width 1); delta = 0.5 pi, r0 = r_f = 6")
    import sys
    eds = (4.0, 8.0, 15.0) if len(sys.argv) < 2 else tuple(float(v) for v in sys.argv[1].split(','))
    pg = Progress(len(eds), 'direct optimisation: decreasing A')
    for ed in eds:
        Afun = lambda e, ed=ed: 0.75 + 0.25*np.tanh(-(e - ed))
        best = None
        for start in (sols[dl], None):
            for T0 in (10.0, 25.0, 50.0):
                x0 = None if start is None else np.concatenate([[T0], start[1:]])
                res, rs, ang, cost = solve(Afun, E, r0, r0, dl, x0=x0, T0=T0)
                if res.success and (best is None or cost < best[3]): best = (res, rs, ang, cost)
        if best is None: print(f"  eta_d={ed:5.1f}: no successful SLSQP run"); continue
        res, rs, ang, cost = best
        print(f"  eta_d={ed:5.1f}: cost {cost:.6f}, T={res.x[0]:.3f}, min r {rs.min():.4f}, angle budget {ang/np.pi:.3f}pi "
              f"(excess {(ang-dl)/np.pi:+.3f}pi), time spent at r < 2.5: {np.mean(rs < 2.5)*res.x[0]:.2f}")
        pg.step()
    pg.close()
    print("Increasing A: A = 1 until eta_d, then 1.3 (tanh), delta = 0.5 pi")
    pg = Progress(len((4.0, 8.0) if len(sys.argv) < 3 else tuple(float(v) for v in sys.argv[2].split(','))), 'direct optimisation: increasing A')
    for ed in ((4.0, 8.0) if len(sys.argv) < 3 else tuple(float(v) for v in sys.argv[2].split(','))):
        Afun = lambda e, ed=ed: 1.15 + 0.15*np.tanh(e - ed)
        best = None
        for T0 in (6.0, 10.0, 16.0):
            res, rs, ang, cost = solve(Afun, E, r0, r0, dl, x0=np.concatenate([[T0], sols[dl][1:]]))
            if res.success and (best is None or cost < best[3]): best = (res, rs, ang, cost)
        if best is None: print(f"  eta_d={ed:5.1f}: no successful SLSQP run"); continue
        res, rs, ang, cost = best
        print(f"  eta_d={ed:5.1f}: cost {cost:.6f}, T={res.x[0]:.3f}, min r {rs.min():.4f}, excess {(ang-dl)/np.pi:+.3f}pi")
        pg.step()
    pg.close()
