"""tau_var_value_function.py -- floating-point diagnostic, not a certificate: global value function W(eta, r, b) of
the reduced problem (a = 0; main text, remark on numerical observations, semi-Lagrangian dynamic programme),
b = remaining angular budget in [0, delta] (b <- max(b, 0): the cap makes the excess angle useless).
  r' = v u,  b' = -w0 sqrt(1 - u^2),  |u| <= 1,  cost omega;  target r = r_f, b = 0; free final time <= T_max;
  state constraint r in [r1, r2] (a domain containing every competitor of cost <= C_cand, by the barriers; the
  domains are computed by tau_var_numbers.py).
Backward semi-Lagrangian dynamic programming, bilinear interpolation (scipy map_coordinates).
W(0, r0, delta) is a NUMERICAL ESTIMATE of the value (neither a proved upper nor a proved lower bound), compared with
the candidate extremal. Control: constant A (shallow dip 5.81303660, global minimum by the theorem on a single maximum).
Bilinear interpolation and the Euler step give errors of uncontrolled sign; W_h is NOT a proved upper bound of the
continuous value; convergence with a point target, state constraint and finite INF is not proved. The effective domain
(r1, r2 aligned with the grid) is printed in the log.
Usage: python3 tau_var_value_function.py MODE [dr,db,de,nu] [C_cand] [r1,r2,Tmax]
  MODE = frozen (constant A), <eta_d> (decreasing A, e.g. 4) or up<eta_d> (increasing A, e.g. up4);
  environment variable DP_INF (default 1e6) sets the value used outside the domain."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'tools'))
from progress import Progress
import numpy as np, sys, time
from scipy.ndimage import map_coordinates

def value(Afun, E, r0, rf, delta, r1, r2, Tmax, dr, db, de, nu):
    f = lambda r: 1 - 2/r
    r1 = rf - dr*np.ceil((rf - r1)/dr - 1e-9)        # align the grid on r_f (domain enlarged downwards: conservative)
    r2 = rf + dr*np.ceil((r2 - rf)/dr - 1e-9)        # align r2 too (domain enlarged upwards); last node = r2 exactly
    rs = rf + dr*np.arange(-int(round((rf - r1)/dr)), int(round((r2 - rf)/dr)) + 1); bs = np.arange(0, delta + db/2, db)
    R, B = rs.size, bs.size
    jf = int(round((rf - r1)/dr)); assert abs(rs[jf] - rf) < 1e-9
    import os; INF = float(os.environ.get('DP_INF', '1e6'))
    W = np.full((R, B), INF); W[jf, 0] = 0.0
    us = np.linspace(-1, 1, nu)
    rr, bb = np.meshgrid(rs, bs, indexing='ij')
    pg = Progress(int(round(Tmax/de)), 'value function: backward time steps')
    for k in range(int(round(Tmax/de)), 0, -1):
        e = (k - 0.5)*de; a = Afun(e)
        N = 1/(1 - a**2*f(rr)/E**2)
        v = f(rr)/np.sqrt(N); w0 = np.sqrt(f(rr))/(rr*np.sqrt(N)); om = a**2*f(rr)/E*de
        best = np.full((R, B), INF)
        for u in us:
            rn = rr + v*u*de; bn = np.maximum(bb - w0*np.sqrt(1 - u*u)*de, 0.0)
            ir = (rn - r1)/dr; ib = bn/db
            val = map_coordinates(W, [ir, ib], order=1, mode='nearest')
            val[(rn < r1) | (rn > r2)] = INF
            best = np.minimum(best, om + val)
        W = np.minimum(best, INF); W[jf, 0] = 0.0
        pg.step()
    pg.close()
    j0 = int(round((r0 - r1)/dr)); kd = int(round(delta/db))
    return W[j0, kd], (rs[0], rs[-1])

if __name__ == "__main__":
    E = np.sqrt(2.0); r0 = rf = 6.0; delta = np.pi/2
    mode = sys.argv[1]
    dr, db, de, nu = [float(x) for x in sys.argv[2].split(',')] if len(sys.argv) > 2 else (0.02, 0.01, 0.02, 21)
    nu = int(nu)
    if mode == 'frozen':
        Afun = lambda e: 1.0; Cc = 5.81303660; r1, r2, Tmax = 2.44, 9.28, 46.0
    elif mode.startswith('up'):
        ed = float(mode[2:]); n = 1.15 + 0.15*np.tanh(-ed)
        Afun = lambda e: (1.15 + 0.15*np.tanh(e - ed))/n
        Cc = float(sys.argv[3]); r1, r2, Tmax = [float(x) for x in sys.argv[4].split(',')]
    else:
        ed = float(mode); n = 0.75 + 0.25*np.tanh(ed)
        Afun = lambda e: (0.75 + 0.25*np.tanh(-(e - ed)))/n
        Cc = float(sys.argv[3]); r1, r2, Tmax = [float(x) for x in sys.argv[4].split(',')]
    db = delta/round(delta/db)                 # make delta a grid point
    t0 = time.time()
    Wv, dom = value(Afun, E, r0, rf, delta, r1, r2, Tmax, dr, db, de, nu)
    print(f"mode {mode}: grid dr={dr}, db={db:.5f}, de={de}, nu={nu}; K input [0,{Tmax}] x [{r1},{r2}], effective r in [{dom[0]:.4f},{dom[1]:.4f}]; "
          f"W(0, r0, delta) = {Wv:.6f}; candidate {Cc:.6f}; candidate - W = {Cc - Wv:+.5f}; {time.time()-t0:.0f}s", flush=True)
