"""vaidya_step_continuation.py -- floating-point diagnostic, not a certificate: proper-time branch of ingoing Vaidya
with an accretion step (main text, remark on an example with three phases: "for a step of duration 1/2 the continued
shallow extremal crosses the escape price near d = 0.6 and coincides with the frozen minimiser for d >= 14").
Data: E = sqrt2, r0 = rf = 6, delta = pi/2, m(v) = 1 for v <= d, cubic C1 step up to m1 on [d, d+eps], m1 afterwards.
Backward shooting from the target (T, psi), cost C = int f/E dv - (rf - r0)/E, escape price B_V(d) of the theorem on the
escape price. Usage: python3 vaidya_step_continuation.py <m1> <eps> > log   (the main text uses m1 = 1.5, eps = 0.5)"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'tools'))
from progress import Progress
import sys, numpy as np, sympy as sp
from scipy.integrate import solve_ivp, quad
from scipy.optimize import least_squares

E = np.sqrt(2.0); R0 = RF = 6.0; DELTA = np.pi/2
M1 = float(sys.argv[1]) if len(sys.argv) > 1 else 1.5
EPS = float(sys.argv[2]) if len(sys.argv) > 2 else 0.5

r, pr, J, m, e = sp.symbols('r p_r J m e', positive=True)
f = 1 - 2*m/r; w = e**2 - f; om = f/e
S = sp.sqrt(e**2*pr**2 + J**2/r**2); h = -w*pr + sp.sqrt(w)*S
L = lambda ex: sp.lambdify((r, pr, J, m, e), ex, 'numpy')
rdot, phidot = L(sp.diff(h, pr)), L(sp.diff(h, J))
prdot = L(-sp.diff(h - om, r)); dm_hom = L(sp.diff(h - om, m)); hfun = L(h)

def mass(v, d):
    z = (v - d)/EPS
    if z <= 0: return 1.0, 0.0
    if z >= 1: return M1, 0.0
    return 1 + (M1 - 1)*(3*z*z - 2*z**3), (M1 - 1)*(6*z - 6*z*z)/EPS

def rhs(v, y, d):
    rr, ph, p, JJ, q = y
    mm, mp = mass(v, d)
    return [rdot(rr, p, JJ, mm, E), phidot(rr, p, JJ, mm, E), prdot(rr, p, JJ, mm, E), 0.0,
            -mp*dm_hom(rr, p, JJ, mm, E)]

def cost_rate(v, y, d):
    return (1 - 2*mass(v, d)[0]/y[0])/E

def shoot(T, psi, d, dense=False):
    """backwards from v=T (target, q=0, h=omega) to v=0"""
    if not (T > 0 and 0 < psi < np.pi): return None
    mm = mass(T, d)[0]; fT = 1 - 2*mm/RF; wT = E*E - fT; omT = fT/E
    Sn = omT/(np.sqrt(wT) - (wT/E)*np.cos(psi))
    if Sn <= 0: return None
    y0 = [RF, DELTA, Sn*np.cos(psi)/E, Sn*RF*np.sin(psi), 0.0, 0.0]
    def rhs6(v, y):
        dy = rhs(v, y[:5], d)
        return dy + [cost_rate(v, y, d)]
    def hit(v, y): return y[0] - 2*mass(v, d)[0] - 1e-7
    hit.terminal = True
    br = sorted({x for x in (d, d + EPS) if 0 < x < T}, reverse=True)
    pts = [T] + br + [0.0]; y = y0
    for a, b in zip(pts[:-1], pts[1:]):
        sol = solve_ivp(rhs6, (a, b), y, method='DOP853', rtol=1e-11, atol=1e-13, events=hit)
        if sol.status != 0: return None
        y = sol.y[:, -1]
    # cost: integrated backwards, hence C = -y[5]
    return dict(r0=y[0], phi0=y[1], pr0=y[2], J=y[3], q0=y[4], C=-y[5])

def F(x, d):
    s = shoot(x[0], x[1], d)
    if s is None: return [1e3, 1e3]
    return [s['r0'] - R0, s['phi0']]

def solve(guess, d):
    try:
        res = least_squares(F, guess, args=(d,), bounds=([1e-3, 1e-4], [80.0, np.pi - 1e-4]),
                            xtol=1e-14, ftol=1e-14, gtol=1e-14)
    except Exception:
        return None
    x = res.x
    if max(abs(np.array(F(x, d)))) > 1e-9: return None
    return x

def jac(x, d, hstep=1e-6):
    c = np.array(F(x, d)); Jm = np.zeros((2, 2))
    for i in range(2):
        xx = x.copy(); xx[i] += hstep; Jm[:, i] = (np.array(F(xx, d)) - c)/hstep
    return np.linalg.det(Jm)

# ---- escape price B_V(d) (main text, theorem on the escape price)
def sp_(v, b, d):
    ww = E*E - 1 + 2*mass(v, d)[0]/b; return -ww + E*np.sqrt(ww)
def sm_(v, b, d):
    ww = E*E - 1 + 2*mass(v, d)[0]/b; return -ww - E*np.sqrt(ww)
def B_V(d):
    vend = d + EPS
    if vend > 0:
        fr = solve_ivp(lambda v, b: [sp_(v, b[0], d)], (vend, 0.0), [2*M1], method='DOP853',
                       rtol=1e-12, atol=1e-14, dense_output=True)
        bfun = lambda v: fr.sol(v)[0] if v < vend else 2*M1
    else:
        bfun = lambda v: 2*M1
    def ev(v, y): return y[0] - bfun(v)
    ev.terminal = True
    fall = solve_ivp(lambda v, y: [sm_(v, y[0], d), E + np.sqrt(E*E - 1 + 2*mass(v, d)[0]/y[0])],
                     (0.0, max(vend, 0) + 50), [R0, 0.0], method='DOP853', rtol=1e-12, atol=1e-14, events=ev)
    vc = fall.t_events[0][0]; Dc = fall.y_events[0][0][1]
    front = 0.0
    if vc < vend:
        front = quad(lambda v: E - np.sqrt(E*E - 1 + 2*mass(v, d)[0]/bfun(v)), vc, vend,
                     points=[x for x in (d,) if vc < x < vend], limit=400, epsabs=1e-13)[0]
    Rinf = quad(lambda s: 1/np.sqrt(E*E - 1 + 2*M1/s), 2*M1, RF, epsabs=1e-13)[0]
    return Dc + front + Rinf, vc, bfun(0.0)

if __name__ == '__main__':
    print(f'# m1={M1} eps={EPS} E=sqrt2 r0=rf=6 delta=pi/2')
    # 1. validation at constant mass m=1 (very large d): frozen reference 5.8130366021, L=2.6258028, T=13.963779
    best = None
    for T in np.linspace(8, 20, 7):
        for psi in np.linspace(0.3, 2.8, 6):
            x = solve(np.array([T, psi]), 1e6)
            if x is not None:
                s = shoot(x[0], x[1], 1e6)
                if best is None or s['C'] < best[1]['C'] - 1e-9: best = (x, s)
    x0, s0 = best
    print(f'VALIDATION m=1: T={x0[0]:.9f} psi={x0[1]:.9f} C={s0["C"]:.10f} J={s0["J"]:.7f}  (ref. 5.8130366021, 2.6258028, 13.963779)')
    # also constant mass m1 (large negative d)
    bestm1 = None
    for T in np.linspace(6, 30, 9):
        for psi in np.linspace(0.3, 2.8, 6):
            x = solve(np.array([T, psi]), -1e6)
            if x is not None:
                s = shoot(x[0], x[1], -1e6)
                if bestm1 is None or s['C'] < bestm1[1]['C'] - 1e-9: bestm1 = (x, s)
    if bestm1:
        print(f'FROZEN m1: T={bestm1[0][0]:.6f} C={bestm1[1]["C"]:.8f}  B_V(-inf)={B_V(-1e6)[0]:.8f}')
    else:
        print('FROZEN m1: no interior solution found')
    print(f'B_V(+inf) (m=1 always before): {B_V(60.0)[0]:.8f}   2R_1(6)={2*quad(lambda s: 1/np.sqrt(1+2/s),2,6)[0]:.8f}')
    # 2. continuation in d, from large values towards 0 and beyond
    print('# d  T  psi  C  B_V  C-B_V  J  q0  det  vc  b(0)')
    x = x0.copy(); ds = list(np.arange(20, -0.001, -0.25))
    pg = Progress(len(ds), 'step continuation: d')
    for d in ds:
        xn = solve(x, d)
        if xn is None:
            print(f'{d:7.3f}  LOST (continuation failed)'); break
        x = xn; s = shoot(x[0], x[1], d)
        B, vc, b0 = B_V(d)
        print(f'{d:7.3f} {x[0]:10.6f} {x[1]:9.6f} {s["C"]:12.8f} {B:12.8f} {s["C"]-B:+.3e} {s["J"]:9.6f} {s["q0"]:+.3e} {jac(x,d):+.4e} {vc:8.4f} {b0:8.5f}')
        sys.stdout.flush()
        pg.step()
    pg.close()
