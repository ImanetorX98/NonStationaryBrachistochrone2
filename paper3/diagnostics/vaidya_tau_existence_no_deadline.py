"""vaidya_tau_existence_no_deadline.py -- floating-point diagnostic, not a certificate: existence without a deadline
below the escape price (main text, theorem on existence with an unbounded arrival time, ingoing Vaidya).
Example: m(v) = 1 + 0.5 tanh(v/5) (m_inf = 1.5), E = sqrt 2, r0 = r_f = 6, launch at v0 = 0.
Lower bound of the escape price: lower = |R_inf(r0)| + R_inf(r_f) = 4.62260 (vaidya_tau_escape_price_bounds.py).
Every rail integrated as an ODE that reaches the target gives an upper bound m*(delta) <= cost; the theorem gives an
attained minimum for every delta with m*(delta) < lower, and by monotonicity of m* on the whole [0, delta].
 (a) homotopy in eps: m_eps = 1 + 0.5 eps tanh(v/5), from eps = 0 (frozen, 5.8130366 at pi/2) to eps = 1;
 (b) continuation in delta of the shallow branch at eps = 1: cost, minimum gap, duration;
 (c) delta_V^- = last angle with cost < lower (numerical lower bound of the existence threshold);
 (d) control that must fail: at constant mass m = 1 the lower bound is 2R(6) = 6.45722 = B exactly; the branch
     continued beyond delta_c must NOT go below B (otherwise the frozen theory would be false).
Importable: vaidya_tau_calibration_identity.py uses mass, solve, info and frozen_start."""
import numpy as np, sympy as sp
from scipy.integrate import solve_ivp, quad
from scipy.optimize import fsolve
r, pr, L, q, E, m, dm = sp.symbols('r p_r L q E m dm', real=True)
f = 1 - 2*m/r; w = E**2 - f; S = sp.sqrt(E**2*pr**2 + L**2/r**2)
h = -w*pr + sp.sqrt(w)*S; om = f/E
args = (r, pr, L, q, E, m, dm)
F = {k: sp.lambdify(args, ex, 'numpy') for k, ex in dict(
    rdot=sp.diff(h, pr), phidot=sp.diff(h, L), prdot=-sp.diff(h - om, r), qdot=-sp.diff(h - om, m)*dm,
    om=om, h=h).items()}
Ev = np.sqrt(2.0); r0 = rf = 6.0
def mass(eps): return (lambda vv: 1 + 0.5*eps*np.tanh(vv/5), lambda vv: 0.1*eps/np.cosh(vv/5)**2)
def flow(x, mf, dense=False):
    mfun, dmfun = mf; pr0, Lv, T = x
    q0 = F['om'](r0, pr0, Lv, 0, Ev, mfun(0.0), 0) - F['h'](r0, pr0, Lv, 0, Ev, mfun(0.0), 0)
    def rhs(vv, y):
        rr, ph, p, qq, c, g = y; mm, d = mfun(vv), dmfun(vv)
        return [F['rdot'](rr, p, Lv, qq, Ev, mm, d), F['phidot'](rr, p, Lv, qq, Ev, mm, d),
                F['prdot'](rr, p, Lv, qq, Ev, mm, d), F['qdot'](rr, p, Lv, qq, Ev, mm, d),
                F['om'](rr, p, Lv, qq, Ev, mm, d), 0.0]
    s = solve_ivp(rhs, (0, T), [r0, 0.0, pr0, q0, 0.0, 0.0], rtol=1e-11, atol=1e-13, dense_output=True)
    return s, q0
def resid(x, mf, delta):
    s, _ = flow(x, mf); return [s.y[0, -1] - rf, s.y[1, -1] - delta, s.y[3, -1]]
def solve(x0, mf, delta):
    x, info, ier, msg = fsolve(resid, x0, args=(mf, delta), full_output=True, xtol=1e-13)
    ok = np.abs(info["fvec"]).max() < 1e-9 and x[2] > 0   # residual test (ier may report slow progress at xtol=1e-13)
    if ok:
        s_, _ = flow(x, mf); ok = s_.success
    return x, ok
def info(x, mf):
    s, q0 = flow(x, mf); vv = np.linspace(0, x[2], 4001); rr = s.sol(vv)[0]
    gap = np.min(rr - 2*mf[0](vv))
    assert s.success and x[2] > 0 and gap > 0, "integration failure, non-positive duration or horizon crossing (sampled gap)"
    return s.y[4, -1], gap, rr.min(), q0
minf = 1.5; hinf = 3.0
Rinf = lambda rr, M=minf: quad(lambda s: 1/np.sqrt(Ev**2 - 1 + 2*M/s), 2*M, rr)[0]
lower = abs(Rinf(r0)) + Rinf(rf)
if __name__ == "__main__":
    print(f"lower bound of escape price: {lower:.6f}")
from scipy.optimize import brentq
def frozen_start():   # inward launch with h = omega at r0, L = 2.6258
    L0 = 2.6258; gg = lambda p: F['h'](r0, p, L0, 0, Ev, 1.0, 0) - F['om'](r0, p, L0, 0, Ev, 1.0, 0)
    ps = np.linspace(-5, 5, 2001); vals = [gg(p) for p in ps]
    roots = [brentq(gg, a, b) for a, b, va, vb in zip(ps, ps[1:], vals, vals[1:]) if va*vb < 0]
    best = None
    for p0 in [p for p in roots if F['rdot'](r0, p, L0, 0, Ev, 1.0, 0) < 0]:
        for T0 in (10.0, 15.0, 20.0):
            xx, ok = solve(np.array([p0, L0, T0]), mass(0.0), np.pi/2)
            if ok:
                c = flow(xx, mass(0.0))[0].y[4, -1]
                if best is None or c < best[1]: best = (xx, c)
    return best[0]
if __name__ == "__main__":
    x0frozen = frozen_start(); x = x0frozen.copy(); ok = True
    c, g, rmin, q0 = info(x, mass(0.0))
    print(f"(a) eps=0: ok={ok}, cost {c:.10f} (frozen 5.8130366021)")
    for eps in np.linspace(0.05, 1.0, 20):
        x, ok = solve(x, mass(eps), np.pi/2); assert ok, f"homotopy failed at eps={eps}"
    c, g, rmin, q0 = info(x, mass(1.0))
    print(f"(a) eps=1, delta=pi/2: cost {c:.8f}, min gap {g:.4f}, r_min {rmin:.4f}, q_in {q0:+.4e}, T_v {x[2]:.4f}")
    mf = mass(1.0); xs = {0.5: x.copy()}
    def cont(direction, stop):
        xx = xs[0.5].copy(); d = 0.5; out = []
        while True:
            dn = d + direction*0.01
            if (direction > 0 and dn > stop + 1e-12) or (direction < 0 and dn < stop - 1e-12): break
            xn, ok = solve(xx, mf, dn*np.pi)
            if not ok:
                ok2 = False
                for sub in (0.005, 0.0025):
                    xn, ok2 = solve(xx, mf, (d + direction*sub)*np.pi)
                    if ok2: dn = d + direction*sub; break
                if not ok2: print(f"    continuation stops after delta = {d:.4f} pi"); break
            xx, d = xn, dn; out.append((d, *info(xx, mf), xx[2]))
        return out
    res = sorted(cont(-1, 0.05) + [(0.5, *info(xs[0.5], mf), xs[0.5][2])] + cont(+1, 0.99))
    print("(b) shallow branch, m = 1 + 0.5 tanh(v/5):  delta/pi  cost  min gap  r_min  q_in  T_v")
    for d, c, g, rmin, q0, T in res:
        if abs(round(d*20) - d*20) < 1e-9 or c > lower - 0.05:
            print(f"    {d:.4f}  {c:.6f}  {g:.4f}  {rmin:.4f}  {q0:+.3e}  {T:8.3f}  {'< lower' if c < lower else '>= lower'}")
    below = [d for d, c, *_ in res if c < lower]
    mono = all(b[1] >= a[1] - 1e-9 for a, b in zip(res, res[1:]))
    print(f"(c) cost nondecreasing along the branch: {mono}; last delta with cost < lower: {max(below):.4f} pi  => "
          f"minimum attained (no deadline) at least for delta in [0, {max(below):.4f} pi]  (not certified)")
    # (d) frozen control: m = 1, B = 2R(6) exact; the criterion m* < B must reproduce the certified delta_c of
    #     q = 2, r0 = 6M, [0.588406, 0.588417] pi (../interval_arithmetic/critical_separation.py), and must FAIL beyond it
    #     (non-attainment, frozen theory).
    mf1 = mass(0.0); B1 = 2*Rinf(r0, 1.0)
    xx = x0frozen.copy(); d = 0.5; prev = (d, info(xx, mf1)[0], xx.copy())
    while True:
        xn, ok = solve(xx, mf1, (d + 0.01)*np.pi)
        assert ok, f"frozen continuation failed at {d + 0.01}"
        xx, d = xn, d + 0.01; c = info(xx, mf1)[0]
        if c >= B1: break
        prev = (d, c, xx.copy())
    lo, hi, xlo = prev[0], d, prev[2]
    for _ in range(40):   # bisection on delta for cost = B, warm start from the lower side
        mid = (lo + hi)/2; xm, ok = solve(xlo, mf1, mid*np.pi); assert ok
        if info(xm, mf1)[0] < B1: lo, xlo = mid, xm
        else: hi = mid
    print(f"(d) frozen m = 1: B = 2R(6) = {B1:.6f}; shallow cost = B at delta = {lo:.6f} pi "
          f"(certified delta_c in [0.588406, 0.588417] pi): {'PASS' if 0.588406 <= lo <= 0.588417 else 'FAIL'}")
    assert 0.588406 <= lo <= 0.588417, "frozen control"
    xa, ok = solve(xlo, mf1, (lo + 0.005)*np.pi); ca = info(xa, mf1)[0]; assert ok and ca >= B1, "must-fail control"
    print(f"    must fail beyond: delta = {lo + 0.005:.4f} pi, cost {ca:.6f} >= B: {ca >= B1}  (criterion gives nothing, as required)")
