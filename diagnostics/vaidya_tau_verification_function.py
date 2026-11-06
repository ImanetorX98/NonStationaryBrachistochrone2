"""vaidya_tau_verification_function.py -- floating-point diagnostic, not a certificate: exact escape price of ingoing
Vaidya through a verification function (calibration; main text, theorem on the escape price: B_V = 5.21295).
Radial problem (projection of K): r' = s in [s-, s+], s+- = -w +- E sqrt(w), cost (f - s)/E.
Field: from the point (v, r), maximal infall s- until contact with b_inf at (v_c, r_c), then chase along b_inf
(s = s+) up to v = inf:  W(v, r) = int_v^{v_c} (E + sqrt w) dv' + Ch(v_c),  Ch(v_c) = int_{v_c}^inf (E - sqrt w)|_{b_inf} dv'.
Verification condition (linear in s, s = s+ suffices): sigma = dW/dr - 1/E >= 0 on the region above b_inf.
If it holds: B_V = W(v0, r0) + R_inf(r_f) for every delta.
 (a) tanh example (m = 1 + 0.5 tanh(v/5), E = sqrt 2, r0 = r_f = 6): W(0, 6), B_V, comparison with
     vaidya_tau_escape_price_bounds.py;
 (b) sigma on a grid of the region;
 (c) families with more switches (rise-fall, fall-rise-fall, wait-fall): none may beat W(0, 6);
 (d) constant mass: W = R(r) - R(2m), sigma = 1/sqrt(w) - 1/E > 0 exactly (control);
 (e) rapid accretion: search for sigma < 0 and, if found, for a curve with more switches that beats the field.
Importable: vaidya_tau_calibration_identity.py uses E, setup and Rinf."""
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import minimize
E = np.sqrt(2.0)
def setup(mfun, minf, Tb=400.0):
    """Frontier bracket: backward from 2m(Tb) gives b_lo <= b_inf, backward from r_f-like level
    above h_inf gives b_hi >= b_inf; the field uses b_lo, the bracket width measures the frontier error. The chase
    integral J(u) = int_u^Tb c_+ is integrated together with b_lo (smooth dense output, no grid interpolation); the
    true tail beyond Tb lies in [0, (h_inf - b_inf(Tb))/c] <= [0, (h_inf - 2m(Tb))/c] and is reported, not added."""
    f = lambda v, r: 1 - 2*mfun(v)/r; w = lambda v, r: E**2 - f(v, r)
    sp_ = lambda v, r: -w(v, r) + E*np.sqrt(w(v, r)); sm = lambda v, r: -w(v, r) - E*np.sqrt(w(v, r))
    lo = solve_ivp(lambda v, y: [sp_(v, y[0]), -(E - np.sqrt(w(v, y[0])))], (Tb, 0.0), [2*mfun(Tb), 0.0], rtol=1e-12,
                   atol=1e-14, dense_output=True, max_step=0.02)
    hi = solve_ivp(lambda v, y: [sp_(v, y[0])], (Tb, 0.0), [2*minf + 3.0], rtol=1e-12, atol=1e-14, dense_output=True,
                   max_step=0.02)
    assert lo.success and hi.success, "frontier integration failed"
    binf = lambda v: lo.sol(v)[0]; bhi = lambda v: hi.sol(v)[0]
    c = np.sqrt(E**2 - 1); tail_max = (2*minf - 2*mfun(Tb))/c
    assert tail_max >= 0, "h_inf must dominate 2m(Tb)"
    Ch = lambda v: lo.sol(v)[1]          # = int_v^Tb c_+ along b_lo; true tail in [0, tail_max] not included
    def fall(v, r, extra=False):
        ev = lambda t, y: y[0] - binf(t); ev.terminal = True; ev.direction = -1
        def rhs(t, y):
            rho, C, xi, I2 = y; mm = mfun(t); ww = w(t, rho)
            return [sm(t, rho), E + np.sqrt(ww), 2*mm/rho**2*(1 + E/(2*np.sqrt(ww)))*xi, 2*mm/rho**2*xi]
        s = solve_ivp(rhs, (v, v + 200), [r, 0.0, 1.0, 0.0], events=ev, rtol=1e-12, atol=1e-14, max_step=0.02)
        assert s.success and s.t_events[0].size == 1 and s.t[-1] >= v, f"contact not found from ({v}, {r})"
        return s
    def W(v, r):
        s = fall(v, r); return s.y[1, -1] + Ch(s.t[-1]), s.t[-1], s.y[0, -1]
    def sigma(v, r):   # analytic sensitivity: (1/E) int (2m/rho^2) xi dt
        s = fall(v, r); return s.y[3, -1]/E
    return dict(f=f, w=w, sp=sp_, sm=sm, binf=binf, bhi=bhi, Ch=Ch, W=W, sigma=sigma, minf=minf, tail=tail_max,
                mfun=mfun, Tb=Tb)
def Rinf(r, minf): return quad(lambda s: 1/np.sqrt(E**2 - 1 + 2*minf/s), 2*minf, r)[0]
def sigma_grid(P, vs, rs_above, hr=1e-4):
    out = []
    for v in vs:
        b = P['binf'](v)
        for dr in rs_above:
            r = b + dr
            fd = (P['W'](v, r + hr)[0] - P['W'](v, r - hr)[0])/(2*hr) - 1/E if dr > 2*hr else np.nan
            out.append((v, r, P['sigma'](v, r), fd))
    return out
def multi(P, v0, r0, seq):
    """seq: list of (mode, duration) with mode in 'rise','fall'; then final fall to contact with b_inf, then chase."""
    v, r, c = v0, r0, 0.0
    for mode, d in seq:
        fn = P['sp'] if mode == 'rise' else P['sm']
        dens = (lambda t, y: E - np.sqrt(P['w'](t, y[0]))) if mode == 'rise' else (lambda t, y: E + np.sqrt(P['w'](t, y[0])))
        ev = lambda t, y: y[0] - P['binf'](t) - 1e-9; ev.terminal = True
        s = solve_ivp(lambda t, y: [fn(t, y[0]), dens(t, y)], (v, v + max(d, 0)), [r, 0.0], events=ev, rtol=1e-11,
                      atol=1e-13, max_step=0.02)
        if s.t_events[0].size: return 1e6   # crossed the frontier: no return (finite penalty, avoids NaN in Nelder-Mead)
        v, r, c = s.t[-1], s.y[0, -1], c + s.y[1, -1]
    Wv = P['W'](v, r)[0]
    return c + Wv
def run(name, mfun, minf, r0, rf, test_multi=True):
    P = setup(mfun, minf)
    W0, vc, rc = P['W'](0.0, r0); Ri = Rinf(rf, minf)
    print(f"[{name}] b_inf(0) = {P['binf'](0.0):.6f}, 2m(0) = {2*mfun(0.0):.4f}, h_inf = {2*minf}; contact v_c = {vc:.4f}, "
          f"r_c = {rc:.6f}; Ch(v_c) = {P['Ch'](vc):.6f} (tail {P['tail']:.1e}); W(0, r0) = {W0:.6f}; R_inf(r_f) = {Ri:.6f}; "
          f"B_V (if sigma >= 0) = {W0 + Ri:.6f}")
    vchk = np.linspace(0, 30, 7)
    width = max(P['bhi'](v) - P['binf'](v) for v in vchk)
    print(f"    frontier bracket b_lo <= b_inf <= b_hi on v in [0, 30]: max width {width:.1e}; tail beyond Tb in [0, {P['tail']:.1e}]")
    assert width < 1e-8
    sg = sigma_grid(P, np.linspace(0, 30, 16), np.concatenate([np.geomspace(1e-3, 1, 8), np.linspace(1.5, 8, 8)]))
    an = np.array([x[2] for x in sg]); fd = np.array([x[3] for x in sg]); ok = np.isfinite(fd)
    imin = int(np.argmin(an))
    print(f"    sigma (analytic sensitivity) on {len(an)} points: min {an.min():.3e} at v = {sg[imin][0]:.2f}, r = {sg[imin][1]:.4f}; "
          f"negative: {(an <= 0).sum()};  FD estimate vs analytic: max |diff| {np.max(np.abs(fd[ok] - an[ok])):.1e}, "
          f"FD negative points: {(fd[ok] < 0).sum()} (FD is a secondary estimate)")
    assert (an > 0).all(), "analytic sigma must be positive"
    if test_multi:
        best = (W0, 'single switch (field)')
        for lab, form in (('rise-fall', lambda z: [('rise', z[0])]), ('fall-rise', lambda z: [('fall', z[0]), ('rise', z[1])]),
                          ('rise-fall-rise', lambda z: [('rise', z[0]), ('fall', z[1]), ('rise', z[2])])):
            for start in ([0.3, 0.3, 0.3], [1.0, 1.0, 1.0], [3.0, 2.0, 2.0]):
                z0 = np.array(start[:{'rise-fall': 1, 'fall-rise': 2, 'rise-fall-rise': 3}[lab]])
                obj = lambda z: multi(P, 0.0, r0, form(np.abs(z)))
                res = minimize(obj, z0, method='Nelder-Mead', options=dict(xatol=1e-6, fatol=1e-9, maxiter=400))
                if res.fun < best[0] - 1e-9: best = (res.fun, f"{lab} {np.abs(res.x)}")
        print(f"    best over multi-switch families: {best[0]:.6f} ({best[1]}); field W(0, r0) = {W0:.6f}; "
              f"field beaten: {best[0] < W0 - 1e-7}")
    return P, W0
if __name__ == "__main__":
    # (a)-(c)
    run("tanh", lambda v: 1 + 0.5*np.tanh(v/5), 1.5, 6.0, 6.0)
    # (d) constant mass: exact W = int_{2m}^{r} 1/sqrt(w); sigma = 1/sqrt(w) - 1/E
    P, W0 = run("const m=1", lambda v: 1.0 + 0*v, 1.0, 6.0, 6.0, test_multi=False)
    print(f"    exact constant-mass W(0, 6) = R(6) = {Rinf(6.0, 1.0):.6f}")
    # (e) rapid accretion
    for vs_, wid, jump, r0 in ((2.0, 0.3, 2.0, 12.0), (4.0, 0.2, 3.0, 14.0), (1.0, 0.5, 1.0, 8.0)):
        run(f"jump v_s={vs_} width={wid} dm={jump} r0=r_f={r0}", lambda v, a=vs_, b=wid, c=jump: 1 + c*(1 + np.tanh((v - a)/b))/2,
            1 + jump, r0, r0)
