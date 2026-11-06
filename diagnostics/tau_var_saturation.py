"""tau_var_saturation.py -- floating-point diagnostic (not a certificate) of the saturation lemma of the
proper-time branch, a = 0, variable A(eta) (Thakurta-Schwarzschild; main text, saturation lemma).

Relaxed curve: spatial path x(s), s in [0, 1], and profile eta(s) with eta' >= G(eta, s) = sqrt(N) |x'|_{h0}
(causal for g_opt, not necessarily null). Cost C = int omega(eta, x) eta' ds, omega = A^2 f / E.
Rail on the same path: eta_r' = G(eta_r, s), eta_r(0) = eta(0); cost C_r = int phi(eta_r, x)|x'| ds,
phi = A^2 f sqrt(N) / E.
Lemma: if A' >= 0, then C >= C_r. Control that must fail: A' < 0 with waiting near 2M."""
import numpy as np
from scipy.integrate import solve_ivp, quad
rng = np.random.default_rng(1)
M, E = 1.0, 1.0

def f(r): return 1 - 2*M/r
def N(A, r): return 1/(1 - A**2*f(r)/E**2)

def path(kind):
    if kind == 'arc':   # r(s), phi(s): dip from 6M through 3M, angle 2
        r = lambda s: 6 - 12*s*(1 - s); ph = lambda s: 2*s
    else:               # near-horizon: down to 2.001M, around, back
        r = lambda s: 2.001 + 3.999*(2*s - 1)**8; ph = lambda s: 3*s
    return r, ph

def speed(r, ph, s, h=1e-6):
    dr = (r(s + h) - r(s - h))/(2*h); dp = (ph(s + h) - ph(s - h))/(2*h)
    return np.sqrt(dr**2/f(r(s))**2 + r(s)**2*dp**2/f(r(s)))      # |x'|_{h0}

def costs(Afun, kind, extra):
    """extra(s) >= 0: relaxed eta' = G + extra (waiting)"""
    r, ph = path(kind)
    G = lambda e, s: np.sqrt(N(Afun(e), r(s)))*speed(r, ph, s)
    om = lambda e, s: Afun(e)**2*f(r(s))/E
    rhs_rel = lambda s, y: [G(y[0], s) + extra(s), om(y[0], s)*(G(y[0], s) + extra(s))]
    rhs_rail = lambda s, y: [G(y[0], s), om(y[0], s)*G(y[0], s)]
    o = dict(rtol=1e-10, atol=1e-12, max_step=1e-3)
    a = solve_ivp(rhs_rel, (1e-6, 1 - 1e-6), [0.0, 0.0], **o).y[:, -1]
    b = solve_ivp(rhs_rail, (1e-6, 1 - 1e-6), [0.0, 0.0], **o).y[:, -1]
    return a[1], b[1], a[0], b[0]

if __name__ == "__main__":
    print("A' >= 0 (lemma): C - C_r must be >= 0")
    worst = np.inf
    for trial in range(40):
        H = rng.uniform(0, 0.05); A0 = rng.uniform(0.2, 0.6)
        Afun = lambda e, H=H, A0=A0: min(A0*np.exp(H*e), 0.95)          # non-decreasing, below E
        c = rng.uniform(0, 5, 3)
        extra = lambda s, c=c: c[0]*np.exp(-((s - c[1]/5)/0.1)**2) + 0.3*c[2]*s
        kind = 'arc' if trial % 2 else 'dive'
        C, Cr, eta, etar = costs(Afun, kind, extra)
        worst = min(worst, C - Cr)
        if trial < 6: print(f"  {kind:4s} H={H:.3f} A0={A0:.2f}: C={C:.6f} C_r={Cr:.6f} diff={C-Cr:+.3e} eta_end {eta:.3f} >= {etar:.3f}")
    print(f"  min over 40 trials of C - C_r = {worst:+.3e}  -> {'PASS' if worst >= -1e-9 else 'FAIL'}")
    print("Control 1 (A' < 0 but drop before the waiting matters: no gain expected, not a test)")
    Afun = lambda e: 0.3 + 0.6/(1 + np.exp(20*(e - 5)))
    r, ph = path('dive')
    # wait where r is closest to 2M (s near 1/2)
    for W in (0.0, 20.0, 60.0, 120.0):
        extra = lambda s, W=W: W*np.exp(-((s - 0.5)/0.02)**2)
        C, Cr, eta, etar = costs(Afun, 'dive', extra)
        print(f"  waiting weight {W:5.1f}: C={C:.6f} C_r={Cr:.6f} diff={C-Cr:+.4e} (eta_end {eta:.2f} vs {etar:.2f})")

def control2():
    """Must fail for A' < 0. Arc from 6M down to r_b and back (angle 2). A = 0.9 until the frozen rail would
    arrive, then drops to 0.3 (tanh). Waiting W at the bottom (eta' = G + W * narrow bump at s = 1/2) costs
    A^2 f W, small when f(r_b) is small, and moves the second half after the drop."""
    for rb in (2.05, 2.2, 3.0):
        r = lambda s: 6 - 4*(6 - rb)*s*(1 - s); ph = lambda s: 2*s
        sp = lambda s: speed(r, ph, s)
        o = dict(rtol=1e-10, atol=1e-12, max_step=1e-3)
        G0 = lambda s, y: [np.sqrt(N(0.9, r(s)))*sp(s)]
        eb = solve_ivp(G0, (1e-6, 0.5), [0.0], **o).y[0, -1]; ee = solve_ivp(G0, (1e-6, 1 - 1e-6), [0.0], **o).y[0, -1]
        Afun = lambda e: 0.3 + 0.3*(1 - np.tanh(40*(e - ee)))
        res = []
        for W in (0.0, 0.5*(ee - eb), ee - eb + 1, 2*(ee - eb)):
            extra = lambda s, W=W: W*np.exp(-((s - 0.5)/0.005)**2)/(0.005*np.sqrt(np.pi))
            G = lambda e, s: np.sqrt(N(Afun(e), r(s)))*sp(s)
            om = lambda e, s: Afun(e)**2*f(r(s))/E
            rhs = lambda s, y: [G(y[0], s) + extra(s), om(y[0], s)*(G(y[0], s) + extra(s))]
            y = solve_ivp(rhs, (1e-6, 1 - 1e-6), [0.0, 0.0], rtol=1e-10, atol=1e-12, max_step=2e-4).y[:, -1]
            res.append((W, y[1]))
        print(f"  r_b={rb}: rail eta bottom {eb:.2f}, end {ee:.2f}; " + ", ".join(
            f"W={W:.1f}: C={c:.4f} (C-C_r {c - res[0][1]:+.4f})" for W, c in res))

if __name__ == "__main__":
    print("Control 2 (A' < 0, must show C < C_r): ")
    control2()
