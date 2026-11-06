"""tau_eq_sector_change_checks.py -- floating-point checks, not a certificate, that COULD contradict the proof of the
theorem on a single transversal sector change in the example on the Kerr equator with spin (main text: "Floating-point
root-finding (not certified) gives tau_+^* = 0.160195 and tau_-^* = 0.170226, with T_s' - 1 = -1.10 and -1.09").
(1) Dilation system (W' = G/tau + D_Y G W, W(T) = 0, at T = tau) against (d_T + d_tau)R by independent centred
    differences.
(2) Roots T_s(tau) = tau in the two senses (Newton on the shooting + bisection in tau) and comparison with the rational
    enclosures.
(3) At the roots: H_diag > 49/500, Gamma_T < 71/500 (and > 91/2500), F' < -49/71; F' also by finite differences of
    T_s(tau) - tau (F' = T_s' - 1).
(4) At the roots, along the whole trip: j > 4, |xi| > 31/8, |p/xi| <= 63/7750, t < 1/500, |r - 6| < 1/300,
    sup|u_d| < 7/100, sup|P_d| < 1/25, sup|u_T| < 3/50, sup|P_T| < 1/25, int|u_T'| < 21/500, |D12|/D22 < 7/500."""
import numpy as np, sympy as sp, importlib.util, os, warnings
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
warnings.filterwarnings('ignore')
here = os.path.dirname(os.path.abspath(__file__))
def load(name, fn):
    spec = importlib.util.spec_from_file_location(name, os.path.join(here, fn))
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m
t = load('t64', 'tau_eq_arrival_shooting.py'); e = load('e66', 'tau_eq_injectivity_scan.py')
ok = lambda lab, c: print(('PASS ' if c else 'FAIL ') + lab, flush=True)
r, pr, J, al = t.r, t.pr, t.pp, t.al; H = t.H
xi_f = sp.lambdify((r, pr, J, al), J - al**2/(sp.sqrt(2)*r), 'numpy')
d_ = t.N - t.bp**2/t.h0pp; aH = 1/(d_*t.h0rr); bH = t.N/(d_**2*t.h0pp)
tt_f = sp.lambdify((r, pr, J, al), aH*pr**2/(bH*(J - al**2/(sp.sqrt(2)*r))**2), 'numpy')

def F(T, th, tau, s):
    return e.reorder_res(t.backward(T, s*np.pi/2 + th, tau, s)['res'], s)
def DF(T, th, tau, s):
    return e.reorder(t.backward(T, s*np.pi/2 + th, tau, s, want_jac=True)['J'], s)
def solve(tau, s, x0):
    x = np.array(x0, float)
    for _ in range(30):
        Rv = F(x[0], x[1], tau, s)
        if np.max(np.abs(Rv)) < 1e-13: break
        x = x - np.linalg.solve(DF(x[0], x[1], tau, s), Rv)
    return x

def dilation(T, th, s):
    """tau = T: integrates W' = G/tau + D G W backwards from W(T) = 0, together with the ordinary T column; returns
    K = (d_T + d_tau)R and phase-space maxima along the trip."""
    tau = T; psi = s*np.pi/2 + th; cr, cp = t.termf(6.0, 1.0)
    p0, JJ = cr*np.cos(psi), cp*np.sin(psi); j = abs(JJ)
    G0 = t.Gf(6.0, p0, JJ, 1.0)
    def rhs(eta, y):
        rr, p, ur, wp, wf, uT, pT, fT = y; A = t.Aprof(eta, tau)
        G = t.Gf(rr, p, JJ, A); M = np.array(t.J3f(rr, p, JJ, A), dtype=float)
        vd = M @ np.array([ur, wp, 0.0]); vT = M @ np.array([uT, pT, 0.0])
        return [G[0], G[2], G[0]/tau + vd[0], G[2]/tau + vd[2], G[1]/tau + vd[1], vT[0], vT[2], vT[1]]
    y0 = [6.0, p0, 0, 0, 0, -G0[0], -G0[2], -G0[1]]          # T column: datum -G(T), l_T = 0 because A'(tau) = 0
    sol = solve_ivp(rhs, [T, 0.0], y0, rtol=1e-11, atol=1e-13, max_step=T/600)
    rr, p, ud, wp, wf, uT, pT, fT = sol.y
    K = np.array([-s*wf[-1], s*ud[-1]])
    A_ = np.array([t.Aprof(x, tau) for x in sol.t])
    xi = xi_f(rr, p, JJ, A_)
    uTp = np.array([np.array(t.J3f(a1, a2, JJ, a3), dtype=float)[0] @ np.array([b1, b2, 0.0])
                    for a1, a2, a3, b1, b2 in zip(rr, p, A_, uT, pT)])
    ext = dict(j=j, xi=np.abs(xi).min(), Z=np.abs(p/xi).max(), t=tt_f(rr, p, JJ, A_).max(), r=np.abs(rr - 6).max(),
               Ud=np.abs(ud).max(), Pd=np.abs(wp/j).max(), UT=np.abs(uT).max(), PT=np.abs(pT/j).max(),
               NT=np.trapezoid(np.abs(uTp[::-1]), sol.t[::-1]))
    return K, ext

for s, lo, hi, enc in [(1, 0.155, 0.168, (96737/606500, 12124010000/75418769133)),
                       (-1, 0.165, 0.178, (10302037/60650000, 1032630000/6033081379))]:
    print(f"sense s={s:+d}", flush=True)
    x = [0.165, 0.0]
    def Fs(tau):
        global x
        x = solve(tau, s, x); return x[0] - tau
    ts = brentq(Fs, lo, hi, xtol=1e-12)
    xs = solve(ts, s, x); T, th = xs
    # (1) dilation against diagonal centred differences
    K, ext = dilation(T, th, s)
    h = 1e-5
    Kfd = (F(T + h, th, ts + h, s) - F(T - h, th, ts - h, s))/(2*h)
    ok(f"   (1) K dilation = {K}, diagonal differences = {Kfd}, deviation {np.max(np.abs(K-Kfd)):.1e}",
       np.max(np.abs(K - Kfd)) < 1e-6)
    # (2) root and enclosures
    ok(f"   (2) tau* = {ts:.10f} in ({enc[0]:.7f}, {enc[1]:.7f}); theta* = {th:.2e}", enc[0] < ts < enc[1])
    # (3) H_diag, Gamma_T, F'
    D = DF(T, th, ts, s)
    Hd = K[0] - D[0, 1]/D[1, 1]*K[1]; GT = D[0, 0] - D[0, 1]*D[1, 0]/D[1, 1]; Fp = -Hd/GT
    hh = 1e-5
    Fp_fd = ((solve(ts + hh, s, xs)[0] - (ts + hh)) - (solve(ts - hh, s, xs)[0] - (ts - hh)))/(2*hh)
    ok(f"   (3) H_diag = {Hd:.5f} > 0.098, Gamma_T = {GT:.5f} in (0.0364, 0.142), F' = {Fp:.5f} < -0.6901; "
       f"F' by finite differences = {Fp_fd:.5f}",
       Hd > 49/500 and 91/2500 < GT < 71/500 and Fp < -49/71 and abs(Fp - Fp_fd) < 1e-4)
    # (4) phase-space extremes and bounds of the columns
    rat = abs(D[0, 1])/D[1, 1]
    ok(f"   (4) j = {ext['j']:.3f} > 4, min|xi| = {ext['xi']:.3f} > 3.875, max|p/xi| = {ext['Z']:.2e} <= 8.13e-3, "
       f"max t = {ext['t']:.1e} < 2e-3, max|r-6| = {ext['r']:.1e} < 3.33e-3",
       ext['j'] > 4 and ext['xi'] > 31/8 and ext['Z'] <= 63/7750 and ext['t'] < 1/500 and ext['r'] < 1/300)
    ok(f"       sup|u_d| = {ext['Ud']:.4f} < 0.07, sup|P_d| = {ext['Pd']:.4f} < 0.04, sup|u_T| = {ext['UT']:.4f} < 0.06, "
       f"sup|P_T| = {ext['PT']:.4f} < 0.04, int|u_T'| = {ext['NT']:.4f} < 0.042, |D12|/D22 = {rat:.5f} < 0.014",
       ext['Ud'] < 7/100 and ext['Pd'] < 1/25 and ext['UT'] < 3/50 and ext['PT'] < 1/25 and ext['NT'] < 21/500
       and rat < 7/500)
