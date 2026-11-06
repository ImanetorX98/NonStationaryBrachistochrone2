"""vaidya_tau_calibration_identity.py -- floating-point check, not a certificate: identity of the verification
function of vaidya_tau_verification_function.py (main text, theorem on the escape price: B_V = 5.21295 and the shallow
branch crosses it at delta = 0.522 pi).
Along the maximal infall rho(t) from (v, r) to contact with b_inf, xi = d rho/d r solves xi' = (2m/rho^2)(1 + E/(2 sqrt w)) xi,
xi(v) = 1. Identity (proved in the main text):
   W_r = xi_c/E - int (m/(rho^2 sqrt w)) xi dt,     sigma = W_r - 1/E = (1/E) int_v^{v_c} (2m/rho^2) xi dt  (> 0).
 (1) sympy: d s_-/dr and d sqrt(w)/dr (exact residuals);
 (2) formula against finite differences of W, tanh example and rapid jump;
 (3) constant mass: sigma = 1/sqrt(w(r)) - 1/E exactly;
 (4) control that must fail: without the contact term xi_c/E (that is, ignoring the shift of v_c) the formula does
     not reproduce the finite differences;
 (5) exact escape price B_V = W(0, 6) + R_inf(6) and angle at which the shallow branch crosses it (tanh example)."""
import numpy as np, sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from vaidya_tau_verification_function import E, setup, Rinf
r_, m_, Es = sp.symbols('r m E', positive=True); f_ = 1 - 2*m_/r_; w_ = Es**2 - f_
sm_ = -w_ - Es*sp.sqrt(w_)
res1 = sp.simplify(sp.diff(sm_, r_) - 2*m_/r_**2*(1 + Es/(2*sp.sqrt(w_))))
res2 = sp.simplify(sp.diff(sp.sqrt(w_), r_) + m_/(r_**2*sp.sqrt(w_)))
print("(1) d s_-/dr - (2m/r^2)(1 + E/(2 sqrt w)) =", res1, "; d sqrt(w)/dr + m/(r^2 sqrt w) =", res2)
assert res1 == 0 and res2 == 0, "symbolic identity failed"
def analytic(P, mfun, v, r):
    ev = lambda t, y: y[0] - P['binf'](t); ev.terminal = True; ev.direction = -1
    def rhs(t, y):
        rho, xi, I1, I2 = y; mm = mfun(t); ww = P['w'](t, rho)
        return [P['sm'](t, rho), 2*mm/rho**2*(1 + E/(2*np.sqrt(ww)))*xi, mm/(rho**2*np.sqrt(ww))*xi, 2*mm/rho**2*xi]
    s = solve_ivp(rhs, (v, v + 200), [r, 1.0, 0.0, 0.0], events=ev, rtol=1e-12, atol=1e-14, max_step=0.01)
    assert s.success and s.t_events[0].size == 1 and s.t[-1] > v, "contact not found"
    xi_c, I1, I2 = s.y[1, -1], s.y[2, -1], s.y[3, -1]
    return xi_c/E - I1, I2/E, 1/E - I1   # W_r, sigma, W_r without contact term
cases = [("tanh", lambda v: 1 + 0.5*np.tanh(v/5), 1.5),
         ("jump", lambda v: 1 + 3*(1 + np.tanh((v - 4)/0.2))/2, 4.0),
         ("const", lambda v: 1.0 + 0*v, 1.0)]
for name, mfun, minf in cases:
    P = setup(mfun, minf); worst = 0; worst_bad = np.inf
    for v, dr in ((0.0, 0.05), (0.0, 3.4), (2.0, 0.5), (3.9, 1.0), (4.1, 0.01), (10.0, 2.0), (0.5, 6.0)):
        r = P['binf'](v) + dr; hr = 1e-4
        fd = (P['W'](v, r + hr)[0] - P['W'](v, r - hr)[0])/(2*hr)
        Wr, sig, Wr_bad = analytic(P, mfun, v, r)
        worst = max(worst, abs(Wr - fd)); worst_bad = max(worst_bad if np.isfinite(worst_bad) else 0, abs(Wr_bad - fd))
        extra = f", exact const 1/sqrt(w) - 1/E = {1/np.sqrt(P['w'](v, r)) - 1/E:.8f}" if name == "const" else ""
        assert sig > 0, "analytic sigma must be positive"
        if name == "const": assert abs(sig - (1/np.sqrt(P['w'](v, r)) - 1/E)) < 1e-9
        print(f"  [{name}] v={v:4.1f} r={r:8.4f}: W_r FD {fd:.8f}, formula {Wr:.8f}; sigma = {sig:.8f} (>0: {sig > 0}){extra}")
    print(f"(2) [{name}] max |formula - FD| = {worst:.1e};  (4) without contact term, max |error| = {worst_bad:.2e} "
          f"(must be large: {'FAILS as required' if worst_bad > 1e2*max(worst, 1e-9) else 'NOT failing'})")
    assert worst < 1e-6, "formula vs FD"
    assert worst_bad > 1e2*max(worst, 1e-9), "must-fail control did not fail"
# (5) exact escape price and crossing angle, tanh example
P = setup(cases[0][1], 1.5); BV = P['W'](0.0, 6.0)[0] + Rinf(6.0, 1.5)
print(f"(5) tanh example: B_V = W(0,6) + R_inf(6) = {BV:.6f}")
from vaidya_tau_existence_no_deadline import mass, solve, info, frozen_start
xx = frozen_start()
for e in np.linspace(0.05, 1.0, 20): xx, ok = solve(xx, mass(e), np.pi/2); assert ok
mf = mass(1.0); d = 0.5; prev = None
while True:
    xn, ok = solve(xx, mf, (d + 0.005)*np.pi); assert ok
    c = info(xn, mf)[0]
    if c >= BV: break
    xx, d = xn, d + 0.005
lo, hi = d, d + 0.005
for _ in range(40):
    mid = (lo + hi)/2; xm, ok = solve(xx, mf, mid*np.pi); assert ok
    if info(xm, mf)[0] < BV: lo, xx = mid, xm
    else: hi = mid
assert 0.5 < lo < 0.55
print(f"    shallow branch cost = B_V at delta = {lo:.6f} pi: minimum attained without deadline for delta < {lo:.4f} pi (not certified)")
