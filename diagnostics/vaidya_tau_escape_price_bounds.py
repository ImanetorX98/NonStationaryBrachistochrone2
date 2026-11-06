"""vaidya_tau_escape_price_bounds.py -- floating-point diagnostic, not a certificate: two-sided estimate of the escape
price on the proper-time branch of ingoing Vaidya (main text, theorem on the escape price: b_inf(0) = 2.6014, lower
bound 4.62260). m(v) = 1 + 0.5 tanh(v/5) (non-decreasing, m_inf = 1.5, h_inf = 3), E = sqrt 2, r0 = r_f = 6.
 lower = |R_inf(r0)| + R_inf(r_f),  R_inf(r) = int_{h_inf}^r ds/sqrt(E^2 - 1 + 2 m_inf/s);
 b_inf: radial return frontier, approximated by b_T with large T (backward integration of r' = a);
 prefix: maximal infall r' = -w - E sqrt(w) from (0, r0) until contact with b_inf at (v_c, r_c); D_c = int (f - r')/E dv;
 upper = D_c + (h_inf - r_c)/sqrt(E^2 - 1) + R_inf(r_f);
 explicit competitor (5 phases): infall down to b_inf + gap, maximal rise up to R_eps = h_inf + eps, angular stay
 (pi/2) at R_eps, maximal rise to r_f: its cost must lie between lower and upper (up to O(sqrt eps)).
The exact escape price B_V = 5.21295 is computed by vaidya_tau_verification_function.py."""
import numpy as np
from scipy.integrate import solve_ivp, quad
E = np.sqrt(2.0); r0 = rf = 6.0
m = lambda v: 1 + 0.5*np.tanh(v/5); minf = 1.5; hinf = 2*minf
f = lambda v, r: 1 - 2*m(v)/r; w = lambda v, r: E**2 - f(v, r)
a = lambda v, r: -w(v, r) + E*np.sqrt(w(v, r))
Rinf = lambda r: quad(lambda s: 1/np.sqrt(E**2 - 1 + 2*minf/s), hinf, r)[0]
lower = abs(Rinf(r0)) + Rinf(rf)
# b_T for large T, backward
def bT(T):
    s = solve_ivp(lambda v, y: [a(v, y[0])], (T, 0.0), [rf], rtol=1e-11, atol=1e-13, dense_output=True, max_step=0.5)
    return s.sol
bs = {T: bT(T) for T in (100.0, 200.0, 400.0)}
for T in bs: print(f"b_T(0) for T={T:.0f}: {bs[T](0.0)[0]:.6f}  (2m(0) = {2*m(0):.3f}, h_inf = {hinf})")
binf = bs[400.0]
# prefix: max infall until contact with b_inf
ev = lambda v, y: y[0] - binf(v)[0]; ev.terminal = True
pre = solve_ivp(lambda v, y: [-w(v, y[0]) - E*np.sqrt(w(v, y[0])), (f(v, y[0]) + w(v, y[0]) + E*np.sqrt(w(v, y[0])))/E],
                (0, 300), [r0, 0.0], events=ev, rtol=1e-11, atol=1e-13, max_step=0.01)
vc, rc, Dc = pre.t[-1], pre.y[0, -1], pre.y[1, -1]
upper = Dc + (hinf - rc)/np.sqrt(E**2 - 1) + Rinf(rf)
print(f"contact with b_inf: v_c = {vc:.4f}, r_c = {rc:.6f} (2m(v_c) = {2*m(vc):.6f}); prefix cost D_c = {Dc:.6f}")
print(f"lower = |R_inf(r0)| + R_inf(r_f) = {lower:.6f};  upper = D_c + (h_inf - r_c)/sqrt(E^2-1) + R_inf(r_f) = {upper:.6f};  gap {upper-lower:.6f}")
# explicit competitor
def competitor(gap, eps, delta=np.pi/2):
    ev1 = lambda v, y: y[0] - binf(v)[0] - gap; ev1.terminal = True
    s1 = solve_ivp(lambda v, y: [-w(v, y[0]) - E*np.sqrt(w(v, y[0])), (f(v, y[0]) + w(v, y[0]) + E*np.sqrt(w(v, y[0])))/E],
                   (0, 300), [r0, 0.0], events=ev1, rtol=1e-11, atol=1e-13, max_step=0.01)
    v1, r1, c1 = s1.t[-1], s1.y[0, -1], s1.y[1, -1]
    Reps = hinf + eps
    ev2 = lambda v, y: y[0] - Reps; ev2.terminal = True
    s2 = solve_ivp(lambda v, y: [a(v, y[0]), (f(v, y[0]) - a(v, y[0]))/E], (v1, v1 + 1e5), [r1, 0.0], events=ev2,
                   rtol=1e-11, atol=1e-13, max_step=0.5)
    if not s2.t_events[0].size: return None
    v2, c2 = s2.t[-1], s2.y[1, -1]
    ev3 = lambda v, y: y[1] - delta; ev3.terminal = True
    s3 = solve_ivp(lambda v, y: [f(v, Reps)/E, np.sqrt(max(f(v, Reps)*(1 - f(v, Reps)/E**2), 0))/Reps], (v2, v2 + 1e7),
                   [0.0, 0.0], events=ev3, rtol=1e-11, atol=1e-13, max_step=5.0)
    v3, c3 = s3.t[-1], s3.y[0, -1]
    ev4 = lambda v, y: y[0] - rf; ev4.terminal = True
    s4 = solve_ivp(lambda v, y: [a(v, y[0]), (f(v, y[0]) - a(v, y[0]))/E], (v3, v3 + 1e5), [Reps, 0.0], events=ev4,
                   rtol=1e-11, atol=1e-13, max_step=0.5)
    return c1 + c2 + c3 + s4.y[1, -1], (v1, v2, v3, s4.t[-1]), (c1, c2, c3, s4.y[1, -1])
for gap, eps in ((1e-2, 1e-2), (1e-4, 1e-3), (1e-6, 1e-4)):
    c = competitor(gap, eps)
    if c is None: print(f"  gap={gap:.0e}, eps={eps:.0e}: R_eps not reached"); continue
    tot, vs, cs = c
    print(f"  competitor gap={gap:.0e} eps={eps:.0e}: cost {tot:.6f} (prefix {cs[0]:.4f}, chase {cs[1]:.4f}, arc {cs[2]:.4f}, rise {cs[3]:.4f}); "
          f"times v1={vs[0]:.1f} v2={vs[1]:.1f} v3={vs[2]:.1f} arrival={vs[3]:.1f}; within [lower, upper + O(sqrt eps)]: {lower - 1e-9 <= tot}")
