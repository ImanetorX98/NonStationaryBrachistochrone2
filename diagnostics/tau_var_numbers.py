"""tau_var_numbers.py -- floating-point diagnostic, not a certificate: the numbers of the remark on numerical
observations for the proper-time branch with a variable factor (main text). M = 1, Ehat = sqrt 2, r0 = r_f = 6,
delta = pi/2, A dropping from 1 to 1/2 (tanh of width 1) around eta_d.
 (1) shooting against direct optimisation, with the SAME raw factor A in both (A = 0.75 + 0.25 tanh(eta_d - eta) and
     A = 1.15 + 0.15 tanh(eta - eta_d)): costs of the shooting extremals (tau_var_shooting.py) against the values of
     the direct optimisation read from ../results/diagnostics_tau_var_direct_optimisation.txt; relative differences;
 (2) crossing of the continued interior extremal with the cost of the radial-descent dive
     U_rad = D_down + R_{A_inf}(6) (descent cut at 2 + 1e-7), A normalised to A(0) = 1: bisection in eta_d ("near
     eta_d = 8.03");
 (3) fold of the family continued in eta_d: arclength continuation (as tau_var_fold.py) with a fixed small step, the
     maximum of eta_d along the branch located by a parabola through the three points around the turn ("folds near
     eta_d = 10.35");
 (4) domains of the semi-Lagrangian programme (tau_var_value_function.py): r1 (a curve reaching r1 costs at least the
     candidate), r2 (excursion beyond r2 costs at least the candidate), duration H = C Ehat/(A_min^2 f(r1)).
Parts (1), (2) and (4) reproduce the independent checks of the auditor (Codex); part (3) is new."""
import os, re, sys
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq
import tau_var_shooting as D
import tau_var_fold as F

HERE = os.path.dirname(os.path.abspath(__file__))
E = np.sqrt(2.0); f = lambda r: 1 - 2/r
ok_all = True

print("(1) shooting vs direct optimisation, same raw A", flush=True)
direct = {}
for name in ('diagnostics_tau_var_direct_optimisation.txt',):
    sec = None
    for line in open(os.path.join(HERE, '..', 'results', name)):
        if line.startswith('Decreasing'): sec = False
        if line.startswith('Increasing'): sec = True
        m = re.match(r'\s*eta_d=\s*([\d.]+): cost ([\d.]+)', line)
        if m and sec is not None: direct[(float(m.group(1)), sec)] = float(m.group(2))
rows = [(4., False, [2.77349, .30496, 13.1503]), (8., False, [3.03709, .17371, 17.7413]), (4., True, [2.47262, -.39151, 15.9364])]
for ed, growth, seed in rows:
    if growth:
        A = lambda t: 1.15 + .15*np.tanh(t - ed); dA = lambda t: .15/np.cosh(t - ed)**2
    else:
        A = lambda t: .75 + .25*np.tanh(ed - t); dA = lambda t: -.25/np.cosh(t - ed)**2
    x, info, ier, msg = D.fsolve(D.residual, seed, args=(A, dA, E, 6., 6., np.pi/2), full_output=True, xtol=1e-11)
    res = np.max(np.abs(info['fvec']))
    sol, L = D.flow(x, A, dA, E, 6., n_out=True)
    tt = np.linspace(0, x[2], 1001); r, ph, p, q, c = sol.sol(tt)
    fr = 1 - 2/r; N = 1/(1 - A(tt)**2*fr/E**2)
    H = q + np.sqrt(fr*fr*p*p + fr*L*L/(r*r))/np.sqrt(N) - A(tt)**2*fr/E
    cs = sol.y[4, -1]; cd = direct[(ed, growth)]; rel = (cd - cs)/cs
    good = sol.success and res < 1e-8 and abs(rel) < 1e-4
    ok_all &= good
    print(f"  {'increasing' if growth else 'decreasing'} A, eta_d = {ed}: A(0) = {A(0.):.10f}; shooting {cs:.10f} "
          f"(T {x[2]:.6f}, residual {res:.1e}, max|H*| {np.max(np.abs(H)):.1e}); direct {cd:.6f}; relative difference {rel:+.1e} "
          f"{'PASS' if good else 'FAIL'}", flush=True)

print("(2) crossing with the radial-descent dive, A normalised to A(0) = 1", flush=True)
seed8 = np.array([3.03709298, .17370671, 17.74132463])
def difference(ed, verbose=True):
    n = .75 + .25*np.tanh(ed)
    A = lambda t: (.75 + .25*np.tanh(ed - t))/n; dA = lambda t: (-.25/np.cosh(ed - t)**2)/n
    x, info, ier, msg = D.fsolve(D.residual, seed8, args=(A, dA, E, 6., 6., np.pi/2), full_output=True, xtol=1e-11)
    assert np.max(np.abs(info['fvec'])) < 1e-8
    sol, L = D.flow(x, A, dA, E, 6.)
    ev = lambda t, y: y[0] - (2 + 1e-7); ev.terminal = True
    down = solve_ivp(lambda t, y: [-f(y[0])*np.sqrt(1 - A(t)**2*f(y[0])/E**2), A(t)**2*f(y[0])/E], (0, 200), [6., 0.],
                     events=ev, rtol=1e-11, atol=1e-13, max_step=.2)
    assert down.success and len(down.t_events[0]) == 1
    ainf = .5/n
    up = ainf**2/E*quad(lambda r: 1/np.sqrt(1 - ainf**2*f(r)/E**2), 2, 6, epsabs=1e-12)[0]
    d = sol.y[4, -1] - down.y[1, -1] - up
    if verbose:
        print(f"  eta_d = {ed}: C_ext {sol.y[4, -1]:.10f}, U_rad {down.y[1, -1] + up:.10f}, difference {d:+.3e}", flush=True)
    return d, x
d8, _ = difference(8.); d825, _ = difference(8.25)
root = brentq(lambda e: difference(e, False)[0], 8., 8.25, xtol=1e-8)
_, xr = difference(root, False)
good = d8 < 0 < d825 and abs(root - 8.03) < 0.005
ok_all &= good
print(f"  crossing eta_d* = {root:.8f} (T {xr[2]:.4f}, q0 {xr[1]:+.5f}) {'PASS' if good else 'FAIL'} (main text: near 8.03)", flush=True)

print("(3) fold of the family continued in eta_d (arclength, fixed step)", flush=True)
za = np.array([3.103601, 0.08014, 24.406, 10.25]); zb = np.array([3.107501, 0.07576, 25.409, 10.34375])
za[:3] = D.fsolve(lambda x: F.F3(np.r_[x, za[3]]), za[:3], xtol=1e-12)
zb[:3] = D.fsolve(lambda x: F.F3(np.r_[x, zb[3]]), zb[:3], xtol=1e-12)
pts = [za, zb]; h = 0.01; turned = False
for it in range(400):
    t = (pts[-1] - pts[-2])*F.SC; t /= np.linalg.norm(t)
    zp = pts[-1] + h*t/F.SC
    z, inf, ier, msg = D.fsolve(lambda z: list(F.F3(z)) + [np.dot((z - zp)*F.SC, t)], zp, full_output=True, xtol=1e-12)
    if ier != 1 or np.abs(inf['fvec']).max() > 1e-8: print(f"  continuation failed at step {it}"); break
    pts.append(z)
    if len(pts) >= 4 and pts[-1][3] < pts[-2][3] and pts[-3][3] < pts[-2][3]:
        turned = True; break
if turned:
    e3 = np.array([p[3] for p in pts[-3:]])                       # equally spaced in arclength (step h)
    a2 = (e3[0] - 2*e3[1] + e3[2])/2; a1 = (e3[2] - e3[0])/2; s = -a1/(2*a2)
    emax = e3[1] + a1*s + a2*s*s
    c, rmin, wq = F.info(pts[-2])
    good = abs(emax - 10.35) < 0.01
    ok_all &= good
    print(f"  eta_d turns back after {len(pts) - 2} steps of {h}: maximum eta_d = {emax:.5f} (T {pts[-2][2]:.3f}, "
          f"min r {rmin:.4f}, cost {c:.6f}, min(omega - q) {wq:+.5f}) {'PASS' if good else 'FAIL'} (main text: near 10.35)", flush=True)
else:
    ok_all = False; print("  no fold found: FAIL", flush=True)

print("(4) domains of the semi-Lagrangian programme", flush=True)
def R(a, r): return a*a/E*quad(lambda s: 1/np.sqrt(1 - a*a*f(s)/E**2), 2, r, epsabs=1e-11)[0]
def descent(A, rb):
    ev = lambda t, y: y[0] - rb; ev.terminal = True
    sol = solve_ivp(lambda t, y: [-f(y[0])*np.sqrt(1 - A(t)**2*f(y[0])/E**2), A(t)**2*f(y[0])/E], (0, 1000), [6., 0.],
                    events=ev, rtol=2e-11, atol=1e-13, max_step=.1)
    assert sol.success and len(sol.t_events[0]) == 1
    return sol.t[-1], sol.y[1, -1]
for mode in ('frozen', 'drop', 'up'):
    if mode == 'frozen': A = lambda t: 1.; amin = 1.; C = 5.81303660
    if mode == 'drop':
        n = .75 + .25*np.tanh(4.); A = lambda t: (.75 + .25*np.tanh(4. - t))/n; amin = .5/n; C = 2.566213
    if mode == 'up':
        n = 1.15 - .15*np.tanh(4.); A = lambda t: (1.15 + .15*np.tanh(t - 4.))/n; amin = 1.; C = 9.821112
    def lower(rb):
        T, Dd = descent(A, rb); af = A(T) if mode == 'up' else amin
        return Dd + R(af, 6.) - R(af, rb)
    r1 = brentq(lambda r: lower(r) - C, 2.0001, 5.999)
    r2 = brentq(lambda r: 2*R(amin, r) - 2*R(amin, 6.) - C, 6., 100.)
    H = C*E/(amin**2*f(r1))
    print(f"  {mode:6s}: candidate {C:.6f}; r1 = {r1:.4f}, r2 = {r2:.4f}, H = {H:.2f}", flush=True)
print('SUMMARY:', 'all checks PASS' if ok_all else 'at least one check FAILED')
sys.exit(0 if ok_all else 1)
