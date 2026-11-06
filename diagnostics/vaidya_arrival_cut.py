"""vaidya_arrival_cut.py -- floating-point diagnostic, not a certificate: the cut on the arrival branch of ingoing
Vaidya with geometric criteria (library of vaidya_ray_count.py).
 (a) In-plane conjugate points: null Hamiltonian flow of G with exact sensitivities J = dy/dpsi (SymPy); along the
     geodesic D = r_psi phi_lambda - r_lambda phi_psi; a zero of D after the launch is an in-plane conjugate point;
     tangencies to the cylinder r = r_t (r_lambda = 0) do NOT make D vanish.
 (b) Possible exchanges of the first arrival: for each angle delta of the target (r_t, delta) the geodesics arriving
     there FOUND BY THE SCAN are collected (bracketing on a finite grid of psi for each crossing index, refined near the
     changes of the number of crossings, then brentq; no interpolation between pieces), mirrored ones (2 pi - delta)
     included. The scan can miss tangential roots, pairs in the same cell, launches near the ends and arrivals beyond
     the numerical stops. The first arrival is the smallest v found. INDICATOR (not a definition of a Maxwell point):
     the family that realised the first arrival still exists at the next delta and another arrives earlier. The tie is
     NOT localized by bisection nor confirmed by a double shot: candidates are only printed. A tangential tie without
     exchange can be missed.
 (c) Integration check: |H| (must stay ~0), dv/dlambda > 0, reason for stopping of every shot.
Controls: constant mass (cut at pi)."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'tools'))
from progress import Progress
import numpy as np, sympy as sp, sys
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
v, r, ph, Pv, Pr, L, m, dm, E = sp.symbols('v r phi P_v P_r L m dm E', real=True)
f = 1 - 2*m/r; w = E**2 - f
H = sp.Rational(1, 2)*(-Pv**2/w + 2*Pv*Pr + f*Pr**2 + L**2/r**2)
state = [v, r, ph, Pv, Pr, L]
rhs = [sp.diff(H, Pv), sp.diff(H, Pr), sp.diff(H, L), -dm*sp.diff(H, m), -sp.diff(H, r), sp.Integer(0)]
# m depends on v: d/dv of rhs through m with dm; second derivative of m needed in d(rhs)/dv -> ddm
ddm = sp.symbols('ddm', real=True)
def dvtot(ex):   # total derivative wrt v of an expression depending on v through m and dm
    return sp.diff(ex, m)*dm + sp.diff(ex, dm)*ddm
Jac = sp.Matrix([[ (dvtot(fi) if s_ is v else sp.diff(fi, s_)) for s_ in state] for fi in rhs])
args = (r, Pv, Pr, L, m, dm, ddm, E)
Frhs = sp.lambdify(args, rhs, 'numpy'); FJ = sp.lambdify(args, Jac, 'numpy'); FH = sp.lambdify(args, H, 'numpy')
Ev = np.sqrt(2.0)
def launch(psi, mfun, r0):
    m0 = mfun(0.0); f0 = 1 - 2*m0/r0; w0 = Ev**2 - f0
    pr = np.cos(psi); Lv = r0*np.sin(psi)/np.sqrt(f0)
    S = np.sqrt(w0**2*pr**2 + w0*(f0*pr**2 + Lv**2/r0**2)); pv = w0*pr - S
    # derivatives wrt psi
    dpr = -np.sin(psi); dL = r0*np.cos(psi)/np.sqrt(f0)
    dS = (w0**2*pr*dpr + w0*(f0*pr*dpr + Lv*dL/r0**2))/S; dpv = w0*dpr - dS
    return np.array([0.0, r0, 0.0, pv, pr, Lv]), np.array([0.0, 0.0, 0.0, dpv, dpr, dL])
def shoot(psi, mfun, dmfun, ddmfun, r0, rt=None, lam_max=400.0, vmax=200.0, jac=False, dense=False):
    y0, J0 = launch(psi, mfun, r0)
    def F(l, Y):
        y = Y[:6]; a = (y[1], y[3], y[4], y[5], mfun(y[0]), dmfun(y[0]), ddmfun(y[0]), Ev)
        out = np.array(Frhs(*a), dtype=float)
        if jac: out = np.concatenate([out, np.array(FJ(*a), dtype=float) @ Y[6:12]])
        return out
    evs = []
    if rt is not None:
        e0 = lambda l, Y: Y[1] - rt; evs.append(e0)
    hor = lambda l, Y: Y[1] - 2*mfun(Y[0]) - 1e-6; hor.terminal = True
    late = lambda l, Y: Y[0] - vmax; late.terminal = True
    wind = lambda l, Y: Y[2] - 2*np.pi; wind.terminal = True
    evs += [hor, late, wind]
    Y0 = np.concatenate([y0, J0]) if jac else y0
    s = solve_ivp(F, (0, lam_max), Y0, events=evs, rtol=1e-10, atol=1e-12, dense_output=dense)
    stop = 'horizon' if s.t_events[-3].size else 'v_max' if s.t_events[-2].size else 'winding 2pi' if s.t_events[-1].size else ('lam_max' if s.success else 'FAIL')
    cr = [(l, Y) for l, Y in zip(s.t_events[0], s.y_events[0]) if l > 1e-6] if rt is not None else []
    return s, cr, stop
def all_arrivals(mfun, dmfun, ddmfun, r0, rt, deltas, n=1601):
    psis = np.linspace(1e-4, np.pi - 1e-4, n)
    data = [shoot(p, mfun, dmfun, ddmfun, r0, rt)[1] for p in psis]
    cnt = [len(c) for c in data]; extra = []
    for i in range(n - 1):
        if cnt[i] != cnt[i + 1]:
            a, b, ca = psis[i], psis[i + 1], cnt[i]
            for _ in range(40):
                mid = (a + b)/2
                if len(shoot(mid, mfun, dmfun, ddmfun, r0, rt)[1]) == ca: a = mid
                else: b = mid
            extra += [a - 10.0**(-j) for j in range(3, 13)] + [b + 10.0**(-j) for j in range(3, 13)] + [a, b]
    if extra:
        psis = np.unique(np.concatenate([psis, np.clip(extra, 1e-4, np.pi - 1e-4)]))
        data = [shoot(p, mfun, dmfun, ddmfun, r0, rt)[1] for p in psis]
    res = {}
    for d in deltas:
        sols = []
        for target, mir in ((d, False), (2*np.pi - d, True)):
            for k in range(max(len(c) for c in data)):
                for i in range(len(psis) - 1):
                    a, b = data[i], data[i + 1]
                    if len(a) > k and len(b) > k and (a[k][1][2] - target)*(b[k][1][2] - target) < 0:
                        def g(p):
                            c = shoot(p, mfun, dmfun, ddmfun, r0, rt)[1]
                            return c[k][1][2] - target if len(c) > k else np.nan
                        try:
                            ps = brentq(g, psis[i], psis[i + 1], xtol=1e-13)
                        except ValueError:
                            continue
                        c = shoot(ps, mfun, dmfun, ddmfun, r0, rt)[1]
                        if len(c) > k: sols.append((c[k][1][0], ps, k, mir, c[k][0]))
        sols.sort(); res[d] = sols
    return res
def D_along(psi, mfun, dmfun, ddmfun, r0, lam_end):
    s, _, stop = shoot(psi, mfun, dmfun, ddmfun, r0, None, lam_max=lam_end, jac=True, dense=True)
    ls = np.linspace(0, s.t[-1], 4001)[1:]; Y = s.sol(ls)
    Dv = []; Hmax = 0.0; vdot_min = np.inf
    for j, l in enumerate(ls):
        y = Y[:6, j]; J = Y[6:12, j]
        a = (y[1], y[3], y[4], y[5], mfun(y[0]), dmfun(y[0]), ddmfun(y[0]), Ev)
        k = np.array(Frhs(*a), dtype=float)
        Dv.append(J[1]*k[2] - k[1]*J[2]); Hmax = max(Hmax, abs(FH(*a))); vdot_min = min(vdot_min, k[0])
    Dv = np.array(Dv); sgn = np.sign(Dv[50:])            # skip the launch, where D = 0
    zeros = np.where(np.diff(sgn) != 0)[0]
    return ls[50:][zeros] if len(zeros) else [], Y[2, :], Hmax, vdot_min, stop
def report(name, mfun, dmfun, ddmfun, r0, rt, ndelta=40):
    deltas = np.linspace(0.05, 1.0, ndelta)*np.pi
    res = all_arrivals(mfun, dmfun, ddmfun, r0, rt, deltas)
    print(f"[{name}] r0 = {r0}, r_t = {rt}")
    prev = None; maxwell = []; conj = []; unreached = []
    for d in deltas:
        sols = res[d]
        if not sols: unreached.append(d/np.pi); prev = None; continue
        best = sols[0]
        if prev is not None and d < np.pi - 1e-9:
            # continuation of the previous realising family: same crossing index and mirror flag, closest launch
            same = [x for x in sols if x[2] == prev[1][2] and x[3] == prev[1][3]]
            cont = min(same, key=lambda x: abs(x[1] - prev[1][1])) if same else None
            if cont is not None and cont is not best and best[0] < cont[0] - 1e-9:
                maxwell.append((prev[0], d, cont, best))
        zs, phis, Hm, vdm, stop = D_along(best[1], mfun, dmfun, ddmfun, r0, best[4])
        if len(zs): conj.append((round(d/np.pi, 4), round(best[1], 6), len(zs)))
        prev = (d, best)
    gap = [(d/np.pi, res[d][1][0] - res[d][0][0]) for d in deltas if len(res[d]) > 1]
    print(f"    Maxwell candidates before pi (previous realising family still exists but another arrives earlier): {len(maxwell)}")
    for d0, d1, a, b in maxwell:
        print(f"      between delta = {d0/np.pi:.4f} and {d1/np.pi:.4f} pi: family psi {a[1]:.6f} (k {a[2]}{', mir' if a[3] else ''}, "
              f"v {a[0]:.5f}) beaten by psi {b[1]:.6f} (k {b[2]}{', mir' if b[3] else ''}, v {b[0]:.5f}): to be confirmed by a double shot")
    print(f"    zeros of D (in-plane conjugate points) before arrival on the realising extremal: "
          + ("none" if not conj else str(conj)))
    print(f"    deltas with no arrival found: {['%.3f' % u for u in unreached] or 'none'}")
    if gap: print(f"    smallest gap between the two earliest arrivals for delta < pi: {min(g[1] for g in gap if g[0] < 0.999):.4f} "
                  f"(at delta/pi = {min((g for g in gap if g[0] < 0.999), key=lambda t: t[1])[0]:.3f}); at pi: "
                  f"{[g[1] for g in gap if g[0] > 0.999]}")
    j = np.argmin(np.abs(deltas - np.pi/2))
    if res[deltas[j]]:
        b = res[deltas[j]][0]; zs, phis, Hm, vdm, stop = D_along(b[1], mfun, dmfun, ddmfun, r0, b[4])
        print(f"    integration check at delta = pi/2: max |H| = {Hm:.1e}, min dv/dlambda = {vdm:.3e}")
    # Lemma K check along the realising extremal at delta = pi: does it reach phi = pi inside D?
    return res
if __name__ == "__main__":
    cases = {
        'frozen': ("frozen m = 1", lambda x: 1.0, lambda x: 0.0, lambda x: 0.0),
        'tanh': ("m = 1 + 0.5 tanh(v/5)", lambda x: 1 + 0.5*np.tanh(x/5), lambda x: 0.1/np.cosh(x/5)**2,
                 lambda x: -0.04*np.tanh(x/5)/np.cosh(x/5)**2),
        'jump': ("rapid jump m = 1 + (1 + tanh((v - 4)/0.3))/2",
                 lambda x: 1 + (1 + np.tanh((x - 4)/0.3))/2,
                 lambda x: (1/0.6)/np.cosh(np.clip((x - 4)/0.3, -300, 300))**2,
                 lambda x: -(1/0.6)*(2/0.3)*np.tanh((x - 4)/0.3)/np.cosh(np.clip((x - 4)/0.3, -300, 300))**2)}
    args = sys.argv[1:] or ['frozen', 'tanh', 'jump']
    rts = [6.0, 10.0]
    pg = Progress(len(rts)*len(args), 'vaidya_arrival_cut: (r_t, profile)')
    for rt in rts:
        for k in args: report(*cases[k], 6.0, rt); pg.step()
    pg.close()
