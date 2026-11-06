"""tau_var_shooting.py -- floating-point diagnostic, not a certificate: shooting on the extremals of the Pontryagin
maximum principle for the proper-time branch with a variable factor A(eta) (a = 0, planar motion; main text, proposition
on the PMP for the proper-time branch and remark on numerical observations). M = 1; normalisation of
critical_separation_estimate.py: Ehat = sqrt(q0), A = 1 at launch.

System (independent variable eta; L constant):
  H*(eta, r, p_r, q) = q + |p|/sqrt N - omega,  |p|^2 = f^2 p_r^2 + f L^2/r^2
  r' = dH*/dp_r,  phi' = dH*/dL,  p_r' = -dH*/dr,  q' = -dH*/deta,  cost' = omega.
Unknowns (psi0, q0, T); initial data from H* = 0: |p| = sqrt N (omega - q0), p_r = |p| cos psi0/f,
L = |p| r0 sin psi0/sqrt f.  Conditions: r(T) = r_f, phi(T) = delta, q(T) = 0.
Control: constant A -> q0 = 0 and cost = length of the shallow dip (legs() of critical_separation_estimate.py).
Then continuation in lambda: A_lambda = 1 + lambda (A_target - 1)."""
import numpy as np, sympy as sp, mpmath as mp
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'interval_arithmetic'))
from critical_separation_estimate import legs

eta, r, pr, L, q, E, a, ad = sp.symbols('eta r p_r L q E a ad', real=True)
f = 1 - 2/r; N = 1/(1 - a**2*f/E**2); om = a**2*f/E
P = sp.sqrt(f**2*pr**2 + f*L**2/r**2)
H = q + P/sp.sqrt(N) - om
# A(eta) enters through a (value) and ad (derivative): d/deta acts as ad * d/da
args = (eta, r, pr, L, q, E, a, ad)
F = {k: sp.lambdify(args, v, 'numpy') for k, v in dict(
    rdot=sp.diff(H, pr), phidot=sp.diff(H, L), prdot=-sp.diff(H, r), qdot=-ad*sp.diff(H, a), om=om,
    absp=sp.sqrt(N)*(om - q), ).items()}

def flow(x, Afun, dAfun, Ev, r0, n_out=False):
    psi0, q0, T = x
    a0 = Afun(0.0)
    absp = F['absp'](0, r0, 0, 0, q0, Ev, a0, 0)
    f0 = 1 - 2/r0
    pr0 = absp*np.cos(psi0)/f0; Lv = absp*r0*np.sin(psi0)/np.sqrt(f0)
    def rhs(e, y):
        rr, ph, p, qq, c = y
        A_, dA = Afun(e), dAfun(e)
        return [F['rdot'](e, rr, p, Lv, qq, Ev, A_, dA), F['phidot'](e, rr, p, Lv, qq, Ev, A_, dA),
                F['prdot'](e, rr, p, Lv, qq, Ev, A_, dA), F['qdot'](e, rr, p, Lv, qq, Ev, A_, dA), F['om'](e, rr, p, Lv, qq, Ev, A_, dA)]
    sol = solve_ivp(rhs, (0, T), [r0, 0.0, pr0, q0, 0.0], rtol=1e-11, atol=1e-13, dense_output=n_out)
    return sol, Lv

def residual(x, Afun, dAfun, Ev, r0, rf, delta):
    sol, _ = flow(x, Afun, dAfun, Ev, r0)
    rr, ph, p, qq, c = sol.y[:, -1]
    return [rr - rf, ph - delta, qq]

def frozen_guess(qv, r0, delta):
    ts = [2 + (r0 - 2)*k/200 for k in range(1, 200)]
    vals = [float(legs(qv, r0, t)[0]) for t in ts]
    out = []
    for t1, t2, v1, v2 in zip(ts, ts[1:], vals, vals[1:]):
        if (v1 - delta)*(v2 - delta) < 0:
            t = mp.findroot(lambda t: legs(qv, r0, t)[0] - delta, (t1, t2), solver='anderson')
            d, ell, Lc = legs(qv, r0, t)
            out.append((float(t), float(ell), float(Lc), 'shallow' if v2 < v1 else 'deep'))
    return out

if __name__ == "__main__":
    qv, r0 = 2.0, 6.0; Ev = np.sqrt(qv); delta = 0.5*np.pi
    rho = lambda rr: rr*np.sqrt((rr - 2)/((qv - 1)*rr + 2))          # rho^2 = V with A = 1 (critical_separation_estimate.py)
    print(f"Control: A = 1 constant, q = 2, r0 = r_f = 6, delta = {delta/np.pi:.2f} pi")
    branches = {}
    for t, ell, Lc, kind in frozen_guess(qv, r0, delta):
        psi0 = np.pi - np.arcsin(Lc/rho(r0))
        # crude T guess by integrating with q0 = 0 until phi = delta
        sol, _ = flow((psi0, 0.0, 60.0), lambda e: 1.0, lambda e: 0.0, Ev, r0, n_out=True)
        ee = np.linspace(0, 60, 60001); ph = sol.sol(ee)[1]
        T0 = ee[np.argmax(ph >= delta)]
        x, info, ier, msg = fsolve(residual, [psi0, 0.0, T0], args=(lambda e: 1.0, lambda e: 0.0, Ev, r0, r0, delta),
                                   full_output=True, xtol=1e-12)
        sol, Lv = flow(x, lambda e: 1.0, lambda e: 0.0, Ev, r0)
        cost = sol.y[4, -1]
        print(f"  {kind:7s} dip (periastron {t:.4f}): psi0={x[0]:.6f}, q0={x[1]:+.2e}, T={x[2]:.4f}, cost {cost:.8f} "
              f"vs legs {ell:.8f} (rel {(cost-ell)/ell:+.1e}); min r {sol.y[0].min():.4f}; residual {np.abs(info['fvec']).max():.1e}")
        branches[kind] = x
    print("Continuation to variable A (A = 1 at launch): lambda in [0, 1], A_lambda = 1 + lambda (A_target - 1)")
    targets = {"decreasing, A -> 0.5 around eta_d = 4": (lambda e: 0.75 + 0.25*np.tanh(-(e - 4)), lambda e: -0.25/np.cosh(e - 4)**2),
               "decreasing, A -> 0.5 around eta_d = 8": (lambda e: 0.75 + 0.25*np.tanh(-(e - 8)), lambda e: -0.25/np.cosh(e - 8)**2),
               "increasing, A -> 1.3 around eta_d = 4": (lambda e: 1.15 + 0.15*np.tanh(e - 4), lambda e: 0.15/np.cosh(e - 4)**2)}
    for name, (At, dAt) in targets.items():
        a_launch = At(0.0)
        for kind, x in branches.items():
            xc = x.copy(); ok = True
            for lam in np.linspace(0, 1, 21)[1:]:
                Al = lambda e, lam=lam: 1 + lam*(At(e)/a_launch - 1)      # A(0) = 1 kept
                dAl = lambda e, lam=lam: lam*dAt(e)/a_launch
                xn, info, ier, msg = fsolve(residual, xc, args=(Al, dAl, Ev, r0, r0, delta), full_output=True, xtol=1e-12)
                if ier != 1 or np.abs(info['fvec']).max() > 1e-8: ok = False; break
                xc = xn
            if not ok:
                print(f"  {name}, {kind}: continuation lost at lambda={lam:.2f}"); continue
            sol, Lv = flow(xc, Al, dAl, Ev, r0, n_out=True)
            ee = np.linspace(0, xc[2], 4001); Y = sol.sol(ee)
            print(f"  {name}, {kind}: psi0={xc[0]:.5f}, q0={xc[1]:+.5f}, T={xc[2]:.4f}, cost {sol.y[4,-1]:.6f}, "
                  f"min r {Y[0].min():.4f}, min(omega - q) {np.min(F['om'](0, Y[0], 0, 0, 0, Ev, Al(ee), 0) - Y[3]):.4f}, L={Lv:.5f}")
