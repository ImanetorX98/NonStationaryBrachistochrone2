# -*- coding: utf-8 -*-
"""
Main text, section "Causal tools and an abstract first-arrival theorem": checks for the global hyperbolicity of
the auxiliary metric.

Symbolic (proof): reduction at a = 0, N >= 1 and the speed bound, radial completeness of h0 at r = 2M for a = 0 and
its loss for a != 0, Finsler degeneracy and N -> 1 at 2M, the freezing front.
Numerical (diagnostic only, not a proof): first arrival at a target near the antipode for an exponential A with
sup A < Ehat; the long arc is beaten by the short one.
Usage: python3 global_hyperbolicity.py
"""
import sys
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

FAIL = []


def check(name, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"   {detail}" if detail else ""))
    if not ok:
        FAIL.append(name)


r, M, a, E, A, eps, vt, vx = sp.symbols('r M a E A epsilon v_t v_x', positive=True)
f = 1 - 2*M/r
Dl = r**2 - 2*M*r + a**2
P0 = r**2/(f*Dl); Q0 = Dl/f**2; beta = -2*M*a/(r*f)
Z = A**2*f/E**2; N = 1/(1 - Z)

print("symbolic")
check("a=0: P0 = 1/f^2", sp.simplify(P0.subs(a, 0) - 1/f**2) == 0)
check("a=0: Q0 = r^2/f", sp.simplify(Q0.subs(a, 0) - r**2/f) == 0)
check("N - 1 = Z/(1-Z) >= 0 where Z<1", sp.simplify(N - 1 - Z/(1 - Z)) == 0)
# causal v: -vt^2 + N vx^2 <= 0  =>  vt >= sqrt(N) vx >= vx
check("null: vt = sqrt(N)|vx|_h0 >= |vx|_h0 for N>=1",
      sp.simplify((sp.sqrt(N)*vx)**2 - N*vx**2) == 0)
Lr = sp.integrate(1/f, (r, 2*M + eps, 3*M))
check("a=0: radial h0 distance to 2M diverges", sp.limit(Lr, eps, 0, '+') == sp.oo)
lead = sp.limit(sp.sqrt(P0).subs(r, 2*M + eps)*sp.sqrt(eps), eps, 0, '+')
check("a!=0: sqrt(P0) ~ C/sqrt(r-2M), C finite and nonzero (integrable)",
      lead.is_finite and lead != 0, f"C = {sp.simplify(lead)}")
check("a!=0: |beta|/sqrt(Q0) -> 1 at 2M",
      sp.limit((2*M*a/(r*f))/sp.sqrt(Q0), r, 2*M, '+') == 1)
check("N -> 1 at 2M", sp.limit(N, r, 2*M, '+') == 1)
# freezing front for A > E: Z<1 <=> r < 2M/(1 - E^2/A^2)
rF = 2*M/(1 - E**2/A**2)
check("freezing front: Z(rF) = 1", sp.simplify(Z.subs(r, rF) - 1) == 0)

print("\nnumerical diagnostic (not a proof)")
Mv, Ev, epsv, r0 = 1.0, 1.4, 0.004, 6.0
Afun = lambda t: np.exp(epsv*t)          # sup on [0,60] = e^0.24 = 1.27 < 1.4
tt, rr_, pp, pr_, pf = sp.symbols('tt rr pp pr pf')
ff = 1 - 2*Mv/rr_
Nn = 1/(1 - sp.exp(2*epsv*tt)*ff/Ev**2)
# K = 1/2 [ -pt^2 + (f^2 pr^2 + f pf^2 / r^2)/N ],  h0^{-1} = diag(f^2, f/r^2)
pt = sp.symbols('pt')
K = sp.Rational(1, 2)*(-pt**2 + (ff**2*pr_**2 + ff*pf**2/rr_**2)/Nn)
vars_ = (tt, rr_, pt, pr_, pf)
dK = {s: sp.lambdify(vars_, sp.diff(K, s), 'numpy') for s in (tt, rr_, pt, pr_, pf)}


def rhs(l, y):
    t, rv, ph, ptv, prv, pfv = y
    args = (t, rv, ptv, prv, pfv)
    return [dK[pt](*args), dK[pr_](*args), dK[pf](*args),
            -dK[tt](*args), -dK[rr_](*args), 0.0]


def shoot(theta, phi_target):
    """null geodesic from (0, r0, 0) with launch angle theta; returns eta at which
    Phi reaches phi_target (sign-aware) and r there."""
    N0 = float(Nn.subs({tt: 0, rr_: r0}))
    f0 = 1 - 2*Mv/r0
    # future null covector with p_t = -1: (f^2 pr^2 + f pf^2/r^2)/N = 1
    prv = np.sqrt(N0)*np.cos(theta)/f0
    pfv = np.sqrt(N0)*np.sin(theta)*r0/np.sqrt(f0)
    ev = lambda l, y: y[2] - phi_target
    ev.terminal = True
    far = lambda l, y: (y[1] - 40.0)*(y[1] - 2.2)
    far.terminal = True
    late = lambda l, y: y[0] - 60.0
    late.terminal = True
    sol = solve_ivp(rhs, [0, 400], [0, r0, 0, -1.0, prv, pfv], events=[ev, far, late],
                    rtol=1e-9, atol=1e-11)
    if sol.t_events[0].size == 0:
        return None
    y = sol.y_events[0][0]
    return y[0], y[1]


target_r, target_phi = r0, np.pi - 0.3
res = {}
for sgn, name in [(+1, 'short (+)'), (-1, 'long (-)')]:
    phi_t = target_phi if sgn > 0 else -(2*np.pi - target_phi)
    g = lambda th: (shoot(th, phi_t) or (np.nan, np.nan))[1] - target_r
    ths = np.linspace(0.05, 3.09, 40) * (1 if sgn > 0 else -1)
    vals = [g(th) for th in ths]
    root = None
    for i in range(len(ths) - 1):
        if np.isfinite(vals[i]) and np.isfinite(vals[i+1]) and vals[i]*vals[i+1] < 0:
            root = brentq(g, ths[i], ths[i+1], xtol=1e-12)
            break
    res[name] = shoot(root, phi_t)[0] if root is not None else None
    print(f"  {name}: arrival eta = {res[name]}")
ok = all(v is not None for v in res.values()) and res['short (+)'] < res['long (-)']
check("first arrival = short arc; the long arc beyond the antipode is beaten", ok)
check("sup A on the interval used < Ehat",
      Afun(max(v for v in res.values() if v)) < Ev)

print()
if FAIL:
    print(f"{len(FAIL)} FAILED:", *FAIL, sep="\n  ")
    sys.exit(1)
print("all checks passed")
