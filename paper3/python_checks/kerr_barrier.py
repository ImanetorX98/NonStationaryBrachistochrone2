# -*- coding: utf-8 -*-
"""
Main text, section "Conformal Kerr: a barrier, and small spin": checks for the barrier theorem with rotation.

Symbolic: |beta|/sqrt(Q0) = 2M|a|/(r sqrt(Delta)); theta increasing and in (0,1); g_opt in lapse-shift form.
Numerical: the margins 2C(2M+, r2) quoted after the first-arrival theorem (16.7M, 12.9M, 12.4M at r2 = 6M for
a = 0.5M, 0.9M, 0.99M); direct test of the excursion-cost inequality on random null vectors.
Usage: python3 kerr_barrier.py
"""
import sys
import numpy as np
import sympy as sp
from scipy.integrate import quad

FAIL = []


def check(name, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"   {detail}" if detail else ""))
    if not ok:
        FAIL.append(name)


r, M, a, N = sp.symbols('r M a N', positive=True)
f = 1 - 2*M/r
D = r**2 - 2*M*r + a**2
P0 = r**2/(f*D); Q0 = D/f**2; beta = -2*M*a/(r*f)
print("symbolic")
ratio = sp.simplify(sp.Abs(beta)/sp.sqrt(Q0))
check("|beta|/sqrt(Q0) = 2Ma/(r sqrt(Delta))",
      sp.simplify(sp.Abs(beta)**2/Q0 - (2*M*a/(r*sp.sqrt(D)))**2) == 0)
g = 2*M*a/(r*sp.sqrt(D))
dg = sp.simplify(sp.diff(g, r))
check("2Ma/(r sqrt Delta) is decreasing in r (so theta increasing)",
      sp.simplify(dg + 2*M*a*(2*r**2 - 3*M*r + a**2)/(r**2*D**sp.Rational(3, 2))) == 0,
      "d/dr = -2Ma(2r^2-3Mr+a^2)/(r^2 Delta^{3/2}), and 2r^2-3Mr+a^2 > 0 for r>2M")
# 2r^2-3Mr+a^2 at r=2M+s: 2(2M+s)^2-3M(2M+s)+a^2 = 2M^2+5Ms+2s^2+a^2 > 0
s = sp.symbols('s', positive=True)
check("2r^2-3Mr+a^2 > 0 for r>2M", sp.expand((2*r**2 - 3*M*r + a**2).subs(r, 2*M + s))
      == 2*M**2 + 5*M*s + 2*s**2 + a**2)
check("theta(2M+) = 0 (Finsler degeneracy)", sp.limit(1 - g, r, 2*M, '+') == 0)
# lapse-shift form of -(d eta - beta dPhi)^2 + N(P0 dr^2 + Q0 dPhi^2)
# g_PhiPhi = N Q0 - beta^2, g_etaPhi = beta, g_etaeta = -1
gPP = N*Q0 - beta**2
lapse2 = sp.simplify(1 + beta**2/gPP)
check("lapse^2 = N Q0/(N Q0 - beta^2)", sp.simplify(lapse2 - N*Q0/gPP) == 0)
check("shift^Phi = beta/(N Q0 - beta^2) reproduces g_etaPhi",
      sp.simplify(gPP*(beta/gPP) - beta) == 0)

print("\nnumerical")
Mv = 1.0


def integrand(rv, av):
    fv = 1 - 2*Mv/rv; Dv = rv*rv - 2*Mv*rv + av*av
    return (1 - 2*Mv*abs(av)/(rv*np.sqrt(Dv)))*np.sqrt(rv*rv/(fv*Dv))


table = {}
for av in (0.5, 0.9, 0.99):
    for r2 in (4, 6, 8, 12):
        table[(av, r2)] = 2*quad(lambda x: integrand(x, av), 2*Mv + 1e-12, r2, limit=400)[0]
    print(f"  a={av}: " + ", ".join(f"r2={r2}: {table[(av, r2)]:.2f}" for r2 in (4, 6, 8, 12)))
check("margins finite and positive for a != 0", all(np.isfinite(v) and v > 0 for v in table.values()))
# a=0: the budget diverges as r1 -> 2M
b0 = [2*quad(lambda x: integrand(x, 0.0), 2 + 10**-k, 6)[0] for k in (2, 4, 6)]
check("a=0: budget grows without bound as r1 -> 2M", b0[0] < b0[1] < b0[2],
      f"{[round(x, 2) for x in b0]}")
# direct test of inequality (1) on random null vectors
rng = np.random.default_rng(1)
worst = np.inf
for _ in range(20000):
    av = rng.uniform(-1, 1); rv = 2 + 10**rng.uniform(-4, 1.5); Nv = 1 + 10**rng.uniform(-3, 2)
    fv = 1 - 2/rv; Dv = rv*rv - 2*rv + av*av
    P0v = rv*rv/(fv*Dv); Q0v = Dv/fv**2; bv = -2*av/(rv*fv)
    vr, vp = rng.normal(size=2)
    h = np.sqrt(P0v*vr**2 + Q0v*vp**2)
    ve = np.sqrt(Nv)*h + bv*vp          # future null
    th = 1 - 2*abs(av)/(rv*np.sqrt(Dv))
    worst = min(worst, ve - th*np.sqrt(P0v)*abs(vr))
check("inequality (1): v^eta >= theta sqrt(P0)|v^r| on 20000 random null vectors",
      worst >= -1e-9, f"min slack {worst:.2e}")

print()
if FAIL:
    print(f"{len(FAIL)} FAILED:", *FAIL, sep="\n  ")
    sys.exit(1)
print("all checks passed")
