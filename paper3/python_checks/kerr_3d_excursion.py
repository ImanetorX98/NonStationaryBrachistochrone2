# -*- coding: utf-8 -*-
"""Main text, section "Conformal Kerr: a barrier, and small spin": 3D Kerr data of the auxiliary metric and the
excursion margins (Lemma 2: the polar minimum of the cost density is at the equator, so the 3D margins equal the
equatorial ones). Usage: python3 kerr_3d_excursion.py"""
import sys
import numpy as np
import sympy as sp
from scipy.integrate import quad

FAIL = []


def check(name, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"   {detail}" if detail else ""))
    if not ok:
        FAIL.append(name)


r, th, M, a, s = sp.symbols('r theta M a s', positive=True)
S = r**2 + a**2*sp.cos(th)**2; D = r**2 - 2*M*r + a**2
gtt = -(1 - 2*M*r/S); gtp = -2*M*a*r*sp.sin(th)**2/S
gpp = (r**2 + a**2 + 2*M*a**2*r*sp.sin(th)**2/S)*sp.sin(th)**2
f = -gtt
print("symbolic")
hpp = gpp + gtp**2/f
check("h_phiphi = Delta sin^2/f", sp.simplify(hpp - D*sp.sin(th)**2/f) == 0)
b2 = (gtp/f)**2/(hpp/f)
check("|beta|_h0^2 = 4M^2a^2r^2 sin^2/(Sigma^2 Delta)",
      sp.simplify(b2 - 4*M**2*a**2*r**2*sp.sin(th)**2/(S**2*D)) == 0)
x = sp.symbols('x', positive=True)          # x = sin(theta) in (0,1]
ratio = x/(r**2 + a**2*(1 - x**2))
check("sin/Sigma increasing in sin(theta): max at equator",
      sp.simplify(sp.diff(ratio, x) - (r**2 + a**2 + a**2*x**2)/(r**2 + a**2*(1 - x**2))**2) == 0)
check("h0_rr = Sigma^2/(Delta (Sigma - 2Mr))",
      sp.simplify((S/D)/f - S**2/(D*(S - 2*M*r))) == 0)
# Lemma 2: G(x) = (Sigma - k x)/sqrt(Sigma - 2Mr) strictly decreasing in x = sin(theta)
X = sp.symbols('X', positive=True)
SX = r**2 + a**2 - a**2*X**2; k = 2*M*a*r/sp.sqrt(D)
G = (SX - k*X)/sp.sqrt(SX - 2*M*r)
check("Lemma 2: (Sigma-2Mr)^{3/2} G' = a^2 X(4Mr-Sigma) - k Delta",
      sp.simplify(sp.diff(G, X)*(SX - 2*M*r)**sp.Rational(3, 2)
                  - (a**2*X*(4*M*r - SX) - k*D)) == 0)
check("Lemma 2: 4M^2 Delta - a^2(4M-r)^2 = s(8M^3+4Ma^2) + s^2(4M^2-a^2) at r=2M+s",
      sp.expand((4*M**2*D - a**2*(4*M - r)**2).subs(r, 2*M + s)
                - (s*(8*M**3 + 4*M*a**2) + s**2*(4*M**2 - a**2))) == 0)

print("\nnumerical")
Mv = 1.0
TH = np.linspace(1e-6, np.pi/2, 4001)


def c3_split(rv, av):
    Dv = rv*rv - 2*rv + av*av
    g = lambda Sv: Sv*Sv/(Sv - 2*Mv*rv)
    cands = [rv*rv, rv*rv + av*av] + ([4*Mv*rv] if rv*rv <= 4*Mv*rv <= rv*rv + av*av else [])
    return (1 - 2*abs(av)/(rv*np.sqrt(Dv)))*np.sqrt(min(g(x) for x in cands)/Dv)


def c3_joint(rv, av):
    Sv = rv*rv + av*av*np.cos(TH)**2; Dv = rv*rv - 2*rv + av*av; fv = 1 - 2*rv/Sv
    bv = 2*abs(av)*rv*np.sin(TH)/(Sv*np.sqrt(Dv))
    val = (1 - bv)*np.sqrt(Sv/(Dv*fv))
    return val.min(), TH[val.argmin()]


eq_min = all(abs(c3_joint(rv, av)[1] - np.pi/2) < 1e-3
             for av in (0.3, 0.6, 0.9, 0.99) for rv in np.linspace(2.001, 30, 200))
check("joint polar minimum at the equator on 800 (a,r) samples (consistent with Lemma 2)", eq_min)
split_below = all(c3_split(rv, av) <= c3_joint(rv, av)[0] + 1e-12
                  for av in (0.3, 0.6, 0.9, 0.99) for rv in np.linspace(2.001, 30, 200))
check("split bound <= joint minimum (consistency)", split_below)
def c_eq(rv, av):
    Dv = rv*rv - 2*rv + av*av; fv = 1 - 2/rv
    return (1 - 2*abs(av)/(rv*np.sqrt(Dv)))*rv/np.sqrt(fv*Dv)


for av in (0.5, 0.9, 0.99):
    print(f"  a={av} (Lemma 2, 3D = equatorial): " + ", ".join(
        f"r2={r2}: {2*quad(lambda x: c_eq(x, av), 2 + 1e-10, r2, limit=300)[0]:.2f}"
        for r2 in (4, 6, 8, 12)))
# direct test of the cost inequality on random 3D null vectors
rng = np.random.default_rng(3); worst = np.inf
for _ in range(20000):
    av = rng.uniform(-1, 1); rv = 2 + 10**rng.uniform(-4, 1.5); t = rng.uniform(0.01, np.pi - 0.01)
    Nv = 1 + 10**rng.uniform(-3, 2)
    Sv = rv*rv + av*av*np.cos(t)**2; Dv = rv*rv - 2*rv + av*av; fv = 1 - 2*rv/Sv
    h0 = np.array([Sv/(Dv*fv), Sv/fv, Dv*np.sin(t)**2/fv**2])
    bph = -2*av*rv*np.sin(t)**2/(Sv*fv)
    v = rng.normal(size=3)
    ve = np.sqrt(Nv)*np.sqrt(h0 @ v**2) + bph*v[2]
    worst = min(worst, ve - c_eq(rv, av)*abs(v[0]))
check("Lemma 2 bound: v^eta >= c_eq(r)|v^r| on 20000 random 3D null vectors",
      worst >= -1e-9, f"min slack {worst:.2e}")

print()
if FAIL:
    print(f"{len(FAIL)} FAILED:", *FAIL, sep="\n  ")
    sys.exit(1)
print("all checks passed")
