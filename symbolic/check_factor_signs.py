"""Sign of the common factor F of the group-4 bridge (bridge_group4_group5.wls).

For each of the 12 group-4 certificates, bridge_group4_group5.wls finds P0_json/A0 = P1_json/A1 = F, where A0 + A1 vh is
the geometric inequality computed from the curvature. The certificate proves P0 + P1 sqrt(g) >= 0, so the inequality
follows if F > 0 on the box. This script proves it: F is a product of a rational constant and powers of polynomial
factors in (X, Rv, U, T); each non-monomial factor is shown to have strictly constant sign on [0,1]^4 by exact
Bernstein coefficients with bisection, the monomial factor X (x = R/r > 0) is positive for X > 0, and the sign of the
product is computed from the signs and multiplicities. Exact rational arithmetic throughout.
Usage: python3 check_factor_signs.py   (reads ../results/bridge_group4_factors.json)"""
import json, os, sys
from fractions import Fraction as Fr
from math import comb
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
data = json.load(open(os.path.join(HERE, '..', 'results', 'bridge_group4_factors.json')))
X, Rv, U, T = sp.symbols('X Rv U T'); V = (X, Rv, U, T)


def bernstein(P, degs):
    cur = {m: Fr(int(c.p), int(c.q)) for m, c in P.terms()}
    for ax, n in enumerate(degs):
        nxt = {}
        for m, c in cur.items():
            j = m[ax]
            for i in range(j, n + 1):
                mm = m[:ax] + (i,) + m[ax + 1:]
                nxt[mm] = nxt.get(mm, Fr(0)) + c*Fr(comb(i, j), comb(n, j))
        cur = nxt
    return cur


def strict_sign(expr, maxdepth=14):
    """+1 or -1 if expr has that strict sign on [0,1]^4 (Bernstein with bisection), else None."""
    P0 = sp.Poly(sp.expand(expr), *V)
    stack = [([(Fr(0), Fr(1))]*4, 0)]; sign = None; leaves = 0
    while stack:
        box, depth = stack.pop()
        sub = {v: lo + (hi - lo)*v for v, (lo, hi) in zip(V, box)}
        P = sp.Poly(sp.expand(P0.as_expr().subs({v: sp.Rational(lo.numerator, lo.denominator) + sp.Rational((hi - lo).numerator, (hi - lo).denominator)*v
                                                 for v, (lo, hi) in zip(V, box)}, simultaneous=True)), *V)
        degs = [P.degree(v) for v in V]
        b = list(bernstein(P, degs).values())
        if all(x > 0 for x in b) or all(x < 0 for x in b):
            s = 1 if b[0] > 0 else -1
            if sign is None: sign = s
            elif sign != s: return None
            leaves += 1; continue
        if depth >= maxdepth: return None
        widths = [hi - lo for lo, hi in box]; ax = widths.index(max(widths))
        lo, hi = box[ax]; mid = (lo + hi)/2
        b1 = list(box); b1[ax] = (lo, mid); b2 = list(box); b2[ax] = (mid, hi)
        stack += [(b1, depth + 1), (b2, depth + 1)]
    return sign


ok_all = True
for entry in data:
    sign_total = 1; details = []; ok = True
    for fs, mult in entry['factors']:
        f = sp.sympify(fs, locals={'X': X, 'Rv': Rv, 'U': U, 'T': T})
        if f.is_number:
            s = 1 if f > 0 else -1; details.append(f'const {s:+d}')
        elif f == X:
            s = 1; details.append(f'X^{mult} (+, X > 0)')
        else:
            s = strict_sign(f)
            if s is None:
                ok = False; details.append(f'UNDECIDED factor of degree {sp.Poly(f, *V).total_degree()}^{mult}')
                continue
            details.append(f'{"+" if s > 0 else "-"}^{mult}')
        sign_total *= s**mult
    ok = ok and sign_total == 1
    ok_all &= ok
    print(f"{'PASS' if ok else 'FAIL'} {entry['certificate']}: F > 0 on the box (X > 0); factors: {', '.join(details)}", flush=True)
print('SUMMARY:', 'all 12 factors F are positive on their boxes' if ok_all else 'at least one sign not proved')
sys.exit(0 if ok_all else 1)
