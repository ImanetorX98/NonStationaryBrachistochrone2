"""Bridge between the closed forms of polynomials.py and the curvature computed from the metric (SymPy).

group 1  exact identities U, V0, J^2 (orthonormal-frame formulas of the reference b = infinity) = closed forms;
group 2  exact identities U, V, J^2 (Ricci of G/N from Christoffel symbols, ZAMO frame) = closed forms, and H > 0 form;
group 3  Ricci of gbar_inf = -(1 - f u(eta)) (deta - beta)^2 + h0 with a time-dependent factor, from Christoffel symbols,
         split into R_stat + D1 u' + D2 u'' + D3 u'^2; checks R_stat = w N/w, D1 = c (al0 + al1 s)/w, D2 = -c^2 be/w,
         D3 = D2 f/(2w) at 40 random rational points, evaluated with 50 significant digits (the square roots make an exact
         SymPy simplification slow; the exact identities are checked in Wolfram, ../symbolic/bridge_group3.wls).
Groups 4 and 5 need no bridge: their polynomials are built directly from the curvature in build_certificates.py.
Exit code 0 if and only if every check passes."""
import sys, os, random
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import curvature, polynomials as PL
checks = []
def check(label, cond):
    checks.append(bool(cond)); print(('PASS ' if cond else 'FAIL ') + label, flush=True)

# ---------------------------------------------------------------- group 1
(U, V0, J), (r, q, a) = curvature.orthonormal_UVJ_b_inf()
g1 = PL.group1(); x, u, z = g1['x'], g1['u'], g1['z']
tox = lambda e: sp.simplify(e.subs(r, 2/x).subs(q, 1/u).subs(a, sp.sqrt(z)))
check('group 1: U = x^3 nU/(32(1-x)w)', sp.simplify(tox(U) - g1['U']) == 0)
check('group 1: V0 = x^3 nV/(64(1-x)w^2)', sp.simplify(tox(V0) - g1['V0']) == 0)
check('group 1: J^2 = z x^8 S L^2/(1024(1-x)^2 w)', sp.simplify(tox(J**2) - g1['J2']) == 0)

# ---------------------------------------------------------------- group 2
(U, V, J), (r, a, u) = curvature.ricci_GN_constant()
g2 = PL.group2(); x, R = g2['x'], g2['R']; a2 = R*(2 - R)
sub = lambda e: e.subs(r, R/x).subs(a, sp.sqrt(a2))
D, wt, H = g2['D'], g2['wt'], g2['H']
check('group 2: U = x^3 D P_U/(R^3 wt H)', sp.simplify(sub(U) - x**3*D*g2['P_U']/(R**3*wt*H)) == 0)
check('group 2: V = -x^3 P_V/(R^4 wt^2 H)', sp.simplify(sub(V) + x**3*g2['P_V']/(R**4*wt**2*H)) == 0)
check('group 2: J^2 = 4 x^8 (2-R) D P_J^2/(R^7 wt^3 H^2)', sp.simplify(sub(J**2) - 4*x**8*(2 - R)*D*g2['P_J']**2/(R**7*wt**3*H**2)) == 0)
check('group 2: H R/x^4 = r^4 + a^2 r^2 + 2 a^2 r + 4 a^2 u (so H > 0)',
      sp.simplify(H*R/x**4 - ((R/x)**4 + a2*(R/x)**2 + 2*a2*(R/x) + 4*a2*u)) == 0)

# ---------------------------------------------------------------- group 3
eta, r, ph, a = sp.symbols('eta r phi a', real=True)
u0, u1, u2, cc, ss = sp.symbols('u0 u1 u2 c s', real=True)
uf = sp.Function('u')(eta); X = [eta, r, ph]
f = 1 - 2/r; Del = r**2 - 2*r + a**2; V2 = 1 - f*uf
hrr = r**2/(f*Del); hpp = Del/f**2; bp = -2*a/(r*f)
g = sp.Matrix([[-V2, 0, V2*bp], [0, hrr, 0], [V2*bp, 0, hpp - V2*bp**2]])
gi = g.inv().applyfunc(sp.cancel); dg = [g.diff(v) for v in X]
Gam = [[[sp.cancel(sum(gi[i, l]*(dg[k][l, j] + dg[j][l, k] - dg[l][j, k]) for l in range(3))/2) for k in range(3)] for j in range(3)] for i in range(3)]
jets = lambda e: e.subs(sp.Derivative(uf, (eta, 2)), u2).subs(sp.Derivative(uf, eta), u1).subs(uf, u0)
Ric = sp.zeros(3)
for j in range(3):
    for k in range(j, 3):
        Ric[j, k] = Ric[k, j] = sp.cancel(jets(sum(sp.diff(Gam[i][j][k], X[i]) - sp.diff(Gam[i][j][i], X[k])
                        + sum(Gam[i][i][l]*Gam[l][j][k] - Gam[i][k][l]*Gam[l][j][i] for l in range(3)) for i in range(3))))
V20 = 1 - f*u0
e0 = sp.Matrix([1/sp.sqrt(V20), 0, 0]); er = sp.Matrix([0, 1/sp.sqrt(hrr), 0]); ep = sp.Matrix([bp, 0, 1])/sp.sqrt(hpp)
K = e0 + cc*er + ss*ep
RicK = sp.expand((K.T*Ric*K)[0])
Rst = RicK.subs({u1: 0, u2: 0}); c2 = RicK.coeff(u2, 1); c11 = RicK.coeff(u1, 2); c1 = RicK.coeff(u1, 1).subs(u2, 0)
check('group 3: Ric(K,K) is R_stat + D1 u\' + D3 u\'^2 + D2 u\'\' exactly (no other jet monomials)',
      sp.expand(RicK - (Rst + c1*u1 + c11*u1**2 + c2*u2)) == 0)
xx = 2/r; zz = a**2; ww = 1 - u0 + u0*xx; SS = xx**2*zz - 4*xx + 4; y2 = 1 - xx; vv = sp.sqrt(ww*SS)
g1s = {g1['x']: xx, g1['u']: u0, g1['z']: zz}
nU, nV, LL = g1['nU'].subs(g1s), g1['nV'].subs(g1s), g1['L'].subs(g1s)
TT = 3*zz*u0*xx**4 - 4*zz*u0*xx**3 + zz*u0*xx**2 + 4*zz*xx**3 - zz*xx**2 - 16*u0*xx**3 + 44*u0*xx**2 - 40*u0*xx + 12*u0 - 20*xx**2 + 32*xx - 12
wN = xx**3*nU/(32*y2)*cc**2 + xx**3*nV/(64*y2*ww)*ss**2 - a*xx**4*LL*vv/(32*y2)*ss
al0 = 3*zz*xx**5*sp.sqrt(y2)*ww/(8*vv); al1 = -a*xx**3*sp.sqrt(y2)*TT/(8*SS*ww); be = zz*xx**4*y2/(2*SS)
ids = {'R_stat = wN/w': Rst - wN/ww, 'D1 = c (al0 + al1 s)/w': c1 - cc*(al0 + al1*ss)/ww,
       'D2 = -c^2 be/w': c2 + cc**2*be/ww, 'D3 = D2 f/(2w)': c11 - c2*f/(2*ww)}
scale = {'R_stat = wN/w': Rst, 'D1 = c (al0 + al1 s)/w': c1, 'D2 = -c^2 be/w': c2, 'D3 = D2 f/(2w)': c11}
rng = random.Random(7); worst = {k: 0 for k in ids}
for _ in range(40):
    rv = 2 + sp.Rational(rng.randint(1, 4000), 100); av = sp.Rational(rng.randint(0, 99), 100); uv = sp.Rational(rng.randint(0, 49), 100)
    tv = sp.Rational(rng.randint(-100, 100), 101)
    pt = {r: rv, a: av, u0: uv, cc: (1 - tv**2)/(1 + tv**2), ss: 2*tv/(1 + tv**2)}
    for k, e in ids.items():
        den = abs(sp.N(scale[k].subs(pt), 50)) + sp.Rational(1, 10**30)
        worst[k] = max(worst[k], abs(sp.N(e.subs(pt), 50))/den)
for k, wv in worst.items():
    check(f'group 3: {k} (40 rational points, max relative residual {float(wv):.1e} < 1e-35)', wv < sp.Rational(1, 10**35))
print(f'SUMMARY: {sum(checks)}/{len(checks)} PASS')
sys.exit(0 if all(checks) else 1)
