"""Search side of the portable certificates: Bernstein subdivision that produces the proof trees.

This is the code that GENERATED the certificates in ../certificates. It is not needed to check them: the verifier
../certificates/verify_certificate.py is independent of this file and uses different algorithms.

Types of lemma:
  'radical'         P0 + P1 sqrt(g) >= 0 (chord and tangent envelopes); cells may be delegated to other certificates;
  'positive'        P > 0;
  'nonneg'          P >= 0;
  'alternative'     on every leaf one of the sufficient sets holds (list of indices of polynomials that are all >= 0);
  'radical_square'  P0 + P1 sqrt(g) >= 0 through Q = P1^2 g - P0^2: sufficient sets {P0,P1}, {P1,Q}, {P0,-Q}.
Tree in preorder: ["s", axis] halving; ["l", ...] leaf; ["d", region] delegated cell (the cell must lie in the region).

Two splitting rules were used, and both are kept so that the trees are reproduced exactly:
  rule 'widest'      split the widest side of the cell (groups 1, 2, 3);
  rule 'round_robin' split the axis that has been halved least often (groups 4, 5)."""
import math
import numpy as np, sympy as sp
from fractions import Fraction as Fr
from math import comb


def poly_json(P, VARS):
    p = sp.Poly(sp.expand(P), *VARS)
    return [[list(m), str(sp.Rational(c))] for m, c in p.terms()]


def bern_int(terms, degs):
    """Bernstein coefficients (tensor products of matrices), scaled to integers by a common denominator D."""
    A = np.empty([d + 1 for d in degs], dtype=object); A[...] = Fr(0)
    for m, c in terms: A[tuple(m)] += Fr(c)
    B = A
    for ax, n in enumerate(degs):
        M = np.empty((n + 1, n + 1), dtype=object); M[...] = Fr(0)
        for i in range(n + 1):
            for j in range(i + 1): M[i, j] = Fr(comb(i, j), comb(n, j))
        B = np.moveaxis(np.tensordot(M, np.moveaxis(B, ax, 0), axes=(1, 0)), 0, ax)
    D = 1
    for v in B.flat: D = D*v.denominator//math.gcd(D, v.denominator)
    return np.vectorize(lambda v: int(v*D), otypes=[object])(B), D


def split_int(N, ax):
    """de Casteljau halving along ax, in integers (both halves multiplied by 2^n)."""
    N = np.moveaxis(N, ax, 0); n = N.shape[0] - 1
    cur = [N[i] for i in range(n + 1)]
    L = [cur[0]*2**n]; R = [cur[n]*2**n]
    for k in range(1, n + 1):
        cur = [cur[i] + cur[i + 1] for i in range(len(cur) - 1)]
        L.append(cur[0]*2**(n - k)); R.append(cur[-1]*2**(n - k))
    return np.moveaxis(np.array(L), 0, ax), np.moveaxis(np.array(R[::-1]), 0, ax)


def sqrt_lo(q, k=60):
    return Fr(0) if q <= 0 else Fr(math.isqrt(q.numerator*4**k//q.denominator), 2**k)


def sqrt_hi(q, k=60):
    return Fr(math.isqrt(-(-q.numerator*4**k//q.denominator)) + 1, 2**k)


def mul_terms(T1, T2):
    out = {}
    for m1, c1 in T1:
        for m2, c2 in T2:
            m = tuple(a + b for a, b in zip(m1, m2)); out[m] = out.get(m, Fr(0)) + Fr(c1)*Fr(c2)
    return [[list(m), str(c)] for m, c in out.items() if c != 0]


def lin(T, coef):
    out = {}
    for t_, a_ in zip(T, coef):
        for m, c in t_: out[tuple(m)] = out.get(tuple(m), Fr(0)) + a_*Fr(c)
    return [[list(m), str(c)] for m, c in out.items() if c != 0]


def _radical_leaf(Ns, Ds, e):
    sc = [Fr(1, D*2**e) for D in Ds]
    g0 = max(Fr(int(Ns[3].min()))*sc[3], Fr(0)); g1 = Fr(int(Ns[3].max()))*sc[3]
    s0, s1 = sqrt_lo(g0), sqrt_lo(g1); sig = sqrt_hi((g0 + g1)/2) if g1 > 0 else Fr(1)
    Bc = (s1 - s0)/(g1 - g0) if g1 > g0 else Fr(0); Ac = s0 - Bc*g0
    for A_, B_ in ((Ac, Bc), (sig/2, 1/(2*sig))):
        c0, c1, c2 = sc[0], A_*sc[1], B_*sc[2]
        L = c0.denominator*c1.denominator*c2.denominator
        cm = Ns[0]*(c0.numerator*(L//c0.denominator)) + Ns[1]*(c1.numerator*(L//c1.denominator)) + Ns[2]*(c2.numerator*(L//c2.denominator))
        if int(cm.min()) < 0: return None
    return ['l', str(s0), str(s1), str(sig)]


def search_widest(typ, polys, d, sets=None, delegation=None, maxdep=40, maxcell=600000):
    """Rule 'widest'. polys: dict name -> monomials. delegation: list of (region, [(lo, hi) per axis]) in unit coordinates."""
    if typ == 'radical':
        keys = ['P0', 'P1', 'P1g', 'g']; T = [polys['P0'], polys['P1'], mul_terms(polys['P1'], polys['g']), polys['g']]
    elif typ == 'radical_square':
        Q = lin([mul_terms(mul_terms(polys['P1'], polys['P1']), polys['g']), mul_terms(polys['P0'], polys['P0'])], [Fr(1), Fr(-1)])
        keys = ['P0', 'P1', 'Q']; T = [polys['P0'], polys['P1'], Q]
        sets = [[0, 1], [1, 2], [0, -2]]          # -2 = -Q
    elif typ in ('positive', 'nonneg'):
        keys = ['P']; T = [polys['P']]
    else:
        keys = list(polys.keys()); T = [polys[k] for k in keys]
    degs = [max(max((m[i] for m, _ in t), default=0) for t in T) for i in range(d)]
    arrs = [bern_int(t, degs) for t in T]; Ds = [a[1] for a in arrs]
    tree = []; stack = [([a[0] for a in arrs], 0, [(Fr(0), Fr(1))]*d)]; ncell = 0; fails = 0
    def nonneg(N, i):
        return int(N[abs(i)].min()) >= 0 if i >= 0 else int(N[abs(i)].max()) <= 0
    while stack:
        Ns, e, cell = stack.pop(); ncell += 1
        if delegation:
            hit = next((nm for nm, reg in delegation if all(reg[k][0] <= cell[k][0] and cell[k][1] <= reg[k][1] for k in range(d))), None)
            if hit: tree.append(['d', hit]); continue
        if typ == 'positive':
            if int(Ns[0].min()) > 0: tree.append(['l']); continue
        elif typ == 'nonneg':
            if int(Ns[0].min()) >= 0: tree.append(['l']); continue
        elif typ in ('alternative', 'radical_square'):
            ok = next((si for si, S in enumerate(sets) if all(nonneg(Ns, i) for i in S)), None)
            if ok is not None: tree.append(['l', ok]); continue
        else:
            leaf = _radical_leaf(Ns, Ds, e)
            if leaf: tree.append(leaf); continue
        widths = [c[1] - c[0] for c in cell]
        if min(widths) < Fr(1, 2**maxdep) or ncell > maxcell: fails += 1; tree.append(['FAIL']); continue
        ax = widths.index(max(widths))
        parts = [split_int(N, ax) for N in Ns]
        mid = (cell[ax][0] + cell[ax][1])/2
        cl = list(cell); cr = list(cell); cl[ax] = (cell[ax][0], mid); cr[ax] = (mid, cell[ax][1])
        tree.append(['s', ax])
        stack.append(([p[1] for p in parts], e + degs[ax], cr)); stack.append(([p[0] for p in parts], e + degs[ax], cl))
    stats = dict(nodes=ncell, leaves=sum(1 for t in tree if t[0] == 'l'), delegations=sum(1 for t in tree if t[0] == 'd'), degrees=degs)
    return tree, stats, fails, keys


def search_round_robin(typ, polys, d, maxdep=40, maxcell=500000):
    """Rule 'round_robin' ('radical' or 'positive' only, no delegation)."""
    if typ == 'radical':
        T = [polys['P0'], polys['P1'], mul_terms(polys['P1'], polys['g']), polys['g']]
    else:
        T = [polys['P']]
    degs = [max(max((m[i] for m, _ in t), default=0) for t in T) for i in range(d)]
    arrs = [bern_int(t, degs) for t in T]; Ds = [a[1] for a in arrs]
    tree = []; stack = [([a[0] for a in arrs], 0, [0]*d)]; ncell = 0; fails = 0
    while stack:
        Ns, e, ks = stack.pop(); ncell += 1
        if typ == 'positive':
            if int(Ns[0].min()) > 0: tree.append(['l']); continue
        else:
            leaf = _radical_leaf(Ns, Ds, e)
            if leaf: tree.append(leaf); continue
        if max(ks) >= maxdep or ncell > maxcell: fails += 1; tree.append(['FAIL']); continue
        ax = ks.index(min(ks))
        parts = [split_int(N, ax) for N in Ns]
        tree.append(['s', ax])
        k2 = ks.copy(); k2[ax] += 1
        stack.append(([p[1] for p in parts], e + degs[ax], k2)); stack.append(([p[0] for p in parts], e + degs[ax], k2))
    stats = dict(nodes=ncell, leaves=sum(1 for t in tree if t[0] == 'l'), degrees=degs)
    return tree, stats, fails, None
