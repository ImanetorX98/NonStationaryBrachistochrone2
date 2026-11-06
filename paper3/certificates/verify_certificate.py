"""Independent verifier of portable positivity certificates (Paper III, supplementary material, section
"Curvature certificates on the equator").

It imports no file of the project and reads only the JSON certificates. For each certificate it
 1. rebuilds the polynomials from their monomials (exact rational coefficients);
 2. for 'radical' (P0 + P1 sqrt(g) >= 0): computes the product P1*g itself and checks g >= 0 (all Bernstein
    coefficients non-negative at the root);
 3. converts to the Bernstein basis with the explicit formula b_I = sum_{J <= I} a_J prod_k C(i_k, j_k)/C(n_k, j_k);
 4. walks the tree: every halving along an axis uses the direct binomial formulas
    left_k = 2^{-k} sum_{j <= k} C(k, j) b_j,  right_k = 2^{-(n-k)} sum_{j >= k} C(n-k, j-k) b_j;
    the tree, read in preorder, covers the box by construction (binary partition);
 5. on every leaf checks:
    - 'positive': all coefficients > 0;  'nonneg': all >= 0;
    - 'alternative': the leaf names a sufficient set (from the list declared in the certificate); every polynomial
      of the set has coefficients >= 0 (a negative index -i means polynomial i with its sign changed);
    - 'radical_square': the verifier computes Q = P1^2 g - P0^2 itself; fixed sufficient sets {P0,P1}, {P1,Q}, {P0,-Q}
      (P1 >= 0 and Q >= 0 give P1 sqrt(g) >= |P0|; P0 >= 0 and -Q >= 0 give P0 >= |P1| sqrt(g));
    - delegated cell ["d", name]: the cell (tracked by the verifier) must lie in the declared region;
    - 'radical': with g0 = max(0, min coeff of g), g1 = max coeff of g: 0 <= s0, s0^2 <= g0, 0 <= s1, s1^2 <= g1, sig > 0;
      chord vlo(g) = s0 + (s1 - s0)(g - g0)/(g1 - g0) (vlo = s0 if g1 = g0) and tangent vhi(g) = sig/2 + g/(2 sig);
      since F = P0 + P1 v is affine in v and vlo <= sqrt(g) <= vhi on [g0, g1], it suffices that the coefficients of
      P0 + A P1 + B P1 g are >= 0 for the (A, B) of both lines.
In addition:
 6. format validation: known type, non-negative integer exponents of the right dimension, no duplicate monomials,
    well-formed nodes, non-empty sufficient sets with valid indices, delegations only to declared regions;
 7. composition of delegations: for every delegated region the delegated certificates must be present in the same run
    and VERIFIED, must declare this main certificate and this region, and must satisfy exactly
      P0(phi) = m P0_del,  lam P1(phi) = m P1_del,  g(phi) = lam^2 g_del
    with m, lam monomials with positive coefficient (so F(phi) = m F_del >= 0 where m > 0; where m = 0, that is on the
    boundary of the box, by continuity). If the coverage is declared 'affine', the verifier checks that phi is a diagonal
    affine bijection of [0,1]^d onto the region; otherwise the coverage is an EXTERNAL OBLIGATION, printed at the end;
 8. unique names, acyclic delegation graph, composition in topological order (a certificate that delegates is accepted
    only if all its delegates are accepted as complete lemmas); field 'coverage' only 'affine' or 'external'; full
    validation of the composition fields (phi of the right length, monomials m and lam with non-negative integer
    exponents and positive coefficient). Per-file messages are LOCAL checks; the final outcome is the global one.
Outcome: VERIFIED means that the algebraic lemmas hold, up to the external obligations listed (logic of the
'alternative' certificates, non-affine coverings, bridge with the geometry).
Arithmetic: integers with a common denominator (exact). Usage: python3 verify_certificate.py file1.json [file2.json ...]"""
# Optional progress bar (../tools/progress.py, standard error only); without it this file is self-contained.
try:
    import os as _os, sys as _sys
    _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'tools'))
    from progress import Progress
except ImportError:
    class Progress:
        def __init__(self, total, label=''): pass
        def step(self, n=1): pass
        def close(self): pass
import sys, json, math, time
from fractions import Fraction as Fr
from math import comb
from itertools import product
def parse(terms, d):
    return {tuple(m): Fr(c) for m, c in terms}
def pmul(A, B):
    out = {}
    for ma, ca in A.items():
        for mb, cb in B.items():
            m = tuple(x + y for x, y in zip(ma, mb)); out[m] = out.get(m, Fr(0)) + ca*cb
    return {m: c for m, c in out.items() if c != 0}
def to_bern(P, degs):
    """dictionary monomial -> coeff  ==>  dictionary multi-index -> Bernstein coefficient (Fr), explicit formula axis by axis."""
    d = len(degs)
    cur = dict(P)
    for ax in range(d):
        n = degs[ax]; nxt = {}
        for m, c in cur.items():
            j = m[ax]
            for i in range(j, n + 1):
                mm = m[:ax] + (i,) + m[ax + 1:]
                nxt[mm] = nxt.get(mm, Fr(0)) + c*Fr(comb(i, j), comb(n, j))
        cur = nxt
    full = {}
    for idx in product(*[range(n + 1) for n in degs]): full[idx] = cur.get(idx, Fr(0))
    return full
def to_int(B):
    D = 1
    for v in B.values(): D = D*v.denominator//math.gcd(D, v.denominator)
    return {k: int(v*D) for k, v in B.items()}, D
def halve(B, degs, ax):
    """halving along ax with the direct binomial formulas; returns (L, R) multiplied by 2^n."""
    n = degs[ax]; L = {}; R = {}
    lines = {}
    for idx, v in B.items():
        key = idx[:ax] + idx[ax + 1:]; lines.setdefault(key, [0]*(n + 1))[idx[ax]] = v
    for key, b in lines.items():
        for k in range(n + 1):
            sl = sum(comb(k, j)*b[j] for j in range(k + 1))*2**(n - k)
            sr = sum(comb(n - k, j - k)*b[j] for j in range(k, n + 1))*2**k
            idx = key[:ax] + (k,) + key[ax:]
            L[idx] = sl; R[idx] = sr
    return L, R
TYPES = ('positive', 'nonneg', 'radical', 'radical_square', 'alternative')
def validate(C):
    typ = C.get('type'); d = len(C.get('variables', []))
    if typ not in TYPES: return 'unknown type'
    for name, T in C['polynomials'].items():
        seen = set()
        for m, c in T:
            if len(m) != d or any((not isinstance(e, int)) or e < 0 for e in m): return f'invalid exponents in {name}'
            if tuple(m) in seen: return f'duplicate monomial in {name}'
            seen.add(tuple(m)); Fr(c)
    nreq = {'positive': ['P'], 'nonneg': ['P'], 'radical': ['P0', 'P1', 'g'], 'radical_square': ['P0', 'P1', 'g']}.get(typ)
    if nreq and any(k not in C['polynomials'] for k in nreq): return 'missing polynomials'
    if typ == 'alternative':
        sets = C.get('sufficient_sets'); keys = C.get('polynomial_order')
        if not sets or any(len(S_) == 0 for S_ in sets): return 'empty sufficient set'
        if any(not isinstance(i, int) or abs(i) >= len(keys) for S_ in sets for i in S_): return 'invalid index'
    regions = [r[0] for r in C.get('delegations', [])]
    nsets = len(C['sufficient_sets']) if typ == 'alternative' else 3
    for tok in C['tree']:
        if tok[0] == 's':
            if len(tok) != 2 or not isinstance(tok[1], int) or not 0 <= tok[1] < d: return 'invalid split node'
        elif tok[0] == 'd':
            if len(tok) != 2 or tok[1] not in regions: return 'delegation to an undeclared region'
        elif tok[0] == 'l':
            if typ in ('positive', 'nonneg') and len(tok) != 1: return 'invalid leaf'
            if typ == 'radical' and len(tok) != 4: return 'invalid leaf'
            if typ in ('alternative', 'radical_square') and (len(tok) != 2 or not isinstance(tok[1], int) or not 0 <= tok[1] < nsets): return 'invalid leaf'
        else: return 'unknown node'
    return None
def verify(path):
    t0 = time.time()
    C = json.load(open(path)); err = validate(C)
    if err: return False, 'invalid format: ' + err
    typ = C['type']; d = len(C['variables']); tree = C['tree']
    if any(t[0] == 'FAIL' for t in tree): return False, 'tree with unclosed leaves'
    sets = None
    if typ == 'radical':
        P0 = parse(C['polynomials']['P0'], d); P1 = parse(C['polynomials']['P1'], d); g = parse(C['polynomials']['g'], d)
        P1g = pmul(P1, g); polys = [P0, P1, P1g, g]
    elif typ == 'radical_square':
        P0 = parse(C['polynomials']['P0'], d); P1 = parse(C['polynomials']['P1'], d); g = parse(C['polynomials']['g'], d)
        Q = pmul(pmul(P1, P1), g)
        for m, c in pmul(P0, P0).items(): Q[m] = Q.get(m, Fr(0)) - c
        polys = [P0, P1, {m: c for m, c in Q.items() if c != 0}]; sets = [[0, 1], [1, 2], [0, -2]]
    elif typ == 'alternative':
        polys = [parse(C['polynomials'][k], d) for k in C['polynomial_order']]; sets = C['sufficient_sets']
    else:
        polys = [parse(C['polynomials']['P'], d)]
    regions = {r[0]: [(Fr(a), Fr(b)) for a, b in r[1]] for r in C.get('delegations', [])}
    degs = [max(max((m[i] for m in P), default=0) for P in polys) for i in range(d)]
    ints = [to_int(to_bern(P, degs)) for P in polys]
    Ns = [a[0] for a in ints]; Ds = [a[1] for a in ints]
    if typ in ('radical', 'radical_square') and min(to_int(to_bern(g, degs))[0].values()) < 0:
        return False, 'g not certified >= 0 at the root'
    pos = 0; nleaves = 0; ndelegated = 0
    def nonneg(Ns, i):
        return min(Ns[i].values()) >= 0 if i >= 0 else max(Ns[-i].values()) <= 0
    def rec(Ns, e, cell):
        nonlocal pos, nleaves, ndelegated
        if pos >= len(tree): raise ValueError('truncated tree')
        tok = tree[pos]; pos += 1
        if tok[0] == 's':
            ax = tok[1]; parts = [halve(N, degs, ax) for N in Ns]
            mid = (cell[ax][0] + cell[ax][1])/2
            cl = list(cell); cr = list(cell); cl[ax] = (cell[ax][0], mid); cr[ax] = (mid, cell[ax][1])
            return rec([p[0] for p in parts], e + degs[ax], cl) and rec([p[1] for p in parts], e + degs[ax], cr)
        if tok[0] == 'd':
            ndelegated += 1; reg = regions[tok[1]]
            return all(reg[k][0] <= cell[k][0] and cell[k][1] <= reg[k][1] for k in range(d))
        nleaves += 1
        if typ == 'positive':
            return min(Ns[0].values()) > 0
        if typ == 'nonneg':
            return min(Ns[0].values()) >= 0
        if typ in ('alternative', 'radical_square'):
            S = sets[tok[1]]
            return all(nonneg(Ns, i) for i in S)
        s0, s1, sig = Fr(tok[1]), Fr(tok[2]), Fr(tok[3])
        sc = [Fr(1, D*2**e) for D in Ds]
        g0 = max(Fr(min(Ns[3].values()))*sc[3], Fr(0)); g1 = Fr(max(Ns[3].values()))*sc[3]
        if not (s0 >= 0 and s0*s0 <= g0 and s1 >= 0 and s1*s1 <= g1 and sig > 0): return False
        Bc = (s1 - s0)/(g1 - g0) if g1 > g0 else Fr(0); Ac = s0 - Bc*g0
        for A_, B_ in ((Ac, Bc), (sig/2, 1/(2*sig))):
            c0, c1, c2 = sc[0], A_*sc[1], B_*sc[2]
            L = c0.denominator*c1.denominator*c2.denominator
            m0, m1, m2 = c0.numerator*(L//c0.denominator), c1.numerator*(L//c1.denominator), c2.numerator*(L//c2.denominator)
            if min(Ns[0][k]*m0 + Ns[1][k]*m1 + Ns[2][k]*m2 for k in Ns[0]) < 0: return False
        return True
    sys.setrecursionlimit(10000)
    try:
        ok = rec(Ns, 0, [(Fr(0), Fr(1))]*d)
    except (ValueError, IndexError, KeyError) as ex:
        return False, f'malformed tree: {ex}'
    if not ok: return False, f'leaf rejected (after {nleaves} leaves, node {pos})'
    if pos != len(tree): return False, f'tree with excess nodes ({pos} used out of {len(tree)})'
    return ok, f"leaves {nleaves}, delegated {ndelegated}, nodes {len(tree)}, degrees {degs}, {time.time() - t0:.1f}s"
def pcompose(P, phi, d_del):
    """P (dictionary, variables of the main certificate) composed with phi (list of dictionaries in the delegate's variables)."""
    cache = {}
    def pw(i, k):
        if (i, k) not in cache:
            cache[(i, k)] = {tuple([0]*d_del): Fr(1)} if k == 0 else pmul(pw(i, k - 1), phi[i])
        return cache[(i, k)]
    out = {}
    for m, c in P.items():
        term = {tuple([0]*d_del): c}
        for i, k in enumerate(m):
            if k: term = pmul(term, pw(i, k))
        for mm, cc in term.items(): out[mm] = out.get(mm, Fr(0)) + cc
    return {m: c for m, c in out.items() if c != 0}
def scale(P, Q):   # product of polynomials (dictionaries)
    return pmul(P, Q)
def equal(A, B):
    keys = set(A) | set(B)
    return all(A.get(k, Fr(0)) == B.get(k, Fr(0)) for k in keys)
def positive_monomial(T):
    return len(T) == 1 and Fr(T[0][1]) > 0
def composition(certs, outcomes):
    """checks the delegations; returns (ok, messages, external obligations)."""
    ok = True; msgs = []; obligations = []
    for name, C in certs.items():
        for reg, box, delegates in C.get('delegations', []):
            boxF = [(Fr(a), Fr(b)) for a, b in box]
            d = len(C['variables'])
            P = {k: parse(C['polynomials'][k], d) for k in ('P0', 'P1', 'g')}
            for dn in delegates:
                if dn not in certs: ok = False; msgs.append(f'{name}/{reg}: delegate {dn} missing'); continue
                if not outcomes.get(dn): ok = False; msgs.append(f'{name}/{reg}: delegate {dn} not verified'); continue
                D = certs[dn]; comp = D.get('composition')
                if not comp or comp['main'] != name or comp['region'] != reg:
                    ok = False; msgs.append(f'{name}/{reg}: {dn} does not declare this composition'); continue
                dd = len(D['variables'])
                if comp.get('coverage') not in ('affine', 'external'):
                    ok = False; msgs.append(f'{dn}: coverage value not allowed: {comp.get("coverage")!r}'); continue
                if not isinstance(comp.get('phi'), list) or len(comp['phi']) != d:
                    ok = False; msgs.append(f'{dn}: phi of wrong length'); continue
                bad = False
                for T in comp['phi'] + [comp['m'], comp['lam']]:
                    for mm, cc in T:
                        if len(mm) != len(D['variables']) or any((not isinstance(e, int)) or e < 0 for e in mm): bad = True
                        Fr(cc)
                if bad: ok = False; msgs.append(f'{dn}: invalid exponents in the composition data'); continue
                if not (positive_monomial(comp['m']) and positive_monomial(comp['lam'])):
                    ok = False; msgs.append(f'{dn}: m or lam not positive monomials'); continue
                phi = [parse(T, dd) for T in comp['phi']]; m = parse(comp['m'], dd); lam = parse(comp['lam'], dd)
                Dp = {k: parse(D['polynomials'][k], dd) for k in ('P0', 'P1', 'g')}
                c0 = equal(pcompose(P['P0'], phi, dd), pmul(m, Dp['P0']))
                c1 = equal(pmul(lam, pcompose(P['P1'], phi, dd)), pmul(m, Dp['P1']))
                c2 = equal(pcompose(P['g'], phi, dd), pmul(pmul(lam, lam), Dp['g']))
                if not (c0 and c1 and c2):
                    ok = False; msgs.append(f'{dn}: composition identities failed (P0 {c0}, P1 {c1}, g {c2})'); continue
                if comp['coverage'] == 'affine':
                    aff = True
                    for i, ph in enumerate(phi):
                        for mm in ph:
                            if sum(mm) > 1 or (sum(mm) == 1 and mm[i] != 1): aff = False
                        lo = ph.get(tuple([0]*dd), Fr(0)); e_i = tuple(1 if j == i else 0 for j in range(dd)); sl = ph.get(e_i, Fr(0))
                        if (min(lo, lo + sl), max(lo, lo + sl)) != boxF[i]: aff = False
                    if not aff: ok = False; msgs.append(f'{dn}: affine coverage declared but not verified')
                    else: msgs.append(f'{name}/{reg}: covered by {dn} (affine bijection verified; exact identities verified)')
                else:
                    msgs.append(f'{name}/{reg}: {dn} composition verified (exact identities)')
            if any(certs.get(dn, {}).get('composition', {}).get('coverage') == 'external' for dn in delegates):
                obligations.append(f'covering of region {reg} of {name} by {delegates} (non-affine maps): mathematical argument')
            if any(any(sum(mm) > 0 for mm, _ in certs.get(dn, {}).get('composition', {}).get('m', [])) for dn in delegates):
                obligations.append(f'{name}/{reg}: at the points with m = 0 (boundary of the box) the inequality follows by continuity')
    return ok, msgs, obligations
def acyclic_graph(certs):
    """DFS on the graph main -> delegates; returns (ok, topological order from delegates to main certificates, message)."""
    state = {}; order = []
    def children(n): return [dn for r in certs[n].get('delegations', []) for dn in (r[2] if len(r) > 2 else [])]
    def dfs(n, stack):
        if state.get(n) == 1: return f'delegation cycle: {" -> ".join(stack + [n])}'
        if state.get(n) == 2: return None
        state[n] = 1
        for f in children(n):
            if f not in certs: continue           # a missing delegate is reported by the composition step
            e = dfs(f, stack + [n])
            if e: return e
        state[n] = 2; order.append(n); return None
    for n in certs:
        e = dfs(n, [])
        if e: return False, [], e
    return True, order, ''
if __name__ == '__main__':
    all_ok = True; certs = {}; outcomes = {}; obligations = []
    pg = Progress(len(sys.argv[1:]), 'certificates verified')
    for p in sys.argv[1:]:
        C = json.load(open(p))
        if C.get('name') in certs:
            all_ok = False; print(f'REJECTED   {p}: duplicate name {C.get("name")!r}'); continue
        ok, msg = verify(p); all_ok &= ok
        certs[C['name']] = C; outcomes[C['name']] = ok
        if C['type'] == 'alternative': obligations.append(f"{C['name']}: the deduction 'sufficient set => claim' (field logic)")
        print(f"{'VERIFIED (local)' if ok else 'REJECTED (local)'} {p}: {msg}", flush=True)
        pg.step()
    pg.close()
    okg, order, msgg = acyclic_graph(certs)
    if not okg:
        print('COMPOSITION:', msgg); print('OUTCOME: at least one certificate rejected or composition failed'); sys.exit(1)
    # composition in topological order: the complete outcome of a certificate requires that of its delegates
    complete = {}
    for n in order:
        dele = [dn for r in certs[n].get('delegations', []) for dn in (r[2] if len(r) > 2 else [])]
        complete[n] = outcomes[n] and all(complete.get(dn, False) for dn in dele)
    outcomes = complete
    for name, C in certs.items():     # every delegated region must declare at least one delegate
        for r in C.get('delegations', []):
            if len(r) < 3 or not r[2]: all_ok = False; print(f'COMPOSITION: {name}/{r[0]}: no delegated certificate declared')
    okc, msgs, obl = composition(certs, outcomes); all_ok &= okc; obligations += obl
    for m_ in msgs: print('COMPOSITION:', m_)
    print('EXTERNAL OBLIGATIONS (not checked by the verifier):')
    print('  - bridge polynomials <-> geometry (construction identities, see README)')
    for o in obligations: print('  -', o)
    print('OUTCOME:', 'all certificates verified' if all_ok else 'at least one certificate rejected or composition failed'); sys.exit(0 if all_ok else 1)
