"""Rebuild the 33 curvature certificates from the metric and compare them with ../certificates.

Default mode: compute the curvature (curvature.py), the polynomials (polynomials.py) and their box versions, and check
that every polynomial, delegation, composition datum and sufficient set equals EXACTLY the one in the published
certificate (monomial by monomial, rational coefficients). About ten minutes.
--search: also rerun the Bernstein subdivision (search.py) and check that every proof tree and its statistics are
reproduced exactly. Most of the time goes to group4_B4_I1 (about ten minutes); about half an hour in all.
--write DIR: with --search, also write the regenerated certificates to DIR.
--only PREFIX[,PREFIX...]: restrict to certificates whose name starts with one of the prefixes (e.g. group1,group4_B1).
Exit code 0 if and only if everything matches."""
import sys, os, json, time, argparse
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'tools'))
from progress import Progress
import sympy as sp
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import curvature, polynomials as PL
from search import poly_json, search_widest, search_round_robin

ap = argparse.ArgumentParser()
ap.add_argument('--search', action='store_true'); ap.add_argument('--write'); ap.add_argument('--only')
args = ap.parse_args()
CERT = os.path.join(HERE, '..', 'certificates')
only = args.only.split(',') if args.only else None
want = lambda name: only is None or any(name.startswith(p) for p in only)
results = []; t_start = time.time()
pg = Progress(sum(want(f[:-5]) for f in os.listdir(CERT) if f.startswith('group') and f.endswith('.json')), 'certificates rebuilt')


def as_dict(terms): return {tuple(m): Fr(c) for m, c in terms}


def compare(name, typ, polys, variables, rule, sets=None, delegation=None, composition=None, delegations_field=None):
    """compare with the published certificate; optionally rerun the search."""
    pub = json.load(open(os.path.join(CERT, name + '.json')))
    ok = pub['type'] == typ and pub['variables'] == variables
    ok &= set(pub['polynomials']) == set(polys) and all(as_dict(pub['polynomials'][k]) == as_dict(polys[k]) for k in polys)
    if sets is not None: ok &= pub.get('sufficient_sets') == sets
    if composition is not None: ok &= pub.get('composition') == composition
    if delegations_field is not None: ok &= pub.get('delegations') == delegations_field
    msg = 'polynomials and data identical' if ok else 'DIFFERENT polynomials or data'
    if args.search:
        t0 = time.time()
        if rule == 'widest':
            tree, stats, fails, keys = search_widest(typ, polys, len(variables), sets=sets, delegation=delegation)
        else:
            tree, stats, fails, keys = search_round_robin(typ, polys, len(variables))
        same_tree = (tree == pub['tree']) and stats == pub['statistics'] and fails == 0
        if typ == 'alternative': same_tree &= keys == pub['polynomial_order']
        ok &= same_tree
        msg += f"; tree {'identical' if same_tree else 'DIFFERENT'} ({stats['leaves']} leaves, {time.time() - t0:.0f}s)"
        if args.write:
            new = dict(pub); new['tree'] = tree; new['statistics'] = stats; new['polynomials'] = polys
            os.makedirs(args.write, exist_ok=True)
            json.dump(new, open(os.path.join(args.write, name + '.json'), 'w'), separators=(',', ':'))
    results.append((name, ok))
    print(f"{'MATCH   ' if ok else 'MISMATCH'} {name}: {msg}", flush=True)
    pg.step()


# ------------------------------------------------------------------------------ group 1: box (X, U, Z)
if want('group1'):
    g = PL.group1()
    Xs, Us, Zs = sp.symbols('X U Z'); V = (Xs, Us, Zs)
    box = {g['x']: Xs, g['u']: Us/2, g['z']: sp.Rational(1089, 2500)*Zs}
    for name, P in [('minus_nU', -g['nU']), ('minus_nV', -g['nV']), ('P3', g['P3']), ('P6', g['P6'])]:
        compare('group1_' + name, 'nonneg', {'P': poly_json(sp.expand(P.subs(box, simultaneous=True)), V)}, ['X', 'U', 'Z'], 'widest')

# ------------------------------------------------------------------------------ group 2: box (X, Rv, U)
if want('group2'):
    g = PL.group2()
    Xs, Rs_, Us = sp.symbols('X Rv U'); V = (Xs, Rs_, Us)
    box = {g['x']: Xs, g['R']: 1 + Rs_, g['u']: sp.Rational(3, 5)*Us}
    polys = {k: poly_json(sp.expand(e.subs(box, simultaneous=True)), V) for k, e in
             [('PV', g['P_V']), ('P3', g['P3']), ('P4', g['P4']), ('P5', g['P5']), ('mPU', -g['P_U']), ('P6', g['P6'])]}
    compare('group2_Kinf', 'alternative', polys, ['X', 'Rv', 'U'], 'widest', sets=[[0, 1, 2], [0, 1, 3], [0, 1, 4, 5]])

# ------------------------------------------------------------------------------ group 3: box (Y, T, U, A), K1 = 5/4, K2 = 1/3
if want('group3'):
    g = PL.group3(); y, t, u, a = g['y'], g['t'], g['u'], g['a']
    K1, K2 = sp.Rational(5, 4), sp.Rational(1, 3); UMAX, AMAX = sp.Rational(2, 5), sp.Rational(33, 50)
    Y, T, U, A = sp.symbols('Y T U A'); V = (Y, T, U, A); VN = ['Y', 'T', 'U', 'A']
    WS = sp.expand(g['w']*g['S'])
    def mono(e, c): return [[e, str(c)]]
    def strip_common(Q0, Q1, lead):
        k = min(min(m[0] for m in sp.Poly(Q0, lead).monoms()), min(m[0] for m in sp.Poly(Q1, lead).monoms()))
        return sp.expand(sp.cancel(Q0/lead**k)), sp.expand(sp.cancel(Q1/lead**k)), k
    for nm, label in [('I1p', 'I1plus'), ('I1m', 'I1minus'), ('I2', 'I2')]:
        if not (want(f'group3_{label}') or want('group3')): continue
        ex = {'I1p': K1*g['mwN'] - (g['cal0'] + g['cal1s']), 'I1m': K1*g['mwN'] + (g['cal0'] + g['cal1s']),
              'I2': K2*g['mwN'] - g['c']**2*g['be']}[nm]
        R0, R1, den = g['lin_in_v'](ex); g_ = (y**2 - 1)**3
        P0 = sp.expand(-sp.cancel(R0/g_)); P1 = sp.expand(-sp.cancel(R1/g_))
        main = f'group3_{label}_main'
        # main box, with delegation of the corners C1 and C2
        box = {y: Y, t: 2*T - 1, u: UMAX*U, a: AMAX*A}
        pol = {'P0': poly_json(P0.subs(box, simultaneous=True), V), 'P1': poly_json(P1.subs(box, simultaneous=True), V),
               'g': poly_json(WS.subs(box, simultaneous=True), V)}
        C1box = [(Fr(0), Fr(1, 8)), (Fr(15, 16), Fr(1)), (Fr(0), Fr(1)), (Fr(0), Fr(1))]
        C2box = [(Fr(0), Fr(1, 8)), (Fr(0), Fr(1)), (Fr(0), Fr(1)), (Fr(0), Fr(1, 8))]
        dfield = [['C1', [[str(lo), str(hi)] for lo, hi in C1box], [f'group3_{label}_C1']],
                  ['C2', [[str(lo), str(hi)] for lo, hi in C2box], [f'group3_{label}_C2A', f'group3_{label}_C2B']]]
        compare(main, 'radical', pol, VN, 'widest', delegation=[('C1', C1box), ('C2', C2box)], delegations_field=dfield)
        # corner C1: affine image
        cb = {y: Y/8, t: 2*(sp.Rational(15, 16) + T/16) - 1, u: UMAX*U, a: AMAX*A}
        pol = {'P0': poly_json(P0.subs(cb, simultaneous=True), V), 'P1': poly_json(P1.subs(cb, simultaneous=True), V),
               'g': poly_json(WS.subs(cb, simultaneous=True), V)}
        comp = dict(main=main, region='C1', phi=[mono([1, 0, 0, 0], Fr(1, 8)), mono([0, 0, 0, 0], Fr(15, 16)) + mono([0, 1, 0, 0], Fr(1, 16)),
                    mono([0, 0, 1, 0], 1), mono([0, 0, 0, 1], 1)], m=mono([0, 0, 0, 0], 1), lam=mono([0, 0, 0, 0], 1), coverage='affine')
        compare(f'group3_{label}_C1', 'radical_square', pol, VN, 'widest', composition=comp)
        # sector C2A: a = 2 y kappa, kappa = A in [0,1], y = Y/8
        yy = Y/8; cb = {y: yy, t: 2*T - 1, u: UMAX*U, a: 2*yy*A}
        Q0 = sp.expand(P0.subs(cb, simultaneous=True)); Q1 = sp.expand((2*yy*P1).subs(cb, simultaneous=True))
        Q0, Q1, k = strip_common(Q0, Q1, Y)
        gh = sp.expand((g['w']*(A**2*(1 - y**2)**2 + 1)).subs({y: yy, u: UMAX*U}, simultaneous=True))
        comp = dict(main=main, region='C2', phi=[mono([1, 0, 0, 0], Fr(1, 8)), mono([0, 1, 0, 0], 1), mono([0, 0, 1, 0], 1),
                    mono([1, 0, 0, 1], Fr(25, 66))], m=mono([k, 0, 0, 0], 1), lam=mono([1, 0, 0, 0], Fr(1, 4)), coverage='external')
        compare(f'group3_{label}_C2A', 'radical_square', {'P0': poly_json(Q0, V), 'P1': poly_json(Q1, V), 'g': poly_json(gh, V)},
                VN, 'widest', composition=comp)
        # sector C2B: a = (AMAX/8) A, y = a eta/2, eta = Y
        aa = AMAX*A/8; cb = {a: aa, y: aa*Y/2, t: 2*T - 1, u: UMAX*U}
        Q0 = sp.expand(P0.subs(cb, simultaneous=True)); Q1 = sp.expand((aa*P1).subs(cb, simultaneous=True))
        Q0, Q1, k = strip_common(Q0, Q1, A)
        gh = sp.expand((g['w']*((1 - y**2)**2 + Y**2)).subs({y: aa*Y/2, u: UMAX*U}, simultaneous=True))
        comp = dict(main=main, region='C2', phi=[mono([1, 0, 0, 1], Fr(33, 800)), mono([0, 1, 0, 0], 1), mono([0, 0, 1, 0], 1),
                    mono([0, 0, 0, 1], Fr(1, 8))], m=mono([0, 0, 0, k], 1), lam=mono([0, 0, 0, 1], Fr(33, 400)), coverage='external')
        compare(f'group3_{label}_C2B', 'radical_square', {'P0': poly_json(Q0, V), 'P1': poly_json(Q1, V), 'g': poly_json(gh, V)},
                VN, 'widest', composition=comp)

# ------------------------------------------------------------------------------ group 4: box (X, Rv, U, T), bands B1-B4
if want('group4'):
    t0 = time.time()
    tables, syms = PL.group4_tables(curvature.ricci_h_jets(), curvature.hessian_r())
    x, Rr, u0 = syms
    print(f'  curvature tables of group 4 computed ({time.time() - t0:.0f}s)', flush=True)
    bands = {'B1': ('8/5', '44/25', '3/2', '7/2'), 'B2': ('143/100', '8/5', '3', '9'),
             'B3': ('131/100', '143/100', '5', '19'), 'B4': ('57/50', '131/100', '14', '110')}
    Xs, Rs_, Us, Ts = sp.symbols('X Rv U T'); V = (Xs, Rs_, Us, Ts)
    for band, vals in bands.items():
        if not want(f'group4_{band}') and not want('group4'): continue
        R_lo, R_hi, K1, K2 = [sp.Rational(v) for v in vals]
        res, t, ghat = PL.group4_assemble(tables, syms, K1, K2)
        xb = 10*R_lo/(11*R_lo - 1)
        box = {x: xb*Xs, Rr: R_lo + (R_hi - R_lo)*Rs_, u0: sp.Rational(2, 5)*Us, t: 2*Ts - 1}
        for nm in ['S', 'I2', 'I1']:
            if not want(f'group4_{band}_{nm}') and not want(f'group4_{band}') and not want('group4'): continue
            P0, P1 = res[nm]
            k = min(min(m[0] for m in sp.Poly(P0, x).monoms()), min(m[0] for m in sp.Poly(P1, x).monoms()))
            P0 = sp.expand(sp.cancel(P0/x**k)); P1 = sp.expand(sp.cancel(P1/x**k))
            polys = {'P0': poly_json(P0.subs(box, simultaneous=True), V), 'P1': poly_json(P1.subs(box, simultaneous=True), V),
                     'g': poly_json(ghat.subs(box, simultaneous=True), V)}
            compare(f'group4_{band}_{nm}', 'radical', polys, ['X', 'Rv', 'U', 'T'], 'round_robin')

# ------------------------------------------------------------------------------ group 5: box (Y, Rv, U)
if want('group5'):
    g = PL.group5()
    Ys, Rs_, Us = sp.symbols('Y Rv U'); V = (Ys, Rs_, Us)
    box = {g['y']: Ys, g['R']: sp.Rational(57, 50) + (sp.Rational(44, 25) - sp.Rational(57, 50))*Rs_, g['u']: sp.Rational(2, 5)*Us}
    sa, sq = g['sa'], g['sq']
    for name, P in [('den_alpha', sa*g['da']), ('minus_alpha', -sa*g['na']), ('den_Q', sq*g['dq']), ('Q', sq*g['nq'])]:
        compare('group5_' + name, 'positive', {'P': poly_json(sp.expand(P.subs(box, simultaneous=True)), V)}, ['Y', 'Rv', 'U'], 'round_robin')

pg.close()
n_ok = sum(ok for _, ok in results)
print(f'SUMMARY: {n_ok}/{len(results)} certificates rebuilt identically ({time.time() - t_start:.0f}s)')
sys.exit(0 if results and n_ok == len(results) else 1)
