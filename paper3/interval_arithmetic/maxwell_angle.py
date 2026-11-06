"""maxwell_angle.py -- interval certificate of an interior Maxwell angle of the proper-time branch at equal radii
(A = 1, a = 0, M = 1), q = 1.001, r0 = 1500 (supplementary material, "The interior Maxwell angle").

Certified statement: there is delta* in (delta_a, delta_b) such that at separation delta* exactly two interior
geodesics, the dives with periastron in N1 and in N2, are global minima in the completion (the vertical dive and
every other dive are strictly longer).

Notation: for a dive of periastron t, L = rho(t), Delta(t), ell(t), S = ell - L Delta, F_delta = S + delta L;
d(delta) = min_{L in [0, rho0]} F_delta(L); dS/dt = -Delta dL/dt, dF_delta/dt = (delta - Delta) dL/dt.

Steps (all in interval arithmetic, mpmath.iv, q an interval that contains 1.001 both decimal and binary):
 A. on N_k = [a_k, b_k]: Delta' < 0 (jets), Delta(a_k) > delta_b, Delta(b_k) < delta_a. For every delta in
    [delta_a, delta_b] there is a single dive L_k(delta) in N_k, and it is the minimum of F_delta on N_k (F'' = -Delta' > 0
    in the variable L).
 B. ell_k(delta) increases with delta (d ell_k/d delta = L_k > 0). ell_k(delta_a), ell_k(delta_b) are enclosed by
    locating the crossing Delta = delta by bisection and enclosing ell on the box.
    G(delta) = ell_2 - ell_1: G(delta_a) < 0 < G(delta_b), G' = L_2 - L_1 > 0 -> a single delta* with G = 0.
 C. U = min_k ell_k(delta_b) (sup of ell_k over the interval) for the edges and the vertical dive. In between, outside
    N1 u N2, for every delta in the interval F_delta(L) > ell_k(delta) for some k, with
    d/d delta [F_delta(L) - ell_k(delta)] = L - L_k(delta): if L >= max L_k the test at delta_a suffices, if
    L <= min L_k the one at delta_b. Edge t in (2, t_lo]: S = 2R0 - int_0^L Delta >= 2R0 - pi L, F >= 2R0 - pi L(t_lo).
    Far, t >= t_far: S >= 0, F >= delta_a L(t_far). In between: bisection with
    F(T) in F(m) + (delta - Delta(T)) L'(T) (T - m), F(m) = ell(m) - L(m) Delta(m) + delta L(m).
    Vertical dive: 2R0 > U.
 Conclusion: at delta*, min F = ell_1 = ell_2 is attained only in N1 and N2.
Usage (from this folder): python3 maxwell_angle.py ; 10 processes, about 75 minutes. Archive: leaves_maxwell/."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'tools'))
from progress import Progress
import sys, time, json, os, hashlib
from multiprocessing import Pool
from mpmath import iv, mp, mpf
iv.dps = 20; mp.dps = 25
from jets_nested import make, integrate, R0_enclosure, qint, down, up, manifest
from jets_t import J, jsqrt

QS, R0 = '1.001', 1500.0
Qi = qint(QS)
GA = make(Qi, R0, 'angle'); GL = make(Qi, R0, 'length')
PI = iv.pi

def mpf_raw(x):
    """exact JSON form of an mpf: (sign, mantissa, exponent, bitcount)"""
    s_, m_, e_, b_ = x._mpf_; return [int(s_), int(m_), int(e_), int(b_)]

def iv_raw(x):
    return [mpf_raw(mp.mpf(x.a)), mpf_raw(mp.mpf(x.b))]

def enc(x):  # exact, picklable
    return None if x is None else x._mpi_
def dec(r):
    return None if r is None else iv.mpf([mp.make_mpf(r[0]), mp.make_mpf(r[1])])

def Ljet(T):
    t = J(T, 1, 0); Q1 = J(Qi - 1)
    return jsqrt(t*t*(t - J(2))/(Q1*t + J(2)))        # .v = L, .d = dL/dt

def point(t, K=512):
    """enclosures of Delta, ell at the point t"""
    T = iv.mpf(t)
    return integrate(GA, T, K)[0], integrate(GL, T, K)[0]

def crossing_ell(a, b, delta, K=2048, tol=1e-7):
    """Delta decreasing on [a, b], Delta(a) > delta > Delta(b): bisect to a box [a, b] containing the crossing,
    return enclosure of ell over that box (ell(T) via mean value: ell(m) + ell'(T)(T - m))."""
    while b - a > tol:
        m = 0.5*(a + b)
        D = integrate(GA, iv.mpf(m), K)[0]
        if (D > delta) is True: a = m
        elif (D < delta) is True: b = m
        else: break
    T = iv.mpf([a, b]); m = iv.mpf(0.5*(a + b))
    lm = integrate(GL, m, K)[0]
    dl = integrate(GL, T, 64)[1]                         # d ell/dt on the box
    return lm + dl*(T - m), (a, b)

def _cross(task):
    k, name, a, b, dl = task
    e, box = crossing_ell(a, b, dec(dl))
    return enc(e), box

def neg_chunk(job):
    """Delta' < 0 on [a, b]: mean-value form Delta'(T) in Delta'(m) + Delta''(T)(T - m), bisection"""
    a, b = job
    stack = [(a, b)]; worst = None; n = 0; leaves = []
    while stack:
        x, y = stack.pop(); n += 1
        T = iv.mpf([x, y]); m = iv.mpf(0.5*(x + y))
        d1m = integrate(GA, m, 256)[1]; d2T = integrate(GA, T, 64)[2]
        d1 = d1m + d2T*(T - m)
        if (d1 < 0) is True:
            worst = d1.b if worst is None else max(worst, d1.b)
            leaves.append([x.hex(), y.hex(), mpf_raw(mp.mpf(d1.b))]); continue
        if y - x < 1e-5: return (a, b, None, n, leaves)
        mm = 0.5*(x + y); stack += [(mm, y), (x, mm)]
    return (a, b, enc(iv.mpf(worst)), n, leaves)

def neg_derivative(a, b, pool, nchunks=40):
    e = [a + (b - a)*k/nchunks for k in range(nchunks)] + [b]
    res = pool.map(neg_chunk, [(e[k], e[k+1]) for k in range(nchunks)])
    leaves = sorted(l for r in res for l in r[4])
    if any(r[2] is None for r in res): return None, leaves
    return iv.mpf(max(dec(r[2]).b for r in res)), leaves

def below_chunk(job):
    """Delta < delta_a on [a, b]: mean value Delta(T) in Delta(m) + Delta'(T)(T - m), bisection"""
    a, b, Da = job; Da = dec(Da)
    stack = [(a, b)]; n = 0
    while stack:
        x, y = stack.pop(); n += 1
        T = iv.mpf([x, y]); m = iv.mpf(0.5*(x + y))
        try:
            D = integrate(GA, m, 256)[0] + integrate(GA, T, 64)[1]*(T - m)
            ok = (D < Da) is True
        except Exception:
            ok = False
        if ok: continue
        if y - x < 1e-7: return (a, b, False, n)
        mm = 0.5*(x + y); stack += [(mm, y), (x, mm)]
    return (a, b, True, n)

def F_box(a, b, Dl, K=256, Kb=32):
    """enclosure of F_delta(L(t)) for t in [a, b] and delta in the interval Dl"""
    T = iv.mpf([a, b]); m = iv.mpf(0.5*(a + b))
    Dm = integrate(GA, m, K)[0]; lm = integrate(GL, m, K)[0]; Lm = Ljet(m).v
    Fm = lm - Lm*Dm + Dl*Lm
    DT = integrate(GA, T, Kb)[0]; LpT = Ljet(T).d
    DT = iv.mpf([max(DT.a, iv.mpf(0).a), min(DT.b, PI.b)])     # 0 < Delta < pi (lemma on the dips)
    return Fm + (Dl - DT)*LpT*(T - m)

def box_ok(x, y, ref, K):
    """ref = list of (Lk_min, Lk_max, ell_k(delta_a) enclosure, ell_k(delta_b) enclosure, Da, Db).
    d/d delta [F_delta(L) - ell_k(delta)] = L - L_k(delta): if L >= max L_k it suffices
    F_{delta_a}(L) > ell_k(delta_a); if L <= min L_k, F_{delta_b}(L) > ell_k(delta_b)."""
    Lb = Ljet(iv.mpf([x, y])).v
    for k, (Lmin, Lmax, ea, eb, Da, Db) in enumerate(ref):
        if (Lb >= Lmax) is True:
            F = F_box(x, y, Da, K=K)
            if (F > ea) is True: return True, (F - ea).a, (k, 'a')
        if (Lb <= Lmin) is True:
            F = F_box(x, y, Db, K=K)
            if (F > eb) is True: return True, (F - eb).a, (k, 'b')
    return False, None, None

def cover(job):
    a, b, ref = job
    ref = [(dec(r[0]), dec(r[1]), dec(r[2]), dec(r[3]), dec(r[4]), dec(r[5])) for r in ref]
    stack = [(a, b)]; n = 0; worst = None; leaves = []
    while stack:
        x, y = stack.pop(); n += 1
        ok = False
        for K in (256, 1024, 2048):
            if K > 256 and (y - x) > {1024: 0.5, 2048: 0.05}[K]: break   # refine K only on narrow boxes
            try:
                ok, marg, which = box_ok(x, y, ref, K)
            except Exception:
                ok = False
            if ok: break
        if ok:
            worst = marg if worst is None else min(worst, marg)
            leaves.append([x.hex(), y.hex(), K, which[0] + 1, which[1], mpf_raw(mp.mpf(marg))])
            continue
        if y - x < 1e-6: return (a, b, False, n, (x, y), leaves)
        m = 0.5*(x + y); stack += [(m, y), (x, m)]
    return (a, b, True, n, enc(iv.mpf(worst)) if worst is not None else None, leaves)

if __name__ == "__main__":
    t0 = time.time()
    print(manifest())
    print(f"sha256[:16] {hashlib.sha256(open(__file__, 'rb').read()).hexdigest()[:16]}  {os.path.basename(__file__)}  (this generator)")
    ARCH = {'q': QS, 'q_interval': iv_raw(Qi), 'r0': R0, 'generator_sha256': hashlib.sha256(open(__file__, 'rb').read()).hexdigest()}
    print(f"q interval {Qi}, r0 = {R0}"); sys.stdout.flush()
    t1, t2 = 26.3298461438805, 119.577471305058
    di = mp.mpf('0.677109581418351')*mp.pi
    eps = 1e-5
    Da, Db = iv.mpf(float(di) - eps), iv.mpf(float(di) + eps)
    Dl = iv.mpf([Da.a, Db.b])
    # A. neighbourhoods
    N = []
    pool = Pool(10)
    pg = Progress(2, 'Maxwell A: neighbourhoods')
    for a, b in [(21.3298461438805, 30.0), (105.0, 140.0)]:
        d1, lvA = neg_derivative(a, b, pool)
        ARCH.setdefault('A', []).append({'N': [a.hex(), b.hex()], 'leaves_t_lo_t_hi_dDelta_dt_upper': lvA})
        Dlo, _ = point(a); Dhi, _ = point(b)
        okA = d1 is not None and (Dlo > Db) is True and (Dhi < Da) is True
        print(f"A: N=[{a}, {b}]: Delta' <= {up(d1.b, 9) if d1 is not None else 'FAIL'}, Delta(a) - delta_b >= {down((Dlo - Db).a, 9)}, "
              f"delta_a - Delta(b) >= {down((Da - Dhi).a, 9)}: {'OK' if okA else 'FAIL'}"); sys.stdout.flush()
        ARCH['A'][-1].update({'Delta_a_enclosure': iv_raw(Dlo), 'Delta_b_enclosure': iv_raw(Dhi), 'ok': okA})
        N.append((a, b, okA))
        pg.step()
    pg.close()
    # B. tie by IVT
    ells = {}
    tasks = [(k, name, a, b, enc(dl)) for k, (a, b, _) in enumerate(N) for name, dl in [('a', Da), ('b', Db)]]
    res = pool.map(_cross, tasks)
    boxes_t = {}
    for (k, name, *_), (e, box) in zip(tasks, res):
        e = dec(e); ells[(k, name)] = e; boxes_t[(k, name)] = box
        print(f"B: ell_{k+1}(delta_{name}) in [{down(e.a, 9)}, {up(e.b, 9)}]  (crossing t in [{box[0]!r}, {box[1]!r}])"); sys.stdout.flush()
    Ga = ells[(1, 'a')] - ells[(0, 'a')]; Gb = ells[(1, 'b')] - ells[(0, 'b')]
    okB = (Ga < 0) is True and (Gb > 0) is True
    print(f"B: G(delta_a) <= {up(Ga.b, 9)}, G(delta_b) >= {down(Gb.a, 9)}: {'OK' if okB else 'FAIL'}"); sys.stdout.flush()
    ARCH['delta_a'] = iv_raw(Da); ARCH['delta_b'] = iv_raw(Db)
    ARCH['B'] = [{'k': k + 1, 'delta': name, 'ell_enclosure': iv_raw(ells[(k, name)]),
                  'crossing_t_box': [boxes_t[(k, name)][0].hex(), boxes_t[(k, name)][1].hex()]} for k, name, *_ in tasks]
    U = iv.mpf(min(ells[(0, 'b')].b, ells[(1, 'b')].b))
    # C. dive and edges
    twoR0 = R0_enclosure(Qi, R0, 4096)
    okDive = (twoR0 > U) is True
    t_lo, t_far = 2.01, 1300.0
    Llo = Ljet(iv.mpf(t_lo)).v; Lfar = Ljet(iv.mpf(t_far)).v
    okLo = (twoR0 - PI*Llo > U) is True
    okFar = (Da*Lfar > U) is True
    print(f"C: 2R0 - U >= {down((twoR0 - U).a, 6)} {'OK' if okDive else 'FAIL'}; "
          f"edge (2, {t_lo}]: 2R0 - pi L - U >= {down((twoR0 - PI*Llo - U).a, 6)} {'OK' if okLo else 'FAIL'}; "
          f"far [{t_far}, r0): delta_a L - U >= {down((Da*Lfar - U).a, 3)} {'OK' if okFar else 'FAIL'}"); sys.stdout.flush()
    # reference data per k: L_k ranges over [delta_a, delta_b] (crossing t moves from the a-box to the b-box)
    REF = []
    for k in (0, 1):
        tlo = min(boxes_t[(k, 'a')][0], boxes_t[(k, 'b')][0]); thi = max(boxes_t[(k, 'a')][1], boxes_t[(k, 'b')][1])
        Lr = Ljet(iv.mpf([tlo, thi])).v
        REF.append((enc(iv.mpf(Lr.a)), enc(iv.mpf(Lr.b)), enc(ells[(k, 'a')]), enc(ells[(k, 'b')]), enc(Da), enc(Db)))
    # C0. near 2M: Delta < pi (lemma on the dips) and S = 2R0 - int_0^L Delta give
    #     F_delta(L) - U >= (2R0 - U) - (pi - delta_a) L, decreasing in L: check it at L(t_e).
    t_e = 2.15
    Le = Ljet(iv.mpf(t_e)).v
    c0 = twoR0 - U - (PI - Da)*Le
    okC0 = (c0 > 0) is True
    ARCH['C_analytic'] = {'twoR0': iv_raw(twoR0), 'U': iv_raw(U), 'dive_ok': okDive, 'far_t': t_far, 'far_ok': okFar,
                          'C0_t_e': t_e, 'C0_margin': iv_raw(c0), 'C0_ok': okC0}
    ARCH['REF_Lk_range'] = [[iv_raw(dec(r[0])), iv_raw(dec(r[1]))] for r in REF]
    print(f"C0: (2, {t_e}]: (2R0 - U) - (pi - delta_a) L(t_e) >= {down(c0.a, 6)}: {'OK' if okC0 else 'FAIL'}"); sys.stdout.flush()
    # C. middle cover, split into chunks for the pool
    (a1, b1, _), (a2, b2, _) = N
    segs = [(t_e, a1), (b1, a2), (b2, t_far)]
    jobs = []
    for x, y in segs:
        n = max(1, int((y - x)/2))
        e = [x + (y - x)*k/n for k in range(n)] + [y]
        jobs += [(e[k], e[k+1], REF) for k in range(n)]
    print(f"C: {len(jobs)} chunks"); sys.stdout.flush()
    allok = True; worst = None; boxes = 0; leavesC = []
    if True:
        pg = Progress(len(jobs), 'Maxwell C: chunks of the middle cover')
        for a, b, ok, n, w, lv in pool.imap(cover, jobs):
            boxes += n; leavesC += lv
            if not ok:
                allok = False; print(f"C: FAIL on [{a}, {b}] near {w}"); sys.stdout.flush()
            else:
                w = dec(w).a if w is not None else None
                worst = w if worst is None else (min(worst, w) if w is not None else worst)
            pg.step()
        pg.close()
    leavesC.sort(key=lambda l: float.fromhex(l[0]))
    gap = {a1: b1, a2: b2}             # N1 = [a1, b1], N2 = [a2, b2] are covered by step A
    tiled = bool(leavesC) and float.fromhex(leavesC[0][0]) == t_e and float.fromhex(leavesC[-1][1]) == t_far \
        and all(float.fromhex(l2[0]) == gap.get(float.fromhex(l1[1]), float.fromhex(l1[1]))
                for l1, l2 in zip(leavesC, leavesC[1:]))
    allok = allok and tiled
    ARCH['C_leaves_t_lo_t_hi_K_ref_end_margin'] = leavesC
    print(f"C: {len(leavesC)} accepted leaves tile [{t_e}, {a1}] u [{b1}, {a2}] u [{b2}, {t_far}]: {'OK' if tiled else 'FAIL'}")
    print(f"C: middle cover {'OK' if allok else 'FAIL'}, {boxes} boxes visited (incl. split), smallest margin F - ell_k >= {down(iv.mpf(worst), 6) if worst is not None else 'n/a'}")
    final = all(x[2] for x in N) and okB and okDive and okFar and okC0 and allok
    ARCH['certified'] = final
    os.makedirs('leaves_maxwell', exist_ok=True)
    with open('leaves_maxwell/maxwell_q1.001_r01500.json', 'w') as fh: json.dump(ARCH, fh)
    print(f"archive: leaves_maxwell/maxwell_q1.001_r01500.json, sha256[:16] "
          f"{hashlib.sha256(open('leaves_maxwell/maxwell_q1.001_r01500.json', 'rb').read()).hexdigest()[:16]}")
    print(f"RESULT: {'CERTIFIED' if final else 'NOT certified'} interior Maxwell angle delta* in "
          f"[{down(Da.a/PI.b if False else (Da/PI).a, 12)}, {up((Db/PI).b, 12)}] pi; {time.time()-t0:.0f}s")
