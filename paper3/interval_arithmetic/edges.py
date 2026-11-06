"""edges.py -- the two edges of the unimodality certificate (unimodality.py covers [t_lo, t_hi] only).

With rho = L x:  Delta(L) = 2 int_1^{rho0/L} w(Lx)/(x sqrt(x^2-1)) dx,  w = dR/drho in (0,1), and
   Delta'(L)/2 = int_0^{arcosh(rho0/L)} w'(L cosh th) dth - w(rho0)/sqrt(rho0^2 - L^2).
Regular variable r = 2 + s^2 (M = 1, A = 1: angles and signs do not depend on the constant A).

Edge L -> 0, L in (0, L_lo]: [0, s0] is split into [0, s_1] and cells [e_k, e_{k+1}], with e_N = s0 EXACT (the
interval S0 = sqrt(r0 - 2) is enclosed), and w' >= m_k on the cell (enclosed on a cell that contains the true one).
Piece k contributes m_k [arcosh(rho_{k+1}/L) - arcosh(rho_k/L)]; the difference increases with L and tends to
log(rho_{k+1}/rho_k) as L -> 0: for m_k > 0 the limit (the minimum) is used, for m_k < 0 the value at L_lo (the
maximum). First piece: m_1 arcosh(rho_1/L) >= m_1 arcosh(rho_1/L_lo).
One checks that the sum minus w(rho0)/sqrt(rho0^2 - L_lo^2) is > 0: everything decreases in L, so it holds for
every L in (0, L_lo].
Edge L -> rho0, L in [L_hi, rho0): Delta'(L)/2 <= B arcosh(rho0/L) - w(rho0)/sqrt(rho0^2 - L^2) with |w'| <= B on
[0, s0]; the right-hand side DECREASES in L, so it suffices to check it < 0 at L_hi (lower end point of the
interval of L_hi). B is enclosed on cells that cover [0, S0.b].
q is an interval (jets_nested.qint): for '1.2' it contains both 6/5 and the binary64 value.
t_lo, t_hi are the same floats as in unimodality.py. Usage (from this folder): python3 edges.py (about a minute)."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'tools'))
from progress import Progress
import sys
from mpmath import iv, mp, mpf
iv.dps = 20; mp.dps = 30
sys.path.insert(0, '.')
from jets_t import J, jsqrt
from jets_nested import qint, down, up, manifest

def ivacosh(y):
    return iv.log(y + iv.sqrt(y*y - 1))

def jets(Qi, S):
    Q = J(Qi); Q1 = J(Qi - 1); one = J(1); two = J(2); C = J(1/iv.sqrt(Qi))
    s = J(S, 1, 0)
    r = two + s*s
    rho = s*r/jsqrt(Q1*r + two)
    f = (s*s)/r
    N = one/(one - f/Q)
    Rp = two*s*C*jsqrt(N)            # dR/ds as a jet: .v = R', .d = R''
    rp, rpp = rho.d, rho.dd
    w = Rp.v/rp
    wprime = (Rp.d*rp - Rp.v*rpp)/(rp*rp*rp)
    return rho.v, w, wprime

def cells(lo, hi, n):
    """n contiguous cells covering [lo, hi]: interior edges are binary points shared by
    neighbours, the outer edges are lo and hi themselves (points or intervals)."""
    lo_f = float(mpf(lo.a)) if isinstance(lo, type(iv.mpf(0))) else float(lo)
    hi_f = float(mpf(hi.a)) if isinstance(hi, type(iv.mpf(0))) else float(hi)
    inner = [lo_f + (hi_f - lo_f)*k/n for k in range(1, n)]
    edges = [lo] + inner + [hi]
    for x, y in zip(inner, inner[1:]): assert x < y
    return edges

def check(qs, r0, t_lo, t_hi, n=400):
    Qi = qint(qs)
    S0 = iv.sqrt(iv.mpf(r0) - 2)                      # exact end point, enclosed
    rho0, w0, _ = jets(Qi, S0)
    # B: |w'| on cells covering [0, S0.b]
    E = cells(0.0, iv.mpf(S0.b), n); B = iv.mpf(0)
    for a, b in zip(E, E[1:]):
        S = iv.mpf([iv.mpf(a).a, iv.mpf(b).b]); wp = jets(Qi, S)[2]
        B = iv.mpf(max(abs(wp.a), abs(wp.b), B.b))
    Llo = jets(Qi, iv.sqrt(iv.mpf(t_lo) - 2))[0]; Lhi = jets(Qi, iv.sqrt(iv.mpf(t_hi) - 2))[0]
    # lower edge, piecewise
    L = iv.mpf(Llo.b)
    s1v = 3*float(mp.sqrt(t_lo - 2))                  # a binary point; only needs rho_1 > L_lo
    rho_1 = jets(Qi, iv.mpf(s1v))[0]
    E1 = cells(0.0, s1v, 50); m1 = None
    for a, b in zip(E1, E1[1:]):
        wp = jets(Qi, iv.mpf([a, b]))[2]; m1 = wp.a if m1 is None else min(m1, wp.a)
    lower_ok = False; total = None
    if m1 > 0 and (L < rho_1.a) is True:
        total = m1*ivacosh(iv.mpf(rho_1.a)/L)
        E = cells(s1v, S0, n)
        assert (iv.mpf(E[-2]) < S0) is True           # last inner edge strictly inside
        rhos = [jets(Qi, iv.mpf(e))[0] for e in E]    # last one is rho(S0) = rho0, an interval
        for k in range(n):
            S = iv.mpf([iv.mpf(E[k]).a, iv.mpf(E[k + 1]).b])
            mi = jets(Qi, S)[2].a
            ra, rb = rhos[k], rhos[k + 1]
            if mi > 0: total += mi*iv.log(iv.mpf(rb.a)/iv.mpf(ra.b))
            else:      total += mi*(ivacosh(iv.mpf(rb.b)/L) - ivacosh(iv.mpf(ra.a)/L))
        total -= w0/iv.sqrt(rho0*rho0 - L*L)
        lower_ok = (total > 0) is True
    # upper edge, at the smallest L of the edge
    L2 = iv.mpf(Lhi.a)
    upper = B*ivacosh(rho0/L2) - w0/iv.sqrt(rho0*rho0 - L2*L2)
    upper_ok = (upper < 0) is True
    return lower_ok, upper_ok, total, upper

if __name__ == "__main__":
    print(manifest()); sys.stdout.flush()
    eps = 1e-4
    pg = Progress(20, 'edges: (q, r0) cases')
    for qs in ['1.2', '1.5', '2', '4', '10']:
        for r0 in [4.0, 6.0, 10.0, 20.0]:
            t_lo = float(2 + eps*(r0 - 2)); t_hi = float(r0 - eps*(r0 - 2))
            lo_ok, up_ok, tot, upv = check(qs, r0, t_lo, t_hi)
            print(f"q={qs:>4} r0={r0:5.1f}: lower edge (2, {t_lo!r}] {'OK' if lo_ok else 'FAIL'} "
                  f"(Delta'/2 >= {down(tot.a, 6) if tot is not None else 'n/a'}), "
                  f"upper edge [{t_hi!r}, r0) {'OK' if up_ok else 'FAIL'} (Delta'/2 <= {up(upv.b, 6)})")
            sys.stdout.flush()
            pg.step()
    pg.close()
