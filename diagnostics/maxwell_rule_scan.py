"""maxwell_rule_scan.py -- floating-point diagnostic (mpmath, not certified): selection of the minimum of the
proper-time branch at equal radii when Delta has several maxima (main text, theorem on the Maxwell rule: "at q = 1.01,
r0 = 10^4 M, where Delta has two maxima, ... the minimizer is unique below delta_c").

Idea: S(L) = ell - L Delta, S' = -Delta (first-variation lemma), S(0) = 2 R0. For an angle delta,
F_delta(L) = S(L) + delta L has F' = delta - Delta(L): the critical points are the dips of angle delta, and
F_delta(dip) = ell(dip). So d(delta) = min_{L in [0, rho0]} F_delta(L): the minima are the contact points of S with its
lower convex hull; the bridges of the hull are the ties (Maxwell's equal-area rule). The bridge that starts at L = 0
gives delta_c.

Checks: (i) min_L F_delta on the grid = min over the dips of angle delta (direct search of the crossings Delta = delta)
and 2R0; (ii) equal areas on the bridges; (iii) case q = 2 (a single maximum, theorem): a single bridge, the one of
delta_c, which must agree with the certified table."""
import sys, time
import mpmath as mp
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'interval_arithmetic'))
from critical_separation_estimate import legs, dive
mp.mp.dps = 20

def profile(q, r0, n=160):
    r0 = mp.mpf(r0)
    # grid in t: dense near 2 and near r0 (log spacing from both ends)
    ts = []
    for k in range(n + 1):
        x = mp.mpf(k)/n                                   # 0..1
        if x <= mp.mpf(1)/2: ts.append(2 + (r0 - 2)*mp.mpf(1)/2*mp.mpf(10)**(-8*(1 - 2*x)))
        else:                ts.append(r0 - (r0 - 2)*mp.mpf(1)/2*mp.mpf(10)**(-8*(2*x - 1)))
    ts = sorted(set(t for t in ts if 2 < t < r0))
    pts = []
    for t in ts:
        D, l, L = legs(q, r0, t)
        pts.append((L, D, l, l - L*D, t))
    return pts

def lower_hull(P):
    H = []
    for p in P:
        while len(H) >= 2 and ((H[-1][0]-H[-2][0])*(p[1]-H[-2][1]) - (H[-1][1]-H[-2][1])*(p[0]-H[-2][0])) <= 0:
            H.pop()
        H.append(p)
    return H

def analyse(q, r0, n=160):
    t0 = time.time()
    pts = profile(q, r0, n)
    R2 = dive(q, r0)
    P = [(mp.mpf(0), R2, -1)] + [(p[0], p[3], i) for i, p in enumerate(pts)]
    H = lower_hull(P)
    # bridges: hull edges skipping at least 3 grid points
    out = []
    for a, b in zip(H, H[1:]):
        if b[2] - a[2] > 3:
            slope = (b[1] - a[1])/(b[0] - a[0]); delta = -slope
            # equal area check: int_{La}^{Lb} (Delta - delta) dL = S(La) - S(Lb) - delta (Lb - La) = 0 by construction;
            # independent: trapezoid of Delta over the grid points between
            seg = [(mp.mpf(0), mp.mpf(0))] if a[2] < 0 else []
            seg += [(pts[i][0], pts[i][1]) for i in range(max(a[2], 0), b[2] + 1)]
            area = sum((seg[j+1][0]-seg[j][0])*((seg[j][1]+seg[j+1][1])/2 - delta) for j in range(len(seg)-1))
            out.append(dict(La=a[0], Lb=b[0], ta=(pts[a[2]][4] if a[2] >= 0 else mp.mpf(2)), tb=pts[b[2]][4],
                            delta=delta, area=area, from_dive=(a[2] < 0)))
    # maxima of Delta on the grid
    Ds = [p[1] for p in pts]
    maxima = [(pts[i][4], Ds[i]) for i in range(1, len(Ds)-1) if Ds[i] > Ds[i-1] and Ds[i] > Ds[i+1]]
    minima = [(pts[i][4], Ds[i]) for i in range(1, len(Ds)-1) if Ds[i] < Ds[i-1] and Ds[i] < Ds[i+1]]
    # check (i): for sample deltas, min over grid of F_delta vs brute force over crossings
    worst = 0
    for frac in [0.1, 0.3, 0.5, 0.7, 0.9, 0.97]:
        delta = frac*max(Ds)
        Fmin = min([R2] + [p[3] + delta*p[0] for p in pts])
        cands = [R2]
        for i in range(len(pts)-1):
            if (Ds[i] - delta)*(Ds[i+1] - delta) < 0:
                w = (delta - Ds[i])/(Ds[i+1] - Ds[i])
                cands.append(pts[i][2] + w*(pts[i+1][2] - pts[i][2]))   # ell interpolated at the crossing
        worst = max(worst, abs(Fmin - min(cands))/min(cands))
    return dict(maxima=maxima, minima=minima, bridges=out, R2=R2, Dmax=max(Ds), worst=worst, secs=time.time()-t0)

if __name__ == "__main__":
    for q, r0 in [(2.0, 6.0), (1.01, 1e4)]:
        res = analyse(q, r0)
        print(f"q={q} r0={r0}: {len(res['maxima'])} maxima of Delta at t = "
              + ", ".join(f"{mp.nstr(t,6)} ({mp.nstr(D/mp.pi,6)} pi)" for t, D in res['maxima'])
              + "; minima at t = " + ", ".join(f"{mp.nstr(t,6)} ({mp.nstr(D/mp.pi,6)} pi)" for t, D in res['minima']))
        for b in res['bridges']:
            print(f"   bridge {'from the dive' if b['from_dive'] else 'interior'}: t in [{mp.nstr(b['ta'],6)}, {mp.nstr(b['tb'],6)}], "
                  f"Maxwell angle {mp.nstr(b['delta']/mp.pi,8)} pi, equal-area residual {mp.nstr(b['area'],3)}")
        print(f"   check min_L F_delta = brute force over crossings: worst rel. diff {mp.nstr(res['worst'],3)} "
              f"(grid interpolation, not exact); {res['secs']:.0f}s")
        sys.stdout.flush()
