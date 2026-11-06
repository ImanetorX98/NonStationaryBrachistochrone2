"""Numerical diagnostic (guidance, not a proof): Thakurta-Kerr on the equator, proper-time branch, constant A
(main text, theorem on the Maxwell rule with a ring; Figure 7).
Metric (c = A^2/E = 1, M = 1): h = N(r)[f r^2/Dl dr^2 + Dl dphi^2], N = q/(q-f), q = E^2/A^2.
rho = sqrt(N Dl), R'(r) = sqrt(N f) r/sqrt(Dl).  Ring r = 2: rho0 = |a|.
Computes: Phi_t, S_t (geodesic tangent to the ring), Delta(L), ell(L) of the symmetric dips,
the number of maxima of Delta on (rho0, rho(r0)), delta_ring, delta_c (Maxwell rule),
and an independent check: distance in the completion by shortest paths on a grid (Dijkstra).
"""
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
import sys

M = 1.0


def fns(a, q):
    f = lambda r: 1 - 2 * M / r
    Dl = lambda r: r * r - 2 * M * r + a * a
    N = lambda r: q / (q - f(r))
    rho = lambda r: np.sqrt(N(r) * Dl(r))
    Rp = lambda r: np.sqrt(N(r) * f(r)) * r / np.sqrt(Dl(r))
    return f, Dl, N, rho, Rp


def tangent(a, q, r):
    """angle and reduced action of the geodesic tangent to the ring (L = rho0 = a)."""
    f, Dl, N, rho, Rp = fns(a, q)
    L = a
    # substitution r = 2 + t^2: regular integrand
    def g_phi(t):
        x = 2 * M + t * t
        return 2 * t * L * Rp(x) / (rho(x) * np.sqrt(rho(x) ** 2 - L * L)) if t > 0 else 0.0
    def g_S(t):
        x = 2 * M + t * t
        return 2 * t * np.sqrt(rho(x) ** 2 - L * L) / rho(x) * Rp(x)
    def g_l(t):
        x = 2 * M + t * t
        return 2 * t * rho(x) * Rp(x) / np.sqrt(rho(x) ** 2 - L * L) if t > 0 else 0.0
    T = np.sqrt(r - 2 * M)
    Phi = quad(g_phi, 0, T, limit=200)[0]
    S = quad(g_S, 0, T, limit=200)[0]
    ell = quad(g_l, 0, T, limit=200)[0]
    return Phi, S, ell


def dip(a, q, r0, rp):
    """symmetric dip with periastron rp between two points at radius r0: (L, Delta, ell)."""
    f, Dl, N, rho, Rp = fns(a, q)
    L = rho(rp)
    def g_phi(t):
        x = rp + t * t
        d = rho(x) ** 2 - L * L
        if t == 0:
            # limit: 2t/sqrt(d) -> 2/sqrt(d'(rp))
            dd = (rho(rp + 1e-9) ** 2 - L * L) / 1e-9
            return 2 * L * Rp(rp) / (rho(rp) * np.sqrt(dd))
        return 2 * t * L * Rp(x) / (rho(x) * np.sqrt(d))
    def g_l(t):
        x = rp + t * t
        d = rho(x) ** 2 - L * L
        if t == 0:
            dd = (rho(rp + 1e-9) ** 2 - L * L) / 1e-9
            return 2 * rho(rp) * Rp(rp) / np.sqrt(dd)
        return 2 * t * rho(x) * Rp(x) / np.sqrt(d)
    T = np.sqrt(r0 - rp)
    Phi = 2 * quad(g_phi, 0, T, limit=200)[0]
    ell = 2 * quad(g_l, 0, T, limit=200)[0]
    return L, Phi, ell


def analyse(a, q, r0, n=400, verbose=True):
    f, Dl, N, rho, Rp = fns(a, q)
    Phit, St, ellt = tangent(a, q, r0)
    dring = 2 * Phit
    # grid in the periastron, denser towards the ring
    s = np.linspace(0, 1, n + 2)[1:-1]
    rps = 2 * M + (r0 - 2 * M) * s ** 2
    data = np.array([dip(a, q, r0, rp) for rp in rps])
    Ls, Ds, ells = data[:, 0], data[:, 1], data[:, 2]
    # number of local maxima of Delta (as a function of the periastron)
    dD = np.diff(Ds)
    sign_changes = np.sum(np.diff(np.sign(dD)) != 0)
    imax = np.argmax(Ds)
    Phimax = Ds[imax]
    # shallow branch: periastra above the maximum; comparison with the wrapping W = 2 S_t + a delta
    W = lambda d: 2 * St + a * d
    shallow = slice(imax, None)
    Dsh, lsh = Ds[shallow], ells[shallow]
    gsh = lsh - W(Dsh)
    # delta_c: zero of g on the shallow branch (g increasing in delta)
    dc = None
    if Phimax >= dring:
        idx = np.where(np.diff(np.sign(gsh)) != 0)[0]
        if len(idx):
            k = idx[0]
            # linear interpolation
            dc = Dsh[k] + (Dsh[k + 1] - Dsh[k]) * (0 - gsh[k]) / (gsh[k + 1] - gsh[k])
    # deep branch above the wrapping?
    deep = slice(0, imax + 1)
    deep_above = np.all(ells[deep] - W(Ds[deep]) > -1e-9)
    # continuity at the ring
    cont = (abs(Ds[0] - dring), abs(ells[0] - (2 * ellt)))
    if verbose:
        print(f"a={a:.3f} q={q:.3f} r0={r0:.2f}: rho0={a:.3f} Phi_t={Phit:.5f} delta_ring={dring:.5f} "
              f"Phi_max={Phimax:.5f} (periastron {rps[imax]:.4f}) sign changes of dDelta={sign_changes} "
              f"deep>=W:{deep_above} delta_c={dc} cont(ring)={cont[0]:.1e},{cont[1]:.1e}")
    return dict(Phit=Phit, St=St, dring=dring, Phimax=Phimax, nchange=sign_changes, dc=dc,
                deep_above=deep_above, Ds=Ds, ells=ells, Ls=Ls, rps=rps)


def dijkstra_check(a, q, r0, delta, nr=160, nphi=360, rmax_fac=1.0):
    """distance in the completion (ring included, r = 2 as a grid row) on a graph with 16 neighbours."""
    import heapq
    f, Dl, N, rho, Rp = fns(a, q)
    # coordinates: R (radial distance from the ring) and phi; metric dR^2 + rho(R)^2 dphi^2
    rs = 2 * M + (r0 * 1.6 - 2 * M) * np.linspace(0, 1, nr) ** 1.5
    Rs = np.array([0.0] + [quad(Rp, 2 * M, x)[0] for x in rs[1:]])
    rhos = np.array([a] + [rho(x) for x in rs[1:]])
    phis = np.linspace(0, 2 * np.pi, nphi, endpoint=False)
    dphi = phis[1]
    steps = [(di, dj) for di in range(-2, 3) for dj in range(-3, 4) if (di, dj) != (0, 0)]
    # source: r0, phi=0; target: r0, phi=delta (approximated on the grid)
    i0 = int(np.argmin(abs(rs - r0)))
    j1 = int(round(delta / dphi)) % nphi
    dist = np.full((nr, nphi), np.inf)
    dist[i0, 0] = 0
    pq = [(0.0, i0, 0)]
    while pq:
        d, i, j = heapq.heappop(pq)
        if d > dist[i, j]:
            continue
        if (i, j) == (i0, j1):
            break
        for di, dj in steps:
            ii, jj = i + di, (j + dj) % nphi
            if ii < 0 or ii >= nr:
                continue
            # length of the segment: midpoint rule in R
            dR = Rs[ii] - Rs[i]
            rm = 0.5 * (rhos[ii] + rhos[i])
            w = np.sqrt(dR * dR + (rm * dj * dphi) ** 2)
            nd = d + w
            if nd < dist[ii, jj]:
                dist[ii, jj] = nd
                heapq.heappush(pq, (nd, ii, jj))
    return dist[i0, j1], rs[i0], j1 * dphi


if __name__ == "__main__":
    print("== shape of Delta(L) and Maxwell rule, equal radii ==")
    for a in [0.1, 0.3, 0.5, 0.7, 0.9, 0.99]:
        for q in [1.2, 1.5, 2.0, 4.0]:
            for r0 in [3.0, 6.0, 12.0]:
                analyse(a, q, r0, n=300)
    print("== independent check (Dijkstra on the completion) ==")
    for (a, q, r0) in [(0.5, 2.0, 6.0), (0.9, 1.5, 4.0), (0.3, 2.0, 6.0)]:
        res = analyse(a, q, r0, n=400, verbose=False)
        for delta in [0.5, 1.0, 2.0, 3.0]:
            dg, rg, dgrid = dijkstra_check(a, q, r0, delta)
            # formula: min(interpolated shallow branch, W)
            Ds, ells = res["Ds"], res["ells"]
            imax = int(np.argmax(Ds))
            Dsh, lsh = Ds[imax:], ells[imax:]
            order = np.argsort(Dsh)
            lint = np.interp(dgrid, Dsh[order], lsh[order]) if dgrid <= Dsh.max() else np.inf
            # the deep branch too (for completeness)
            Dd, ld = Ds[:imax + 1], ells[:imax + 1]
            o2 = np.argsort(Dd)
            ldeep = np.interp(dgrid, Dd[o2], ld[o2]) if Dd.min() <= dgrid <= Dd.max() else np.inf
            W = 2 * res["St"] + a * dgrid if dgrid >= res["dring"] else np.inf
            # grid radius rg ~ r0: no recomputation needed if rg == r0
            print(f"a={a} q={q} r0={r0} (grid r={rg:.3f}) delta={dgrid:.4f}: Dijkstra={dg:.5f} "
                  f"formula min(shallow={lint:.5f}, deep={ldeep:.5f}, W={W:.5f}) = {min(lint, ldeep, W):.5f}")
