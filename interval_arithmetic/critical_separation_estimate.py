"""critical_separation_estimate.py -- floating-point estimate of the critical separation (mpmath, 30 digits).
NOT a certificate: critical_separation.py uses legs() and dive() only to locate the root t_c of ell - 2R0, and then
proves the enclosure in interval arithmetic around it.

a = 0, A = 1, M = 1, q = Ehat^2 > 1. Metric dR^2 + rho^2 dphi^2, dR = c sqrt(N) dr, c = 1/sqrt(q),
rho^2(r) = r^2 (r - 2)/((q - 1) r + 2).
Regular kernel: for r > t > 2
  rho^2(r) - rho^2(t) = (r - t) W(r,t)/[((q-1)r+2)((q-1)t+2)],
  W = (q-1) r t (r + t - 2) + 2 (r^2 + r t + t^2 - 2r - 2t) > 0,
and with r = t + (r0 - t) u^2 the factor u of dr cancels the one of the square root: smooth real integrands on
[0,1]. Delta(t) = 2 int L c sqrt(N) dr/(rho sqrt(rho^2-L^2)), ell(t) = 2 int rho c sqrt(N) dr/sqrt(...), L = rho(t).
The standalone run checks Delta < 2 arccos(L/rho0) and ell'(L) = L Delta'(L) (finite differences) and prints a grid
of delta_c and of the maximum of Delta; the search for maxima and roots uses a grid (it does not exclude further
maxima or roots)."""
import mpmath as mp, hashlib, sys
mp.mp.dps = 30
SRC = hashlib.sha256(open(__file__, 'rb').read()).hexdigest()[:16]

def model(q):
    q = mp.mpf(q); c = 1/mp.sqrt(q)
    f = lambda r: 1 - 2/r
    N = lambda r: 1/(1 - f(r)/q)
    rho2 = lambda r: r*r*(r - 2)/((q - 1)*r + 2)
    W = lambda r, t: (q - 1)*r*t*(r + t - 2) + 2*(r*r + r*t + t*t - 2*r - 2*t)
    den = lambda r, t: ((q - 1)*r + 2)*((q - 1)*t + 2)
    return q, c, N, rho2, W, den

def legs(q, r0, t):
    q, c, N, rho2, W, den = model(q); r0 = mp.mpf(r0); t = mp.mpf(t)
    L = mp.sqrt(rho2(t)); a = r0 - t
    def common(u):
        r = t + a*u*u
        s = mp.sqrt(W(r, t)/den(r, t))          # sqrt((rho^2 - L^2)/(r - t))
        return r, 2*mp.sqrt(a)*c*mp.sqrt(N(r))/s   # = c sqrt(N) dr / sqrt(rho^2 - L^2) per du
    dphi = 2*mp.quad(lambda u: L*common(u)[1]/mp.sqrt(rho2(common(u)[0])), [0, 1])
    ell = 2*mp.quad(lambda u: mp.sqrt(rho2(common(u)[0]))*common(u)[1], [0, 1])
    return dphi, ell, L

def dive(q, r0):
    q, c, N, rho2, W, den = model(q)
    return 2*c*mp.quad(lambda r: mp.sqrt(N(r)), [2, r0])

if __name__ == "__main__":
    print(f"source sha256[:16] = {SRC}")
    # checks on q = 2, r0 = 6
    q, r0 = 2.0, 6.0; rho0 = mp.sqrt(model(q)[3](mp.mpf(r0)))
    ok = True
    for t in [2.001, 2.1, 2.5, 3.0, 4.0, 5.5]:
        d, l, L = legs(q, r0, t)
        ok &= d < 2*mp.acos(L/rho0)
    print(f"check Delta < 2 arccos(L/rho0) on samples: {ok}")
    t0 = mp.mpf(3); h = mp.mpf('1e-8')
    (d1, l1, L1), (d2, l2, L2) = legs(q, r0, t0 - h), legs(q, r0, t0 + h)
    dl = (l2 - l1)/(L2 - L1); dd = (d2 - d1)/(L2 - L1); Lm = (L1 + L2)/2
    print(f"check ell'(L) = L Delta'(L) at t=3: {mp.nstr(dl,12)} vs {mp.nstr(Lm*dd,12)}  rel.diff {mp.nstr(abs(dl-Lm*dd)/abs(dl),3)}")
    print("q, r0 : delta_c/pi, Delta_max/pi   (grid search, not certified)")
    for q in [1.2, 1.5, 2.0, 4.0, 10.0]:
        row = []
        for r0 in [4.0, 6.0, 10.0, 20.0]:
            D = dive(q, r0)
            ts = [2 + (r0 - 2)*mp.mpf(10)**(-7 + 7*k/40) for k in range(41)]
            ts[-1] = r0 - (r0 - 2)*mp.mpf('1e-7')
            vals = [legs(q, r0, t) for t in ts]
            dphis = [v[0] for v in vals]; gaps = [v[1] - D for v in vals]
            nmax = sum(1 for k in range(1, 40) if dphis[k] > dphis[k-1] and dphis[k] > dphis[k+1])
            ncross = sum(1 for k in range(40) if gaps[k]*gaps[k+1] < 0)
            k = max(range(41), key=lambda i: dphis[i])
            lo, hi = ts[max(k-1, 0)], ts[min(k+1, 40)]; gr = (mp.sqrt(5) - 1)/2
            for _ in range(70):
                a1 = hi - gr*(hi - lo); b1 = lo + gr*(hi - lo)
                if legs(q, r0, a1)[0] > legs(q, r0, b1)[0]: hi = b1
                else: lo = a1
            dmax = legs(q, r0, (lo + hi)/2)[0]
            kc = next(i for i in range(40) if gaps[i]*gaps[i+1] < 0)
            tc = mp.findroot(lambda t: legs(q, r0, t)[1] - D, (ts[kc], ts[kc+1]), solver='anderson')
            dc = legs(q, r0, tc)[0]
            row.append(f"r0={r0:4.0f}: {mp.nstr(dc/mp.pi, 7)}, {mp.nstr(dmax/mp.pi, 7)} [maxima {nmax}, crossings {ncross}]")
        print(f"q={q:5.2f}  " + " | ".join(row)); sys.stdout.flush()
