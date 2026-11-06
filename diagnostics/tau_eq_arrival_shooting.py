"""tau_eq_arrival_shooting.py -- floating-point guide, not a certificate: shooting from the arrival for the example on
the Kerr equator with spin (library of tau_eq_injectivity_scan.py, tau_eq_uniqueness_checks.py and
tau_eq_sector_change_checks.py).
Example: M = 1, a = 1/2, Ehat = sqrt 2, x0 = (6, 0), target (6, s/50) with s = +1 (prograde) or -1 (retrograde),
A = 1 - (1/2) sin^2(pi eta/tau) on [0, tau], then 1 (main text, example on the Kerr equator with spin).
(1) Two-parameter shooting from the arrival (T, psi): p(T) on the ellipse |p|_{h0^-1} = omega sqrt(N) (H(T) = 0
    exactly), p_r = omega sqrt(N h0rr) cos psi, p_phi = omega sqrt(N h0pp) sin psi; backward integration down to
    eta = 0; residual E(T, psi) = X(0) - x0. Grid over ALL directions psi and T in (0, 35/44); zeros refined by Newton
    with the Jacobian from the variational equations (not finite differences). Failures are recorded too.
    Competitive filter (main text, proposition on localisation): cost <= 7/80, T < 35/44, 23/4 < r < 25/4.
    Diagnostics: chi_h = min |p + omega beta|_{h0^-1} (intrinsic norm), sigma_min and condition number of the 2x2
    Jacobian (scale: E in coordinates (r, phi), parameters (T, psi)), T - tau.
(2) Short-pulse limit for tau >= 1/5 (exploratory): E/delta against the limit map (0, s) - t v0(psi), t = T/delta."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'tools'))
from progress import Progress
import numpy as np, sympy as sp, time, sys, warnings
from scipy.integrate import solve_ivp
warnings.filterwarnings('ignore')

M, a, E_ = 1.0, 0.5, np.sqrt(2.0)
r, pr, pp, al = sp.symbols('r pr pp al', real=True)
f = 1 - 2*M/r; Dl = r**2 - 2*M*r + a**2
h0rr = r**2/(f*Dl); h0pp = Dl/f**2; bp = -2*M*a/(r - 2*M)
om = al**2*f/E_; N = 1/(1 - al**2*f/E_**2)
xr, xp = pr, pp + om*bp
xb = xp*bp/h0pp; b2 = bp**2/h0pp; xx = xr**2/h0rr + xp**2/h0pp
H = -om + (sp.sqrt(xb**2 + (N - b2)*xx) - xb)/(N - b2)
Z = [r, pr]                                    # phi cyclic, p_phi constant
G = [sp.diff(H, pr), sp.diff(H, pp), -sp.diff(H, r)]     # (rdot, phidot, prdot)
J3 = sp.Matrix(G).jacobian([r, pr, pp])        # derivatives with respect to (r, p_r, p_phi)
Gf = sp.lambdify((r, pr, pp, al), G, 'numpy')
J3f = sp.lambdify((r, pr, pp, al), J3, 'numpy')
Hf = sp.lambdify((r, pr, pp, al), H, 'numpy')
costf = sp.lambdify((r, pr, pp, al), om*(1 - bp*sp.diff(H, pp)), 'numpy')
chif = sp.lambdify((r, pr, pp, al), sp.sqrt(xr**2/h0rr + xp**2/h0pp), 'numpy')
termf = sp.lambdify((r, al), [om*sp.sqrt(N*h0rr), om*sp.sqrt(N*h0pp)], 'numpy')
dterm_dal = sp.lambdify((r, al), [sp.diff(om*sp.sqrt(N*h0rr), al), sp.diff(om*sp.sqrt(N*h0pp), al)], 'numpy')

def Aprof(eta, tau):
    return 1.0 - 0.5*np.sin(np.pi*eta/tau)**2 if 0 <= eta <= tau else 1.0
def dAprof(eta, tau):
    return -0.5*np.pi/tau*np.sin(2*np.pi*eta/tau) if 0 <= eta <= tau else 0.0

def backward(T, psi, tau, s, want_jac=False):
    """integrates from eta=T to 0. State y = (r, phi, p_r, p_phi, cost) + variational columns."""
    al_T = Aprof(T, tau)
    cr, cp = termf(6.0, al_T)
    p0r, p0p = cr*np.cos(psi), cp*np.sin(psi)
    # derivatives of the terminal datum: with respect to psi, and with respect to T (through A(T); x_f fixed)
    dpsi = np.array([0.0, 0.0, -cr*np.sin(psi), cp*np.cos(psi)])
    dcr, dcp = dterm_dal(6.0, al_T); dA = dAprof(T, tau)
    dTterm = np.array([0.0, 0.0, dcr*dA*np.cos(psi), dcp*dA*np.sin(psi)])
    def rhs(eta, y):
        rr, ph, prr, ppp = y[:4]; alp = Aprof(eta, tau)
        g = Gf(rr, prr, ppp, alp)
        out = [g[0], g[1], g[2], 0.0, costf(rr, prr, ppp, alp)]
        if want_jac:
            J = np.array(J3f(rr, prr, ppp, alp), dtype=float)      # d(rdot,phidot,prdot)/d(r,pr,pp)
            for k in range(2):
                v = y[5 + 4*k: 9 + 4*k]                          # (dr, dphi, dpr, dpp)
                dv = J @ np.array([v[0], v[2], v[3]])
                out += [dv[0], dv[1], dv[2], 0.0]
        return out
    y0 = [6.0, s/50, p0r, p0p, 0.0]
    if want_jac: y0 += list(dpsi) + list(dTterm)
    # two pieces: [T, tau] with constant A (free step) and [min(T,tau), 0] inside the window (step tau/150)
    class _Cat:
        pass
    segs = []
    if T > tau:
        s1 = solve_ivp(rhs, [T, tau], y0, rtol=1e-10, atol=1e-12, dense_output=want_jac)
        if not s1.success: return None
        segs.append(s1); y0 = s1.y[:, -1]
    s2 = solve_ivp(rhs, [min(T, tau), 0.0], y0, rtol=1e-10, atol=1e-12, max_step=min(T, tau)/150,
                   dense_output=want_jac)
    if not s2.success: return None
    segs.append(s2)
    sol = _Cat(); sol.y = np.concatenate([g.y for g in segs], axis=1); sol.success = True
    if want_jac:
        def _sol(eta, segs=segs, tau=tau, T=T):
            return segs[0].sol(eta) if (len(segs) == 2 and eta >= tau) else segs[-1].sol(eta)
        sol.sol = lambda e: np.array([_sol(x) for x in np.atleast_1d(e)]).T if np.ndim(e) else _sol(e)
    if np.any(sol.y[0] <= 2.0 + 1e-6):
        return None
    yE = sol.y[:, -1]
    res = np.array([yE[0] - 6.0, yE[1] - 0.0])
    cost = -yE[4]                                               # integrated backwards
    out = dict(res=res, cost=cost, rmin=sol.y[0].min(), rmax=sol.y[0].max(), pphi=p0p)
    if want_jac:
        # dX(0)/dpsi = column 1; dX(0)/dT = Phi(0,T)[dZ_T/dT] - G(T) transported: the second variational column carries
        # dZ_T/dT; the term -G(Z_T) is transported by the same flow: it is added by integrating a third column
        Jpsi = yE[5:7]
        gT = Gf(6.0, p0r, p0p, al_T)
        # transport of -G(Z_T): linear flow, known combination -> separate integration
        def rhs_g(eta, v):
            yy = sol.sol(eta); J = np.array(J3f(yy[0], yy[2], yy[3], Aprof(eta, tau)), dtype=float)
            dv = J @ np.array([v[0], v[2], v[3]]); return [dv[0], dv[1], dv[2], 0.0]
        vg = [-gT[0], -gT[1], -gT[2], 0.0]
        if T > tau:
            sg1 = solve_ivp(rhs_g, [T, tau], vg, rtol=1e-10, atol=1e-12); vg = sg1.y[:, -1]
        sg = solve_ivp(rhs_g, [min(T, tau), 0.0], vg, rtol=1e-10, atol=1e-12, max_step=min(T, tau)/150)
        JT = yE[9:11] + sg.y[:2, -1]
        out['J'] = np.column_stack([JT, Jpsi])
        etas = np.linspace(0, T, 300); Y = sol.sol(etas)
        out['chi'] = min(chif(Y[0, k], Y[2, k], Y[3, k], Aprof(e, tau)) for k, e in enumerate(etas))
    return out

def search(tau, s, nT=24, npsi=48, Tmax=35/44):
    Ts = np.linspace(0.05, Tmax*0.999, nT); psis = np.linspace(0, 2*np.pi, npsi, endpoint=False)
    grid = np.full((nT, npsi), np.inf); fails = 0
    for i, T in enumerate(Ts):
        for j, ps in enumerate(psis):
            o = backward(T, ps, tau, s)
            if o is None: fails += 1; continue
            grid[i, j] = np.linalg.norm(o['res']*np.array([1.0, 6.0]))   # scale: phi weighted by the radius
    # seeds: local minima of the grid
    seeds = []
    for i in range(nT):
        for j in range(npsi):
            v = grid[i, j]
            if not np.isfinite(v): continue
            nb = [grid[ii % nT if 0 <= ii < nT else i, jj % npsi] for ii in (i-1, i, i+1) for jj in (j-1, j, j+1)
                  if 0 <= ii < nT]
            if v <= min(nb) and v < 0.5: seeds.append((v, Ts[i], psis[j]))
    seeds.sort()
    sols = []
    for v, T, ps in seeds[:25]:
        u = np.array([T, ps]); ok = False
        for it in range(40):
            o = backward(u[0], u[1], tau, s, want_jac=True)
            if o is None: break
            if np.linalg.norm(o['res']) < 1e-12: ok = True; break
            try: du = np.linalg.solve(o['J'], -o['res'])
            except np.linalg.LinAlgError: break
            lam = 1.0
            while lam > 1e-3:
                un = u + lam*du
                if un[0] > 0 and un[0] < Tmax: break
                lam /= 2
            u = un
        if not ok: continue
        u[1] = u[1] % (2*np.pi)
        if any(abs(u[0] - q['T']) < 1e-7 and abs(np.angle(np.exp(1j*(u[1] - q['psi'])))) < 1e-6 for q in sols): continue
        sv = np.linalg.svd(o['J'], compute_uv=False)
        sols.append(dict(T=u[0], psi=u[1], cost=o['cost'], res=np.linalg.norm(o['res']), chi=o['chi'], rmin=o['rmin'],
                         rmax=o['rmax'], smin=sv[-1], cond=sv[0]/sv[-1],
                         comp=(o['cost'] <= 7/80 and u[0] < 35/44 and 23/4 < o['rmin'] and o['rmax'] < 25/4)))
    return sorted(sols, key=lambda q: q['cost']), fails, len(seeds)

def constants_6216(tau, s, delta=1/50, nt=25, npsi=72):
    """comparison of E/delta with the limit map (0,s) - t v0(psi), t=T/delta in [c,C]."""
    hrr, hpp, N0, bphi = 216/97, 873/16, 1.5, -0.25
    def v0(ps):
        n = np.array([np.cos(ps)/np.sqrt(hrr), np.sin(ps)/np.sqrt(hpp)])
        return n/(np.sqrt(N0) + bphi*n[1])
    ts = np.linspace(5.0, 16.0, nt); psis = np.linspace(0, 2*np.pi, npsi, endpoint=False)
    dev = 0.0; minlim = np.inf
    for t in ts:
        for ps in psis:
            o = backward(delta*t, ps, tau, s)
            if o is None: continue
            Ered = np.array([o['res'][0], o['res'][1] - 0.0])/delta
            # residual defined as X(0) - x0 = (r(0)-6, phi(0)); limit: x0 - (xf - T v0) -> E/delta ~ (0,s) - t v0 ... with sign:
            Elim = -(np.array([0.0, s]) - t*v0(ps))          # X(0) = xf - T v0  =>  X(0)-x0 = (0, s delta) - T v0
            Elim = -Elim
            dev = max(dev, np.linalg.norm(Ered - Elim))
            minlim = min(minlim, np.linalg.norm(Elim))
    return dev, minlim

if __name__ == "__main__":
    t0 = time.time()
    taus = [0.01, 0.05, 0.1, 0.15, 0.155, 0.16, 0.165, 0.17, 0.175, 0.18, 0.2, 0.5, 1.0]
    print("(1) Shooting from the arrival (T, psi), all terminal directions; competitive filter C<=7/80, T<35/44, r in (23/4,25/4)", flush=True)
    pg = Progress(2*len(taus), 'arrival shooting (1): (sense, tau)')
    for s, lab in [(+1, 'prograde'), (-1, 'retrograde')]:
        for tau in taus:
            sols, fails, nseeds = search(tau, s)
            comp = [q for q in sols if q['comp']]
            print(f"{lab:10s} tau={tau:<6} zeros found={len(sols)} (competitive {len(comp)}), seeds={nseeds}, "
                  f"failed integrations={fails} [{time.time()-t0:.0f}s]", flush=True)
            for q in sols:
                print(f"     T={q['T']:.6f} T-tau={q['T']-tau:+.6f} psi={q['psi']:.5f} cost={q['cost']:.8f} "
                      f"{'COMP' if q['comp'] else 'not comp'} chi_h={q['chi']:.4f} smin={q['smin']:.3e} cond={q['cond']:.2e} "
                      f"r=[{q['rmin']:.4f},{q['rmax']:.4f}] res={q['res']:.1e}", flush=True)
            pg.step()
    pg.close()
    print("(2) Short-pulse limit (tau >= 1/5): max |E/delta - E_lim| and min |E_lim| over t in [5,16], psi in S^1", flush=True)
    pg = Progress(8, 'arrival shooting (2): (sense, tau)')
    for s in (+1, -1):
        for tau in [0.2, 0.5, 1.0, 5.0]:
            dev, minlim = constants_6216(tau, s)
            print(f"   s={s:+d} tau={tau}: max deviation = {dev:.4e}   (min |E_lim| on the grid = {minlim:.3e})", flush=True)
            pg.step()
    pg.close()
    print(f"total time {time.time()-t0:.0f}s")
