"""tau_eq_free_arrival_shooting.py -- floating-point diagnostic, not a certificate: free-arrival shooting for the
proper-time branch on the equator of slowly rotating Kerr (main text, remark after the existence theorem for the
example: "multi-start free-arrival shooting returned one accepted solution at each of twelve sampled values of
tau in [0.12, 0.20] in each sense").
Example: M = 1, a = 1/2, Ehat = sqrt 2, x0 = (6, 0), target (6, s/50) with s = +1 (prograde) or -1 (retrograde),
A = 1 - (1/2) sin^2(pi eta/tau) on [0, tau], then 1 (main text, example on the Kerr equator with spin).
Unknowns: p(0) = (p_r, p_phi) and T. Conditions: X(T) = x_f (2), H_{A(T)}(X(T), P(T)) = 0 (free arrival).
Canonical flow Xdot = H_p, Pdot = -H_x forwards; cost c = omega(1 - beta(xdot)).
Checks: final residual, convergence, rail (beta(v) + sqrt(N)|v| = 1), p + omega beta away from 0, T > 0.
Multi-start search of branches: many starting points; distinct branches are grouped by (T, p(0)).
Sector: early if T < tau, crossing if T >= tau. Negative test: no early branch for tau <= 0.14."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'tools'))
from progress import Progress
import numpy as np, sympy as sp, time, itertools, warnings
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
warnings.filterwarnings('ignore')

M, a, E = 1.0, 0.5, np.sqrt(2.0)
r, pr, pp, al = sp.symbols('r pr pp al', real=True)
f = 1 - 2*M/r; Dl = r**2 - 2*M*r + a**2
h0rr = r**2/(f*Dl); h0pp = Dl/f**2; bp = -2*M*a/(r - 2*M)
om = al**2*f/E; N = 1/(1 - al**2*f/E**2)
xr, xp = pr, pp + om*bp
xb = xp*bp/h0pp; b2 = bp**2/h0pp; xx = xr**2/h0rr + xp**2/h0pp
H = -om + (sp.sqrt(xb**2 + (N - b2)*xx) - xb)/(N - b2)
Hp = [sp.diff(H, pr), sp.diff(H, pp)]
Hf = sp.lambdify((r, pr, pp, al), H, 'numpy')
Gf = sp.lambdify((r, pr, pp, al), [Hp[0], Hp[1], -sp.diff(H, r)], 'numpy')
costf = sp.lambdify((r, pr, pp, al), om*(1 - bp*Hp[1]), 'numpy')
railf = sp.lambdify((r, pr, pp, al), bp*Hp[1] + sp.sqrt(N)*sp.sqrt(h0rr*Hp[0]**2 + h0pp*Hp[1]**2), 'numpy')
chif = sp.lambdify((r, pr, pp, al), sp.sqrt(xr**2 + xp**2), 'numpy')

def Aprof(eta, tau):
    return 1.0 - 0.5*np.sin(np.pi*eta/tau)**2 if 0 <= eta <= tau else 1.0

def flow(p0, T, tau, dense=False):
    def rhs(eta, y):
        g = Gf(y[0], y[2], y[3], Aprof(eta, tau))
        return [g[0], g[1], g[2], 0.0, costf(y[0], y[2], y[3], Aprof(eta, tau))]
    # p_phi conserved (phi cyclic); its derivative is -H_phi = 0
    sol = solve_ivp(rhs, [0, T], [6.0, 0.0, p0[0], p0[1], 0.0], rtol=1e-11, atol=1e-13, max_step=min(tau, T)/200,
                    dense_output=dense)
    return sol

def residual(u, tau, xf):
    p0, T = u[:2], u[2]
    if T <= 1e-4: return [1e3, 1e3, 1e3]
    sol = flow(p0, T, tau)
    if not sol.success: return [1e3, 1e3, 1e3]
    y = sol.y[:, -1]
    return [y[0] - xf[0], (y[1] - xf[1])*10, Hf(y[0], y[2], y[3], Aprof(T, tau))]

def solve_all(tau, sign):
    xf = np.array([6.0, sign/50])
    found = []
    # starting points: p_phi ~ sign*4.26 scaled, small p_r, T on a grid
    for scale, prr, T0 in itertools.product([0.6, 0.8, 1.0, 1.2], [-0.3, 0.0, 0.3], [0.12, 0.15, 0.18, 0.22, 0.28, 0.4]):
        u0 = [prr, sign*4.26*scale, T0]
        u, info, ier, msg = fsolve(residual, u0, args=(tau, xf), full_output=True, xtol=1e-12)
        res = np.linalg.norm(residual(u, tau, xf))
        if ier != 1 or res > 1e-9 or u[2] <= 0: continue
        if any(abs(u[2] - g['T']) < 1e-6 and np.linalg.norm(u[:2] - g['p0']) < 1e-5 for g in found): continue
        sol = flow(u[:2], u[2], tau, dense=True)
        etas = np.linspace(0, u[2], 400)
        Y = sol.sol(etas)
        rail = max(abs(railf(Y[0, k], Y[2, k], Y[3, k], Aprof(e, tau)) - 1) for k, e in enumerate(etas))
        chi = min(chif(Y[0, k], Y[2, k], Y[3, k], Aprof(e, tau)) for k, e in enumerate(etas))
        rmin, dphi = Y[0].min(), Y[1, -1]
        found.append(dict(T=u[2], p0=u[:2].copy(), cost=sol.y[4, -1], res=res, rail=rail, chi=chi, rmin=rmin,
                          sector='early' if u[2] < tau else 'crossing', dphi=dphi))
    return sorted(found, key=lambda g: g['cost'])

if __name__ == "__main__":
    t0 = time.time()
    print("Negative test and scan 0.14 < tau < 0.20 (free-arrival shooting, all branches found):", flush=True)
    pg = Progress(24, 'free-arrival shooting: (sense, tau)')
    for sign, lab in [(+1, 'prograde'), (-1, 'retrograde')]:
        for tau in [0.12, 0.14, 0.145, 0.15, 0.155, 0.16, 0.165, 0.17, 0.175, 0.18, 0.19, 0.20]:
            br = solve_all(tau, sign)
            nearly = [g for g in br if g['sector'] == 'early']
            flag = "  !! NEGATIVE TEST FAILED" if (tau <= 0.14 and nearly) else ""
            print(f"{lab:10s} tau={tau:.3f}: {len(br)} branches{flag}", flush=True)
            for g in br:
                print(f"     T={g['T']:.6f} [{g['sector']}] cost={g['cost']:.8f} p0=({g['p0'][0]:+.5f},{g['p0'][1]:+.5f}) "
                      f"res={g['res']:.1e} rail={g['rail']:.1e} chi_min={g['chi']:.3f} r_min={g['rmin']:.4f}", flush=True)
            pg.step()
    pg.close()
    print(f"total time {time.time()-t0:.0f}s")
