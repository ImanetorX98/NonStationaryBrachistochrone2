"""tau_eq_injectivity_scan.py -- floating-point guide, not a certificate: exploration of the injectivity of the
shooting map on the rectangle Q_s (main text, triangular criterion and theorem on global uniqueness in the example).
Reordered shooting map R_s = (-s E_phi, s E_r) on Q_s: T in [3/25, 21/100], psi = s pi/2 + theta, |theta| <= 1/5.
For each (tau, s) and each grid node: the actual DR (variational equations of tau_eq_arrival_shooting.py), the DR of the
reference (geometry frozen at r = 6, exact A_tau, constant p), the deviation ||DR - DR_ref||_2.
Formulations compared:
  (a) parameters (T, theta), Sym(DR) >= kappa I
  (b) parameters (T, z), z = T sin theta: convex domain (trapezium), Sym(DR~)
  (c) as (a)/(b) with a constant diagonal weight W = diag(1, w): Sym(W DR) >= kappa I (w optimised on the sample)
Recorded: minimum eigenvalue and its position, margin of the reference, maximum deviation, failed integrations, and the
quantities of the triangular criterion. Writes ../results/tau_eq_injectivity_data.npy."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'tools'))
from progress import Progress
import numpy as np, sympy as sp, time, importlib.util, os, warnings
from scipy.integrate import solve_ivp
warnings.filterwarnings('ignore')
here = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('t64', os.path.join(here, 'tau_eq_arrival_shooting.py'))
t64 = importlib.util.module_from_spec(spec); spec.loader.exec_module(t64)
Gf, J3f, termf, dterm_dal, Aprof, dAprof = t64.Gf, t64.J3f, t64.termf, t64.dterm_dal, t64.Aprof, t64.dAprof

def ref_jac(T, psi, tau, s):
    """reference: r frozen at 6, constant p. X(0) = x_f - int_0^T v(p(T); A(eta)) d eta.
    DR_ref: psi column = -int dv/dp . dp/dpsi; T column = -v(T) - int dv/dp . dp/dT (through A(T))."""
    al_T = Aprof(T, tau); cr, cp = termf(6.0, al_T)
    p = (cr*np.cos(psi), cp*np.sin(psi))
    dps = np.array([-cr*np.sin(psi), cp*np.cos(psi)])
    dcr, dcp = dterm_dal(6.0, al_T); dA = dAprof(T, tau)
    dT = np.array([dcr*dA*np.cos(psi), dcp*dA*np.sin(psi)])
    def rhs(eta, y):
        alp = Aprof(eta, tau)
        g = Gf(6.0, p[0], p[1], alp); J = np.array(J3f(6.0, p[0], p[1], alp), dtype=float)
        Jv = J[:2, 1:3]                                  # d(rdot,phidot)/d(pr,pphi)
        return [g[0], g[1], *(Jv @ dps), *(Jv @ dT)]
    ts = [T, tau, 0.0] if T > tau else [T, 0.0]
    y = np.zeros(6)
    for k in range(len(ts) - 1):
        sol = solve_ivp(rhs, [ts[k], ts[k+1]], y, rtol=1e-11, atol=1e-13,
                        max_step=(min(T, tau)/150 if ts[k+1] == 0.0 else np.inf))
        if not sol.success: return None
        y = sol.y[:, -1]
    # backward integrals: y = -int_0^T (...)
    res = np.array([y[0], s/50 + y[1]])            # X(0) - x0 = x_f + (backward integral) - x0
    gT = Gf(6.0, p[0], p[1], al_T)
    JT = np.array([y[4], y[5]]) - np.array([gT[0], gT[1]])
    Jp = np.array([y[2], y[3]])
    return dict(res=res, J=np.column_stack([JT, Jp]))

def reorder(J, s):
    return np.array([[0, -s], [s, 0]]) @ J

def to_z(DR, T, th):
    """parameters (T, z), z = T sin theta: d/dT|_z = d_T - (tan th / T) d_th ; d/dz = d_th/(T cos th)"""
    C = np.array([[1.0, 0.0], [-np.tan(th)/T, 1.0/(T*np.cos(th))]])
    return DR @ C

def lam(M):
    return np.linalg.eigvalsh((M + M.T)/2)[0]

def reorder_res(res, s):
    return np.array([-s*res[1], s*res[0]])

if __name__ == "__main__":
    # triangular criterion (alternative to (a)): (i) d_theta R2 > 0 on Q_s; (ii) R2(T,-1/5) < 0 < R2(T,1/5) for every T;
    # (iii) det DR > 0 on the curve {R2 = 0}. Then theta*(T) is C^1 on the whole [3/25,21/100] and g(T)=R1(T,theta*(T))
    # has g' = det DR / d_theta R2 > 0: at most one zero. The minimum of det over the whole Q_s is recorded too.
    t0 = time.time()
    taus = [0.02, 0.06, 0.1, 0.14, 0.155, 0.16, 0.165, 0.17, 0.175, 0.185, 0.2, 0.3, 0.5, 1.0, 2.0, 5.0]
    Ts = np.linspace(3/25, 21/100, 9); ths = np.linspace(-0.2, 0.2, 9)
    ws = np.exp(np.linspace(np.log(0.05), np.log(20), 61))
    store = {}
    print("Exploration on Q_s; 9x9 grid in (T,theta)", flush=True)
    pg = Progress(2*len(taus), 'injectivity scan: (sense, tau)')
    for s in (+1, -1):
        for tau in taus:
            fails = 0; recs = []
            for T in Ts:
                for th in ths:
                    psi = s*np.pi/2 + th
                    o = t64.backward(T, psi, tau, s, want_jac=True)
                    q = ref_jac(T, psi, tau, s)
                    if o is None or q is None or not np.all(np.isfinite(o['J'])): fails += 1; continue
                    recs.append((T, th, reorder(o['J'], s), reorder(q['J'], s), reorder_res(o['res'], s),
                                 reorder_res(q['res'], s)))
            store[(s, tau)] = recs
            la = min((lam(D), T, th) for T, th, D, Dr, R, Rr in recs)
            laref = min(lam(Dr) for T, th, D, Dr, R, Rr in recs)
            dev = max(np.linalg.norm(D - Dr, 2) for T, th, D, Dr, R, Rr in recs)
            bw = max((min(lam(np.diag([1, w]) @ D) for T, th, D, Dr, R, Rr in recs), w) for w in ws)
            d22 = min(D[1, 1] for T, th, D, Dr, R, Rr in recs)
            d11 = min(D[0, 0] for T, th, D, Dr, R, Rr in recs)
            det = min(np.linalg.det(D) for T, th, D, Dr, R, Rr in recs)
            bp = min(R[1] for T, th, D, Dr, R, Rr in recs if th > 0.199)
            bm = max(R[1] for T, th, D, Dr, R, Rr in recs if th < -0.199)
            r20 = max(abs(R[1]) for T, th, D, Dr, R, Rr in recs if abs(th) < 1e-9)
            print(f"s={s:+d} tau={tau:<5} fail={fails} | (a) lam={la[0]:.4f}@(T={la[1]:.3f},th={la[2]:+.2f}) ref={laref:.4f} "
                  f"weighted={bw[0]:.4f}(w={bw[1]:.2f}) dev={dev:.4f} | triang: min d_thR2={d22:.4f} min d_TR1={d11:.4f} "
                  f"min det={det:.5f} R2(+1/5)>={bp:.5f} R2(-1/5)<={bm:.5f} max|R2(T,0)|={r20:.2e} [{time.time()-t0:.0f}s]",
                  flush=True)
            pg.step()
    pg.close()
    np.save(os.path.join(here, '..', 'results', 'tau_eq_injectivity_data.npy'), np.array([(k, v) for k, v in store.items()], dtype=object),
            allow_pickle=True)
    print(f"total time {time.time()-t0:.0f}s")
