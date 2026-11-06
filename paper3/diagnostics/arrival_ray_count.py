"""arrival_ray_count.py -- floating-point diagnostic, not a certificate. Arrival branch of conformal Schwarzschild,
variable A, q(eta) = Ehat^2/A^2 >= 3/2 (main text, theorem on the swept angle for variable factors). Null geodesics of
the representative Gb = -V^2 deta^2 + h0 (the same curves as for g_opt), launched from (0, r0) with direction alpha0
(cos alpha = r_eta/(f V)). For each delta: all same-sense geodesics found by the scan that cross r = r_T at angle delta
(bracketing on a fine grid of alpha0 for each crossing index + brentq); prediction of the theorem: exactly one for delta
in (0, pi). Sign of D_J = r_a phi_l - r_l phi_a along the realising ray (exact sensitivities via SymPy): no zero before pi.
Profiles (Ehat^2 = 2, M = 1): constant A = 1 (q = 2); A from 1 to 1.15, fast (tanh of width 0.5 at eta = 5: q from 2 to
1.512); A from 1.15 to 1 (q from 1.512 to 2); oscillating A = 1.075 + 0.075 sin(eta/2) (q between 1.512 and 2)."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'tools'))
from progress import Progress
import numpy as np, sympy as sp, sys
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
eta, r, Pe, Pr, L, q, qp, qpp = sp.symbols('eta r P_eta P_r L q qp qpp', real=True)
M = 1; f = 1 - 2*M/r; V2 = 1 - f/q
H = sp.Rational(1, 2)*(-Pe**2/V2 + f**2*Pr**2 + f*L**2/r**2)
state = [eta, r, sp.Symbol('phi'), Pe, Pr, L]
rhs = [sp.diff(H, Pe), sp.diff(H, Pr), sp.diff(H, L), -qp*sp.diff(H, q), -sp.diff(H, r), sp.Integer(0)]
def dtot(ex): return sp.diff(ex, q)*qp + sp.diff(ex, qp)*qpp
Jac = sp.Matrix([[(dtot(fi) if s_ is eta else sp.diff(fi, s_)) for s_ in state] for fi in rhs])
args = (r, Pe, Pr, L, q, qp, qpp)
Frhs = sp.lambdify(args, rhs, 'numpy'); FJ = sp.lambdify(args, Jac, 'numpy')
E2 = 2.0
def prof(kind):
    if kind == 'const': A = (lambda e: 1.0, lambda e: 0.0, lambda e: 0.0)
    elif kind == 'up':   A = (lambda e: 1.075 + 0.075*np.tanh((e - 5)/0.5), lambda e: 0.15/np.cosh((e - 5)/0.5)**2, lambda e: -0.6*np.tanh((e - 5)/0.5)/np.cosh((e - 5)/0.5)**2)
    elif kind == 'down': A = (lambda e: 1.075 - 0.075*np.tanh((e - 5)/0.5), lambda e: -0.15/np.cosh((e - 5)/0.5)**2, lambda e: 0.6*np.tanh((e - 5)/0.5)/np.cosh((e - 5)/0.5)**2)
    else:                A = (lambda e: 1.075 + 0.075*np.sin(e/2), lambda e: 0.0375*np.cos(e/2), lambda e: -0.01875*np.sin(e/2))
    a, da, dda = A
    qf = lambda e: E2/a(e)**2; qd = lambda e: -2*E2*da(e)/a(e)**3; qdd = lambda e: E2*(6*da(e)**2/a(e)**4 - 2*dda(e)/a(e)**3)
    return qf, qd, qdd
def launch(al, r0, qf):
    f0 = 1 - 2/r0; q0 = qf(0.0); V0 = np.sqrt(1 - f0/q0)
    pr = np.cos(al)/f0; Lv = r0*np.sin(al)/np.sqrt(f0)                 # spatial covector of unit h0-norm
    pe = -V0*np.sqrt(f0**2*pr**2 + f0*Lv**2/r0**2)
    dpr = -np.sin(al)/f0; dL = r0*np.cos(al)/np.sqrt(f0); dpe = 0.0       # |p|_h0 = 1 for all alpha
    return np.array([0.0, r0, 0.0, pe, pr, Lv]), np.array([0, 0, 0, dpe, dpr, dL], float)
def shoot(al, r0, Q, rt=None, jac=False, dense=False, lam_max=400.0, emax=300.0):
    qf, qd, qdd = Q; y0, J0 = launch(al, r0, qf)
    def Fn(l, Y):
        a = (Y[1], Y[3], Y[4], Y[5], qf(Y[0]), qd(Y[0]), qdd(Y[0]))
        out = np.array(Frhs(*a), float)
        if jac: out = np.concatenate([out, np.array(FJ(*a), float) @ Y[6:12]])
        return out
    evs = []
    if rt is not None: evs.append(lambda l, Y: Y[1] - rt)
    hor = lambda l, Y: Y[1] - 2.0 - 1e-7; hor.terminal = True
    late = lambda l, Y: Y[0] - emax; late.terminal = True
    wind = lambda l, Y: Y[2] - 2*np.pi; wind.terminal = True
    evs += [hor, late, wind]
    Y0 = np.concatenate([y0, J0]) if jac else y0
    s = solve_ivp(Fn, (0, lam_max), Y0, events=evs, rtol=1e-10, atol=1e-12, dense_output=dense)
    cr = [(l, Y) for l, Y in zip(s.t_events[0], s.y_events[0]) if l > 1e-6] if rt is not None else []
    return s, cr
def run(kind, r0, rt, nd=40, n=1601):
    Q = prof(kind)
    als = np.linspace(1e-4, np.pi - 1e-4, n)
    data = [shoot(a, r0, Q, rt)[1] for a in als]
    cnt = [len(c) for c in data]; extra = []
    for i in range(n - 1):
        if cnt[i] != cnt[i + 1]:
            a_, b_, ca = als[i], als[i + 1], cnt[i]
            for _ in range(40):
                mid = (a_ + b_)/2
                if len(shoot(mid, r0, Q, rt)[1]) == ca: a_ = mid
                else: b_ = mid
            extra += [a_ - 10.0**(-j) for j in range(3, 13)] + [b_ + 10.0**(-j) for j in range(3, 13)] + [a_, b_]
    if extra:
        als = np.unique(np.concatenate([als, np.clip(extra, 1e-4, np.pi - 1e-4)])); data = [shoot(a, r0, Q, rt)[1] for a in als]
    deltas = np.linspace(0.02, 0.98, nd)*np.pi; counts = []; dzeros = 0
    for d in deltas:
        sols = []
        for k in range(max(len(c) for c in data)):
            for i in range(len(als) - 1):
                a, b = data[i], data[i + 1]
                if len(a) > k and len(b) > k and (a[k][1][2] - d)*(b[k][1][2] - d) < 0:
                    g = lambda x: (lambda c: c[k][1][2] - d if len(c) > k else np.nan)(shoot(x, r0, Q, rt)[1])
                    try: x = brentq(g, als[i], als[i + 1], xtol=1e-13)
                    except ValueError: continue
                    c = shoot(x, r0, Q, rt)[1]
                    if len(c) > k and not any(abs(x - z[1]) < 1e-8 for z in sols): sols.append((c[k][1][0], x, c[k][0]))
        counts.append(len(sols))
        if sols:
            v, x, lam = min(sols)
            s, _ = shoot(x, r0, Q, None, jac=True, dense=True, lam_max=lam)
            ls = np.linspace(0, s.t[-1], 3001)[30:]; Y = s.sol(ls); Dv = []
            for j in range(len(ls)):
                a = (Y[1, j], Y[3, j], Y[4, j], Y[5, j], Q[0](Y[0, j]), Q[1](Y[0, j]), Q[2](Y[0, j]))
                k_ = np.array(Frhs(*a), float); J = Y[6:12, j]; Dv.append(J[1]*k_[2] - k_[1]*J[2])
            if np.any(np.diff(np.sign(Dv)) != 0): dzeros += 1
    qf = Q[0]; qs = [qf(e) for e in np.linspace(0, 60, 601)]
    print(f"[{kind}] q in [{min(qs):.3f}, {max(qs):.3f}], r0 = {r0}, r_T = {rt}: counts of same-sense rays over {nd} deltas in "
          f"[0.02, 0.98] pi: {sorted(set(counts))}{' (at ' + ', '.join(f'{d/np.pi:.3f}' for d, c in zip(deltas, counts) if c != 1) + ')' if any(c != 1 for c in counts) else ''}; "
          f"realising rays with a zero of D_J: {dzeros}", flush=True)
if __name__ == "__main__":
    pg = Progress(3*len(sys.argv[1:] or ['const', 'up', 'down', 'osc']), 'arrival_ray_count: (r_T, profile)')
    for rt in (6.0, 10.0, 4.0):
        for kind in (sys.argv[1:] or ['const', 'up', 'down', 'osc']): run(kind, 6.0, rt); pg.step()
    pg.close()
