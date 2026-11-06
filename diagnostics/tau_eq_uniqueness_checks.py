"""tau_eq_uniqueness_checks.py -- floating-point checks, not a certificate, that COULD contradict the proof of the
theorem on global uniqueness in the example on the Kerr equator with spin (main text: "Floating-point checks, not
certified, give D_{2theta}/T >= 0.52 on the sampled shots in Q_s, and Gamma' >= 0.10 on the separately sampled
radial-return curves F_2 = 0").
(A) table of coefficients sampled on the two domains: the whole Q_s (shell, |p|/j < 77/1000, j in (8/9, sqrt291/4])
    and return (|r - 6| < 1/50, |p| <= 63/2000); A in [1/2, 1]; both signs.
(B) rank-one angular identity: dphi' = (H_rJ + Z H_rp) u - Z u' along actual shots.
(C) theta column on Q_s: sup|U| < 0.213, sup|V| < 0.223, min V > 0.189, dU/dt_b > 0.23, D22 > 0.23 T.
(D) on the curve R2 = 0: |D12|/D22 < 6/125; T column: sup|u| < 0.7, sup|P| < 0.47, int|u'| < 0.61,
    D11 > 0.07, |D21| < 0.7; Gamma' = D11 - D12 D21/D22 > 7/200."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'tools'))
from progress import Progress
import numpy as np, sympy as sp, importlib.util, os, time, warnings
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
warnings.filterwarnings('ignore')
here = os.path.dirname(os.path.abspath(__file__))
def load(name, fn):
    spec = importlib.util.spec_from_file_location(name, os.path.join(here, fn))
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m
t = load('t64', 'tau_eq_arrival_shooting.py')
ok = lambda lab, c: print(('PASS ' if c else 'FAIL ') + lab, flush=True)
r, pr, J, al = t.r, t.pr, t.pp, t.al; H = t.H
g = al**2/(sp.sqrt(2)*r); xi = J - g
Hp = sp.diff(H, pr); HJ = sp.diff(H, J)
ex = dict(Hpr=sp.diff(H, pr, r), Hpp=sp.diff(H, pr, 2), HpJ=sp.diff(H, pr, J), HJJ=sp.diff(H, J, 2),
          Hrr=sp.diff(H, r, 2), HrJ=sp.diff(H, r, J), Hr=sp.diff(H, r), phid=HJ, rdot=Hp,
          omr=sp.diff(t.om, r), xir=sp.diff(xi, r), g=g, xi=xi)
fn = {k: sp.lambdify((r, pr, J, al), v, 'numpy') for k, v in ex.items()}
def coeffs(R, P, JJ, A):
    q = {k: f(R, P, JJ, A)*np.ones(np.broadcast(R, P, JJ, A).shape) for k, f in fn.items()}
    j = np.abs(JJ); Z = P/q['xi']; eps = q['g']/j
    return dict(A=j*q['Hpp'], B=q['Hpr'], L=q['Hrr']/j, HrJ=q['HrJ'],
                f=(q['omr'] - q['xir']*q['phid'] + q['g']*q['HrJ'])/j, q=j*q['Hpp']*Z*eps, Hrj=q['Hr']/j, Z=Z)

# ---- (A)
print("(A) table of coefficients, sampled maxima/minima", flush=True)
tab = {'glob': dict(Alo=29/20, Ahi=48/11, B=7/50, L=1/25, HrJ=4/125, f=57/1000, q=11/200, Hrj=2/25),
       'ret': dict(Alo=37/20, Ahi=48/11, B=8/125, L=1/25, HrJ=13/500, f=7/125, q=13/500, Hrj=2/25)}
for dom in ('glob', 'ret'):
    for s in (1, -1):
        if dom == 'glob':
            R, X, JJ, A = np.meshgrid(np.linspace(23/4, 25/4, 31), np.linspace(-1, 1, 31),
                                      np.linspace(8/9 + 1e-9, np.sqrt(291)/4, 41), np.linspace(.5, 1, 11), indexing='ij')
            P = X*0.077*JJ                                       # |p|/j < 77/1000
        else:
            R, P, JJ, A = np.meshgrid(np.linspace(6 - 1/50, 6 + 1/50, 21), np.linspace(-63/2000, 63/2000, 31),
                                      np.linspace(8/9 + 1e-9, np.sqrt(291)/4, 41), np.linspace(.5, 1, 11), indexing='ij')
        c = coeffs(R, P, s*JJ, A); b = tab[dom]
        res = dict(Alo=c['A'].min(), Ahi=c['A'].max(), B=abs(c['B']).max(), L=abs(c['L']).max(),
                   HrJ=abs(c['HrJ']).max(), f=abs(c['f']).max(), q=abs(c['q']).max(), Hrj=abs(c['Hrj']).max())
        good = res['Alo'] > b['Alo'] and all(res[k] < b[k] for k in b if k != 'Alo')
        ok(f"   {dom} s={s:+d}: " + ", ".join(f"{k}={res[k]:.4f}({'>' if k=='Alo' else '<'}{b[k]:.4f})" for k in b), good)

# ---- integration with projective and angular variational equations
def shoot(T, th, tau, s, col, n=600):
    psi = s*np.pi/2 + th; alT = t.Aprof(T, tau); cr, cp = t.termf(6.0, alT)
    p0, JJ = cr*np.cos(psi), cp*np.sin(psi); j = abs(JJ)
    G = t.Gf(6.0, p0, JJ, alT)
    if col == 'theta':
        u0, dp0, dJ, dphi0 = 0.0, -cr*np.sin(psi), cp*np.cos(psi), 0.0
    else:
        dcr, dcp = t.dterm_dal(6.0, alT); dA = t.dAprof(T, tau)
        u0, dp0, dJ, dphi0 = -G[0], -G[2] + dcr*dA*np.cos(psi), dcp*dA*np.sin(psi), -G[1]
    h = dJ/JJ
    def rhs(eta, y):
        rr, p, u, dp, dphi = y; A = t.Aprof(eta, tau)
        M = np.array(t.J3f(rr, p, JJ, A), dtype=float); v = M @ np.array([u, dp, dJ])
        gg = t.Gf(rr, p, JJ, A)
        return [gg[0], gg[2], v[0], v[2], v[1]]
    sol = solve_ivp(rhs, [T, 0.0], [6.0, p0, u0, dp0, dphi0], rtol=1e-11, atol=1e-13, max_step=T/n)
    rr, p, u, dp, dphi = sol.y
    P = (dp - p/JJ*dJ)/j
    A_ = np.array([t.Aprof(e, tau) for e in sol.t])
    c = coeffs(rr, p, JJ*np.ones_like(rr), A_)
    up = c['B']*u + c['A']*P - c['q']*h                              # u' (in eta) from the projective system
    q = {k: fn[k](rr, p, JJ, A_) for k in ('HrJ', 'Hpr')}
    dphi_id = (q['HrJ'] + c['Z']*q['Hpr'])*u - c['Z']*up              # angular identity
    M_all = [np.array(t.J3f(a1, a2, JJ, a3), dtype=float) for a1, a2, a3 in zip(rr, p, A_)]
    dphi_orig = np.array([M[1] @ np.array([a, b_, dJ]) for M, a, b_ in zip(M_all, u, dp)])
    uprime_orig = np.array([M[0] @ np.array([a, b_, dJ]) for M, a, b_ in zip(M_all, u, dp)])
    tv = np.trapezoid(np.abs(uprime_orig[::-1]), sol.t[::-1])
    return dict(t=sol.t, u=u, P=P, s=s, h=h, dphi0=dphi[-1], u0=u[-1], j=j, iderr=np.max(np.abs(dphi_id - dphi_orig)),
                uperr=np.max(np.abs(up - uprime_orig)), up=uprime_orig, tv=tv, phid=G[1])

# ---- (B)+(C)
print("(B)+(C) theta column on Q_s", flush=True)
t0 = time.time(); W = dict(U=0, V=0, Vmin=9, slope=9, D22T=9, iderr=0, uperr=0)
pg = Progress(12, 'uniqueness (B)+(C): (sense, tau)')
for s in (1, -1):
    for tau in [0.02, 0.16, 0.17, 0.5, 1.0, 5.0]:
        for T in np.linspace(3/25, 21/100, 4):
            for th in np.linspace(-0.2, 0.2, 5):
                o = shoot(T, th, tau, s, 'theta')
                U, V = s*o['u'], -s*o['P']
                W['U'] = max(W['U'], abs(U).max()); W['V'] = max(W['V'], abs(V).max()); W['Vmin'] = min(W['Vmin'], V.min())
                W['slope'] = min(W['slope'], (s*(-o['up'])).min())     # dU/dt_b = -s u'
                W['D22T'] = min(W['D22T'], s*o['u0']/T)
                W['iderr'] = max(W['iderr'], o['iderr']); W['uperr'] = max(W['uperr'], o['uperr'])
        pg.step()
pg.close()
print(f"   [{time.time()-t0:.0f}s]")
ok(f"   projective u' equals the variational one (max err {W['uperr']:.1e}); angular identity (max err {W['iderr']:.1e})",
   W['uperr'] < 1e-10 and W['iderr'] < 1e-10)
ok(f"   sup|U| = {W['U']:.4f} < 0.213, sup|V| = {W['V']:.4f} < 0.223", W['U'] < 0.213 and W['V'] < 0.223)
ok(f"   min V = {W['Vmin']:.4f} > 0.189, min dU/dt_b = {W['slope']:.4f} > 0.23, min D22/T = {W['D22T']:.4f} > 0.23",
   W['Vmin'] > 0.189 and W['slope'] > 0.23 and W['D22T'] > 0.23)

# ---- (D)
print("(D) curve R2 = 0: ratio D12/D22 and T column", flush=True)
W = dict(rat=0, UT=0, VT=0, NT=0, D11=9, D21=0, gam=9)
pg = Progress(16, 'uniqueness (D): (sense, tau)')
for s in (1, -1):
    for tau in [0.02, 0.1, 0.16, 0.17, 0.2, 0.5, 1.0, 5.0]:
        for T in np.linspace(3/25, 21/100, 5):
            R2 = lambda th: s*t.backward(T, s*np.pi/2 + th, tau, s)['res'][0]
            ths = brentq(R2, -0.2, 0.2, xtol=1e-14)
            a_ = shoot(T, ths, tau, s, 'theta'); b_ = shoot(T, ths, tau, s, 'T')
            D22 = s*a_['u0']; D12 = -s*a_['dphi0']; D21 = s*b_['u0']; D11 = -s*b_['dphi0']
            W['rat'] = max(W['rat'], abs(D12)/D22); W['UT'] = max(W['UT'], abs(b_['u']).max())
            W['VT'] = max(W['VT'], abs(b_['P']).max()); W['NT'] = max(W['NT'], b_['tv'])
            W['D11'] = min(W['D11'], D11); W['D21'] = max(W['D21'], abs(D21)); W['gam'] = min(W['gam'], D11 - D12*D21/D22)
        pg.step()
    print(f"   s={s:+d} done [{time.time()-t0:.0f}s]", flush=True)
pg.close()
ok(f"   |D12|/D22 = {W['rat']:.5f} < 6/125 = 0.048", W['rat'] < 0.048)
ok(f"   T column: sup|u| = {W['UT']:.4f} < 0.7, sup|P| = {W['VT']:.4f} < 0.47, int|u'| = {W['NT']:.4f} < 0.61",
   W['UT'] < 0.7 and W['VT'] < 0.47 and W['NT'] < 0.61)
ok(f"   D11 = {W['D11']:.4f} > 0.07, |D21| = {W['D21']:.4f} < 0.7, Gamma' = {W['gam']:.4f} > 0.035",
   W['D11'] > 0.07 and W['D21'] < 0.7 and W['gam'] > 0.035)
