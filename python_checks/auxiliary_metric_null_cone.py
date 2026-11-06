# -*- coding: utf-8 -*-
"""Main text, section "The rail as the null cone of an auxiliary metric": the rail is the null cone of g_opt
(3+1 Kerr), and the Hamiltonian H_eta of Paper II is the support function of that cone (Proposition on extremals).
Usage: python3 auxiliary_metric_null_cone.py"""
import sys
import numpy as np
import sympy as sp

FAIL = []


def check(name, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"   {detail}" if detail else ""))
    if not ok:
        FAIL.append(name)


print("symbolic: 3+1 Kerr threading")
r, th, M, a, A, E = sp.symbols('r theta M a A Ehat', positive=True)
S = r**2 + a**2*sp.cos(th)**2
D = r**2 - 2*M*r + a**2
gtt = -(1 - 2*M*r/S)
gtp = -2*M*a*r*sp.sin(th)**2/S
gpp = (r**2 + a**2 + 2*M*a**2*r*sp.sin(th)**2/S)*sp.sin(th)**2
f = -gtt
bph = gtp/f                                   # beta = (g_tphi/f) dphi
h = sp.diag(S/D, S, gpp + gtp**2/f)           # threading spatial metric (r, theta, phi)
ue, ur, ut, up = sp.symbols('u_eta u_r u_th u_ph', real=True)
ux = sp.Matrix([ur, ut, up])
# Kerr in threading form: g_K = -f (d eta - beta)^2 + h
gK = sp.zeros(4, 4)
gK[0, 0] = gtt; gK[0, 3] = gK[3, 0] = gtp
gK[1, 1] = S/D; gK[2, 2] = S; gK[3, 3] = gpp
U = sp.Matrix([ue, ur, ut, up])
thr = -f*(ue - bph*up)**2 + (ux.T*h*ux)[0]
check("g_Kerr = -f(d eta - beta)^2 + h", sp.simplify((U.T*gK*U)[0] - thr) == 0)
check("h_phiphi = Delta sin^2/f", sp.simplify(h[3-1, 3-1] - D*sp.sin(th)**2/f) == 0)

# rail: -g(u, d_eta) = Ehat and g(u,u) = -1, with g = A^2 g_K
w = sp.symbols('w', positive=True)            # w = u^eta - beta(u)
hu = sp.symbols('hu', nonnegative=True)       # h(u,u)
Z = A**2*f/E**2
N = 1/(1 - Z)
# -g(u,d_eta) = A^2 f w  ;  g(u,u) = A^2(-f w^2 + hu)
wsol = E/(A**2*f)
husol = sp.solve(sp.Eq(A**2*(-f*wsol**2 + hu), -1), hu)[0]
check("rail + normalisation => w^2 = N h(u,u)/f  (null for g_opt)",
      sp.simplify(wsol**2 - N*husol/f) == 0)
# converse: v future null for g_opt, w_v = sqrt(N h0(v,v)); rescale lambda so that rail holds
wv, h0v = sp.symbols('w_v h0v', positive=True)
lam = E/(A**2*f*wv)
gnorm = A**2*lam**2*(-f*wv**2 + f*h0v)
check("converse: rescaled g_opt-null vector is g-unit",
      sp.simplify(gnorm.subs(h0v, wv**2/N) + 1) == 0)
check("n^2 = N (Paper II Randers index)",
      sp.simplify((E/sp.sqrt(E**2 - A**2*f))**2 - N) == 0)

b2 = bph**2/(h[2, 2]/f)                      # |beta|^2_{h0}
check("f g_phiphi + g_tphi^2 = Delta sin^2", sp.simplify(f*gpp + gtp**2 - D*sp.sin(th)**2) == 0)
check("1 - |beta|^2_h0 = f g_phiphi/(Delta sin^2)  (so |beta|_h0<1 iff f>0)",
      sp.simplify(1 - b2 - f*gpp/(D*sp.sin(th)**2)) == 0)
# inverse of g_opt restricted to (eta, phi) block (r, theta decouple): g^{eta eta}
Nn = sp.symbols('N', positive=True)
Gb = sp.Matrix([[-1, bph], [bph, Nn*h[2, 2]/f - bph**2]])
check("g_opt^{eta eta} = -1 + |beta|^2_h0/N  (d eta timelike since N>=1>|beta|^2)",
      sp.simplify(Gb.inv()[0, 0] - (-1 + b2/Nn)) == 0)

print("\nsymbolic: equatorial data agree with Paper II eq. (randers-data)")
fe = 1 - 2*M/r; De = r**2 - 2*M*r + a**2
ev = {th: sp.pi/2}
check("h0_rr = r^2/(f Delta)", sp.simplify((h[0, 0]/f).subs(ev) - r**2/(fe*De)) == 0)
check("h0_phiphi = Delta/f^2", sp.simplify((h[2, 2]/f).subs(ev) - De/fe**2) == 0)
check("beta = -2Ma/(r f) dphi", sp.simplify(bph.subs(ev) + 2*M*a/(r*fe)) == 0)
check("|beta|_h0 = 2M|a|/(r sqrt Delta) at equator",
      sp.simplify((bph**2/(h[2, 2]/f)).subs(ev) - 4*M**2*a**2/(r**2*De)) == 0)

print("\nnumeric: Paper II H_eta + 1 is the support function of the g_opt cone")
rng = np.random.default_rng(7)
worst = 0.0
for _ in range(5000):
    Mv = 1.0; av = rng.uniform(-1, 1); rv = 2 + 10**rng.uniform(-3, 1.5)
    Av = rng.uniform(0.3, 1.0); Ev = Av*rng.uniform(1.01, 3.0)
    fv = 1 - 2*Mv/rv; Dv = rv*rv - 2*Mv*rv + av*av
    vb2 = 1 - Av**2*fv/Ev**2
    P = rv*rv + av*av + 2*Mv*av*av/rv
    Pb = P + Av**2*(2*Mv*av/rv)**2/Ev**2
    ph0 = (2*Mv*av/rv)*vb2/Pb
    R = np.sqrt(Dv*vb2/Pb)
    pr, pp = rng.normal(size=2)*rng.uniform(0.1, 10)
    hsup = pp*ph0 + R*np.sqrt((Dv/rv**2)*pr**2 + pp**2/Pb)     # = H_eta + 1
    # g_opt (equatorial) and its inverse
    Nv = 1/vb2; P0 = rv*rv/(fv*Dv); Q0 = Dv/fv**2; bv = -2*Mv*av/(rv*fv)
    G = np.array([[-1, 0, bv], [0, Nv*P0, 0], [bv, 0, Nv*Q0 - bv*bv]])
    Pcov = np.array([-hsup, pr, pp])
    q = Pcov @ np.linalg.solve(G, Pcov)
    scale = hsup**2 + pr**2/(Nv*P0) + pp**2/(Nv*Q0)
    worst = max(worst, abs(q)/scale)
check("g_opt^{-1}(P,P) = 0 for P = (-(H_eta+1), p_r, p_phi), 5000 random samples",
      worst < 1e-10, f"max rel {worst:.1e}")

print()
if FAIL:
    print(f"{len(FAIL)} FAILED:", *FAIL, sep="\n  ")
    sys.exit(1)
print("all checks passed")
