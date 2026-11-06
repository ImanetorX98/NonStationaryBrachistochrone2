"""Floating-point diagnostic, not a certificate: the symmetric case r0 = r_T = d of conjugate points before pi for a
constant q < 3/2 (main text, proposition on sharpness at q = 3/2 and its remark; "q = 1.2, r0 = r_T = 1280M,
Phi = 0.949 pi" is line (3) at d = 1280). M = 1.
Coordinate rho = r^2 sqrt(q)/sqrt((r-2)((q-1)r+2)) (outer branch r > r_c), periastron at rho = b.
Phi(b; d) = angle swept from r = d (incoming) back to r = d.
Checks:
 (1) dphi/dr = b sqrt((q-1)r+2)/(sqrt(r) sqrt(Q4)), Q4 = q r^4 - b^2 (r-2)((q-1)r+2), against the equations of
     conjugate_points_static.py;
 (2) exact decomposition Phi = 2 arccos(b/rho_d) + 2 A1 sqrt(rho_d^2-b^2)/(b rho_d) + 2 int b G/(rho sqrt(rho^2-b^2));
 (3) zero of d_b Phi against the prediction b ~ sqrt(|A1| rho_d), Phi ~ pi - 4 sqrt(|A1|/rho_d);
 (4) zero of d_b Phi = zero of the in-plane mode (ODE of conjugate_points_static.py) on return to r = d;
 (5) the main part of Phi_b and the size of the remainder.
Usage: python3 conjugate_symmetric_asymptotics.py"""
import numpy as np, sys
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq
sys.path.insert(0, '.')
from conjugate_points_static import rhs, Fr
def rho(r, q): return r**2*np.sqrt(q)/np.sqrt((r-2)*((q-1)*r+2))
def A1(q): return np.sqrt(q)*(2*q-3)/(q-1)**1.5
def rp_of_b(b, q):  # periastron on the outer branch: rho(r) = b
    rc = brentq(lambda r: (q-1)*r**2 + (6-3*q)*r - 8, 2.0001, 3.0)
    return brentq(lambda r: rho(r, q) - b, rc, 1e12)
def dphidr(r, b, q):
    Q4 = q*r**4 - b**2*(r-2)*((q-1)*r+2)
    return b*np.sqrt((q-1)*r+2)/(np.sqrt(r)*np.sqrt(Q4))
def Phi(b, d, q):
    rp = rp_of_b(b, q)
    # Q4 = (r - rp) Q3 exactly (synthetic division); r = rp + s^2 removes the square-root singularity
    c4 = np.array([q, 0., 0., 0., 0.]) - b*b*np.array([0., 0., q-1, 2-2*(q-1), -4.])
    Q3, _ = np.polydiv(c4, np.array([1., -rp]))
    def g(s):
        r = rp + s*s
        return 2*b*np.sqrt((q-1)*r+2)/(np.sqrt(r)*np.sqrt(np.polyval(Q3, r)))
    val, err = quad(g, 0, np.sqrt(d - rp), limit=400, epsabs=1e-13, epsrel=1e-13)
    return 2*val
def G_of_rho(x, q):
    r = brentq(lambda r: rho(r, q) - x, 3.0, 1e14)
    # F = (dr/drho)/(fV); dr/drho = 1/rho_r
    h = 1e-6*r
    rr = (rho(r+h, q) - rho(r-h, q))/(2*h)
    f = 1 - 2/r; V = np.sqrt(1 - f/q)
    return 1/(f*V*rr) - 1 - A1(q)/x
if __name__ == '__main__':
    q = 1.2
    # (1)
    for r, b in [(50., 30.), (400., 100.), (3000., 250.)]:
        a = np.arcsin(b/rho(r, q)); y = rhs(0, [r, 0, a, 0, 0], q)
        print(f"(1) r={r} b={b}: dphi/dr formula {dphidr(r,b,q):.12e} vs ODE {y[1]/y[0]:.12e}")
    # (2)
    for d, b in [(1280., 150.), (10000., 400.)]:
        rd = rho(d, q)
        main = 2*np.arccos(b/rd) + 2*A1(q)*np.sqrt(rd**2-b**2)/(b*rd)
        rest, _ = quad(lambda t: 2*G_of_rho(b*t, q)/(t*np.sqrt(t*t-1)), 1, rd/b, limit=400)
        print(f"(2) d={d} b={b}: Phi direct {Phi(b,d,q):.10f}  decomposition {main+rest:.10f}  (remainder {rest:.3e})")
    # (3) and (4)
    for d in [1280., 5000., 20000., 1e5]:
        rd = rho(d, q); h = 1e-3
        dPhi = lambda b: (Phi(b+h, d, q) - Phi(b-h, d, q))/(2*h)
        bs = np.geomspace(20, 0.5*rd, 60); vals = [dPhi(b) for b in bs]
        zs = [brentq(dPhi, bs[i], bs[i+1], xtol=1e-9) for i in range(len(bs)-1) if np.sign(vals[i]) != np.sign(vals[i+1])]
        bz = zs[-1]  # large root
        pred_b = np.sqrt(abs(A1(q))*rd); pred_P = np.pi - 4*np.sqrt(abs(A1(q))/rd)
        # (4) ODE: incoming launch from d with sin a0 = bz/rho_d, zero of u on return to r = d
        a0 = np.pi - np.arcsin(bz/rd)
        ev = lambda s, y, q: y[0] - d; ev.direction = 1; ev.terminal = True
        sol = solve_ivp(rhs, [1e-9, 1e12], [d, 0, a0, 1e-9, 1.], args=(q,), events=[ev], rtol=1e-12, atol=1e-12, method='DOP853')
        uret = sol.y_events[0][0][3]; urange = np.max(np.abs(sol.y[3]))
        print(f"(3) d={d:.0f} rho_d={rd:.1f}: zeros of d_bPhi {['%.3f'%z for z in zs]}; b={bz:.3f} (predicted {pred_b:.3f}), Phi/pi={Phi(bz,d,q)/np.pi:.6f} (predicted {pred_P/np.pi:.6f})")
        print(f"(4)   u on return to r=d: {uret:.3e} (scale of u {urange:.3e}), phi/pi on return {sol.y_events[0][0][1]/np.pi:.6f}")

def F_of_r(r, q):
    h = 1e-6*r
    rr = (rho(r+h, q) - rho(r-h, q))/(2*h)
    f = 1 - 2/r; V = np.sqrt(1 - f/q)
    return 1/(f*V*rr)
def check_derivative(q=1.2):
    """(5) Phi_b = -2F(rho_d)/sqrt(rho_d^2-b^2) + 2|A1| sqrt(rho_d^2-b^2)/(b^2 rho_d) + E, |E| <= (pi/2) K' / b^3,
    K' = sup_{rho >= b} rho^3 |G'(rho)| (estimated on a grid)."""
    a1 = abs(A1(q))
    rs = np.geomspace(10, 1e8, 400); rh = rho(rs, q)
    G = np.array([F_of_r(r, q) - 1 + a1/x for r, x in zip(rs, rh)])
    Gp = np.gradient(G, rh)
    for d, b in [(20000., 549.16), (20000., 200.), (1e5, 1291.48)]:
        rd = rho(d, q); hb = 1e-3*np.sqrt(rd)
        fd = (Phi(b+hb, d, q) - Phi(b-hb, d, q))/(2*hb)
        s = np.sqrt(rd*rd - b*b)
        main = -2*F_of_r(d, q)/s + 2*a1*s/(b*b*rd)
        Kp = np.max(np.abs(Gp[rh >= b])*rh[rh >= b]**3)
        print(f"(5) d={d:.0f} b={b}: Phi_b finite differences {fd:.6e}, main part {main:.6e}, E*b^3 = {(fd-main)*b**3:.2f}, (pi/2)K' = {np.pi/2*Kp:.2f}")
if __name__ == '__main__':
    check_derivative()
