"""tau_var_fold.py -- floating-point diagnostic, not a certificate: arclength continuation of the family of extremals
in the parameter eta_d (instant of the drop of A: 1 -> 0.5, tanh of width 1). q = 2, r0 = r_f = 6, delta = pi/2.
Unknowns z = (psi0, q0, T, eta_d); equations: r(T) = 6, phi(T) = pi/2, q(T) = 0, (z - z_pred).t = 0.
Distinguishes a fold (eta_d turns back) from an escape (T -> infinity, r_min -> 2M). The position of the fold
(eta_d = 10.35, main text, remark on numerical observations) is computed by tau_var_numbers.py."""
import numpy as np
import tau_var_shooting as D
qv, r0 = 2.0, 6.0; Ev = np.sqrt(qv); delta = 0.5*np.pi
def mk(ed):
    n = 0.75 + 0.25*np.tanh(ed)
    return (lambda e: (0.75 + 0.25*np.tanh(-(e - ed)))/n), (lambda e: (-0.25/np.cosh(e - ed)**2)/n)
SC = np.array([1.0, 1.0, 0.05, 0.2])          # scaling of (psi0, q0, T, eta_d) in the arclength metric
def F3(z):
    Al, dAl = mk(z[3]); return D.residual(z[:3], Al, dAl, Ev, r0, r0, delta)
def info(z):
    Al, dAl = mk(z[3]); sol, Lv = D.flow(z[:3], Al, dAl, Ev, r0, n_out=True)
    ee = np.linspace(0, z[2], 4001); Y = sol.sol(ee)
    return sol.y[4, -1], Y[0].min(), np.min(D.F['om'](0, Y[0], 0, 0, 0, Ev, Al(ee), 0) - Y[3])
if __name__ == "__main__":
    za = np.array([3.103601, 0.08014, 24.406, 10.25]); zb = np.array([3.107501, 0.07576, 25.409, 10.34375])
    za[:3] = D.fsolve(lambda x: F3(np.r_[x, za[3]]), za[:3], xtol=1e-12)
    zb[:3] = D.fsolve(lambda x: F3(np.r_[x, zb[3]]), zb[:3], xtol=1e-12)
    h = 0.05
    for it in range(400):
        t = (zb - za)*SC; t /= np.linalg.norm(t)
        zp = zb + h*t/SC
        g = lambda z: list(F3(z)) + [np.dot((z - zp)*SC, t)]
        z, inf, ier, msg = D.fsolve(g, zp, full_output=True, xtol=1e-12)
        res = np.abs(inf['fvec']).max()
        if ier != 1 or res > 1e-8 or z[2] <= 0:
            h /= 2
            if h < 1e-5: print("step too small; stop"); break
            continue
        za, zb = zb, z
        if it % 3 == 0:
            c, rmin, wq = info(z)
            print(f"eta_d={z[3]:8.4f} psi0={z[0]:.6f} q0={z[1]:+.5f} T={z[2]:8.3f} cost {c:.6f} min r {rmin:.5f} "
                  f"min(omega-q) {wq:+.5f}", flush=True)
        h = min(h*1.3, 0.3)
        if z[3] > 16 or z[3] < 6 or z[2] > 300: print("left the window; stop"); break
