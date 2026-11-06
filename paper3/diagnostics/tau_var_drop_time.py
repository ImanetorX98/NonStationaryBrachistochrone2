"""tau_var_drop_time.py -- floating-point diagnostic, not a certificate: critical instant of the drop of the factor
(main text, remark on numerical observations: "the continued interior extremal crosses the escape price ... near
eta_d = 8.03"). A: 1 -> 0.5 (tanh, width 1, centre eta_d), q = 2 (Ehat = sqrt 2), r0 = r_f = 6, delta = pi/2.
Comparison:
 V(eta_d) = D_down(eta_d) + R_{0.5}(6)  (limit dive: radial descent on a rail with its own clock, waiting, ascent
            with A_inf)
 c_int(eta_d) = cost of the interior extremal along the family continued from eta_d = 8 (lower branch) and, backwards,
                down to eta_d = 2 (shooting, tau_var_shooting.py).
The crossing itself (bisection) is computed by tau_var_numbers.py."""
import numpy as np
from scipy.integrate import solve_ivp, quad
import tau_var_shooting as D
E = np.sqrt(2.0); f = lambda r: 1 - 2/r; N = lambda a, r: 1/(1 - a**2*f(r)/E**2)
up = 0.25/E*quad(lambda r: np.sqrt(N(0.5, r)), 2, 6)[0]
def V(ed):
    """dive with radial descent, with the SAME normalisation of A as the extremal (A(0) = 1, Ehat unchanged; dividing
    by A(0)^2 instead would not be homogeneous, because N changes)"""
    Al, _ = mk(ed); ainf = Al(1e9)
    ev = lambda e, y: y[0] - (2 + 1e-7); ev.terminal = True
    sol = solve_ivp(lambda e, y: [-f(y[0])/np.sqrt(N(Al(e), y[0])), Al(e)**2*f(y[0])/E], (0, 1e5), [6.0, 0.0],
                    events=ev, rtol=1e-11, atol=1e-13, max_step=0.2)
    return sol.y[1, -1] + ainf**2/E*quad(lambda r: np.sqrt(N(ainf, r)), 2, 6)[0]
def mk(ed):
    n = 0.75 + 0.25*np.tanh(ed)
    return (lambda e: (0.75 + 0.25*np.tanh(-(e - ed)))/n), (lambda e: (-0.25/np.cosh(e - ed)**2)/n)
qv, r0 = 2.0, 6.0; Ev = np.sqrt(qv); delta = 0.5*np.pi
def cont(x, eds):
    out = []
    for ed in eds:
        Al, dAl = mk(ed)
        xn, info, ier, msg = D.fsolve(D.residual, x, args=(Al, dAl, Ev, r0, r0, delta), full_output=True, xtol=1e-12)
        if ier != 1 or np.abs(info['fvec']).max() > 1e-9: out.append((ed, None)); break
        x = xn; sol, Lv = D.flow(x, Al, dAl, Ev, r0, n_out=True)
        ee = np.linspace(0, x[2], 3001)
        out.append((ed, (sol.y[4, -1], sol.sol(ee)[0].min(), x[2])))
    return out
if __name__ == "__main__":
    x8 = np.array([3.03709298, 0.17370671, 17.74132463])
    rows = cont(x8, np.arange(8.0, 1.99, -0.5))[::-1] + cont(x8, list(np.arange(8.25, 10.30, 0.25)) + [10.34])[0:]
    print(" eta_d   c_int      V        c_int - V   min r    T")
    for ed, c in sorted(rows, key=lambda r: r[0]):
        v = V(ed)
        if c is None: print(f"{ed:6.2f}   continuation lost"); continue
        print(f"{ed:6.2f}  {c[0]:.6f}  {v:.6f}  {c[0]-v:+.6f}  {c[1]:.4f}  {c[2]:7.3f}", flush=True)
