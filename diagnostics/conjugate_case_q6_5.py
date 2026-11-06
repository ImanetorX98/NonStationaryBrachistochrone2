"""Floating-point diagnostic, not a certificate: the in-plane conjugate point quoted in the remark after the proposition
on sharpness at q = 3/2. Constant factor A = 1, Ehat^2 = q = 6/5, launch radius r0 = 10^4 M, impact parameter 10 b_c
(b_c = min rho, the critical impact parameter of the Fermat surface; Clairaut: rho sin a = b). The zero of the in-plane
Jacobi field is located with two integrator tolerances; the remark quotes Delta phi = 0.9451 pi, r = 530.09 M.
Usage: python3 conjugate_case_q6_5.py"""
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.integrate import solve_ivp
from conjugate_points_static import rhs, Fr

q, r0 = 1.2, 1.0e4
rho = lambda r: Fr(r, q)[1]
res = minimize_scalar(rho, bounds=(2.0001, 50.0), method='bounded', options={'xatol': 1e-12})
bc = res.fun
a0 = np.pi - np.arcsin(10*bc/rho(r0))            # incoming ray
print(f"b_c = {bc:.10f} (at r = {res.x:.8f}); launch angle pi - a0 = {np.pi - a0:.6e}")
for tol in (1e-10, 1e-12):
    ev_u = lambda s, y, qq: y[3]; ev_u.direction = -1; ev_u.terminal = True
    ev_pi = lambda s, y, qq: y[1] - np.pi; ev_pi.terminal = True
    sol = solve_ivp(rhs, [1e-9, 1e9], [r0, 0., a0, 1e-9, 1.], args=(q,), events=[ev_u, ev_pi],
                    rtol=tol, atol=tol*1e-2, method='DOP853')
    if len(sol.t_events[0]):
        r_c, ph_c = sol.y_events[0][0][0], sol.y_events[0][0][1]
        print(f"tol {tol:.0e}: in-plane zero of u at phi/pi = {ph_c/np.pi:.10f}, r = {r_c:.6f} M")
    else:
        print(f"tol {tol:.0e}: no zero of u before pi")
