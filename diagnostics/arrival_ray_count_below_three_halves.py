"""arrival_ray_count_below_three_halves.py -- floating-point diagnostic, not a certificate: descent below q = 3/2
(a = 0), the "earlier searches at r0 = r_T = 20M" of the remark after the proposition on sharpness at q = 3/2.
Profiles with q_min = 1.2 (R(1.2) = 12.995 M): (i) constant q = 1.2 (A = 1, Ehat^2 = 1.2); (ii) increasing A with
Ehat^2 = 2, q from 2 to 1.2 (A from 1 to sqrt(2/1.2)); (iii) oscillating A, q in [1.2, 2]. For each delta: same-sense
rays found by the scan and zeros of D_J on the ray that realises the arrival (reuses arrival_ray_count.py). Cases:
r0 = r_T = 6 (inside R: the theorem predicts one ray and no zero); r0 = 6, r_T = 20 and r0 = r_T = 20 (outside R: the
theorem is silent)."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'tools'))
from progress import Progress
import numpy as np, sys
import arrival_ray_count as Z
Ahi = np.sqrt(2/1.2)
def myprof(kind):
    if kind == 'c12':
        E2 = 1.2; a, da, dda = (lambda e: 1.0, lambda e: 0.0, lambda e: 0.0)
    elif kind == 'up12':
        E2 = 2.0; m, h = (1 + Ahi)/2, (Ahi - 1)/2
        a, da, dda = (lambda e: m + h*np.tanh((e - 5)/0.5), lambda e: (h/0.5)/np.cosh(np.clip((e - 5)/0.5, -300, 300))**2,
                      lambda e: -2*(h/0.25)*np.tanh((e - 5)/0.5)/np.cosh(np.clip((e - 5)/0.5, -300, 300))**2)
    else:
        E2 = 2.0; m, h = (1 + Ahi)/2, (Ahi - 1)/2
        a, da, dda = (lambda e: m + h*np.sin(e/2), lambda e: (h/2)*np.cos(e/2), lambda e: -(h/4)*np.sin(e/2))
    qf = lambda e: E2/a(e)**2; qd = lambda e: -2*E2*da(e)/a(e)**3; qdd = lambda e: E2*(6*da(e)**2/a(e)**4 - 2*dda(e)/a(e)**3)
    return qf, qd, qdd
Z.prof = myprof
Pq = lambda q, r: 3*(2*q - 3)*(q - 1)*r**2 + (33*q - 38)*r + 40
if __name__ == "__main__":
    cases = [(6.0, 6.0), (6.0, 20.0), (20.0, 20.0)]
    pg = Progress(len(cases)*len(sys.argv[1:] or ['c12', 'up12', 'osc12']), 'below 3/2: (r0, r_T, profile)')
    for r0, rt in cases:
        for kind in (sys.argv[1:] or ['c12', 'up12', 'osc12']):
            Z.run(kind, r0, rt, nd=30)
            pg.step()
    pg.close()
