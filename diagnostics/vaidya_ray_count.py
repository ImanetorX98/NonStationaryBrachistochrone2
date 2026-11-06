"""vaidya_ray_count.py -- floating-point diagnostic, not a certificate: the Morse count in the universal cover of the
plane for ingoing Vaidya (main text, theorem on the antipodal cut: "non-certified root searches with Ehat^2 = 2, r0 = 6
and r_T in {6, 10} found one same-sense ray for each sampled lifted angle in [0.02 pi, 1.9 pi] for constant mass,
m = 1 + tanh(v/5)/2 and a jump of m from 1 to 2 within Delta v = 0.6").
If the Morse relations for Fermat's principle hold in (D~, G_3) (equatorial plane, universal cover, contractible and
globally hyperbolic) and all null geodesics to the lifted receiver {r = r_T, phi~ = delta} have index 0 and are
non-degenerate (total angle delta < pi), then there is EXACTLY ONE. Test: for delta on a fine grid, the same-sense (not
mirrored) geodesics found by the scan are counted (vaidya_arrival_cut.all_arrivals, bracketing + brentq for each
crossing index). Prediction: 1 for each delta in (0, pi). The grid is extended beyond pi (lifted angle delta in
(pi, 2pi)) as an exploration: there the index is not guaranteed.
Control that must fail ("low" argument): Ehat^2 < 3/2 at large radii, where in-plane conjugate points can appear:
count 3 (two minima and a saddle) on an interval of delta."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'tools'))
from progress import Progress
import numpy as np, sys
import vaidya_arrival_cut as Dn
def count(name, mfun, dmfun, ddmfun, r0, rt, dmax=1.0, nd=60, n=1601):
    deltas = np.linspace(0.02, dmax, nd)*np.pi
    res = Dn.all_arrivals(mfun, dmfun, ddmfun, r0, rt, deltas, n=n)
    cnt = []
    for d in deltas:
        same = [x for x in res[d] if not x[3]]
        # deduplicate (same launch found from two brackets)
        uniq = []
        for x in sorted(same, key=lambda t: t[1]):
            if not uniq or abs(x[1] - uniq[-1][1]) > 1e-8: uniq.append(x)
        cnt.append((d/np.pi, len(uniq), [round(u[0], 4) for u in uniq]))
    bad = [c for c in cnt if c[1] != 1]
    print(f"[{name}] r0 = {r0}, r_T = {rt}, Ehat^2 = {Dn.Ev**2:.3f}: delta/pi in [0.02, {dmax}], {nd} values")
    print(f"    same-sense geodesics found: counts {sorted(set(c[1] for c in cnt))}; deltas with count != 1: "
          + ("none" if not bad else ", ".join(f"{b[0]:.3f} ({b[1]}: v = {b[2]})" for b in bad[:12]) + (" ..." if len(bad) > 12 else "")))
    return cnt
if __name__ == "__main__":
    one = lambda x: 1.0; zero = lambda x: 0.0
    tanh = (lambda x: 1 + 0.5*np.tanh(x/5), lambda x: 0.1/np.cosh(x/5)**2, lambda x: -0.04*np.tanh(x/5)/np.cosh(x/5)**2)
    jump = (lambda x: 1 + (1 + np.tanh((x - 4)/0.3))/2, lambda x: (1/0.6)/np.cosh(np.clip((x - 4)/0.3, -300, 300))**2,
            lambda x: -(1/0.6)*(2/0.3)*np.tanh((x - 4)/0.3)/np.cosh(np.clip((x - 4)/0.3, -300, 300))**2)
    which = sys.argv[1] if len(sys.argv) > 1 else 'main'
    if which == 'main':
        Dn.Ev = np.sqrt(2.0)
        pg = Progress(3, 'vaidya_ray_count: r_T')
        for rt in (6.0, 10.0, 4.0):
            count("frozen", one, zero, zero, 6.0, rt, dmax=1.9)
            count("tanh", *tanh, 6.0, rt, dmax=1.9)
            count("jump", *jump, 6.0, rt, dmax=1.9)
            pg.step()
        pg.close()
    else:
        pg = Progress(3, 'vaidya_ray_count low: (E2, r0)')
        for E2, r0 in ((1.05, 20.0), (1.05, 50.0), (1.2, 50.0)):
            Dn.Ev = np.sqrt(E2)
            count(f"frozen low q", one, zero, zero, r0, r0, dmax=1.95, nd=79)
            pg.step()
        pg.close()
