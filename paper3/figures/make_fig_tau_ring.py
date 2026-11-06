"""Figure 7 of the main text (numerical illustration, not certified): proper-time branch on the equator, constant factor.
a = M/2, Q = Ehat^2/A^2 = 2, r0 = rf = 6M (theorem on the Maxwell rule with a ring).
(a) equatorial plane: the ring r = 2M and the minimisers for three separations: two shallow dips (delta < delta_c) and
    the tangent--arc--tangent curve of the completion (delta > delta_c), which no rail attains.
(b) length (units c M, c = A^2/Ehat) against delta: shallow and deep dips, the escape line W(delta), delta_ring, delta_c.
Quadratures from ../diagnostics/equator_ring_frozen.py. Usage: python3 make_fig_tau_ring.py"""
import sys, os, importlib.util, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.optimize import brentq
here = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('ring', os.path.join(here, '..', 'diagnostics', 'equator_ring_frozen.py'))
ring = importlib.util.module_from_spec(spec); spec.loader.exec_module(ring)
plt.rcParams.update({'font.size': 10, 'font.family': 'serif', 'mathtext.fontset': 'cm'})

a, q, r0, M = 0.5, 2.0, 6.0, 1.0
rho, Rp = ring.fns(a, q)[3], ring.fns(a, q)[4]
res = ring.analyse(a, q, r0, n=500, verbose=False)
Ds, ells, Ls, rps = res['Ds'], res['ells'], res['Ls'], res['rps']
im = int(np.argmax(Ds)); dring, dmax, St = res['dring'], res['Phimax'], res['St']
W = lambda d: 2*St + a*d
# delta_c: zero of ell_shallow - W on the shallow branch
sh_D, sh_l, sh_rp = Ds[im:], ells[im:], rps[im:]
g = sh_l - W(sh_D); k = np.where(np.diff(np.sign(g)) != 0)[0][0]
dc = sh_D[k] + (sh_D[k+1] - sh_D[k])*(-g[k])/(g[k+1] - g[k])
print(f"delta_ring={dring:.4f}  Phi_max={dmax:.4f}  delta_c={dc:.4f}")

def phi_profile(L, r_from, r_to, n=300):
    """angle swept by the Clairaut-L geodesic from r_from down to r_to (r_to >= periastron)."""
    rr = r_to + (r_from - r_to)*np.linspace(0, 1, n)**2
    out = []
    for x in rr:
        T = np.sqrt(x - r_to)
        Tf = np.sqrt(r_from - r_to)
        f_ = lambda t: 2*t*L*Rp(r_to + t*t)/(rho(r_to + t*t)*np.sqrt(max(rho(r_to + t*t)**2 - L*L, 1e-300))) if t > 0 else \
            (2*L*Rp(r_to)/(rho(r_to)*np.sqrt((rho(r_to + 1e-9)**2 - L*L)/1e-9)) if rho(r_to) > L*(1 + 1e-12) or abs(rho(r_to) - L) < 1e-12 else 0.0)
        out.append(quad(f_, T, Tf, limit=200)[0])
    return rr, np.array(out)

def shallow_curve(delta):
    j = int(np.argmin(abs(sh_D - delta)))
    rp = brentq(lambda t: ring.dip(a, q, r0, t)[1] - delta, sh_rp[max(j-2, 0)], sh_rp[min(j+2, len(sh_rp)-1)])
    L = rho(rp)
    rr, ph = phi_profile(L, r0, rp)                     # ph: angle from x0 to radius rr (decreasing r)
    half = delta/2
    # incoming leg: angle grows from 0 at r0 to half at rp; outgoing mirrored
    rr, ph = rr[::-1], ph[::-1]                         # from x0 (r0, angle 0) down to the periastron
    x1, y1 = rr*np.cos(ph), rr*np.sin(ph)
    x2, y2 = rr*np.cos(delta - ph), rr*np.sin(delta - ph)   # outgoing leg, from the periastron back to r0
    return np.concatenate([x1, x2[::-1]]), np.concatenate([y1, y2[::-1]]), ring.dip(a, q, r0, rp)[2]

def wrap_curve(delta):
    rr, ph = phi_profile(a, r0, 2*M)                   # tangent from r0 to the ring: L = rho_h = a (c = 1)
    Pt = ph[0]
    rr, ph = rr[::-1], ph[::-1]                         # from x0 down to the ring
    arc = np.linspace(Pt, delta - Pt, 200)
    xs = np.concatenate([rr*np.cos(ph), 2*np.cos(arc), (rr*np.cos(delta - ph))[::-1]])
    ys = np.concatenate([rr*np.sin(ph), 2*np.sin(arc), (rr*np.sin(delta - ph))[::-1]])
    return xs, ys, W(delta)

fig, axs = plt.subplots(1, 2, figsize=(8.6, 4.0), gridspec_kw={'width_ratios': [1, 1.15]})
ax = axs[0]
th = np.linspace(0, 2*np.pi, 400)
ax.fill(2*np.cos(th), 2*np.sin(th), color='0.85', zorder=0)
ax.plot(2*np.cos(th), 2*np.sin(th), 'k', lw=0.8)
ax.plot(6*np.cos(th), 6*np.sin(th), color='0.7', lw=0.5, ls=':')
cases = [(0.9, 'C0', r'$\delta=0.9$'), (1.75, 'C2', r'$\delta=1.75$'), (2.6, 'C3', r'$\delta=2.6$ (not attained)')]
for d, col, lab in cases:
    if d < dc:
        xs, ys, l = shallow_curve(d); ax.plot(xs, ys, color=col, lw=1.4, label=lab)
    else:
        xs, ys, l = wrap_curve(d); ax.plot(xs, ys, color=col, lw=1.4, ls='--', label=lab)
    ax.plot(6*np.cos(d), 6*np.sin(d), 'o', color=col, ms=4)
    print(f"  delta={d}: length={l:.4f}")
ax.plot(6, 0, 'ko', ms=5)
ax.annotate('ring $r=2M$', xy=(2*np.cos(-0.6), 2*np.sin(-0.6)), xytext=(1.0, -4.8), fontsize=8,
            arrowprops=dict(arrowstyle='->', lw=0.6))
ax.text(6.25, -0.55, '$x_0$', fontsize=9)
ax.set_aspect('equal'); ax.set_xlim(-7, 7.3); ax.set_ylim(-5.6, 7.3)
ax.set_xlabel('$r\\cos\\varphi/M$'); ax.set_ylabel('$r\\sin\\varphi/M$')
ax.legend(fontsize=7.5, loc='upper left', frameon=False)
ax.set_title('(a) minimisers and completion escape curve', fontsize=10)

ax = axs[1]
order_sh = np.argsort(Ds[im:]); order_dp = np.argsort(Ds[:im+1])
shD = np.concatenate([[0.0], Ds[im:][order_sh]]); shL = np.concatenate([[0.0], ells[im:][order_sh]])
ax.plot(shD, shL, 'k', lw=1.2, label='shallow dips')
ax.plot(Ds[:im+1][order_dp], ells[:im+1][order_dp], 'k--', lw=1.0, label='deep dips')
dd = np.linspace(dring, np.pi, 100)
ax.plot(dd, W(dd), color='C3', lw=1.2, label=r'ring curves $W(\delta)$')
dx = np.linspace(0, dc, 200); env_sh = np.interp(dx, shD, shL)
ax.plot(dx, env_sh, color='C0', lw=3.2, alpha=0.35, zorder=0)
ax.plot(np.linspace(dc, np.pi, 50), W(np.linspace(dc, np.pi, 50)), color='C0', lw=3.2, alpha=0.35, zorder=0,
        label='infimum')
for x, lab, ha in [(dring, r'$\delta_{\rm ring}$', 'center'), (dc, r'$\delta_c$', 'right'), (dmax, r'$\Phi_{\max}$', 'left')]:
    ax.axvline(x, color='0.6', lw=0.6, ls=':')
    ax.text(x, 0.35, ' ' + lab + ' ', ha=ha, fontsize=8.5, color='0.25')
ax.set_xlim(0, np.pi); ax.set_ylim(0, 1.1*W(np.pi))
ax.set_xlabel('angular separation $\\delta$'); ax.set_ylabel('length $/(cM)$, $c=A^2/\\hat E$')
ax.legend(fontsize=7.5, loc='center right', frameon=False)
ax.set_title('(b) Maxwell rule with a ring', fontsize=10)
fig.suptitle('Proper-time branch on the equator, constant factor: $a=M/2$, $\\hat E^2/A^2=2$ (numerical)', fontsize=10)
fig.tight_layout()
fig.savefig(os.path.join(here, 'fig7_tau_ring.pdf')); fig.savefig(os.path.join(here, 'fig7_tau_ring.png'), dpi=200)
print('fig7 done')
