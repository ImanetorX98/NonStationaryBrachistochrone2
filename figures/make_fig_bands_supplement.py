"""Figure of the supplementary material (section "The whole family of static frames"): the bands in (q_min, q_max)
covered by the swept-angle theorem and by the convex-front argument. Usage: python3 make_fig_bands_supplement.py"""
import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 10, 'font.family': 'serif', 'mathtext.fontset': 'cm'})
plt.rcParams.update({'font.size': 9})
fig, ax = plt.subplots(figsize=(3.3, 3.0))
qmin = np.linspace(1.0, 3.0, 2001)
qmax_top = 8.0
ax.plot([1, 3], [1, 3], 'k', lw=0.8)                                   # diagonal q_min = q_max
m32 = qmin >= 1.5
ax.fill_between(qmin[m32], qmin[m32], qmax_top, color='0.78', lw=0, label='Thm angle, any factor')
ax.plot([1.5, 1.5], [1.5, qmax_top], 'k-', lw=1.3)
m53 = qmin >= 5/3
ax.fill_between(qmin[m53], qmin[m53], qmax_top, facecolor='none', edgecolor='0.35', hatch='///', lw=0,
                label='convex fronts')
mc = (qmin > 1.5) & (qmin < 5/3)
qc = np.minimum((9 - 5*qmin[mc])/(10 - 6*qmin[mc]), qmax_top)
ax.fill_between(qmin[mc], qmin[mc], qc, facecolor='none', edgecolor='0.35', hatch='///', lw=0)
ax.plot(qmin[mc], qc, color='0.35', lw=0.9)
ax.text(1.04, 6.5, 'radius\nbound in\nthe main\ntext', fontsize=7.5, color='0.25')
ax.set_xlim(1, 3); ax.set_ylim(1, qmax_top)
ax.set_xlabel('$q_{\\min}$'); ax.set_ylabel('$q_{\\max}$')
ax.legend(loc='lower right', fontsize=8.5, frameon=True)
fig.tight_layout(pad=0.3); import os; here = os.path.dirname(os.path.abspath(__file__))
fig.savefig(os.path.join(here, 'fig2_bands_supplement.pdf')); fig.savefig(os.path.join(here, 'fig2_bands_supplement.png'), dpi=300)
