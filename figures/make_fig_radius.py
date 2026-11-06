"""Figure 2 of the main text: radius bound R(q_min) of the theorem on the swept angle with a radius bound (a = 0, exact
closed form), log scale; markers: numerical in-plane conjugate points before pi for constant A (floating-point
searches, ../diagnostics/), plotted at max(r0, r_conj). Usage: python3 make_fig_radius.py (writes fig2_radius.*)"""
import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 9, 'font.family': 'serif', 'mathtext.fontset': 'cm'})
q = np.linspace(1.0005, 1.4995, 4000)
a = 3*(3 - 2*q)*(q - 1); b = 33*q - 38
R = (b + np.sqrt(b**2 + 160*a))/(2*a)
top = 1e9
fig, ax = plt.subplots(figsize=(3.3, 3.0))
ax.fill_between(q, 2, np.minimum(R, top), color='0.78', lw=0, label='covered: $r_0,r_T<R(q_{\\min})$')
ax.plot(q, np.minimum(R, top), 'k', lw=1.3)
ax.fill_between([1.5, 1.6], 2, top, color='0.78', lw=0); ax.axvline(1.5, color='k', lw=1.3)
ax.axhline(8, color='0.3', lw=0.8, ls='--'); ax.text(1.01, 9.5, '$8M$', fontsize=8.5, color='0.2')
# numerical conjugate points before pi, constant A: (q, max(d, r_conj)) from the conjugate-point searches in ../diagnostics/
pts = [(1.05, 431.4), (1.1, 720.1), (1.2, 2338.0), (1.3, 5959.5), (1.4, 210836.), (1.45, 2164428.)]
ax.plot([p[0] for p in pts], [p[1] for p in pts], 'kv', ms=4.5, mfc='w', label='conjugate before $\\pi$ (numerical)')
ax.set_yscale('log'); ax.set_xlim(1.0, 1.6); ax.set_ylim(2, top)
ax.set_xlabel('$q_{\\min}$'); ax.set_ylabel('$\\max(r_0,r_T)/M$')
ax.legend(loc='upper left', fontsize=8, frameon=True)
fig.tight_layout(pad=0.3); import os; here = os.path.dirname(os.path.abspath(__file__))
fig.savefig(os.path.join(here, 'fig2_radius.pdf')); fig.savefig(os.path.join(here, 'fig2_radius.png'), dpi=300)
print('fig2_radius done')
