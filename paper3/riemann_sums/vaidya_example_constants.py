"""Constants of the Vaidya example (supplementary material, section "Constants of the Vaidya example").

Data: Ehat = sqrt(2), r0 = rf = 6, delta = pi/2, m0 = 1, m1 = 3/2. Every bound below is proved in outward-rounded
interval arithmetic (mpmath.iv); no floating-point value is used in a decision.

Method for the integrals without elementary form (the dips at mass 3/2). After scaling to M = 1 and r0 = 4,
    rho(r)^2 = r^2 (r - 2)/(r + 2),   R'(r) = sqrt(r/(r + 2)),   L = rho(t)  (Clairaut constant, periastron t),
    rho^2 - L^2 = (r - t) B_t(r),   B_t(r) = [(t + 2) r^2 + (t^2 - 4) r + 2t^2 - 4t] / ((r + 2)(t + 2)),
    angle  Delta(t) = 2 int_t^4 L R' / (rho sqrt(rho^2 - L^2)) dr,
    length ell(t)   = 2 int_t^4 R' rho / sqrt(rho^2 - L^2) dr,
    dive   2R(4)    = 2 int_2^4 R' dr.
The substitution r = t + (4 - t) sin^2(theta), theta in [0, pi/2], removes the inverse square root at the periastron:
dr / sqrt(r - t) = 2 sqrt(4 - t) cos(theta) dtheta. The interval [0, pi/2] is cut into n cells; on each cell the
integrand is enclosed by evaluating it in interval arithmetic on the whole cell, so
    integral in  sum_cells |cell| * f(cell)
with no monotonicity assumption. This is a different route from the monotone Riemann sums with integer square roots
quoted in the supplementary material, and serves as a second, independent verification of the same enclosures.

The elementary quantities of the late-accretion paragraph (F_a, Lambda, Phi, the piecewise-constant competitor) and
the rational inequalities of the early-window paragraph are evaluated in the same interval arithmetic.
Usage: python3 vaidya_example_constants.py [n_cells]   (default 6000; writes ../results/vaidya_example_constants.json)"""
import sys, json, os, time
from fractions import Fraction as Fr
import mpmath
from mpmath import iv, mpf

iv.prec = 120
N = int(sys.argv[1]) if len(sys.argv) > 1 else 6000
checks = []
def check(label, cond):
    checks.append({'check': label, 'pass': bool(cond)})
    print(('PASS ' if cond else 'FAIL ') + label, flush=True)
def q(fr):                      # exact rational -> thin interval
    fr = Fr(fr); return iv.mpf(fr.numerator)/fr.denominator
def lt(a, b): return a.b < b.a            # certainly a < b
def asinh(x): return iv.log(x + iv.sqrt(x*x + 1))
def acosh(x): return iv.log(x + iv.sqrt(x*x - 1))
def show(x):   # end points as 40-digit decimal strings, for reporting only (decisions use the binary end points)
    return [mpmath.libmp.to_str(x._mpi_[0], 40), mpmath.libmp.to_str(x._mpi_[1], 40)]

# ---------------------------------------------------------------- late accretion, mass 1 (M = 1), D = 91
def F(a, r): return iv.sqrt(r*(r + a)) - a*asinh(iv.sqrt(r/a))
def Phi(r): return 2*acosh(iv.sqrt(r/2))
two, three, six = q(2), q(3), q(6)
Lam = F(two, six) - F(two, q(Fr(11, 5))) + F(three, six) - F(three, q(Fr(11, 5)))
twoR6 = 2*(F(two, six) - F(two, two))
bands = [Fr(6), Fr(527, 100), Fr(19, 4), Fr(443, 100), Fr(433, 100)]
tans = [Fr(26, 25), Fr(38, 25), Fr(127, 50), Fr(289, 50)]
angle = iv.mpf(0); length = iv.mpf(0)
for i, tn in enumerate(tans):
    hi, lo = q(bands[i]), q(bands[i + 1])
    angle += 2*q(tn)*(Phi(hi) - Phi(lo))
    length += 2*iv.sqrt(1 + q(tn)**2)*(F(two, hi) - F(two, lo))
bound91 = q(91)/(11*iv.sqrt(two))
check('competitor sweeps at least pi/2: angle > pi/2', lt(iv.pi/2, angle))
check('competitor length < 5.8390', lt(length, q(Fr(58390, 10000))))
check('5.8390 < 91/(11 sqrt 2)', lt(q(Fr(58390, 10000)), bound91))
check('91/(11 sqrt 2) < Lambda(11/5)', lt(bound91, Lam))
check('Lambda(11/5) < 2R(6)', lt(Lam, twoR6))
check('T0 < 11 sqrt(2) * 5.8390 < 91', lt(11*iv.sqrt(two)*q(Fr(58390, 10000)), q(91)))

# ---------------------------------------------------------------- critical separation at mass 3/2 (scaled: M = 1, r0 = 4)
def rho2(r): return r*r*(r - 2)/(r + 2)
def Rp(r): return iv.sqrt(r/(r + 2))
def B(r, t): return ((t + 2)*r*r + (t*t - 4)*r + 2*t*t - 4*t)/((r + 2)*(t + 2))
def dip(t, n):
    """enclosures of Delta(t) and ell(t) by cellwise interval evaluation in theta."""
    L = iv.sqrt(rho2(t)); c = 2*iv.sqrt(4 - t); h = (iv.pi/2)/n
    A = iv.mpf(0); S = iv.mpf(0)
    for k in range(n):
        th = iv.mpf([(k*h).a, ((k + 1)*h).b])        # the cell, outward
        s = iv.sin(th); r = t + (4 - t)*s*s
        r = iv.mpf([max(r.a, t.a), min(r.b, mpf(4))])  # r lies in [t, 4] by construction
        Bt = B(r, t); rho = iv.sqrt(rho2(r)); cs = iv.cos(th)
        cs = iv.mpf([max(cs.a, 0), cs.b])
        common = Rp(r)*c*cs/iv.sqrt(Bt)
        A += L/rho*common
        S += rho*common
    width = iv.mpf([h.a, h.b])
    return 2*A*width, 2*S*width
ta, tb = q(Fr(29, 10)), q(Fr(297, 100))
check('B_t > 0 on [t, 4] at both periastra (interval evaluation on the whole range)',
      all(B(iv.mpf([t.a, mpf(4)]), t).a > 0 for t in (ta, tb)))
t0 = time.time()
Da, _ = dip(ta, N)
Db, ellb = dip(tb, N)
dive = iv.mpf(0); m = N
hd = q(2)/m
for k in range(m):           # 2R(4) = 2 int_2^4 sqrt(r/(r+2)) dr, cellwise enclosure
    r = iv.mpf([(2 + k*hd).a, (2 + (k + 1)*hd).b]); dive += Rp(r)*iv.mpf([hd.a, hd.b])
dive = 2*dive
dive_closed = 2*(F(two, q(4)) - F(two, two))
print(f'  Delta(29/10)  in {show(Da)}\n  Delta(297/100) in {show(Db)}\n  ell(297/100)  in {show(ellb)}\n'
      f'  2R(4)         in {show(dive)} (cells), {show(dive_closed)} (closed form); {time.time() - t0:.0f}s', flush=True)
delta_bar = q(Fr(39, 25))
check('Delta(29/10) > 39/25', lt(delta_bar, Da))
check('Delta(297/100) < 39/25', lt(Db, delta_bar))
check('ell(297/100) > 2R(4)', lt(dive_closed, ellb) and lt(dive, ellb))
check('39/25 < pi/2', lt(delta_bar, iv.pi/2))
check('cell enclosure of 2R(4) contains the closed form', dive.a <= dive_closed.a and dive_closed.b <= dive.b)
sup = {'Delta_ta': (Fr(1605672, 10**6), Fr(1605738, 10**6)), 'Delta_tb': (Fr(1558655, 10**6), Fr(1558714, 10**6)),
       'ell_tb': (Fr(3091180, 10**6), Fr(3091274, 10**6)), 'dive': (Fr(3081722, 10**6), Fr(3081750, 10**6))}
for name, X in (('Delta_ta', Da), ('Delta_tb', Db), ('ell_tb', ellb), ('dive', dive_closed)):
    lo, hi = sup[name]
    check(f'{name}: enclosure meets the interval printed in the supplement', not (X.b < q(lo).a or q(hi).b < X.a))

# ---------------------------------------------------------------- early window, h* = 1/1000 (mass 3/2, r0 = rf = 6)
H = Fr(1, 1000)
check('layer keeps r in [5, 7]: |rdot| < 4', 6 - 4*H > 5 and 6 + 4*H < 7)
check('angular speed sqrt(2)/r <= sqrt(2)/5 < 3/10', lt(iv.sqrt(two)/5, q(Fr(3, 10))))
check('monotone frozen geodesic between radii in [5,7] sweeps < pi/3: sqrt(125/392) > 1/2', Fr(125, 392) > Fr(1, 4))
check('rho(5)^2/rho(7)^2 = 125/392 at mass 3/2', Fr(5**2*(5 - 3), 5 + 3)/Fr(7**2*(7 - 3), 7 + 3) == Fr(125, 392))
check('1/(a sqrt(f(a))) < 1/3 on [5,7] (worst case a = 5, f = 1 - 3/a)', Fr(5)**2*(1 - Fr(3, 5)) > 9)
check('157/100 < pi/2', lt(q(Fr(157, 100)), iv.pi/2))
gap = Fr(157, 100) - Fr(3, 10)*H - (Fr(39, 25) + Fr(4, 3)*H)
check('gap at h = 1/1000 equals 251/30000 > 0', gap == Fr(251, 30000))
check('translated step: d <= 1/2000 and epsilon = 1/2000 end the accretion by 1/1000', Fr(1, 2000) + Fr(1, 2000) <= H)

out = {'method': 'outward-rounded interval arithmetic (mpmath.iv, 120 bits); cellwise range enclosure in theta, '
                  'r = t + (4 - t) sin^2(theta)', 'cells': N,
       'Delta_29_10': show(Da), 'Delta_297_100': show(Db), 'ell_297_100': show(ellb),
       'two_R4_cells': show(dive), 'two_R4_closed_form': show(dive_closed),
       'competitor_angle': show(angle), 'competitor_length': show(length), 'Lambda_11_5': show(Lam), 'two_R6': show(twoR6),
       'checks': checks, 'passed': sum(c['pass'] for c in checks), 'total': len(checks)}
res = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
os.makedirs(res, exist_ok=True)
json.dump(out, open(os.path.join(res, 'vaidya_example_constants.json'), 'w'), indent=1)
print(f"SUMMARY: {out['passed']}/{out['total']} PASS")
sys.exit(0 if out['passed'] == out['total'] else 1)
