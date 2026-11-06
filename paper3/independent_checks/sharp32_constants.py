#!/usr/bin/env python3
"""Independent check of the constants of the sharpness example below q = 3/2 (supplementary material, "Constants in
the sharpness example"; main text, proposition on sharpness at q = 3/2 and lemma on the radius bound), plus the
elementary reasons written by hand in the main text in place of decision procedures.

Rational intervals and exact symbolic identities only: no numerical optimisation; SymPy is used for identities.
Written through OpenAI Codex as a second, independent implementation (see the paragraph on AI tools in the paper);
it imports no other file of this repository. Usage: python3 sharp32_constants.py
"""
from fractions import Fraction as Q
from math import isqrt
import sympy as s

count = 0
def check(name, condition):
    global count
    assert bool(condition), name
    count += 1
    print('PASS', name)

def sqrt_interval(x, digits=30):
    x = Q(x)
    scale = 10**digits
    n = isqrt(x.numerator * scale**2 // x.denominator)
    lo, hi = Q(n, scale), Q(n+1, scale)
    assert lo**2 <= x < hi**2
    return lo, hi

def atan_interval(x, n=24):
    x = Q(x)
    partial = sum(((-1)**j * x**(2*j+1) / (2*j+1) for j in range(n)), Q(0))
    following = partial + (-1)**n * x**(2*n+1) / (2*n+1)
    return min(partial, following), max(partial, following)

# Machin: pi = 16 atan(1/5) - 4 atan(1/239); alternating-series bounds.
a, b = atan_interval(Q(1,5)), atan_interval(Q(1,239))
pi_lo, pi_hi = 16*a[0]-4*b[1], 16*a[1]-4*b[0]
c_lo, c_hi = (3*v for v in sqrt_interval(6))
em, ep = Q(9,40000), Q(1,100)
sm_lo, sm_hi = sqrt_interval(1-em)
K, D0, bm, bp = 465, 20000, 300, 2000
check('epsilon minus', Q(bm**2,D0**2)==em)
check('epsilon plus', Q(bp**2,D0**2)==ep)
check('D0 > K/c', D0*c_lo > K)
fmin_lo = 1-c_hi/D0-Q(K,D0**2)
check('F_min > 0', fmin_lo > 0)
low_lo = -2/sm_lo+2*c_lo*D0*sm_lo/bm**2-pi_hi*K*D0/bm**3
low_hi = -2/sm_hi+2*c_hi*D0*sm_hi/bm**2-pi_lo*K*D0/bm**3
up_hi = -2*fmin_lo+2*c_hi*D0/bp**2+pi_hi*K*D0/bp**3
up_lo = -2*(1-c_lo/D0-Q(K,D0**2))+2*c_lo*D0/bp**2+pi_lo*K*D0/bp**3
check('lower bound > 0.1832 (exact)', low_lo > Q(1832,10000))
check('upper bound < -1.922 (exact)', up_hi < -Q(1922,1000))
check('negative control: lower bound NOT > 0.19', low_hi < Q(19,100))
check('negative control: upper bound NOT < -1.923', up_lo > -Q(1923,1000))
check('b_plus <= D/2: beta_plus/D0 <= 1/2', Q(bp,D0) <= Q(1,2))
rt3_lo, _ = sqrt_interval(3)
check('Phi < pi: 300 > pi K/(2 sqrt(3)c)', bm*2*rt3_lo*c_lo > pi_hi*K)

r, t, e, z, q = s.symbols('r t e z q', positive=True)
def identity(name, expression):
    check(name, s.simplify(expression)==0)

rho2 = 6*r**4/((r-2)*(r+10))
identity('rho(10)^2 = 375', rho2.subs(r,10)-375)
check('rho_star < 20 < 300', Q(375)<20**2<bm**2)
rc = 2*s.sqrt(19)-6
identity('r_c: r_c^2+12r_c-40 = 0', rc**2+12*rc-40)
check('2 < r_c < 3', 2 < rc and rc < 3)
identity('derivative of rho^2 with the sign of the critical polynomial',
         s.diff(rho2,r)-12*r**3*(r**2+12*r-40)/((r-2)**2*(r+10)**2))
identity('rho(3)^2 = 486/13', rho2.subs(r,3)-s.Rational(486,13))
check('b_c <= rho(3) < 7 < 300', Q(486,13)<7**2<bm**2)
check('rho(8200)^2 > D0^2', Q(6*8200**4,8198*8210)>D0**2)

# On t>=1 and 0<e<1: all displayed derivative factors have known sign.
identity('derivative of -2/sqrt(1-e/t) positive',
         s.diff(-2/s.sqrt(1-e/t),t)-e/(t**2*(1-e/t)**s.Rational(3,2)))
identity('derivative of sqrt(1-e/t) positive',
         s.diff(s.sqrt(1-e/t),t)-e/(2*t**2*s.sqrt(1-e/t)))
identity('derivative of -1/sqrt(t) positive',
         s.diff(-1/s.sqrt(t),t)-1/(2*t**s.Rational(3,2)))
identity('derivative of +1/sqrt(t) negative',
         s.diff(1/s.sqrt(t),t)+1/(2*t**s.Rational(3,2)))
check('negative control: the positive upper term increases between t=1 and 4',
      sqrt_interval(1-ep/4)[0] > sqrt_interval(1-ep)[1])
identity('D R1 bound scales as t^(-1/2)',
         D0*t/(bm*s.sqrt(t))**3-D0/(bm**3*s.sqrt(t)))

# Optional lightweight checks of the added elementary reasons.
c0, c1, c2 = 3*(2*q-3)*(q-1), 33*q-38, 3*(q-1)
num, den = z*(c0+c1*z+40*z*z), c2+4*z
identity('numerator of H_q prime', s.diff(num,z)*den-num*s.diff(den,z)
         -(c0*c2+2*c1*c2*z+(4*c1+120*c2)*z*z+320*z**3))
check('c1 positive at q=3/2 and increasing', c1.subs(q,s.Rational(3,2))>0 and s.diff(c1,q)>0)
f, u = s.symbols('f u', real=True)
w = 1-f*u
identity('margin for -f/(2w) < 1/2', w+f-(1+f*(1-u)))
delta = r*r-2*r+s.Rational(1,4)
hr, hp = r**3/((r-2)*delta), delta/(1-2/r)**2
identity('log derivative of h_r', s.diff(hr,r)/hr-(3/r-1/(r-2)-2*(r-1)/delta))
identity('log derivative of h_phi', s.diff(hp,r)/hp-(2*(r-1)/delta-4/(r*(r-2))))
identity('margin Delta < (r-1)^2', (r-1)**2-delta-s.Rational(3,4))
identity('upper bound of the log derivative of h_r', 3/r-1/(r-2)-2/(r-1)-(6-4*r)/(r*(r-1)*(r-2)))
identity('lower bound of the log derivative of h_phi', 2/(r-1)-4/(r*(r-2))-2*(r*r-4*r+2)/(r*(r-1)*(r-2)))
check('shell beyond 2+sqrt(2)', s.Rational(23,4)>2+s.sqrt(2))
identity('zeros of the radius polynomial',
         86*(q-(36+s.sqrt(1210))/86)*(q-(36-s.sqrt(1210))/86)-(86*q*q-72*q+1))
print('CERTIFIED low in', float(low_lo), float(low_hi))
print('CERTIFIED upper in', float(up_lo), float(up_hi))
print('rho_star/M =', s.sqrt(375).evalf(15), '; b_c/M =', s.sqrt(rho2.subs(r,rc)).evalf(15))
print(f'PASS {count}; FAIL 0 (negative controls included)')
