"""Independent check of the constants of the Vaidya example (supplementary material, "Constants of the Vaidya
example"), in exact rational arithmetic.

Part 1, late accretion (D = 91): primitives F_a and Phi, the genus-two discriminant of the dip differential, and
exact elementary enclosures of the swept angle and length of the piecewise competitor, of Lambda(11/5), of the escape
price 2R(6) and of 91/(11 sqrt 2), with rational square-root bounds, the atanh series for the logarithm with an
explicit positive tail, and pi from Machin's formula.
Part 2, critical separation at mass 3/2 and early window: monotone left/right Riemann sums with 16384 intervals and
square roots enclosed between consecutive multiples of 10^-30 for Delta(29/10), Delta(297/100), ell(297/100) and
2R(4); the rational inequalities of the early window h* = 1/1000.

Written through OpenAI Codex as a second, independent implementation (see the paragraph on AI tools in the paper);
the computational parts are reproduced unchanged, without the parts that read private files. It imports no other
file of this repository; ../riemann_sums/vaidya_example_constants.py checks the same bounds by a different method.
Usage: python3 vaidya_example_constants_exact.py
"""
from fractions import Fraction as F
from functools import lru_cache
from math import isqrt
import json, os
import sympy as sp

checks = []
def check(name, c):
    c = bool(c); checks.append({'name': name, 'pass': c})
    print(('PASS ' if c else 'FAIL ') + name, flush=True)
    if not c: raise AssertionError(name)


def part1_late_accretion():
    r,a,t=sp.symbols('r a t',positive=True)
    Fa=sp.sqrt(r*(r+a))-a*sp.asinh(sp.sqrt(r/a))
    check('primitive F_a derivative for r,a>0',sp.simplify(sp.diff(Fa,r)-sp.sqrt(r/(r+a)))==0)
    # Phi derivative is equivalently derivative of 2 log(sqrt(r)+sqrt(r-2)).
    u=sp.symbols('u',positive=True)
    rr=2+u
    Ph=2*sp.log(sp.sqrt(rr)+sp.sqrt(u))
    check('Phi derivative for r>2',sp.simplify(sp.diff(Ph,u)-1/sp.sqrt(rr*u))==0)
    check('R-prime/rho at mass one',sp.factor((r/(r+2))/(r*r*(r-2)/(r+2))-1/(r*(r-2)))==0)
    B=(t+2)*r*r+(t*t-4)*r+2*t*t-4*t
    H=r*(r-2)*(r+2)*(r-t)*B
    disc=16777216*t**10*(t-2)**5*(t+2)**3*(t*t-8*t-4)*(t*t+2*t-4)**2
    check('sextic discriminant factorization',sp.factor(sp.discriminant(H,r)-disc)==0)
    check('quadratic discriminant factorization',sp.factor(sp.discriminant(B,r)-(t-2)*(t+2)*(t*t-8*t-4))==0)
    check('entire periastro interval lies in (2,3)',2<F(29,10)<F(297,100)<3)
    # On 2<t<3, t^2-8t-4 < 9-16-4 <0; t^2+2t-4 >4>0.
    check('discriminant signs bounded analytically on (2,3)',9-16-4<0 and 4+4-4>0)
    rho2=r*r*(r-2)/(r+2);L2=t*t*(t-2)/(t+2)
    check('dip angle differential squared uses the sextic',
          sp.factor((L2*r/(r+2)/(rho2*(rho2-L2)))*H - t*t*(t-2)*(r+2)**2)==0)

    # Exact elementary enclosures. Root scale 10^-40, log series N=40.
    SCALE=10**40; ROUND=10**30; N=40
    def root(v):
        v=F(v);k=isqrt(v.numerator*SCALE*SCALE//v.denominator)
        upper=k if k*k*v.denominator==v.numerator*SCALE*SCALE else k+1
        return F(k,SCALE),F(upper,SCALE)
    def add(x,y):return x[0]+y[0],x[1]+y[1]
    def neg(x):return -x[1],-x[0]
    def sub(x,y):return add(x,neg(y))
    def mul(x,y):
        vals=[a*b for a in x for b in y];return min(vals),max(vals)
    def scale(x,c):return mul(x,(F(c),F(c)))
    def div(x,y):
        assert y[0]>0
        return mul(x,(1/y[1],1/y[0]))
    def rounded(x):
        return F((x[0]*ROUND).__floor__(),ROUND),F((x[1]*ROUND).__ceil__(),ROUND)
    @lru_cache(None)
    def log_unit(z):
        assert 1<=z<=2
        y=(z-1)/(z+1)
        partial=2*sum((y**(2*k+1)/F(2*k+1) for k in range(N)),F(0))
        tail=2*y**(2*N+1)/(F(2*N+1)*(1-y*y))
        return partial,partial+tail
    @lru_cache(None)
    def log_exact(z):
        assert z>0
        if z<1:return neg(log_exact(1/z))
        k=0
        while z>=2:z/=2;k+=1
        return rounded(add(log_unit(z),scale(log_unit(F(2)),k)))
    def log_interval(z):return log_exact(z[0])[0],log_exact(z[1])[1]
    def phi_diff(hi,lo):
        num=add(root(hi),root(hi-2));den=add(root(lo),root(lo-2))
        return scale(log_interval(div(num,den)),2)
    def Fdiff(aa,hi,lo):
        first=sub(root(hi*(hi+aa)),root(lo*(lo+aa)))
        ratio=div(add(root(hi),root(hi+aa)),add(root(lo),root(lo+aa)))
        return sub(first,scale(log_interval(ratio),aa))
    b=[F(6),F(527,100),F(19,4),F(443,100),F(433,100)]
    ks=[F(26,25),F(38,25),F(127,50),F(289,50)]
    check('rational radii strictly decrease and stay above 2',all(x>y for x,y in zip(b,b[1:])) and b[-1]>2)
    angle=(F(0),F(0));length=(F(0),F(0))
    for hi,lo,k in zip(b,b[1:],ks):
        angle=add(angle,scale(phi_diff(hi,lo),2*k))
        length=add(length,scale(mul(root(1+k*k),Fdiff(F(2),hi,lo)),2))
    angle=rounded(angle);length=rounded(length)
    lam=rounded(add(Fdiff(F(2),F(6),F(11,5)),Fdiff(F(3),F(6),F(11,5))))
    dive=rounded(scale(Fdiff(F(2),F(6),F(2)),2))
    threshold=rounded(div((F(91),F(91)),scale(root(F(2)),11)))
    # Pi from Machin alternating arctangent series, exact rational tails.
    def atan_bounds(inv,n):
        a=sum((F((-1)**k,(2*k+1)*inv**(2*k+1)) for k in range(n)),F(0))
        b=a+F((-1)**n,(2*n+1)*inv**(2*n+1))
        return min(a,b),max(a,b)
    pi=rounded(sub(scale(atan_bounds(5,60),16),scale(atan_bounds(239,15),4)))
    check('Machin tangent and positive branch',
          (F(120,119)-F(1,239))/(1+F(120,119*239))==1 and 4*(F(1,5)-F(1,375))-F(1,239)>0)
    check('angle strictly exceeds pi/2',angle[0]>pi[1]/2)
    check('printed rational angle lower bound 1.57395 valid',angle[0]>F('1.57395'))
    check('length strictly below 5.8390',length[1]<F('5.8390'))
    check('5.8390 < 91/(11 sqrt2)',F('5.8390')<threshold[0])
    check('91/(11 sqrt2) < Lambda',threshold[1]<lam[0])
    check('Lambda < escape price',lam[1]<dive[0])
    check('Lambda > 5.94',lam[0]>F('5.94'))
    check('corner radius safely exterior',b[-1]>F(4))
    check('exact squared duration inequality',2*(11*F('5.8390'))**2<91**2)
    return {'angle': angle, 'length': length, 'Lambda': lam, 'escape_2R6': dive, 'threshold_91': threshold}


def part2_mass_three_halves():
    r,t,u,J,E,m=sp.symbols('r t u J E m',positive=True)
    rho2=r*r*(r-2)/(r+2)
    B=((t+2)*r*r+(t*t-4)*r+2*t*t-4*t)/((r+2)*(t+2))
    check('desingularization exact divided difference',
          sp.factor(rho2-rho2.subs(r,t)-(r-t)*B)==0)
    check('rho squared convex for r>2: exact second derivative',
          sp.factor(sp.diff(rho2,r,2)-(2-32/(r+2)**3))==0)
    check('divided difference increasing in r: positive shifted numerator',
          all(c>0 for c in sp.Poly(sp.fraction(sp.factor(sp.diff(B,r)))[0].subs({r:2+u,t:2+J}).expand(),u,J).coeffs()))

    # Analytic prefix estimates (no intervals here).
    p=sp.symbols('p',real=True)
    w=E*E-1+2*m/r
    S=sp.sqrt(E*E*p*p+J*J/(r*r))
    h=-w*p+sp.sqrt(w)*S
    hr=sp.diff(h,r)
    hr_formula=2*m*p/(r*r)-m*S/(r*r*sp.sqrt(w))-sp.sqrt(w)*J*J/(r**3*S)
    check('prefix radial derivative of support, no m prime',sp.simplify(hr-hr_formula)==0)
    check('angular speed outside small-costate strip identity',
          sp.simplify(sp.diff(h,J)-sp.sqrt(w)*J/(r*r*S))==0)

    SCALE=10**30
    def sqrt_lo(v):return isqrt(v.numerator*SCALE*SCALE//v.denominator)
    def sqrt_hi(v):
        lo=sqrt_lo(v)
        return lo if lo*lo*v.denominator==v.numerator*SCALE*SCALE else lo+1
    def sqsum(fun,n,increasing):
        left=sum(sqrt_lo(fun(F(i,n))) for i in range(n))
        right=sum(sqrt_hi(fun(F(i,n))) for i in range(1,n+1))
        if increasing:return F(left,n*SCALE),F(right,n*SCALE)
        # For decreasing functions, right is a lower sum, left an upper sum.
        low=sum(sqrt_lo(fun(F(i,n))) for i in range(1,n+1))
        high=sum(sqrt_hi(fun(F(i,n))) for i in range(n))
        return F(low,n*SCALE),F(high,n*SCALE)
    def BB(rr,tt):return ((tt+2)*rr*rr+(tt*tt-4)*rr+2*tt*tt-4*tt)/((rr+2)*(tt+2))
    def ang2(x,tt):
        rr=tt+(4-tt)*x*x
        LL=tt*tt*(tt-2)/(tt+2)
        return 16*(4-tt)*LL/(rr*(rr-2)*BB(rr,tt))
    def length2(x,tt):
        rr=tt+(4-tt)*x*x
        return 16*(4-tt)*rr**3*(rr-2)/((rr+2)**2*BB(rr,tt))
    n=16384
    ta,tb=F(29,10),F(297,100)
    for tt in (ta,tb):
        ll=16*(4-t)*r**3*(r-2)/((r+2)**2*B)
        numerator=sp.fraction(sp.factor(sp.diff(ll.subs(t,sp.Rational(tt.numerator,tt.denominator)),r)))[0]
        shifted=sp.Poly(sp.expand(numerator.subs(r,sp.Rational(tt.numerator,tt.denominator)+u)),u)
        check(f'length integrand increasing at t={tt}: positive derivative coefficients',all(c>0 for c in shifted.all_coeffs()))
    aa=sqsum(lambda x:ang2(x,ta),n,False)
    ab=sqsum(lambda x:ang2(x,tb),n,False)
    lb=sqsum(lambda x:length2(x,tb),n,True)
    # 2 R(4), with M=1. All costs scale by 3/2 for physical M=3/2.
    dive=sqsum(lambda x:4*(2+2*x)/(4+2*x),n,True)
    dive=(2*dive[0],2*dive[1])
    delta_bar=F(39,25)
    check('Delta(29/10) > 39/25',aa[0]>delta_bar)
    check('Delta(297/100) < 39/25',ab[1]<delta_bar)
    check('ell(297/100) > 2R(4)',lb[0]>dive[1])
    check('rational angle bound 39/25 < pi/2 using pi>3.14',delta_bar<F(157,100))

    # Explicit early layer at physical M=3/2, r0=rf=6, E^2=2.
    H=F(1,1000)
    eps=eta=F(1,2000)
    check('epsilon+eta1 = 1/1000',eps+eta==H)
    check('prefix radial shell [5,7] from |r prime|<4',6-4*H>5 and 6+4*H<7)
    check('prefix angular speed <3/10 from sqrt2<3/2 and r>=5',F(3,2)/5==F(3,10))
    check('endpoint action Lipschitz bound: 1/(r sqrt(f))<1/3 on [5,7]',F(5)**2*(1-F(3,5))>9)
    check('monotone frozen tail angle <pi/3: rho(5)/rho(7)>1/2',F(125,392)>F(1,4))
    check('residual angle exceeds pi/3 using pi>3 and pi<22/7',
          F(3,2)-F(3,10)*H>F(22,21))
    check('residual angle exceeds perturbed critical threshold',
          F(157,100)-F(3,10)*H > delta_bar+F(4,3)*H)
    check('strict early gap at h=1/1000 has rational margin',
          F(157,100)-delta_bar-(F(3,10)+F(4,3))*H == F(251,30000))
    return {'Delta_29_10': aa, 'Delta_297_100': ab, 'ell_297_100': lb, 'two_R4': dive}


if __name__ == '__main__':
    res1 = part1_late_accretion()
    res2 = part2_mass_three_halves()
    enc = lambda v: {'lower': str(v[0]), 'upper': str(v[1]), 'decimal_diagnostic': [float(v[0]), float(v[1])]}
    out = {k: enc(v) for k, v in {**res1, **res2}.items()}
    out['checks'] = checks
    res = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
    os.makedirs(res, exist_ok=True)
    json.dump(out, open(os.path.join(res, 'independent_vaidya_example_constants_exact.json'), 'w'), indent=1)
    for k, v in out.items():
        if k != 'checks': print(k, v['decimal_diagnostic'])
    print(f'SUMMARY: {sum(c["pass"] for c in checks)}/{len(checks)} PASS')
