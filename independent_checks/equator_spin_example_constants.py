"""Independent check of the rational constants of the example on the equator with spin (supplementary material,
"The example on the equator with spin": localization and uniqueness). Exact arithmetic and symbolic identities, no
interval certificate; it does not integrate the real flow.
Written through OpenAI Codex as a second, independent implementation (see the paragraph on AI tools in the paper);
it imports no other file of this repository. Usage: python3 equator_spin_example_constants.py
"""
from fractions import Fraction as F


def check(name, condition):
    if not condition:
        raise AssertionError(name)
    print("PASS", name)


hr, hp = F(216, 97), F(873, 16)
nmin, nmax, b0 = F(12, 11), F(3, 2), F(1, 873)
dmin = nmin - b0
check("B2 exact boundary sign", F(36, 1067)*F(149, 750)**2 > F(63, 2000)**2)
check("reference v squared", F(44, 2619) < F(13, 100)**2)
check("K lower", F(97, 324) > F(109, 200)**2)
check("K upper", F(1067, 2592) < F(13, 20)**2)
check("kappa K upper", F(44, 2619) < F(13, 100)**2)
check("sqrt(b) lower", F(32, 2619) > F(11, 100)**2)
check("sqrt(b) upper", nmin/(dmin**2*hp) < F(67, 500)**2)
check("abs(k) upper", 1/(4*dmin*hp) < F(43, 10000))
check("shape square-root bound", F(305, 289)**3 < F(87, 80)**2)
q = F(3, 25)*F(109, 200)*F(49, 50)/F(113, 100)*F(80, 87)
p = F(11, 100)*F(289, 305)-F(43, 10000)-F(67, 500)*F(1, 25)*F(13, 100)*F(46, 7)/F(17, 20)**3
e12 = F(21, 100)*F(13, 100)*F(1, 5)*F(89, 100)/F(17, 20)**3
e21 = F(13, 20)*F(1, 5)/F(17, 20)+F(13, 20)*F(1, 5)*F(46, 7)*F(13, 100)/F(17, 20)**2
check("D22 exact fraction", q == F(21364, 409625))
check("D11 exact fraction", p == F(1939954523, 20978510000))
check("D12 exact fraction", e12 == F(24297, 3070625))
check("D21 exact fraction", e21 == F(15509, 50575))
check("reference E1 stronger", q > F(13, 250))
check("reference E2 stronger", p > F(9, 100))
check("reference E3 stronger", e12 < F(1, 125))
check("reference E4 stronger", e21 < F(8, 25))
check("reference determinant margin", F(9, 100)*F(13, 250) > F(1, 125)*F(8, 25))
check("conditional transfer E1", q-F(1, 500) > F(1, 20))
check("conditional transfer E2", p-F(1, 500) > F(1, 12))
check("conditional transfer E3", e12+F(1, 500) < F(1, 100))
check("conditional transfer E4", e21+F(9, 100) < F(2, 5))
check("requested determinant margin", F(1, 240) > F(1, 250))
check("kappa upper", hr/hp < F(203, 1000)**2)
check("normalized p bound", F(203, 1000)*F(10, 49)+F(567, 16000) < F(77, 1000))
check("real a/sqrt(b) upper", 58/(4*F(184, 169)) < F(15, 4)**2)
speed = F(15, 4)*F(77, 1000)*F(64, 55)
check("real speed fraction", speed == F(84, 250))
check("all Q radial shell", F(21, 100)*speed < F(3, 40))
rlo = F(23, 4)
deltalo = rlo*rlo-2*rlo+F(1, 4)
hrlo = rlo*rlo/((1-2/rlo)*deltalo)
check("hr shell upper endpoint", hrlo == F(12167, 5235) and hrlo < F(12, 5))
check("hp shell lower", F(184621, 3600) >49)
check("N shell upper", F(50, 33) < F(25, 16))
check("radial shape bound", 29*(F(63, 2000)*F(72, 55))**2 < F(1, 20))
check("radial square-root bound", F(21, 20) < F(41, 40)**2)
m = F(7, 3)*F(64, 73)*F(40, 41)
M = F(15, 4)*F(64, 55)
check("weight fractions", m == F(17920, 8979) and M == F(48, 11))
check("weight ratio", M/m == F(26937, 12320) and M/m < F(9, 4))
c = F(3, 5)
check("weighted mean inequality", M*(1-c)**2-m*c*c <0)
check("terminal p improved", F(3, 20)*F(3, 5)*F(21, 100) == F(189, 10000))
heading = F(11, 100)
check("return heading 0.11", F(36, 1067)*(heading-heading**3/6)**2 > F(189, 10000)**2)
excursion = F(21, 200)*F(63, 2000)*F(36, 55)
check("return shell 0.02 squared", excursion**2*58 < F(1, 50)**2)

# Symbolic checks are identities, not parameter-domain bounds.
import sympy as sp

theta = sp.symbols("theta", real=True)
v, gamma, K, kappa, ell, rootb = sp.symbols("v gamma K kappa ell rootb", positive=True)
for sense in (-1, 1):
    X = sp.cos(theta)-sense*v
    L = sp.sqrt(X*X+gamma*sp.sin(theta)**2)
    radial = K*sp.sin(theta)/L
    angular = rootb*X/L
    check(f"ref D22 identity s={sense}", sp.simplify(sp.diff(radial, theta)-K*X*(1-sense*v*sp.cos(theta))/L**3) ==0)
    # rootb*gamma=kappa*K follows from the frozen coefficient definitions.
    check(f"ref D12 identity s={sense}", sp.simplify(sp.diff(angular, theta)+rootb*gamma*sp.sin(theta)*(1-sense*v*sp.cos(theta))/L**3) ==0)
    check(f"ref D11 integrand s={sense}", sp.simplify(sp.diff(angular, v)*(-v*ell)-sense*rootb*gamma*sp.sin(theta)**2*v*ell/L**3) ==0)
    check(f"ref D21 integrand s={sense}", sp.simplify(sp.diff(radial, v)*(-v*ell)+sense*K*sp.sin(theta)*X*v*ell/L**3) ==0)
A = sp.symbols("A", positive=True)
C = sp.sqrt(2)*A*A/sp.Integer(3)*sp.sqrt(hp.numerator/(sp.Integer(hp.denominator)*(1-A*A/3)))
logfactor = 2/A+A/(3-A*A)
check("terminal log derivative", sp.simplify(sp.diff(C,A)/C-logfactor) ==0)
check("terminal log factor at 1/2", logfactor.subs(A,sp.Rational(1,2)) ==sp.Rational(46,11))
# Monotonicity: derivative = -2/A²+(3+A²)/(3-A²)².
# On [1/2,1], the first summand <=-2, the second <=1.
check("terminal log factor derivative identity", sp.simplify(sp.diff(logfactor,A)-(-2/A**2+(3+A*A)/(3-A*A)**2)) ==0)

r = sp.symbols("r", positive=True)
pr, J = sp.symbols("pr J", real=True)
a, b, k, g, omega = (sp.Function(n)(r) for n in ("a","b","k","g","omega"))
xi = J-g
S = sp.sqrt(a*pr*pr+b*xi*xi)
H = -omega+S-k*xi
Hpp = sp.diff(H,pr,pr)
HpJ = sp.diff(H,pr,J)
HJJ = sp.diff(H,J,J)
check("real Hpp", sp.simplify(Hpp-a*b*xi*xi/S**3) ==0)
check("real HpJ", sp.simplify(HpJ+a*b*pr*xi/S**3) ==0)
check("real HJJ", sp.simplify(HJJ-a*b*pr*pr/S**3) ==0)
check("rank one cancellation", sp.simplify(Hpp*pr+HpJ*xi) ==0)
Hr, Hrp, HrJ, velphi = sp.diff(H,r), sp.diff(H,r,pr), sp.diff(H,r,J), sp.diff(H,J)
identity = pr*Hrp+J*HrJ-Hr-sp.diff(omega,r)+sp.diff(xi,r)*velphi-g*HrJ
check("real projective radial force identity", sp.simplify(identity) ==0)
check("real projective u source", sp.simplify(Hpp*pr/J+HpJ+Hpp*pr*g/(J*xi)) ==0)

print("\nExact analytic constants:")
for name, value in (("q_ref",q),("p_ref",p),("e12_ref",e12),("e21_ref",e21),("weight ratio",M/m)):
    print(name, value, float(value))
print("No real-model derivative budget established by this script.")
