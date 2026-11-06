"""jets_nested.py -- nested jets (t outer, u inner) and midpoint quadrature with a rigorous remainder, for the
interval certificates of the proper-time branch at equal radii (supplementary material, section "The certified
critical separations and the Maxwell example").

For F smooth on a cell [u0, u1] of width w and centre m (Taylor with Lagrange remainder):
    int_{u0}^{u1} F du  in  F(m) w + [inf F'', sup F''](cell) w^3/24,
because the linear term integrates to zero and int (u-m)^2/2 du = w^3/24. The total error is O(1/K^2) instead of
the O(1/K) of the sum |U_k| F([U_k]) used in jets_t.py.

G is built with jets of order 2 in t whose coefficients are jets of order 2 in u with interval coefficients
(mpmath.iv): g.d.v = d/dt G at the point u, g.d.dd = d2/du2 d/dt G, and so on. No finite differences: all
derivatives are exact and enclosed with outward rounding."""
from mpmath import iv
iv.dps = 20
IVT = type(iv.mpf(0))

class JG:
    """second-order jet with coefficients in any ring (iv intervals or other JG)"""
    __slots__ = ("v", "d", "dd")
    def __init__(s, v, d=0, dd=0): s.v, s.d, s.dd = v, d, dd
    def __add__(s, o):
        if not isinstance(o, JG): return JG(s.v + o, s.d, s.dd)
        return JG(s.v + o.v, s.d + o.d, s.dd + o.dd)
    __radd__ = __add__
    def __neg__(s): return JG(-s.v, -s.d, -s.dd)
    def __sub__(s, o): return s + (-o)
    def __rsub__(s, o): return (-s) + o
    def __mul__(s, o):
        if not isinstance(o, JG): return JG(s.v*o, s.d*o, s.dd*o)
        return JG(s.v*o.v, s.d*o.v + s.v*o.d, s.dd*o.v + 2*(s.d*o.d) + s.v*o.dd)
    __rmul__ = __mul__
    def __truediv__(s, o): return s*inv(o)
    def __rtruediv__(s, o): return inv(s)*o

def inv(x):
    if not isinstance(x, JG): return 1/x
    i0 = inv(x.v); d = -(x.d*i0*i0)
    dd = (2*(x.d*x.d)*i0 - x.dd)*i0*i0
    return JG(i0, d, dd)

def gsqrt(x):
    if not isinstance(x, JG): return iv.sqrt(x)
    r = gsqrt(x.v); h = inv(2*r); d = x.d*h; dd = (x.dd - 2*(d*d))*h
    return JG(r, d, dd)

def lift(c):
    """interval constant as a (t,u)-jet"""
    return JG(JG(iv.mpf(c), iv.mpf(0), iv.mpf(0)), JG(iv.mpf(0), iv.mpf(0), iv.mpf(0)),
              JG(iv.mpf(0), iv.mpf(0), iv.mpf(0)))

def tvar(T):
    z = JG(iv.mpf(0), iv.mpf(0), iv.mpf(0))
    return JG(JG(iv.mpf(T), iv.mpf(0), iv.mpf(0)), JG(iv.mpf(1), iv.mpf(0), iv.mpf(0)), z)

def uvar(U):
    z = JG(iv.mpf(0), iv.mpf(0), iv.mpf(0))
    return JG(JG(iv.mpf(U), iv.mpf(1), iv.mpf(0)), z, z)

def make(q, r0, which):
    """which = 'angle': integrand of Delta; 'length': integrand of ell. Same regular kernel as
    critical_separation_estimate.py: r = t + a u^2, rho^2(r)-rho^2(t) = (r-t) W/den, dr = 2 a u du."""
    qi = iv.mpf(q); r0i = iv.mpf(r0)
    Q = lift(qi); Q1 = lift(qi - 1); R0 = lift(r0i); C = lift(1/iv.sqrt(qi))
    def rho2(x): return x*x*(x - 2)/(Q1*x + 2)
    def G(U, T):
        t = tvar(T); u = uvar(U)
        a = R0 - t; r = t + a*(u*u)
        f = 1 - 2*inv(r); N = inv(1 - f/Q)
        rho = gsqrt(rho2(r))
        W = Q1*r*t*(r + t - 2) + 2*(r*r + r*t + t*t - 2*r - 2*t)
        den = (Q1*r + 2)*(Q1*t + 2)
        common = 4*gsqrt(a)*C*gsqrt(N)/gsqrt(W/den)
        if which == 'angle': return common*gsqrt(rho2(t))/rho
        return common*rho
    return G

def integrate(G, T, K):
    """enclosures of int_0^1 of (G, dG/dt, d2G/dt2) du, midpoint rule + Lagrange remainder"""
    out = [iv.mpf(0), iv.mpf(0), iv.mpf(0)]
    w = iv.mpf(1)/K; w3 = w*w*w/24
    for k in range(K):
        gm = G((iv.mpf(2*k + 1))/(2*K), T)
        gc = G(iv.mpf([k, k + 1])/K, T)
        for i, comp in enumerate(('v', 'd', 'dd')):
            out[i] += getattr(gm, comp).v*w + getattr(gc, comp).dd*w3
    return out

def R0_enclosure(q, r0, K):
    """2 R0 = 2 int_2^{r0} c sqrt(N) dr = 2 int_2^{r0} sqrt(r/((q-1) r + 2)) dr, same rule in r"""
    qi = iv.mpf(q); s = iv.mpf(0); r0i = iv.mpf(r0)
    w = (r0i - 2)/K; w3 = w*w*w/24
    def F(R):
        r = JG(iv.mpf(R), iv.mpf(1), iv.mpf(0))
        return gsqrt(r/(r*(qi - 1) + 2))
    for k in range(K):
        m = 2 + (r0i - 2)*(2*k + 1)/(2*K)
        cell = iv.mpf([(2 + (r0i - 2)*k/K).a, (2 + (r0i - 2)*(k + 1)/K).b])
        s += F(m).v*w + F(cell).dd*w3
    return 2*s

# --- parameters, directed conversions, provenance ---
from fractions import Fraction
import hashlib, sys, platform, math

def qint(qs):
    """q as an interval that contains both the decimal rational (e.g. 6/5 for '1.2') and its binary64 value:
    hull of iv.mpf(string) and iv.mpf(float(string))."""
    a = iv.mpf(str(qs)); b = iv.mpf(float(qs))
    return iv.mpf([min(a.a, b.a), max(a.b, b.b)])

def exact(x):
    """end point of an mpmath interval (or a point) as an exact Fraction"""
    from mpmath.libmp import to_rational
    if isinstance(x, IVT):
        lo, hi = x._mpi_
        assert lo == hi, "exact() needs a point (an interval end point)"
        p, q = to_rational(lo)
    else:
        p, q = to_rational(x._mpf_)
    return Fraction(p, q)

def down(x, d):
    """decimal string <= x (directed rounding downwards, d decimals)"""
    fr = exact(x); k = math.floor(fr * 10**d)
    return _dec(k, d)

def up(x, d):
    fr = exact(x); k = math.ceil(fr * 10**d)
    return _dec(k, d)

def _dec(k, d):
    s = '-' if k < 0 else ''; k = abs(k)
    return f"{s}{k // 10**d}.{k % 10**d:0{d}d}"

DEPS = ['jets_nested.py', 'jets_t.py', 'unimodality.py', 'edges.py',
        'critical_separation.py', 'critical_separation_estimate.py']

def manifest():
    import mpmath, os
    here = os.path.dirname(os.path.abspath(__file__))
    lines = [f"python {platform.python_version()}  mpmath {mpmath.__version__}  iv.dps {iv.dps}"]
    for fn in DEPS:
        p = os.path.join(here, fn)
        if os.path.exists(p):
            lines.append(f"sha256[:16] {hashlib.sha256(open(p, 'rb').read()).hexdigest()[:16]}  {fn}")
    return '\n'.join(lines)
