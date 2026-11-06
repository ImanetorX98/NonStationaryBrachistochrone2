"""jets_t.py -- first version of the interval certification of the unimodality of Delta(t) (angle of the symmetric
dives of the proper-time branch; a = 0, A = 1, M = 1, q = Ehat^2, r0 = r_f). unimodality.py uses its jets J, jsqrt
and its zeroth-order quadrature (make, enclose) as a fast first stage; edges.py and maxwell_angle.py use J, jsqrt.

Delta(t) = int_0^1 G(u,t) du,  G = 4 L sqrt(a) c sqrt(N(r)) / (rho(r) sqrt(W/den)),
r = t + a u^2, a = r0 - t, L = rho(t), rho^2(x) = x^2 (x-2)/((q-1)x+2),
W = (q-1) r t (r+t-2) + 2(r^2 + r t + t^2 - 2r - 2t), den = ((q-1)r+2)((q-1)t+2), c = 1/sqrt(q).

Method (rigorous on [t_lo, t_hi] = [2 + eps0, r0 - eps1]):
  - interval arithmetic mpmath.iv (outward rounding);
  - Taylor jets in t up to order 2 (value, d/dt, d2/dt2) with interval coefficients: exact derivatives of the
    integrand, no finite differences;
  - integral enclosed by sum_k |U_k| * G([U_k], T): it contains the true value because G is continuous;
  - bisection of T until on every leaf Delta' > 0 (P), Delta' < 0 (N) or Delta'' < 0 (C).
The criterion used by the standalone run below (P* C+ N*) is superseded by the one in unimodality.py
([PC]*[CN]* with at least one C), which is the one used in the paper. The edges (2, t_lo) and (t_hi, r0), where the
integrand degenerates (L -> 0, a -> 0), are not covered here: see edges.py."""
import sys, time
from mpmath import iv, mpf
iv.dps = 20

class J:
    """second-order jet in t with interval coefficients"""
    __slots__ = ("v", "d", "dd")
    def __init__(s, v, d=0, dd=0):
        s.v = v if isinstance(v, type(iv.mpf(0))) else iv.mpf(v)
        s.d = d if isinstance(d, type(iv.mpf(0))) else iv.mpf(d)
        s.dd = dd if isinstance(dd, type(iv.mpf(0))) else iv.mpf(dd)
    @staticmethod
    def c(x): return x if isinstance(x, J) else J(x)
    def __add__(s, o): o = J.c(o); return J(s.v + o.v, s.d + o.d, s.dd + o.dd)
    __radd__ = __add__
    def __sub__(s, o): o = J.c(o); return J(s.v - o.v, s.d - o.d, s.dd - o.dd)
    def __rsub__(s, o): return J.c(o) - s
    def __neg__(s): return J(-s.v, -s.d, -s.dd)
    def __mul__(s, o):
        o = J.c(o); return J(s.v*o.v, s.d*o.v + s.v*o.d, s.dd*o.v + 2*s.d*o.d + s.v*o.dd)
    __rmul__ = __mul__
    def recip(s):
        iv0 = 1/s.v; d = -s.d*iv0*iv0
        dd = (2*s.d*s.d*iv0 - s.dd)*iv0*iv0
        return J(iv0, d, dd)
    def __truediv__(s, o): return s*J.c(o).recip()
    def __rtruediv__(s, o): return J.c(o)*s.recip()
def jsqrt(x):
    r = iv.sqrt(x.v); d = x.d/(2*r); dd = (x.dd - 2*d*d)/(2*r)
    return J(r, d, dd)

def make(q, r0):
    qi = iv.mpf(q); r0i = iv.mpf(r0); ci = 1/iv.sqrt(qi)
    one = J(1); two = J(2); Q = J(qi); Q1 = J(qi - 1); R0 = J(r0i); C = J(ci)
    def rho2(x): return x*x*(x - two)/(Q1*x + two)
    def G(U, T):
        t = J(T, 1, 0)
        a = R0 - t
        r = t + a*J(U*U)
        f = one - two/r
        N = one/(one - f/Q)
        L = jsqrt(rho2(t))
        rho = jsqrt(rho2(r))
        W = Q1*r*t*(r + t - two) + two*(r*r + r*t + t*t - two*r - two*t)
        den = (Q1*r + two)*(Q1*t + two)
        return J(4)*L*jsqrt(a)*C*jsqrt(N)/(rho*jsqrt(W/den))
    return G

def enclose(G, T, K):
    s1 = iv.mpf(0); s2 = iv.mpf(0)
    for k in range(K):
        U = iv.mpf([k, k + 1])/K
        g = G(U, T)
        s1 += g.d/K; s2 += g.dd/K
    return s1, s2

def certify(q, r0, eps=1e-4, Kmax=512, minw=1e-9):
    G = make(q, r0)
    tlo = float(2 + eps*(r0 - 2)); thi = float(r0 - eps*(r0 - 2))
    stack = [(tlo, thi)]; leaves = []; evals = 0
    while stack:
        a, b = stack.pop(); T = iv.mpf([a, b]); label = None   # a, b binary floats: exact endpoints
        K = 64
        while K <= Kmax:
            try:
                d1, d2 = enclose(G, T, K); evals += K
            except Exception:          # sqrt/division not certified on this (too wide) box
                break
            if d1 > 0: label = 'P'; break
            if d1 < 0: label = 'N'; break
            if d2 < 0: label = 'C'; break
            if K * 4 > Kmax: break
            K *= 4
        if label: leaves.append((a, b, label)); continue
        if b - a < minw: return None, [['?', 0]], evals, (a, b)
        m = 0.5*(a + b); stack.append((m, b)); stack.append((a, m))
    leaves.sort()
    seq = ''.join(l for _, _, l in leaves)
    import re
    ok = re.fullmatch(r'P*C+N*|P+|N+', seq) is not None
    # compress for report
    comp = []
    for _, _, l in leaves:
        if comp and comp[-1][0] == l: comp[-1][1] += 1
        else: comp.append([l, 1])
    return ok, comp, evals, (tlo, thi)

if __name__ == "__main__":
    cases = [(q, r0) for q in [1.2, 1.5, 2.0, 4.0, 10.0] for r0 in [4.0, 6.0, 10.0, 20.0]]
    if len(sys.argv) > 1: cases = [tuple(map(float, sys.argv[1].split(',')))]
    for q, r0 in cases:
        t0 = time.time()
        ok, comp, evals, rng = certify(q, r0)
        tag = "CERTIFIED unimodal" if ok else ("NOT certified" if ok is False else f"FAILED to resolve at width limit near t in [{rng[0]:.9f}, {rng[1]:.9f}]")
        print(f"q={q:5.2f} r0={r0:5.1f}: {tag} on t in [{float(rng[0]):.6f}, {float(rng[1]):.6f}]; "
              f"pattern {''.join(f'{l}x{n} ' for l, n in comp)}; {evals} integrand evals, {time.time()-t0:.0f}s")
        sys.stdout.flush()
