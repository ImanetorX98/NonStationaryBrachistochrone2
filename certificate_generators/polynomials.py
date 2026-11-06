"""Polynomials of the curvature certificates, in the physical variables (M = 1).

group1()  constant factor, reference b = infinity: x = 2/r, u = 1/q, z = a^2, w = 1 - u + u x, S = x^2 z - 4x + 4,
          U = x^3 nU/(32(1-x)w), V0 = x^3 nV/(64(1-x)w^2), J^2 = z x^8 S L^2/(1024(1-x)^2 w);
          P3 = nV^2 - 4 z x^2 S L^2 w^3 (V0^2 - J^2 = x^6 P3/(4096 (1-x)^2 w^4)), P6 = 4 nU^2 - z x^2 S L^2 w.
group3()  variable factor, same reference: y = sqrt(1 - x), t = tan(alpha/2), v = sqrt(w S); R_stat = (w/d^2) N(s),
          N = U c^2 + V0 s^2 + J s, D1 = c (al0 + al1 s)/d^2, D2 = -c^2 be/d^2.
group2()  representative G/N, constant factor, every |a| < M: R = r_+ (a^2 = R(2 - R)), x = R/r,
          U = x^3 D P_U/(R^3 wt H), V = -x^3 P_V/(R^4 wt^2 H), J^2 = 4 x^8 (2-R) D P_J^2/(R^7 wt^3 H^2),
          wt = R(1-u) + 2 u x > 0, D = (1-x)(R(1+x) - 2x) >= 0, H > 0.
group4()  weighted representative e^{2 sigma} G/N, sigma = -p log(1 - R/r), variable factor (built from curvature.py).
group5()  radial confinement: d_r G(X+-, X+-) = alpha +- beta sqrt(N Delta) for the null tangents X+- of G.
The closed forms of groups 1-3 are identified with the curvature in bridge.py."""
import sympy as sp

# ------------------------------------------------------------------------------------------------ group 1
x1, u1, z1 = sp.symbols('x u z', positive=True)


def group1():
    x, u, z = x1, u1, z1
    w = 1 - u + u*x
    S = x**2*z - 4*x + 4
    L = 8*u*x**2 - 14*u*x + 6*u + 5*x - 6
    nU = (2*u**2*x**5*z - 4*u**2*x**4*z + 2*u**2*x**3*z + 6*u*x**4*z - 12*u*x**3*z + 6*u*x**2*z + 4*u*x**2 - 8*u*x + 4*u
          + 5*x**3*z - 6*x**2*z - 6*x**2 + 14*x - 8)
    nV = (4*u**3*x**6*z - 12*u**3*x**5*z + 12*u**3*x**4*z - 4*u**3*x**3*z + 24*u**2*x**5*z - 66*u**2*x**4*z - 32*u**2*x**4
          + 60*u**2*x**3*z + 128*u**2*x**3 - 18*u**2*x**2*z - 192*u**2*x**2 + 128*u**2*x - 32*u**2 + 31*u*x**4*z - 61*u*x**3*z
          - 48*u*x**3 + 30*u*x**2*z + 144*u*x**2 - 144*u*x + 48*u + 10*x**3*z - 12*x**2*z - 12*x**2 + 28*x - 16)
    K = z*x**2*S*L**2
    return dict(x=x, u=u, z=z, w=w, S=S, L=L, nU=nU, nV=nV,
                P3=sp.expand(nV**2 - 4*K*w**3), P6=sp.expand(4*nU**2 - K*w),
                U=x**3*nU/(32*(1 - x)*w), V0=x**3*nV/(64*(1 - x)*w**2), J2=z*x**8*S*L**2/(1024*(1 - x)**2*w))


# ------------------------------------------------------------------------------------------------ group 3
def group3():
    """Rationalized quantities for the jets. Returns a dict with symbols y, t, u, a, v and
    mwN = w(-N), cal0 = c al0, cal1s = c s al1, c, be, w, S, and the reduction lin_in_v."""
    g1 = group1()
    y, t, u, a, v = sp.symbols('y t u a v', real=True)
    x = 1 - y**2; z = a**2
    sub = {g1['x']: x, g1['u']: u, g1['z']: z}
    nU = g1['nU'].subs(sub); nV = g1['nV'].subs(sub)
    w = 1 - u + u*x; S = x**2*z - 4*x + 4
    L = g1['L'].subs(sub)
    T = 3*z*u*x**4 - 4*z*u*x**3 + z*u*x**2 + 4*z*x**3 - z*x**2 - 16*u*x**3 + 44*u*x**2 - 40*u*x + 12*u - 20*x**2 + 32*x - 12
    s = 2*t/(1 + t**2); c = (1 - t**2)/(1 + t**2)
    wU = x**3*nU/(32*y**2); wV0 = x**3*nV/(64*y**2*w); wJ = -a*x**4*L*v/(32*y**2)
    mwN = -(wU*c**2 + wV0*s**2 + wJ*s)                     # w (-N) > 0 by the constant-factor result
    cal0 = c*3*z*x**5*y*w/(8*v)
    cal1s = -c*s*a*x**3*y*T/(8*S*w)
    be = z*x**4*y**2/(2*S)
    def lin_in_v(expr):
        """rational expression in v -> (R0, R1, den) with expr * den = R0 + R1 v, den > 0, after v^2 = w S."""
        num, den = sp.fraction(sp.together(expr))
        num = sp.expand(num)
        pv = sp.Poly(num, v)
        R0 = 0; R1 = 0
        for (k,), co in pv.terms():
            red = (w*S)**(k//2)
            if k % 2 == 0: R0 += co*red
            else: R1 += co*red
        return sp.expand(R0), sp.expand(R1), den
    return dict(y=y, t=t, u=u, a=a, v=v, mwN=mwN, cal0=cal0, cal1s=cal1s, c=c, be=be, w=w, S=S, T=T, lin_in_v=lin_in_v)


# ------------------------------------------------------------------------------------------------ group 2
def group2():
    x, R = sp.symbols('x R', positive=True)
    u = sp.Symbol('u', real=True)
    wt = R*(1 - u) + 2*u*x
    D = (1 - x)*(R*(1 + x) - 2*x)
    H = -(R**3*x**2 - R**3 + 2*R**2*x**3 - 2*R**2*x**2 + 4*R*u*x**4 - 4*R*x**3 - 8*u*x**4)
    P_U = (-3*R**3*u*x**2 - R**3*u + 3*R**3*x**2 - 16*R**2*u**2*x**3 + 21*R**2*u*x**3 + 6*R**2*u*x**2 - 6*R**2*x**2
           + 28*R*u**2*x**4 + 32*R*u**2*x**3 - 42*R*u*x**3 - 56*u**2*x**4)
    P_V = (9*R**6*u**2*x**4 - 14*R**6*u**2*x**2 + 5*R**6*u**2 - 15*R**6*u*x**4 + 23*R**6*u*x**2 - 8*R**6*u + 6*R**6*x**4 - 9*R**6*x**2 + 3*R**6
           + 28*R**5*u**3*x**5 - 28*R**5*u**3*x**3 - 81*R**5*u**2*x**5 - 36*R**5*u**2*x**4 + 111*R**5*u**2*x**3 + 28*R**5*u**2*x**2
           - 22*R**5*u**2*x + 48*R**5*u*x**5 + 60*R**5*u*x**4 - 81*R**5*u*x**3 - 46*R**5*u*x**2 + 19*R**5*u*x - 24*R**5*x**4 + 6*R**5*x**3
           + 18*R**5*x**2 - 108*R**4*u**3*x**6 - 112*R**4*u**3*x**5 + 148*R**4*u**3*x**4 + 56*R**4*u**3*x**3 + 126*R**4*u**2*x**6
           + 324*R**4*u**2*x**5 - 208*R**4*u**2*x**4 - 222*R**4*u**2*x**3 + 24*R**4*u**2*x**2 - 192*R**4*u*x**5 - 10*R**4*u*x**4
           + 162*R**4*u*x**3 + 24*R**4*x**4 - 12*R**4*x**3 - 16*R**3*u**4*x**7 + 124*R**3*u**3*x**7 + 432*R**3*u**3*x**6
           - 120*R**3*u**3*x**5 - 296*R**3*u**3*x**4 - 504*R**3*u**2*x**6 - 200*R**3*u**2*x**5 + 488*R**3*u**2*x**4 + 192*R**3*u*x**5
           - 100*R**3*u*x**4 + 32*R**2*u**4*x**8 + 64*R**2*u**4*x**7 - 496*R**2*u**3*x**7 - 336*R**2*u**3*x**6 + 464*R**2*u**3*x**5
           + 504*R**2*u**2*x**6 - 248*R**2*u**2*x**5 - 128*R*u**4*x**8 - 64*R*u**4*x**7 + 496*R*u**3*x**7 - 192*R*u**3*x**6 + 128*u**4*x**8)
    P_J = (3*R**4*u**2*x**2 - 3*R**4*u**2 - 6*R**4*u*x**2 + 6*R**4*u + 3*R**4*x**2 - 3*R**4 - 21*R**3*u**2*x**3 - 6*R**3*u**2*x**2
           + 19*R**3*u**2*x + 21*R**3*u*x**3 + 12*R**3*u*x**2 - 19*R**3*u*x - 6*R**3*x**2 - 20*R**2*u**3*x**4 + 48*R**2*u**2*x**4
           + 42*R**2*u**2*x**3 - 24*R**2*u**2*x**2 - 42*R**2*u*x**3 + 36*R*u**3*x**5 + 40*R*u**3*x**4 - 96*R*u**2*x**4 - 72*u**3*x**5)
    # (C1) P_U < 0, (C2) P_V > 0, (C3) P3 = P_V^2 - 4 R (2-R) x^2 D wt P_J^2 > 0 (V^2 > J^2),
    # (C4) P4 = -P_U (R wt D P_U + P_V) - x^2 (2-R) P_J^2 > 0 (value at the vertex negative: J^2 < -4 U (U - V))
    P3 = sp.expand(P_V**2 - 4*R*(2 - R)*x**2*D*wt*P_J**2)
    P4 = sp.expand(-P_U*(R*wt*D*P_U + P_V) - x**2*(2 - R)*P_J**2)
    P5 = sp.expand(-(P_V + R*wt*D*P_U)); P6 = sp.expand(R*wt*D*P_U**2 - 4*x**2*(2 - R)*P_J**2)
    return dict(x=x, R=R, u=u, wt=wt, D=D, H=H, P_U=P_U, P_V=P_V, P_J=P_J, P3=P3, P4=P4, P5=P5, P6=P6)


# ------------------------------------------------------------------------------------------------ group 4
def group4_tables(jets, hess, p=sp.Rational(1, 2)):
    """Coefficient tables in (x, R, u0) of -R_sigma, D1^2 and C c^2 as polynomials in (c, s, v), v = a sqrt(Delta w),
    for the weight sigma_p = -p log(1 - R/r). jets = curvature.ricci_h_jets(), hess = curvature.hessian_r().
      R_sig = R_st - sigma' Hess r(k,k) + (sigma'^2 - sigma'') c^2/hrr
      D1 = 2c [R01 sqrt(H')/(r^2 sqrt w) + s R12 sqrt(Delta)/sqrt(H')]   (coefficient of u'; unchanged by the weight)
      Cc^2 = -D2                                                          (coefficient of -u''; unchanged by the weight)
    with H' = r^4 + a^2 r^2 + 2a^2 r + 4a^2 u."""
    out, lap, hrr, hpp, (r, a, u0) = jets
    Hd, (rh, ah, uh) = hess
    Hd = {k: v.subs({rh: r, ah: a, uh: u0}, simultaneous=True) for k, v in Hd.items()}
    x, Rr = sp.symbols('x R', real=True)
    Dl = r**2 - 2*r + a**2; f = 1 - 2/r; w = 1 - f*u0; Hp = r**4 + a**2*r**2 + 2*a**2*r + 4*a**2*u0
    R = Rr
    sig1 = -p*R/(r*(r - R)); sig2 = p*R*(2*r - R)/(r**2*(r - R)**2)
    R00, R11, R22, R02 = out['R00'][(0, 0)], out['R11'][(0, 0)], out['R22'][(0, 0)], out['R02'][(0, 0)]
    Rsig = {(0, 0, 0): R00/lap - sig1*Hd['Hr00']/lap,
            (2, 0, 0): R11/hrr - sig1*Hd['Hr11']/hrr + (sig1**2 - sig2)/hrr,
            (0, 2, 0): R22/hpp - sig1*Hd['Hr22']/hpp,
            (0, 1, 1): 2*(R02 - sig1*Hd['Hr02'])/(a*Dl*w)}
    R01, R12 = out['R01'][(1, 0)], out['R12'][(1, 0)]
    D1sq = {(2, 0, 0): 4*R01**2*Hp/(r**4*w), (2, 1, 1): 8*R01*R12/(a*r**2*w), (2, 2, 0): 4*R12**2*Dl/Hp}
    Cc2 = {(2, 0, 0): -out['R00'][(0, 1)]/lap}
    def to_x(co):
        e = co.subs(r, Rr/x)
        e = sp.together(e)
        n, d = sp.fraction(e)
        n = sp.expand(n).subs(a**2, Rr*(2 - Rr)); d = sp.expand(d).subs(a**2, Rr*(2 - Rr))
        if n.has(a) or d.has(a):
            n = sp.expand(n.subs(a, sp.sqrt(Rr*(2 - Rr)))); d = sp.expand(d.subs(a, sp.sqrt(Rr*(2 - Rr))))
        return sp.factor(sp.cancel(n/d))
    tables = {'Rsig': {m: to_x(c) for m, c in Rsig.items()}, 'D1sq': {m: to_x(c) for m, c in D1sq.items()},
              'Cc2': {m: to_x(c) for m, c in Cc2.items()}}
    return tables, (x, Rr, u0)


def group4_assemble(tables, syms, K1, K2):
    """Inequalities S: -R_sig >= 0, I2: K2(-R_sig) - C c^2 >= 0, I1: K1^2 R_sig^2 - D1^2 >= 0 as P0 + P1 vh >= 0,
    vh^2 = ghat > 0 (vh = x a sqrt(Delta w)), with the sign of each common denominator taken at a sample point; the
    denominators have constant sign on the domain (products of factors of known sign, see the supplementary material).
    Returns (dict name -> (P0, P1), t, ghat)."""
    x, Rr, u0 = syms
    t = sp.Symbol('t', real=True); vh = sp.Symbol('vh', positive=True)
    ghat = sp.expand(Rr*(2 - Rr)*(1 - x)*(Rr*(1 + x) - 2*x)*(Rr*(1 - u0) + 2*u0*x))
    def poly_parts(dct, scale=1):
        dens = [sp.fraction(co)[1] for co in dct.values()]
        L = sp.lcm_list(dens)
        dcs = max(m[0] + m[1] for m in dct); dv = max(m[2] for m in dct)
        num = 0
        for (ic, is_, iv), co in dct.items():
            num += sp.expand(sp.cancel(scale*co*L)*(1 - t**2)**ic*(2*t)**is_*(1 + t**2)**(dcs - ic - is_)*vh**iv*x**(dv - iv))
        return sp.expand(num), sp.factor(L), dcs, dv
    def red(num):
        pv = sp.Poly(num, vh); P0 = 0; P1 = 0
        for (k,), co in pv.terms():
            if k % 2 == 0: P0 += co*ghat**(k//2)
            else: P1 += co*ghat**(k//2)
        return sp.expand(P0), sp.expand(P1)
    pt = {x: sp.Rational(1, 3), Rr: sp.Rational(17, 10), u0: sp.Rational(1, 5), t: sp.Rational(1, 7)}
    sign_at = lambda e: sp.sign(e.subs(pt))
    mR = {m: -co for m, co in tables['Rsig'].items()}
    numR, LR, dcsR, dvR = poly_parts(mR)
    numC, LC, dcsC, dvC = poly_parts(tables['Cc2'])
    numD, LD, dcsD, dvD = poly_parts(tables['D1sq'])
    res = {'S': red(sign_at(LR)*numR)}
    L2 = sp.lcm(LR, LC)
    num2 = sp.expand(K2*numR*sp.cancel(L2/LR)*(1 + t**2)**(dcsC - dcsR if dcsC > dcsR else 0)*x**(max(dvR, dvC) - dvR)
                     - numC*sp.cancel(L2/LC)*(1 + t**2)**(dcsR - dcsC if dcsR > dcsC else 0)*x**(max(dvR, dvC) - dvC))
    res['I2'] = red(sign_at(L2)*num2)
    L1 = sp.lcm(LR**2, LD)
    eR = 2*dcsR; eD = dcsD; ev = max(2*dvR, dvD)
    num1 = sp.expand(K1**2*numR**2*sp.cancel(L1/LR**2)*(1 + t**2)**(max(eR, eD) - eR)*x**(ev - 2*dvR)
                     - numD*sp.cancel(L1/LD)*(1 + t**2)**(max(eR, eD) - eD)*x**(ev - dvD))
    res['I1'] = red(sign_at(L1)*num1)
    return res, t, ghat


# ------------------------------------------------------------------------------------------------ group 5
def group5():
    """alpha and Q = alpha^2 - beta^2 N Delta as numerator/denominator in (y, R, u), with r = R + y (R - 1)/10,
    a^2 = R(2 - R). G_ee = -(r-2)/r, G_ep = -2a/r, G_pp = H'/(r W), W = r - u(r - 2) = r w > 0, N Delta = r Delta/W."""
    r, a2, u, y, R = sp.symbols('r a2 u y R', real=True)
    W = r - u*(r - 2); Hp = r**4 + a2*r**2 + 2*a2*r + 4*a2*u; Dl = r**2 - 2*r + a2
    Gee = -(r - 2)/r; Gpp = Hp/(r*W)
    dGee = sp.diff(Gee, r); dGpp = sp.diff(Gpp, r)
    #   alpha = dGee - 2 dGep Gep/Gpp + dGpp (Gep^2 + N Delta)/Gpp^2,  beta^2 = 4 (dGep - dGpp Gep/Gpp)^2/Gpp^2
    alpha = dGee - 2*(2/r**2)*(-2/r)*a2/Gpp + dGpp*(4*a2/r**2 + r*Dl/W)/Gpp**2
    beta2 = 4*a2*(2/r**2 + dGpp*(2/r)/Gpp)**2/Gpp**2
    Q = alpha**2 - beta2*(r*Dl/W)
    sub = {r: R + y*(R - 1)/10, a2: R*(2 - R)}
    def numden(e):
        e = sp.together(e.subs(sub))
        n, d = sp.fraction(e)
        return sp.expand(n), sp.expand(d)
    na, da = numden(alpha); nq, dq = numden(Q)
    pt = {R: sp.Rational(3, 2), y: sp.Rational(1, 2), u: sp.Rational(1, 5)}
    sa = int(sp.sign(da.subs(pt))); sq = int(sp.sign(dq.subs(pt)))
    return dict(y=y, R=R, u=u, na=na, da=da, nq=nq, dq=dq, sa=sa, sq=sq)
