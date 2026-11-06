"""Curvature on the equator of conformal Kerr, computed from the metric (SymPy, exact). M = 1.

Coordinates (eta, r, phi) on the equatorial plane; f = 1 - 2/r, Delta = r^2 - 2r + a^2, u = 1/q.
  G (extended auxiliary metric, constant factor):
      G = [[-f, 0, -2a/r], [0, N r^2/Delta, 0], [-2a/r, 0, (r^2 + a^2(1 + 2/r) + 4a^2 u/r^2)/(1 - f u)]],  N = 1/(1 - f u).
  h = G/N with a time-dependent factor u = u(eta), w = 1 - f u:
      h = [[-f w, 0, -2a w/r], [0, r^2/Delta, 0], [-2a w/r, 0, r^2 + a^2(1 + 2/r) + 4a^2 u/r^2]].
Functions:
  ricci_GN_constant()      Ric(k,k) = U c^2 + V s^2 + J s for G/N (constant factor), ZAMO frame       (group 2)
  ricci_h_jets()           Ric of h with the jets (u0, u1, u2) = (u, u', u''), split by monomials in (u1, u2),
                           contracted with E0 = d_eta - (h_ep/h_pp) d_phi, E1 = d_r, E2 = d_phi       (groups 3 and 4)
  hessian_r()              Hess_h r contracted with E0, E1, E2 (for the weighted representative)       (group 4)
  orthonormal_UVJ_b_inf()  U, V0, J of the reference b = infinity, from the orthonormal-frame formulas    (group 1)
Each takes a few seconds; ricci_h_jets about half a minute."""
import sympy as sp


def _christoffel_ricci(G, X, simplify):
    Gi = simplify(G.inv())
    Gam = [[[simplify(sum(Gi[i, l]*(sp.diff(G[l, j], X[k]) + sp.diff(G[l, k], X[j]) - sp.diff(G[j, k], X[l])) for l in range(3))/2)
             for k in range(3)] for j in range(3)] for i in range(3)]
    return Gam


def ricci_GN_constant():
    """Group 2: representative G/N, constant factor. Returns (U, V, J) in (r, a, u) with Ric(k,k) = U c^2 + V s^2 + J s
    for k = e0 + c e1 + s e2 in the orthonormal ZAMO frame."""
    eta, r, phi, a, u = sp.symbols('eta r phi a u', real=True)
    X = [eta, r, phi]
    f = 1 - 2/r; Dl = r**2 - 2*r + a**2; N = 1/(1 - f*u)
    G = sp.Matrix([[-f, 0, -2*a/r], [0, N*r**2/Dl, 0], [-2*a/r, 0, (r**2 + a**2*(1 + 2/r) + 4*a**2*u/r**2)/(1 - f*u)]])/N
    G = sp.simplify(G)
    Gam = _christoffel_ricci(G, X, sp.simplify)
    Ric = sp.zeros(3)
    for j in range(3):
        for k in range(3):
            Ric[j, k] = sp.simplify(sum(sp.diff(Gam[i][j][k], X[i]) - sp.diff(Gam[i][j][i], X[k])
                            + sum(Gam[i][i][l]*Gam[l][j][k] - Gam[i][k][l]*Gam[l][j][i] for l in range(3)) for i in range(3)))
    Gee, Gep, Gpp, Grr = G[0, 0], G[0, 2], G[2, 2], G[1, 1]
    lap = sp.simplify(-(Gee - Gep**2/Gpp))
    e0 = sp.Matrix([1, 0, -Gep/Gpp])/sp.sqrt(lap); e1 = sp.Matrix([0, 1, 0])/sp.sqrt(Grr); e2 = sp.Matrix([0, 0, 1])/sp.sqrt(Gpp)
    E = [e0, e1, e2]
    Rab = sp.Matrix(3, 3, lambda i, j: sp.simplify((E[i].T*Ric*E[j])[0]))
    # Ric(k,k) = R00 + 2 R01 c + 2 R02 s + R11 c^2 + 2 R12 c s + R22 s^2; R01 = R12 = 0 (stationary, axisymmetric)
    assert sp.simplify(Rab[0, 1]) == 0 and sp.simplify(Rab[1, 2]) == 0
    U = sp.simplify(Rab[0, 0] + Rab[1, 1]); V = sp.simplify(Rab[0, 0] + Rab[2, 2]); J = sp.simplify(2*Rab[0, 2])
    return (U, V, J), (r, a, u)


def _h_metric(r, a, uf):
    f = 1 - 2/r; Dl = r**2 - 2*r + a**2; w = 1 - f*uf
    return sp.Matrix([[-f*w, 0, -2*a*w/r], [0, r**2/Dl, 0], [-2*a*w/r, 0, r**2 + a**2*(1 + 2/r) + 4*a**2*uf/r**2]])


def ricci_h_jets():
    """Groups 3 and 4: Ricci of h = G/N with u = u(eta), as rational functions of (r, a) and the jets (u0, u1, u2).
    Returns (split, lap, hrr, hpp, symbols): split[name][(i1, i2)] is the coefficient of u1^i1 u2^i2 in
    R00 = Ric(E0,E0), R01, R02, R11, R12, R22, with E0 = d_eta - (h_ep/h_pp) d_phi (norm^2 = -lap), E1 = d_r, E2 = d_phi."""
    eta, r, phi, a = sp.symbols('eta r phi a', real=True)
    u0, u1, u2 = sp.symbols('u0 u1 u2', real=True)
    uf = sp.Function('u')(eta)
    X = [eta, r, phi]
    h = _h_metric(r, a, uf)
    hi = h.inv().applyfunc(sp.cancel)
    dh = [h.diff(v) for v in X]
    Gam = [[[sp.cancel(sum(hi[i, l]*(dh[k][l, j] + dh[j][l, k] - dh[l][j, k]) for l in range(3))/2) for k in range(3)] for j in range(3)] for i in range(3)]
    def jets(e):
        return e.subs(sp.Derivative(uf, (eta, 2)), u2).subs(sp.Derivative(uf, eta), u1).subs(uf, u0)
    Ric = sp.zeros(3)
    for j in range(3):
        for k in range(j, 3):
            e = sum(sp.diff(Gam[i][j][k], X[i]) - sp.diff(Gam[i][j][i], X[k]) + sum(Gam[i][i][l]*Gam[l][j][k] - Gam[i][k][l]*Gam[l][j][i] for l in range(3)) for i in range(3))
            Ric[j, k] = Ric[k, j] = sp.cancel(jets(e))
    r_, a_ = r, a
    f = 1 - 2/r_; Dl = r_**2 - 2*r_ + a_**2; w = 1 - f*u0
    hee = -f*w; hep = -2*a_*w/r_; hrr = r_**2/Dl; hpp = r_**2 + a_**2*(1 + 2/r_) + 4*a_**2*u0/r_**2
    lap = sp.cancel(-(hee - hep**2/hpp))
    E0 = sp.Matrix([1, 0, -hep/hpp]); E1 = sp.Matrix([0, 1, 0]); E2 = sp.Matrix([0, 0, 1])
    q = lambda A, B: sp.cancel((A.T*Ric*B)[0])
    comps = {'R00': q(E0, E0), 'R01': q(E0, E1), 'R02': q(E0, E2), 'R11': q(E1, E1), 'R12': q(E1, E2), 'R22': q(E2, E2)}
    def split(e):
        n, d = sp.fraction(sp.cancel(e)); P = sp.Poly(sp.expand(n), u1, u2)
        return {m: sp.factor(cf/d) for m, cf in P.terms()}
    out = {nm: split(e) for nm, e in comps.items()}
    return out, lap, hrr, hpp, (r_, a_, u0)


def hessian_r():
    """Group 4: Hess_h r(X, Y) = -Gamma^r_XY for the constant-factor h, contracted with E0, E2 (block) and E1."""
    r, a, u0 = sp.symbols('r a u0', real=True)
    f = 1 - 2/r; Dl = r**2 - 2*r + a**2; w = 1 - f*u0
    h = _h_metric(r, a, u0)
    hrr_inv = Dl/r**2
    Gr = sp.zeros(3)
    for i in (0, 2):
        for j in (0, 2):
            Gr[i, j] = -sp.Rational(1, 2)*hrr_inv*sp.diff(h[i, j], r)
    Gr[1, 1] = sp.Rational(1, 2)*hrr_inv*sp.diff(h[1, 1], r)
    hee, hep, hpp, hrr = h[0, 0], h[0, 2], h[2, 2], h[1, 1]
    lap = sp.cancel(-(hee - hep**2/hpp))
    E0 = sp.Matrix([1, 0, -hep/hpp]); E1 = sp.Matrix([0, 1, 0]); E2 = sp.Matrix([0, 0, 1])
    q = lambda A, B: sp.factor(sp.cancel(-(A.T*Gr*B)[0]))
    Hd = {'Hr00': q(E0, E0), 'Hr02': q(E0, E2), 'Hr22': q(E2, E2), 'Hr11': q(E1, E1), 'lap': lap, 'hpp': hpp, 'hrr': hrr}
    # check: Hess(e0 + e1) = Delta/(2 r^2) [h_pp'/h_pp - w'/w - 2/r]
    assert sp.simplify(Hd['Hr00']/lap + Hd['Hr11']/hrr - Dl/(2*r**2)*(sp.diff(hpp, r)/hpp - sp.diff(w, r)/w - 2/r)) == 0
    return Hd, (r, a, u0)


def orthonormal_UVJ_b_inf():
    """Group 1: reference b = infinity, constant factor. Ric(K,K) = U cos^2 + V0 sin^2 + J sin in the radial length l of
    h0, with A = W_l/W, B = rho_l/rho, C = W B_l/rho (B = wind), W^2 = 1 - f/q, rho = sqrt(Delta)/f, dl = sqrt(Hr) dr."""
    r, q, a = sp.symbols('r q a', positive=True)
    f = 1 - 2/r; Dl = r**2 - 2*r + a**2
    W = sp.sqrt(1 - f/q); rho = sp.sqrt(Dl)/f; Hr = r**2/(f*Dl); Bw = -2*a/(r*f)
    dl = lambda e: sp.diff(e, r)/sp.sqrt(Hr)
    A_ = dl(W)/W; B_ = dl(rho)/rho; C_ = W*dl(Bw)/rho
    U = A_*B_ - dl(B_) - B_**2 + C_**2
    V0 = dl(A_) + A_**2 - dl(B_) - B_**2 + C_**2
    J = -(dl(C_) + 2*A_*C_)
    return (U, V0, J), (r, q, a)
