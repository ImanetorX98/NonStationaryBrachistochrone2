# Portable positivity certificates (curvature on the equator of conformal Kerr)

These 33 files are the certificates of the supplementary material, section *Curvature certificates on the equator*
(table *The curvature certificates on the equator*). They prove the polynomial inequalities behind the propositions
on the absence of conjugate points on the equator of conformal Kerr (`K_inf`, constant factor; `K_var`, variable
factor) and the radial confinement used in case (2) of the latter.

Each file is a closed algebraic lemma: a sign condition for a polynomial with rational coefficients on a box
`[0,1]^d`, together with its proof. The verifier `verify_certificate.py` checks the proofs in exact rational
arithmetic. It is a single self-contained file (it shows the progress bar of `../tools/` only if that file is present)
and does not use the code that generated the certificates
(`../certificate_generators/`); it uses different algorithms (see below).

## Quick start

```sh
python3 verify_certificate.py *.json     # about 20 minutes; most of it on group4_B4_I1.json
python3 tamper_tests.py                  # 15 deliberate alterations, all must be rejected (a few seconds)
shasum -a 256 -c SHA256SUMS              # integrity of the files
```

Only the Python standard library is needed (`fractions`, `math`, `itertools`, `json`). Expected output:
`../results/certificates_verify_all.txt` and `../results/certificates_tamper_tests.txt`; the last line of the first is
`OUTCOME: all certificates verified`, and the exit code is 0.

## The certificates

The groups follow the supplementary material. Box variables are affine images of the physical ones unless stated:
`u = 1/q`, `t = tan(alpha/2)` for the direction of the null vector, and `R = r_+/M` with `a^2 = R(2 - R)`. The exact
map of each file is in its field `map`.

| Group | Files | Statement | Type | Leaves |
|---|---|---|---|---|
| (1) static part, constant factor, `q >= 2`, `|a| <= 33M/50` | `group1_minus_nU`, `group1_minus_nV`, `group1_P3`, `group1_P6` | `U <= 0`, `V <= 0`, `V^2 >= J^2`, `J^2 <= 4U^2` for `Ric = U(1 - s^2) + V s^2 + J s` | nonneg | 1 each |
| (2) `K_inf`, constant factor, `q >= 5/3`, every `|a| <= M` | `group2_Kinf` | one of three sufficient sets on each leaf | alternative | 19 |
| (3) `K_var` (1), jets, `q >= 5/2`, `|a| <= 33M/50` | `group3_{I1plus,I1minus,I2}_{main,C1,C2A,C2B}` | `|D_1| <= (5/4)(-R_stat)` with both signs (I1plus, I1minus), `|D_2| <= (1/3)(-R_stat)` (I2) | radical (main), radical_square (corners) | 2007, 2024, 2229 in the main boxes |
| (4) `K_var` (2), spin bands B1–B4, weight `sigma = -(1/2) log(1 - R/r)` | `group4_B{1,2,3,4}_{S,I2,I1}` | `-R_sigma >= 0` (S), `K_2(-R_sigma) - C c^2 >= 0` (I2), `K_1^2 R_sigma^2 - D_1^2 >= 0` (I1) | radical | B4, I1: 1414 |
| (5) radial confinement in case (2) | `group5_den_alpha`, `group5_den_Q`, `group5_minus_alpha`, `group5_Q` | signs of the two denominators, `-alpha > 0`, `alpha^2 - beta^2 N Delta > 0` | positive | 1 each |

Spin bands of group (4), with `R = r_+/M`:

| Band | R | K_1 | K_2 |
|---|---|---|---|
| B1 | [8/5, 44/25] | 3/2 | 7/2 |
| B2 | [143/100, 8/5] | 3 | 9 |
| B3 | [131/100, 143/100] | 5 | 19 |
| B4 | [57/50, 131/100] | 14 | 110 |

In group (3) the main certificate of each inequality delegates two corners where both `y` and either `1 - T` or `a`
are small. The corner C1 (`y <= 1/8`, `T >= 15/16`) is the affine image of a separate box. The corner C2
(`y <= 1/8`, `A <= 1/8`) is covered by two sectors: C2A (`a <= 2y`, with `a = 2yA'`) and C2B (`2y <= a`, with
`y = aY'/2`), in which the common power of the small parameter is divided out of `P0` and `P1`.

## File format

- `name`, `type`, `description`, `variables`, `map` (documentation of the link between box and physical variables).
- `polynomials`: lists of monomials `[[exponents], "p/q"]`, rational numbers in exact textual form.
- `tree`: the proof, read in preorder.
  - `["s", k]` halves the current cell along axis `k`.
  - `["l"]` is a leaf of a `positive` or `nonneg` lemma: all Bernstein coefficients of the polynomial on the cell are
    `> 0` (or `>= 0`).
  - `["l", s0, s1, sigma]` is a leaf of a `radical` lemma (below).
  - `["l", i]` is a leaf of an `alternative` or `radical_square` lemma: sufficient set number `i` holds on the cell.
  - `["d", region]` is a delegated cell: it is covered by other certificates, and the verifier checks that the cell
    lies in the declared region.
- `delegations`: `[region, box, [delegated certificates]]`.
- `composition` (in delegated certificates): `main`, `region`, the substitution `phi`, the positive monomials `m` and
  `lam`, and `coverage` (`affine` or `external`).
- `sufficient_sets`, `polynomial_order`, `logic` (for `alternative`): the sets and the reason why each implies the
  claim.
- `statistics`: nodes, leaves, delegations and degrees, for information only.

**Types of lemma.**
- `positive`: `P > 0` on the box. `nonneg`: `P >= 0`.
- `radical`: `P0 + P1 sqrt(g) >= 0`, with `g >= 0`, through chord and tangent envelopes. On a leaf let
  `g0 = max(0, min of the Bernstein coefficients of g)` and `g1` their maximum. The rationals `s0, s1 >= 0` satisfy
  `s0^2 <= g0`, `s1^2 <= g1`, so the chord through `(g0, s0)` and `(g1, s1)` lies below `sqrt(g)` by concavity; for
  `sigma > 0`, `sqrt(g) <= sigma/2 + g/(2 sigma)`. Since `F = P0 + P1 v` is affine in `v`, it suffices that the
  coefficients of `P0 + A P1 + B P1 g` are `>= 0` for the two pairs `(A, B)` given by the two lines.
- `radical_square`: `P0 + P1 sqrt(g) >= 0` through `Q = P1^2 g - P0^2`, computed by the verifier itself. On each leaf
  one of `{P0 >= 0, P1 >= 0}`, `{P1 >= 0, Q >= 0}` (then `P1 sqrt(g) >= |P0|`), `{P0 >= 0, -Q >= 0}` (then
  `P0 >= |P1| sqrt(g)`) holds.
- `alternative`: on each leaf one of the declared sufficient sets holds, that is every polynomial of the set is
  `>= 0`. The deduction from the sets to the claim is mathematics, stated in the field `logic`, and is not checked
  by the verifier.

**Covering and positivity.** The tree is a binary partition of the box, so the leaves cover it by construction. On
every cell, `P = sum_I b_I prod_j C(n_j, i_j) z_j^{i_j} (1 - z_j)^{n_j - i_j}` with the variables rescaled to the
cell; if `b_I >= 0` for every `I` (or `> 0`), then `P >= 0` (or `> 0`) on the cell.

## How the verifier differs from the generators

| Step | Generators | Verifier |
|---|---|---|
| Bernstein coefficients | tensor products of matrices (numpy) | explicit formula `b_I = sum_{J <= I} a_J prod_k C(i_k, j_k)/C(n_k, j_k)`, on dictionaries |
| Halving a cell | successive averages (de Casteljau) | direct binomial sums, `left_k = 2^{-k} sum_{j <= k} C(k, j) b_j` |
| Product `P1 g` | computed | recomputed by the verifier |
| `g >= 0` | assumed | checked (Bernstein coefficients at the root `>= 0`) |

The verifier also checks the format (known type, exponents, dimension, duplicate monomials, nodes, non-empty
sufficient sets, indices, delegations only to declared regions), the shape of the tree (truncated, excess nodes), the
envelope data of every leaf (`s0^2 <= g0`, `s1^2 <= g1`, `sigma > 0`), unique names, an acyclic delegation graph, and,
for the corners of group (3), the exact identities `P0(phi) = m P0_del`, `lam P1(phi) = m P1_del`,
`g(phi) = lam^2 g_del` with positive monomials `m`, `lam`. If a coverage is declared `affine`, it checks that `phi` is
a diagonal affine bijection of `[0,1]^d` onto the region (corner C1).

Per-file messages are marked `(local)`; the outcome that counts is the global one on the last line.

## What the verifier does not check

`VERIFIED` means that the algebraic lemmas hold, up to four steps that are mathematical and are proved separately
(supplementary material, paragraph *What the verifier does not check*). The verifier prints them at the end as
external obligations.

1. **The bridge with the geometry.** Each polynomial equals the geometric inequality of the main text written in the
   box variables, up to a positive factor whose sign is checked factor by factor. This is an exact symbolic identity,
   computed by the generators in `../certificate_generators/` and reconstructed independently from the curvature
   (`../symbolic/`, `../independent_checks/`).
2. **The implications** from the sufficient sets of group (2) to `max_s Ric <= 0`, including degenerate cases, and
   the analogous deduction for group (1).
3. **The covering of C2 by the two sectors.** A point with `y <= 1/8` and `a <= 33/400` has either `a <= 2y` (take
   `Y = 8y`, `A' = a/(2y)` in C2A) or `2y <= a` (take `A = 400a/33`, `Y' = 2y/a` in C2B); the case `y = a = 0` follows
   by continuity.
4. **Continuity at the boundary points where the monomial `m` vanishes.** The composed expression is polynomial in the
   box variables and in the square root, hence continuous there.

## Tamper tests

`tamper_tests.py` applies fifteen deliberate alterations and checks that each is rejected for the expected reason:
a perturbed `P0` (radical and radical_square), a false envelope, a truncated tree, a split replaced by a leaf, a sign
change of `g`, a delegated cell outside its region, a wrong sufficient set, a negated `nonneg` polynomial, a
delegation without delegate, an empty sufficient set, an altered composition map, a missing delegate, circular
delegation, and an unknown coverage value. A tampering only counts if it makes the claim false or the proof
incomplete; the script exits with code 0 only if all fifteen are rejected.
