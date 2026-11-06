# Certificate generators

The code that produced the 33 curvature certificates in `../certificates`, starting from the metric. Checking the
certificates does not need any of it: `../certificates/verify_certificate.py` is independent of these files. What this
folder adds is the derivation of every certificate polynomial from the curvature, so that the first external
obligation printed by the verifier (the bridge between polynomials and geometry) can be rechecked.

| File | Content |
|---|---|
| `curvature.py` | Christoffel symbols and Ricci tensors on the equator, computed from the metric: the representative `G/N` with constant factor (group 2); `h = G/N` with a time-dependent factor, as rational functions of the jets `(u, u', u'')` (groups 3 and 4); the Hessian of `r` (group 4); the orthonormal-frame formulas of the reference `b = infinity` (group 1). |
| `polynomials.py` | The polynomials in the physical variables: closed forms for groups 1 to 3, the weighted representative of group 4 assembled from `curvature.py`, and `alpha`, `Q` for group 5. |
| `bridge.py` | Identities between the closed forms and the curvature: exact for groups 1 and 2; for group 3, the Ricci tensor of the time-dependent metric recomputed from scratch and compared at 40 rational points with 50 digits (the exact identities are in `../symbolic/bridge_group3.wls`). |
| `build_certificates.py` | Rebuilds every certificate polynomial in its box variables and checks that it equals the published one exactly, together with delegations, composition data and sufficient sets. With `--search` it also reruns the subdivision and checks that every proof tree is reproduced exactly. |
| `search.py` | The Bernstein subdivision that produced the proof trees (two splitting rules, both kept so that the trees are reproduced exactly). |

## Usage

```sh
python3 bridge.py                         # about 4 minutes
python3 build_certificates.py             # polynomials and data, about 10 minutes
python3 build_certificates.py --search    # also the proof trees, about half an hour
python3 build_certificates.py --only group1,group2 --search   # a quick subset
```

Requirements: Python 3 with SymPy and NumPy (versions in `../requirements.txt`). Expected outputs are in
`../results/`.

## Sign of the denominators

In group 4 the inequalities are multiplied by a common denominator whose sign is taken at one sample point
(`group4_assemble`). The denominators are products of factors of known constant sign on the box (`r`, `Delta`, `w`,
`r - R`, `H'`, powers of `1 + t^2`); this is part of the mathematical argument in the supplementary material, not
something the sample point proves. In group 5 the sign of each denominator is itself certified: `group5_den_alpha`
and `group5_den_Q` prove that the signed denominators are positive on the whole box.
