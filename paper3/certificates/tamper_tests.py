"""Checks that must fail: the verifier must reject tampered certificates.

Fifteen deliberate alterations of genuine certificates (or small hand-made ones); each must be rejected, for the
reason printed next to it. Exit code 0 if and only if all fifteen are rejected.
Usage: python3 tamper_tests.py   (from any directory; the certificates are read next to this file)"""
import json, copy, subprocess, sys, os, tempfile, shutil
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
VERIFIER = os.path.join(HERE, 'verify_certificate.py')
def load(name): return json.load(open(os.path.join(HERE, name)))
base = load('group4_B2_I2.json')
def run(C, label):
    d = tempfile.mkdtemp(); p = os.path.join(d, f'tampered_{label}.json'); json.dump(C, open(p, 'w'))
    r = subprocess.run([sys.executable, VERIFIER, p], capture_output=True, text=True); shutil.rmtree(d)
    r.stdout = r.stdout.replace(d + os.sep, '')
    outcome = 'REJECTED (expected)' if r.returncode != 0 else 'ACCEPTED (ERROR!)'
    print(f"{label}: {outcome} | {r.stdout.splitlines()[0][:110]}")
    return r.returncode != 0
ok = []
# 1. P0 perturbed: a large positive constant is subtracted (the polynomial becomes negative somewhere)
C = copy.deepcopy(base); M = max(abs(Fr(c)) for _, c in C['polynomials']['P0'])
k0 = next((i for i, (m, _) in enumerate(C['polynomials']['P0']) if m == [0, 0, 0, 0]), None)
if k0 is None: C['polynomials']['P0'].append([[0, 0, 0, 0], str(-M)])
else: C['polynomials']['P0'][k0][1] = str(Fr(C['polynomials']['P0'][k0][1]) - M)
ok.append(run(C, '1_P0_perturbed'))
# 2. false envelope: s0 of one leaf replaced by a value that is too large
C = copy.deepcopy(base); i = next(k for k, t in enumerate(C['tree']) if t[0] == 'l')
C['tree'][i][1] = str(Fr(C['tree'][i][2])*2 + 1); ok.append(run(C, '2_s0_too_large'))
# 3. truncated tree: the last leaf is removed
C = copy.deepcopy(base); C['tree'] = C['tree'][:-1]; ok.append(run(C, '3_truncated_tree'))
# 4. a split turned into a leaf (covering with cells that are too large)
C = copy.deepcopy(base); j = next(k for k, t in enumerate(C['tree']) if t[0] == 's')
leaf = next(t for t in C['tree'] if t[0] == 'l'); C['tree'][j] = leaf; ok.append(run(C, '4_split_into_leaf'))
# 5. sign of g changed
C = copy.deepcopy(base); C['polynomials']['g'] = [[m, str(-Fr(c))] for m, c in C['polynomials']['g']]; ok.append(run(C, '5_g_negated'))
print('ALL TAMPERINGS REJECTED' if all(ok) else 'AT LEAST ONE TAMPERING ACCEPTED')
# ---- other certificate types ----
ok2 = []
# 6. delegation: the declared region C2 is moved (the delegated cell is no longer contained in it)
C = load('group3_I2_main.json'); C['delegations'][1][1][0] = ['1/2', '1']; ok2.append(run(C, '6_delegation_outside_region'))
# 7. alternative: the leaf names a set that does not hold (set 1 = P5 is forced where P4 or P6 was needed)
C = load('group2_Kinf.json')
for tok in C['tree']:
    if tok[0] == 'l' and tok[1] != 1: tok[1] = 1
ok2.append(run(C, '7_wrong_alternative'))
# 8. radical_square: P0 altered (a negative constant is added)
C = load('group3_I2_C1.json'); M = max(abs(Fr(c)) for _, c in C['polynomials']['P0'])
k0 = next((i for i, (m, _) in enumerate(C['polynomials']['P0']) if m == [0, 0, 0, 0]), None)
if k0 is None: C['polynomials']['P0'].append([[0, 0, 0, 0], str(-M)])
else: C['polynomials']['P0'][k0][1] = str(Fr(C['polynomials']['P0'][k0][1]) - M)
ok2.append(run(C, '8_square_P0_perturbed'))
# 9. nonneg (group 1): sign changed
C = load('group1_P3.json'); C['polynomials']['P'] = [[m, str(-Fr(c))] for m, c in C['polynomials']['P']]; ok2.append(run(C, '9_nonneg_negated'))
print('OTHER TYPES: ALL TAMPERINGS REJECTED' if all(ok2) else 'OTHER TYPES: AT LEAST ONE TAMPERING ACCEPTED')
# ---- hand-made counterexamples and composition ----
ok3 = []
def run_multi(files_dict, label):
    d = tempfile.mkdtemp(); paths = []
    for fn, C in files_dict.items():
        pth = os.path.join(d, fn); json.dump(C, open(pth, 'w')); paths.append(pth)
    r = subprocess.run([sys.executable, VERIFIER] + paths, capture_output=True, text=True)
    shutil.rmtree(d); r.stdout = r.stdout.replace(d + os.sep, '')
    outcome = 'REJECTED (expected)' if r.returncode != 0 else 'ACCEPTED (ERROR!)'
    lines = [l for l in r.stdout.splitlines() if 'REJECTED' in l or ('COMPOSITION' in l and ('missing' in l or 'failed' in l or 'not ' in l))]
    print(f"{label}: {outcome} | {(lines[0] if lines else r.stdout.splitlines()[-1])[:120]}")
    return r.returncode != 0
# 10. nonneg with P = -1 and everything delegated to a region [0,1] without a delegated certificate
C = {'name': 'false', 'type': 'nonneg', 'description': '', 'variables': ['X'], 'map': '', 'polynomials': {'P': [[[0], '-1']]},
     'tree': [['d', 'R']], 'delegations': [['R', [['0', '1']], ['nonexistent']]]}
ok3.append(run_multi({'false.json': C}, '10_delegation_without_delegate'))
# 11. alternative with an empty set
C = {'name': 'empty', 'type': 'alternative', 'description': '', 'variables': ['X'], 'map': '', 'polynomials': {'P': [[[0], '-1']]},
     'polynomial_order': ['P'], 'sufficient_sets': [[]], 'tree': [['l', 0]]}
ok3.append(run_multi({'empty.json': C}, '11_empty_set'))
# 12. composition: phi altered in the delegate C2A
fs = {f: load(f) for f in ['group3_I2_main.json', 'group3_I2_C1.json', 'group3_I2_C2A.json', 'group3_I2_C2B.json']}
fs['group3_I2_C2A.json']['composition']['phi'][3] = [[[1, 0, 0, 1], '1/3']]
ok3.append(run_multi(fs, '12_composition_phi_altered'))
# 13. composition: missing delegate
fs = {f: load(f) for f in ['group3_I2_main.json', 'group3_I2_C1.json', 'group3_I2_C2A.json']}
ok3.append(run_multi(fs, '13_missing_delegate'))
print('COUNTEREXAMPLES AND COMPOSITION: ALL REJECTED' if all(ok3) else 'AT LEAST ONE ACCEPTED')
# ---- circular delegation and unknown coverage ----
ok4 = []
def certA(name, other):
    return {'name': name, 'type': 'radical', 'description': '', 'variables': ['X'], 'map': '',
            'polynomials': {'P0': [[[0], '-1']], 'P1': [[[0], '0']][:0], 'g': [[[0], '1']]},
            'tree': [['d', 'R']], 'delegations': [['R', [['0', '1']], [other]]],
            'composition': {'main': other, 'region': 'R', 'phi': [[[[1], '1']]], 'm': [[[0], '1']], 'lam': [[[0], '1']], 'coverage': 'affine'}}
ok4.append(run_multi({'A.json': certA('A', 'B'), 'B.json': certA('B', 'A')}, '14_circular_delegation'))
main = {'name': 'P', 'type': 'radical', 'description': '', 'variables': ['X'], 'map': '',
        'polynomials': {'P0': [[[1], '1'], [[0], '-1/2']], 'P1': [], 'g': [[[0], '1']]},
        'tree': [['d', 'R']], 'delegations': [['R', [['0', '1']], ['D']]]}
dele = {'name': 'D', 'type': 'radical', 'description': '', 'variables': ['X'], 'map': '',
        'polynomials': {'P0': [[[0], '1/2']], 'P1': [], 'g': [[[0], '1']]}, 'tree': [['l', '0', '1', '1']],
        'composition': {'main': 'P', 'region': 'R', 'phi': [[[[0], '1']]], 'm': [[[0], '1']], 'lam': [[[0], '1']], 'coverage': 'affine_typo'}}
ok4.append(run_multi({'P.json': main, 'D.json': dele}, '15_unknown_coverage'))
print('CIRCULAR DELEGATION AND UNKNOWN COVERAGE: ALL REJECTED' if all(ok4) else 'AT LEAST ONE ACCEPTED')
total = ok + ok2 + ok3 + ok4
print(f'SUMMARY: {sum(total)}/{len(total)} tamperings rejected')
sys.exit(0 if all(total) else 1)
