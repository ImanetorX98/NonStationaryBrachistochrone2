# Global minimality of controlled-rail brachistochrones: computer-assisted steps

This folder of the repository NonStationaryBrachistochrone2 is the reproducibility package for the paper

> I. Rosignoli, *Global minimality of controlled-rail brachistochrones in conformal Kerr and Vaidya spacetimes*,
> with its supplementary material *Alternative proofs and computer-assisted steps*

(Paper III of the series on controlled-rail brachistochrones; Paper I, concept DOI
[10.5281/zenodo.21739997](https://doi.org/10.5281/zenodo.21739997), and Paper II, concept DOI
[10.5281/zenodo.21739999](https://doi.org/10.5281/zenodo.21739999)).

The paper proves global minimality of the controlled-rail brachistochrones by reading the rail as the null cone of an
auxiliary metric. Where a statement concerns specific transcendental quantities (the sign of curvature polynomials on
the equator of Kerr, the critical separation for given data, the constants of the examples), the text reduces it to
finitely many explicit inequalities. This repository contains the programs that certify those inequalities, the
symbolic checks of the derivations, independent second implementations, and the floating-point diagnostics quoted in
the text, together with the logs of our runs.

## Where to start

- [`RESULTS_MAP.md`](RESULTS_MAP.md): every computer-assisted step of the paper, in the order of the text, with the
  script that checks it, the log and the kind of check (exact, certificate, interval arithmetic, numerical).
- [`certificates/README.md`](certificates/README.md): the 33 portable curvature certificates and their independent
  verifier, the shortest way to check the main computer-assisted step.

## Layout

| Folder | Content | Kind |
|---|---|---|
| [`certificates/`](certificates/) | 33 curvature certificates (JSON), the verifier `verify_certificate.py` (Python standard library only), 15 tampering tests, `SHA256SUMS` | exact certificates |
| [`certificate_generators/`](certificate_generators/) | the code that derives the certificate polynomials from the geometry and searches the proof trees; rebuilds all 33 certificates identically | exact |
| [`interval_arithmetic/`](interval_arithmetic/) | the certified grid of critical separations and the interior Maxwell angle, with the archived leaves | interval arithmetic |
| [`riemann_sums/`](riemann_sums/) | the constants of the Vaidya example | interval arithmetic |
| [`symbolic/`](symbolic/) | 40 Wolfram Language scripts: identities and decisions behind the statements of the main text and of the supplement | exact (a few numerical checks marked `[N]`) |
| [`python_checks/`](python_checks/) | symbolic and numerical checks of the auxiliary metric, of global hyperbolicity and of the Kerr barrier | exact and numerical |
| [`independent_checks/`](independent_checks/) | second implementations: Python scripts written independently of the main ones, and Wolfram Language reconstructions of the certificate bridges and of the leaves | exact |
| [`diagnostics/`](diagnostics/) | the floating-point numbers that the paper labels as diagnostics, not certified | numerical |
| [`paper/`](paper/) | LaTeX sources of the paper and of the supplementary material, with a build script | — |
| [`figures/`](figures/) | the scripts that draw the figures, and the figures | — |
| [`results/`](results/) | the logs of our runs, one per script | — |
| [`tools/`](tools/) | a progress bar shared by the slow scripts (it never writes to the logs) | — |
| `run_all.sh` | reruns the checks and applies an acceptance rule to every log | — |

The folders `certificates/`, `certificate_generators/`, `interval_arithmetic/`, `symbolic/`, `independent_checks/`,
`diagnostics/` and `paper/` have a README with their own map, usage and expected output; the scripts of `python_checks/`,
`riemann_sums/` and `figures/` describe themselves in their docstrings.

## What is a proof and what is not

- **Certificates** (`certificates/`) are closed algebraic lemmas: a sign condition for a polynomial with rational
  coefficients on a box, with a Bernstein subdivision tree as proof. The verifier checks them in exact rational
  arithmetic; it is a single self-contained file (it uses the progress bar of `tools/` only if present) and does not
  use the generating code. The identification of
  each polynomial with the geometric inequality of the text is a separate symbolic computation
  (`certificate_generators/bridge.py`, `symbolic/bridge_group4_group5.wls`, `independent_checks/wolfram_bridge_*.wls`),
  and some steps are proved in the supplementary material (paragraph *What the verifier does not check*).
- **Interval arithmetic** (`interval_arithmetic/`, `riemann_sums/`): every certified sign is decided on enclosures
  computed with outward rounding (`mpmath.iv`). Floating-point estimates are used only to choose where the boxes and
  the evaluation points go, never to decide a sign.
- **Symbolic checks** (`symbolic/` in Wolfram Language, `python_checks/` in SymPy) are exact identities, decisions by
  `Reduce` and `Resolve`, and rational arithmetic. Checks marked `[N]` in their labels, and the numerical parts that
  the docstrings of `python_checks/` declare as diagnostics, are numerical confirmations, not certificates.
- **Diagnostics** (`diagnostics/`) are floating-point computations without error bounds, on finitely many sampled
  cases; the paper labels them as such and does not use them in any proof. `RESULTS_MAP.md` gives the sampled ranges.

## Requirements

- Python 3 with NumPy, SciPy, SymPy, mpmath and Matplotlib. Our runs used Python 3.13.2 on macOS (Apple silicon) with
  the versions pinned in [`requirements.txt`](requirements.txt). The certificate verifier needs only the standard
  library.
- Wolfram Language (`wolframscript` with Wolfram Engine or Mathematica) for `symbolic/` and the `wolfram_*.wls` files
  of `independent_checks/`; our runs used Mathematica 13.3.1. Without it the Python checks still run, and the Wolfram
  jobs are reported as skipped.

## Setting up

From this folder, create an environment with the pinned packages. Our logs were produced with Python 3.13; to
recreate that environment, call that version explicitly:

```
python3.13 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
```

`run_all.sh` uses the first `python3` on the `PATH` (the one of the active environment) and the first `wolframscript`.
Two environment variables select other ones:

```
PYTHON=/path/to/python3 ./run_all.sh
WOLFRAMSCRIPT=/path/to/wolframscript ./run_all.sh
```

Before running anything, `run_all.sh` prints the Python and the `wolframscript` it will use and checks that all
packages of `requirements.txt` are discoverable by this Python (it does not compare their versions with the pinned
ones). If one is missing it stops at once, names the missing packages and the command that installs them, and leaves
`results/` untouched.

## Running

Everything at once, with acceptance checks on the recorded logs:

```
./run_all.sh            # the standard set, about one hour
./run_all.sh --full     # also the long runs (proof-tree search, unimodality, Maxwell angle, full leaf rechecks,
                        # figures, diagnostics): many hours
./run_all.sh --check    # only check the logs already in results/
```

It exits with status 0 only if every job passed (2 if some Wolfram jobs were skipped because `wolframscript` is
missing). With `--check` the status is 2 even when every log passes, because four jobs of the full set write no
single log and are reported as skipped. A failed job is reported as `FAIL`, with the exit status or the acceptance rule that was not met; its log in
`results/` contains the output of the failing script, and the summary lists the failed jobs. An interrupted run can be
resumed with `--from K`; `./run_all.sh --list` prints the numbered jobs. Unknown or conflicting options (for example
`--list` together with `--full`) are rejected.

Each script can also be run from its own folder (with the same Python); the commands and the expected durations are
in the READMEs and docstrings. The quickest checks:

```
cd certificates && python3 verify_certificate.py *.json   # about 20 minutes; last line: OUTCOME: all certificates verified
cd certificates && python3 tamper_tests.py                # 15 alterations, all rejected
cd symbolic && ./run_symbolic.sh                          # about 10 minutes (the certificate bridge alone up to 25)
```

The longest runs are `interval_arithmetic/unimodality.py` (about two hours on 10 processes), the full leaf recheck in
Wolfram Language (one to four hours per certificate) and the diagnostics (about seven hours on one core).

**Capturing Wolfram output.** With `wolframscript -file x.wls > log`, printed lines were occasionally lost in our runs;
through a pipe (`wolframscript -file x.wls | cat > log`) the output was complete. The launcher of `symbolic/` uses the
pipe and rejects incomplete logs.

## Licence and citation

Everything in this folder, the code, the LaTeX sources in `paper/` and the figures, is released under the MIT
licence ([`../LICENSE`](../LICENSE)). The repository is archived on Zenodo with concept DOI
[10.5281/zenodo.24039819](https://doi.org/10.5281/zenodo.24039819), which resolves to the latest version. If you use
it, please cite it as described in [`../CITATION.cff`](../CITATION.cff), together with the paper.
