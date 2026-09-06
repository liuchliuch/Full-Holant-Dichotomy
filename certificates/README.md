# Exact finite verification

Verification code for **A Full Complexity Dichotomy for Complex-Valued
Boolean Holant Problems**, by Chenghua Liu, Boning Meng, and Juqiu Wang.
The [theorem-to-script map](theorem-to-script-map.md) identifies the finite
claims checked by each program.

All companion verification code was written by OpenAI Codex and subsequently
reviewed by the authors.

## Requirements

- Python 3.10 or later, with assertions enabled (do not use `python -O`).
- For the six-port high-charge calculation: a GCC- or Clang-compatible C++17
  compiler and **Z3 4.14.0**. The other entries use only Python's standard library.

## Run

From the repository root:

```sh
python3 certificates/run_all_exact.py
```

This runs all four entries and returns a nonzero exit code if any fails.
Success ends with `PAPER 2 EXACT-VERIFIER SUITE: PASS`.

To install the required Z3 executable in a virtual environment on macOS or Linux:

```sh
python3 -m venv .venv
.venv/bin/python -m pip install z3-solver==4.14.0.0
.venv/bin/python certificates/run_all_exact.py --q6-z3 "$PWD/.venv/bin/z3"
```

On Windows, use the corresponding `Scripts` paths and `z3.exe`.
Use `--q6-compiler PATH` to select a compiler outside `PATH`.
The standalone six-port wrapper uses `--compiler PATH` and `--z3 PATH`.

| Entry | Calculation |
| --- | --- |
| [`verify_odd_dihedral.py`](verify_odd_dihedral.py) | Order-three projective quaternary set, six-port exclusions, and support classification. |
| [`verify_high_charge_q6.py`](verify_high_charge_q6.py) | Six-port loop/two-copy implication: 24 orbit families and exact polynomial checks. |
| [`verify_high_charge_q8_atlas.py`](verify_high_charge_q8_atlas.py) | Eight-port matching families and support checks. |
| [`verify_high_charge_q8_normalized.py`](verify_high_charge_q8_normalized.py) | Supplementary eight-port linear-system and relation-lattice checks. |

Each entry can also be run as `python3 certificates/FILENAME.py`.
The six-port wrapper builds `verify_high_charge_q6_orbits.cpp` in a temporary
folder, then invokes `verify_high_charge_q6_algebra.py`. All finite objects
and solver queries are generated on each run; no stored run data are required.
The odd-dihedral `--skip-support` and `--skip-switching` options give only
`PARTIAL PASS`; use the default invocation for full coverage.

## Optional output

To retain a transcript and the generated six-port queries and solver responses:

```sh
python3 certificates/run_all_exact.py \
  --transcript verification.log --q6-audit-dir q6-run
```

Use a new output directory. These paths are relative to the invocation directory.
The standalone six-port wrapper uses `--audit-dir` for the same export.

## Scope

The calculations use exact arithmetic, without random sampling or numerical
tolerances. The six-port calculation relies on Z3 4.14.0 UNSAT answers;
no solver proof objects are requested or independently checked. The wrapper
checks the reported version, not a fixed binary hash.

The programs check the finite statements in the map. Algebraic reductions,
gadget realizability, normalization, and arbitrary-arity arguments remain in
the paper. The two eight-port programs reuse matching-generation logic and
are not independent implementations.
