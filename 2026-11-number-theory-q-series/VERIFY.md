# How to verify

## Lean proofs (machine-checked)

Requirements: [elan](https://github.com/leanprover/elan). Each project pins its Lean toolchain and
Mathlib commit; `lake exe cache get` downloads prebuilt Mathlib (do not build Mathlib from source).

### Kanade–Russell–Tsuchioka–Warnaar identity (i), every k

```bash
cd lean/lean-krtw
lake exe cache get
lake build
lake env lean scripts/PrintAxioms.lean      # every statement: [propext, Classical.choice, Quot.sound]
./scripts/compare_statements.sh             # statement file (Mathlib only) matches the proved file
lake env lean scripts/Sanity.lean           # Lean definitions reproduce the first 31 coefficients
python3 scripts/sanity_independent.py       # independent Python values for the same coefficients
```

The headline theorem is `KRTW.krtw_identity` in `KRTWStatement.lean`; read the definitions
there against the printed identity in arXiv:2404.03851.

### E8, E7 and Andrews–Gordon

```bash
cd lean/lean-nahm         # uses ../lean-krtw and ../rr-identities/lean
lake exe cache get
lake build
lake env lean scripts/PrintAxioms.lean
./scripts/compare_statements.sh
lake env lean scripts/Sanity.lean          # first 41 coefficients of both sides
python3 scripts/sanity_independent.py      # the same coefficients, computed independently
```

### τ(n) ≠ ±461, ±599 (Diophantine core)

```bash
cd lean/lean-tau
lake exe cache get
lake build
lake env lean scripts/PrintAxioms.lean
```

### Chebyshev's 1853 assertion at every x > 0 (GRH for χ₋₄ as the only hypothesis)

```bash
cd lean/lean-cheb-allx     # uses ../lean-allx and ../lean-chebyshev
lake exe cache get
lake build
lake env lean scripts/PrintAxioms.lean
```

## Certified computations

- ζ(5) constants: `03-ramanujan-machine-zeta5/U1/verify/`, `U2/verify/`, `U3/verify/` (exact
  rational and high-precision checks; `python3 <script>`).
- Chebyshev-bias bounds under GRH alone: `04-chebyshev-bias-grh-only/density/reflection_certified.py`
  (ball arithmetic, needs `pip install python-flint`; output in `reflection_certified_results.txt`).
- τ(n) ≠ ±461, ±599: `06-ramanujan-tau/tau-p4/check_certificate.py` re-verifies every line of
  `certificate_tau_p4.txt`. τ(p²) ≠ 601: `06-ramanujan-tau/tau-p2/local_k_check.py`.
- Nahm-sum identities: `02-exceptional-nahm-sums/code/` and `certificates/` (exact matrix
  certificates and power-series checks).

Python scripts need `numpy`, `mpmath`, `sympy`, and for ball arithmetic `python-flint`.
