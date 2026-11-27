# lean-nahm: the Andrews–Gordon identities and the E8 and E7 Nahm sums, in Lean 4

Sources of the mathematics: `../rr-identities/proofs/E8-E6.md` (Theorem 1, Corollary 1) and
`../rr-identities/proofs/E7.md` (Theorem 4), both refereed; the Andrews–Gordon identities (Andrews, PNAS 71
(1974), Theorem 1 = *The Theory of Partitions*, Theorem 7.8), which the E8 proof needs at modulus 11 and
which Mathlib does not have.  The library of `../lean-krtw` (q-Pochhammer symbols, Gaussian binomials,
coefficientwise sums, the finite triple product) and the certificates of `../rr-identities/lean`
(RRFusion: Durfee identity) are reused as Lake path dependencies.

## Status (2026-11-27)

**Machine-checked, with no `sorry`, no `axiom`, no `native_decide`** (`#print axioms` gives
`[propext, Classical.choice, Quot.sound]` for every statement, `logs/print_axioms.txt`):

1. **The Andrews–Gordon identities for every modulus `2k+1` and every `1 ≤ i ≤ k`**, in the printed
   product form (`andrews_gordon`, and `andrews_gordon_field` over any field):

       Σ_{n_1,…,n_{k-1} ≥ 0} q^{N_1²+⋯+N_{k-1}² + N_i+⋯+N_{k-1}} / ((q)_{n_1}⋯(q)_{n_{k-1}})
           = ∏_{n ≥ 1, n ≢ 0, ±i (mod 2k+1)} (1 - q^n)⁻¹ .

2. **The E8 vacuum identity** (`e8_vacuum`), with `C(E8)` = Mathlib's `CartanMatrix.E₈` verbatim:

       Σ_{n ∈ ℕ^8} q^{½ nᵀ C(E8)⁻¹ n} / ∏_{i=1}^{8} (q)_{n_i} = 1/(q, q², q³, q⁴, q⁷, q⁸, q⁹, q¹⁰; q¹¹)_∞ ,

   and **its four companions** (`e8_companion_4`, `_3`, `_2`, `_1`; Corollary 1 of `E8-E6.md`):
   `Σ q^{½ nᵀC⁻¹n + B^{(i)}·n}/∏(q)_{n_j} = ∏_{n ≢ 0,±i (11)} (1-q^n)⁻¹` for `i = 1, …, 4`.
   Theorem 1 of `E8-E6.md` (for every linear term `b`) is proved as the internal theorem
   `Nahm.e8_theorem1` (it is stated for the library's tuple sums, so it is not in the statement file).

3. **The E7 identity** (Theorem 4 of `E7.md`), as power series in `X = q^{1/4}`, with `C(E7)` = Mathlib's
   `CartanMatrix.E₇` verbatim, in two forms:

       (one-sided)  Σ_{n ∈ ℕ^7} q^{½ nᵀ C(E7)⁻¹ n} / ∏ (q)_{n_i} = (1/(q;q)_∞) (Σ_{k≥0} q^{3k²/4} + Σ_{k≥1} q^{7k²/4})
       (bilateral)  2 · Σ_{n ∈ ℕ^7} q^{½ nᵀ C(E7)⁻¹ n} / ∏ (q)_{n_i} = (1/(q;q)_∞) Σ_{k ∈ ℤ} (q^{3k²/4} + q^{7k²/4})

   (`e7_identity`, `e7_identity_bilateral`; the printed `1/(2(q)_∞)` is stated without division by
   multiplying by 2; the sum over `k ∈ ℤ` is a genuinely bilateral sum, `Finset.Icc (-N) N` in the
   coefficient of `X^N`; the identity `Σ_ℤ = 2 · (one-sided)` is proved, `e7Theta_eq_two_mul`).

4. The auxiliary statements that make the definitions auditable: the exponents are integers
   (`e8Form_eq_e8Exp`: `½ nᵀ C(E8)⁻¹ n + B·n` is a natural number; `e7Form_eq_e7Exp`: `4·½ nᵀ C(E7)⁻¹ n` is
   a natural number), the inverse matrices are the explicit ones (`invCartanE8_eq_explicit`,
   `invCartanE7_eq_explicit`), the finite sums in the coefficient definitions lose nothing
   (`coeff_agTerm_eq_zero`, `coeff_e8Term_eq_zero`, `coeff_e7Term_eq_zero`, `coeff_e7ThetaTerm_eq_zero`),
   and the truncated products of `qPochInf` are stable (`coeff_qPochInf`).

The statements (21 theorems, 18 definitions) are in `NahmStatement.lean` (imports only Mathlib; proofs
`sorry` on purpose) and are proved verbatim in `NahmSolution.lean` (`logs/compare_statements.txt`).
All statements over `ℚ`; the Andrews–Gordon, E8 vacuum and E7 statements also over an arbitrary field.

## The definitions (`Nahm/Defs.lean`, inlined verbatim in `NahmStatement.lean`)

All series are formal power series in `X` (`K⟦X⟧`); `q = X`, except for E7 where `q = X^4`.
`⁻¹` is Mathlib's inverse in `K⟦X⟧` (all series inverted have constant term 1).

| definition | meaning |
|---|---|
| `qPoch n` | `(q;q)_n = ∏_{i<n} (1 - X^{i+1})` |
| `qPoch4 n` | `(q;q)_n` for `q = X^4`: `∏_{i<n} (1 - X^{4(i+1)})` |
| `qPochInf a d` | `(q^a;q^d)_∞`, by coefficients: coefficient of `X^N` = that of `∏_{i ≤ N} (1 - X^{a+di})` |
| `box k N` | the tuples `Fin k → ℕ` with entries `≤ N` |
| `psum n j` | `N_j = n_j + n_{j+1} + ⋯` (0-based) |
| `agExp k i n` | `N_1² + ⋯ + N_{k-1}² + N_i + ⋯ + N_{k-1}` |
| `agSum k i` | coefficient of `X^N` = `Σ_{n ∈ box (k-1) N} coeff_N (X^{agExp k i n} ∏_j (qPoch n_j)⁻¹)` |
| `agProduct k i` | `(∏_{1 ≤ a ≤ 2k, a ≠ i, a ≠ 2k+1-i} qPochInf a (2k+1))⁻¹` |
| `invCartanE8` | `(CartanMatrix.E₈ cast to ℚ)⁻¹` (Mathlib's matrix inverse) |
| `e8Form B n` | `½ nᵀ (invCartanE8) n + B·n ∈ ℚ` |
| `e8Exp B n` | `⌊e8Form B n⌋₊` (exact: `e8Form_eq_e8Exp`) |
| `e8Sum B` | coefficient of `X^N` = `Σ_{n ∈ box 8 N} coeff_N (X^{e8Exp B n} ∏_i (qPoch n_i)⁻¹)` |
| `invCartanE7`, `e7Form n` | `(CartanMatrix.E₇ cast to ℚ)⁻¹`, `½ nᵀ (invCartanE7) n` |
| `e7Exp n` | `⌊4 · e7Form n⌋₊` (exact: `e7Form_eq_e7Exp`) |
| `e7Sum` | coefficient of `X^N` = `Σ_{n ∈ box 7 N} coeff_N (X^{e7Exp n} ∏_i (qPoch4 n_i)⁻¹)` |
| `e7Theta` | coefficient of `X^N` = `Σ_{k ∈ Icc (-N) N} coeff_N (X^{3k²} + X^{7k²})` |
| `e7ThetaOneSided` | coefficient of `X^N` = `Σ_{k ≤ N} coeff_N (X^{3k²} + [k ≠ 0] X^{7k²})` |

Reading guide.  Each multiple sum is defined by its coefficients as a finite sum over a box; the
theorems `coeff_*Term_eq_zero` say that the summands outside the box have zero coefficient, so the
definition is the coefficient of the infinite sum.  **Labelling.**  The E8 and E7 statements use
Mathlib's `CartanMatrix.E₈`, `CartanMatrix.E₇` (Bourbaki labelling, 0-based in Lean); the vacuum sums do
not depend on the labelling.  The companions' linear terms are the paper's `B^{(i)}` (paper labelling:
chain 1—2—…—7, node 8 on node 5) transported to Bourbaki labelling by paper node `p ↦` Bourbaki node
`(8,7,6,5,4,3,1,2)_p` (1-based; `RRFusion.C8_eq_mathlib`); both vectors are in the docstrings.  The sanity
test checks the transported vectors against the paper's vectors with the paper's matrix (and checks that
the untransported vector gives a different series).

## The proof

| file | content |
|---|---|
| `Nahm/Bailey.lean` | over an arbitrary base `Q` (zero constant term), relative to an arbitrary `A`: the inner sum `Σ_i [N;i] B^i Q^{i²} (BQ^{i+1};Q)_{N-i} = 1` (by a two-term recurrence), **Bailey's lemma** `ρ₁, ρ₂ → ∞` (`IsBP.step`), the unit Bailey pair (`ubp`), the chain from any pair |
| `Nahm/Lattice.lean` | the **change of `a`** of the Bailey lattice (`a = q → a = 1`, `IsBP.toOne`) |
| `Nahm/Cauchy.lean` | `(Q;Q)_∞` (`pinf`), **Cauchy's identity** `Σ_m Q^{m(m+c)}/((Q)_m(Q)_{m+c}) = 1/(Q;Q)_∞`, from the Durfee identity (`RRFusion.durfee_inv`) by a limit |
| `Nahm/Multisum.lean` | tuple sums `Σ_{m ∈ ℕ^r} Q^{Σ(M_l² + b_l M_l)}/∏(Q)_{m_l}` = nested chain sums (`tsum_mterm_eq`) |
| `Nahm/AGTheta.lean` | **Andrews–Gordon in theta form**, all `k`, `i` (`ag_theta_one`, `ag_theta_two`) |
| `Nahm/JTP.lean` | **Jacobi triple product** `(q^c, q^i, q^{c+i}; q^{c+i})_∞ = Σ_{n ≥ 0} (-1)^n q^{c n(n+1)/2 + i n(n-1)/2}(1 - q^{i(2n+1)})` (from `KRTW.limit_tp`); `(q;q)_∞ = ∏_{a=1}^{D} (q^a;q^D)_∞` |
| `Nahm/AG.lean` | the product form (`andrews_gordon_top`) |
| `Nahm/Fusion.lean` | **Lemma 2 of `E8-E6.md` (fusion) for infinite sums** over `ι → ℕ`, any finite `ι` (`fusion`) |
| `Nahm/E8.lean` | `C(E8)⁻¹` (by `decide`), the four fusions `4:=7+8, 6:=1+8, 5:=1+4, 3:=1+2`, the exponent identity (by `ring`), `e8_theorem1`, vacuum and companions |
| `Nahm/E7Lemmas.lean` | Lemma 4 of `E7.md` (`lemma4`), the chain relative to `q^a`, Proposition 5 (`prop5`) |
| `Nahm/Prop6.lean` | Proposition 6 of `E7.md` (`prop6`: the finite identity `c_n = 1 + q^{n²}`) |
| `Nahm/E7.lean`, `Nahm/E7Final.lean` | `C(E7)⁻¹`, the three fusions `7:=1+6, 4:=5+6, 3:=1+2`, peeling `a`, collecting `n = a + 2r`, the final forms |
| `Nahm/Defs.lean`, `Nahm/Final.lean` | statement-level definitions; truncation lemmas, explicit inverses, field versions |

Differences from the paper proofs (none affects the statements):

* Bailey's lemma (ρ₁, ρ₂ → ∞) is **derived**, not quoted: from the finite identity
  `Σ_i [N;i] B^i q^{i²} (Bq^{i+1};q)_{N-i} = 1` (the `b → ∞, d → 0` limit of q-Chu–Vandermonde, proved
  by a recurrence in `N`); Bailey pairs are cleared of the denominators `(Aq;q)_{2n}` as in lean-krtw.
* **Andrews–Gordon for general `i`** is proved by the Bailey lattice (Agarwal–Andrews–Bressoud, 1987):
  `j = k+1-i` steps relative to `q` from the unit pair, the change of `a` (`q → 1`), `i-2` steps relative
  to `1`, and Cauchy's identity for the outermost sum (for `i = 1`, `k-1` steps relative to `q` and
  Cauchy relative to `q`).  Andrews' own 1974 proof (functional equations of `J_{k,i}`) is not used.
  The case needed for E8 is `k = 5`; the companions use `i = 1, …, 4`.
* **Proposition 6 of `E7.md`**: instead of the `₂φ₁` q-Chu–Vandermonde evaluation with negative
  powers of `q`, the vanishing `S_u = 0` (`1 ≤ u ≤ n-1`) is proved with the homogeneous q-binomial
  theorem `∏_{i<u}(w - q^{i+1}) = Σ_r [u;r](-1)^r q^{r(r+1)/2} w^{u-r}` at `w = q^t`, `1 ≤ t ≤ u`, after
  expanding `∏_s (1 - c_s y)` (`Finset.prod_add`), using the symmetry `u ↔ n-u` to keep all exponents
  nonnegative (proof written by a sub-agent against a fixed statement).
* The rearrangements of the formal sums (fibres of the fusion map, `n = a + 2r`, the triangular shifts)
  are done with Mathlib's `tsum` for the product topology of the discrete topology on the coefficient
  field (`K⟦X⟧`, `PowerSeries.WithPiTopology`), with summability from divisibility by powers of `X`.

**No gap or error was found** in `E8-E6.md` §§1–3 or in `E7.md` §§1–3 (all exponent identities of the
fusion chains were confirmed by `ring`; the paper's `B^{(i)}` are confirmed with the transported labels).

## What is NOT machine-checked

* `E8-E6.md` §4 (E6: Theorem 2, Corollary 2) and §6, `E7.md` §6 (Theorem 5, the family `k ≥ 1`;
  it would follow from the same files with `![0, a, a]` replaced by a general chain, not done), the
  refined forms with indeterminates `y_j`, and all modularity statements.
* Theorem 1 of `E8-E6.md` is proved (`e8_theorem1`) for the linear terms `b' = (b_1, b_1+b_2, b_1+b_2+b_3,
  b_1+⋯+b_4)` on the `m`-variables with `b ∈ ℕ⁴` (i.e. `Σ_l b_l M_l` on the `M`-variables), which covers
  the five cases of Corollary 1; the paper's general `b ∈ ℚ⁴` (e.g. negative or decreasing entries) is
  not covered, and the statement file contains only the five instances of Corollary 1.
* The identification of the Lean definitions with the printed identities: a matter of reading
  `Nahm/Defs.lean` (above); the sanity test supports it numerically.  The analytic interpretation
  (`|q| < 1`) is not formalized; everything is an identity of formal power series.
* Novelty/priority claims; the identification with Sun–Wang's (3.20) and with Virasoro characters.

## Sanity test of the definitions

`scripts/Sanity.lean` (`lake env lean scripts/Sanity.lean`, about 20 s) is a computable mirror of the
statement definitions (it calls the very `Nahm.psum`, `Nahm.agExp` of the statement file and copies the
other formulas symbol by symbol; the noncomputable exponents are computed from the explicit inverse
matrices of `invCartanE8_eq_explicit`, `invCartanE7_eq_explicit`).  Its `#guard`s check that the first
**41 coefficients** of both sides of: Andrews–Gordon for `(k,i) ∈ {(2,1),(2,2),(3,1),(3,2),(3,3),(4,2),
(5,1),…,(5,5)}`, the E8 vacuum (sum and the product `1/(q,…,q¹⁰;q¹¹)_∞`), the four E8 companions, and E7
(sum, one-sided and bilateral right sides, the factor 2) equal the integers computed independently by
`scripts/sanity_independent.py` (plain Python: Cartan matrices built from the paper's diagrams and
inverted with exact fractions, the paper's `B^{(i)}` in the paper's labelling, depth-first enumeration;
it also checks each identity itself to `X^40`).  E8 vacuum: `1, 1, 2, 3, 5, 6, 9, 12, 17, 22, 30, 38, …`;
E7 (in `X = q^{1/4}`): `1, 0, 0, 1, 1, 0, 0, 2, 2, 0, 0, 3, 4, …`.  Five negative controls are detected
(wrong residue `i`, wrong `B` vector, the paper-labelled `B^{(4)}` used with the Bourbaki matrix, the
missing factor 2 in the bilateral E7 form).  Outputs: `logs/sanity.txt`, `logs/sanity_independent.txt`.
(The mirror is linked to the noncomputable definitions by inspection, not by a proof.)

## Pins and how to check

Lean `leanprover/lean4:v4.33.0-rc2`, Mathlib `51e6992efd06126df61a496bebf8f49482a4e129` (the pins of
`../lean-krtw` and `../rr-identities/lean`, which are path dependencies).

    curl -sSfL https://raw.githubusercontent.com/leanprover/elan/master/elan-init.sh | sh -s -- -y --default-toolchain none
    export PATH="$HOME/.elan/bin:$PATH"
    cd lean/lean-krtw && lake exe cache get     # Mathlib binaries; never build Mathlib from source
    cd ../lean-nahm && mkdir -p .lake && ln -sfn ../../lean-krtw/.lake/packages .lake/packages
    lake build                                         # "declaration uses sorry" only from NahmStatement.lean (deliberate)
    lake env lean scripts/PrintAxioms.lean
    ./scripts/compare_statements.sh
    lake env lean scripts/Sanity.lean
    python3 scripts/sanity_independent.py

`scripts/compare_statements.sh` checks that `NahmStatement.lean` imports only Mathlib, that its
definition block is byte for byte the block of `Nahm/Defs.lean`, that each statement (docstring
included) occurs byte for byte in `NahmSolution.lean`, that the solution side and the libraries it
imports (`Nahm/`, `../lean-krtw/KRTW/`, `../rr-identities/lean/RRFusion/`) contain no `sorry` / `axiom` /
`native_decide`, and that `#check` / `#print` give the same output from both modules; `comparator.json`
is the configuration for `leanprover/comparator` (not installed here).  `NahmStatement.lean` is
generated by `scripts/gen_statement.py`.  The library is about 2700 lines (plus lean-krtw's).
