# RRFusion: Lean 4 certificates for `proofs/E8-E6.md`

A small Lean 4 + Mathlib project (Lake package `RRFusion`) that machine-checks the finite and
purely algebraic ingredients of `../proofs/E8-E6.md` (Nahm sums of `C(E8)⁻¹` and `C(E6)⁻¹`).

* Toolchain: `leanprover/lean4:v4.33.0-rc2` (see `lean-toolchain`)
* Mathlib: revision `51e6992efd06126df61a496bebf8f49482a4e129` (see `lakefile.toml`, `lake-manifest.json`)
* No `sorry`, no additional axioms, no `native_decide`.

Build (Mathlib must already be fetched/built in `.lake`, e.g. via `lake exe cache get`):

    export PATH=$HOME/.elan/bin:$PATH
    lake build

`RRFusion.lean` imports the four files `RRFusion/Cartan.lean`, `RRFusion/Fusion.lean`,
`RRFusion/Chains.lean`, `RRFusion/Durfee.lean`. Everything lives in the namespace `RRFusion`.

## Conventions

Node labelling as in the document (1-based there, 0-based `Fin` indices in Lean):

    E8:  1 — 2 — 3 — 4 — 5 — 6 — 7          E6:  1 — 2 — 3 — 4 — 5
                         |                                |
                         8                                6

Matrices of Nahm data are functions `ι → ι → ℚ` on a finite index type `ι`; fusing two species
enlarges the index type to `Option ι`, the new species being `none`.

## Main theorems

### `RRFusion/Cartan.lean` — the Cartan matrices and their inverses

Definitions: `C8 : Matrix (Fin 8) (Fin 8) ℤ` (Cartan matrix of E8), `A8` (the integer matrix printed
in §2 of the document), `C6 : Matrix (Fin 6) (Fin 6) ℚ` (Cartan matrix of E6), `A6` (the rational
matrix printed in §4).

| Name | Statement | Certified by |
|---|---|---|
| `C8_mul_A8` | `C8 * A8 = 1`: the printed 8×8 matrix `A8` is a right inverse of the E8 Cartan matrix (over ℤ). | `decide` (kernel evaluation) |
| `A8_mul_C8` | `A8 * C8 = 1`: it is also a left inverse; so `A8 = C(E8)⁻¹`. | `decide` |
| `C6_mul_A6` | `C6 * A6 = 1`: the printed 6×6 rational matrix `A6` is a right inverse of the E6 Cartan matrix (over ℚ). | `fin_cases` on both indices (36 entries), then `simp` + `norm_num` |
| `A6_mul_C6` | `A6 * C6 = 1`: it is also a left inverse; so `A6 = C(E6)⁻¹`. | `fin_cases`, `simp`, `norm_num` |

Auxiliary sanity checks in the same file: `C8_eq_diagram`, `C6_eq_diagram` (the matrices are "2 on the
diagonal, −1 for adjacent nodes" of the diagrams above), `C8_eq_mathlib`, `C6_eq_mathlib` (they are
Mathlib's `CartanMatrix.E₈`, `CartanMatrix.E₆` after an explicit relabelling of the nodes),
`A8_symm`, `A6_symm` — the first five by `decide`, `A6_symm` by `fin_cases` + `simp`.

### `RRFusion/Fusion.lean` — the fusion operation of Lemma 2 and identity (1.2)

Definitions: `Q A B n = ½ Σ_{i,j} A i j n_i n_j + Σ_i B i n_i`; `fuse A a b` (the matrix `A°` of Lemma 2
on `Option ι`: `A°_{ij} = A_{ij} + [{i,j} = {a,b}]`, `A°_{i0} = A°_{0i} = A_{ia} + A_{ib}`,
`A°_{00} = A_{aa} + A_{bb} + 2A_{ab}`); `fuseB B a b` (`B°`: `B°_i = B_i`, `B°_0 = B_a + B_b`);
`extend m c` (the vector `(m, c)`); `pairVec a b` (the vector `v = e_a + e_b`).

| Name | Statement | Certified by |
|---|---|---|
| `fuse_symm` | If `A` is symmetric then `fuse A a b` is symmetric (any index type, any `a`, `b`). | case split on `Option`, `rfl` / `simp only` |
| `Q_fuse` | **Identity (1.2).** For a finite index type `ι`, symmetric `A : ι → ι → ℚ`, any `B`, indices `a ≠ b`, any `m : ι → ℚ` and `c : ℚ`: `Q (fuse A a b) (fuseB B a b) (m, c) = Q A B (m + c·(e_a + e_b)) + m_a · m_b`. | general proof (finite-sum manipulations, `ring`); no case enumeration |

### `RRFusion/Chains.lean` — the fusion chains of §2 (E8) and §4 (E6)

Definitions: `AG4 i j = 2·min(i,j)` for `1 ≤ i,j ≤ 4` (written `2 * (min i j + 1)` with 0-based
indices), the Andrews–Gordon matrix of rank 4; `E8step1 … E8step4`, the four successive fusions
(node 4 := 7 + 8, node 6 := 1 + 8, node 5 := 1 + 4, node 3 := 1 + 2); `nodeE8 : Fin 8 → Option⁴ (Fin 4)`,
the position of each E8 node in the fused index set. Similarly `E6start = (4/3 2/3 0; 2/3 4/3 0; 0 0 2)`,
`E6step1 … E6step3` (node 2 := 5 + 6, node 4 := 1 + 6, node 3 := 1 + 2), `nodeE6 : Fin 6 → Option³ (Fin 3)`.

| Name | Statement | Certified by |
|---|---|---|
| `nodeE8_bijective` | `nodeE8` is a bijection from the 8 nodes of E8 onto the 4-fold fused index set `Option⁴ (Fin 4)` (so the relabelling loses and duplicates nothing). | `decide` |
| `E8_chain` | For all `i j : Fin 8`: `E8step4 (nodeE8 i) (nodeE8 j) = A8 i j`. I.e. the 4-fold fusion of `2 min(i,j)` is, after relabelling, exactly `A8 = C(E8)⁻¹`. | `fin_cases` on both indices (64 entries), `simp`, `norm_num` |
| `E8_chain_linear` | For every `b : Fin 4 → ℚ`, the 4-fold `fuseB` of `b` is, after relabelling, the vector `B(b) = (b₁, b₃, b₁+b₃, b₂+b₄, b₁+b₂+b₄, b₁+b₄, b₂, b₄)` of Theorem 1 (1-based `b`). | `fin_cases`, `simp` |
| `nodeE6_bijective` | `nodeE6` is a bijection from the 6 nodes of E6 onto `Option³ (Fin 3)`. | `decide` |
| `E6_chain` | For all `i j : Fin 6`: `E6step3 (nodeE6 i) (nodeE6 j) = A6 i j`. I.e. the 3-fold fusion of `E6start` is, after relabelling, exactly `A6 = C(E6)⁻¹`. | `fin_cases` (36 entries), `simp`, `norm_num` |
| `E6_chain_linear` | For every `b : Fin 3 → ℚ`, the 3-fold `fuseB` of `b` is, after relabelling, `B = (b₁, b₂+b₃, b₁+b₂+b₃, b₁+b₃, b₂, b₃)` of Theorem 2. | `fin_cases`, `simp` |

### `RRFusion/Durfee.lean` — Lemma 1 (Durfee rectangle identity)

Definitions, for a commutative ring `R` and `q : R`:
`qPoch q n = ∏_{j<n} (1 − q^{j+1}) = (q;q)_n`;
`gauss q n k`, the Gaussian binomial `[n;k]_q`, defined by `[n;0] = 1`, `[0;k+1] = 0`,
`[n+1;k+1] = [n;k] + q^{k+1}[n;k+1]` (q-Pascal);
`qFalling q a c = ∏_{j<c} (1 − q^{a−j})`;
`durfeeSum q a b = Σ_{c=0}^{a} q^{(a−c)(b−c)} · qFalling q a c · gauss q b c`.

| Name | Statement | Certified by |
|---|---|---|
| `durfee_gauss` | **Polynomial form.** In any commutative ring `R`, for every `q : R` and all `a b : ℕ`: `Σ_{c=0}^{min(a,b)} q^{(a−c)(b−c)} · [a;c]_q · [b;c]_q · (q;q)_c = 1`. (In particular in `ℤ[X]` with `q = X`.) | induction, see below |
| `durfee_inv` | **Inverse form.** In any commutative ring, if `u : ℕ → R` satisfies `u n · (q;q)_n = 1` for all `n`, then for all `a b`: `u_a · u_b = Σ_{c=0}^{min(a,b)} q^{(a−c)(b−c)} · u_{a−c} · u_{b−c} · u_c`. | from `durfee_gauss` and the product formula, `linear_combination` |
| `durfee_field` | **Field form, identity (1.1) as printed.** In a field `K`, if `q^n ≠ 1` for all `n > 0`, then `1/((q)_a (q)_b) = Σ_{c=0}^{min(a,b)} q^{(a−c)(b−c)} / ((q)_{a−c} (q)_{b−c} (q)_c)`. | from `durfee_inv` with `u n = ((q)_n)⁻¹` |
| `durfee_ratFunc` | The same identity in the rational function field `RatFunc ℚ = ℚ(q)` with `q = X`. | `durfee_field`; `X^n ≠ 1` by comparing degrees |

Supporting lemmas: `gauss_mul_qPoch` (`[a;c]_q (q)_c = qFalling q a c`, for all `a`, `c`),
`gauss_mul_qPoch_mul_qPoch` (`[n;k]_q (q)_k (q)_{n−k} = (q)_n` for `k ≤ n`, which identifies `gauss` with
the usual Gaussian binomial), `durfeeSum_succ_succ` (the recurrence
`D(a+1,b+1) = (1 − q^{a+1}) D(a,b) + q^{a+1} D(a+1,b)` for `D = durfeeSum`), `durfeeSum_eq_one`
(`D(a,b) = 1`, by induction on `b` from the recurrence and `D(a,0) = D(0,b) = 1`).

The Lean proof is thus a direct induction on the q-Pascal recurrence; it does **not** follow the
generating-function proof of the document (Euler + q-binomial theorem), but proves the same statement.

## What the finite checks rest on

* `decide`: `C8_mul_A8`, `A8_mul_C8`, `nodeE8_bijective`, `nodeE6_bijective` (and the auxiliary
  `C8_eq_diagram`, `C8_eq_mathlib`, `C6_eq_diagram`, `C6_eq_mathlib`, `A8_symm`). These are evaluated by
  the Lean kernel; no `native_decide` (no trust in the compiler) is used anywhere.
* `fin_cases` + `simp` + `norm_num`: `C6_mul_A6`, `A6_mul_C6`, `E8_chain`, `E6_chain`,
  `E8_chain_linear`, `E6_chain_linear` (and `A6_symm`): the finitely many entries are enumerated and each
  is a rational-arithmetic identity, proved by proof-producing tactics (checked by the kernel).
* `Q_fuse`, `fuse_symm` and everything in `Durfee.lean` are general theorems, not finite checks.

## `#print axioms`

Output of `lake env lean` on a scratch file containing `import RRFusion`, `open RRFusion` and one
`#print axioms` line per theorem (the scratch file is not part of the project):

    'RRFusion.C8_mul_A8' depends on axioms: [propext, Classical.choice, Quot.sound]
    'RRFusion.A8_mul_C8' depends on axioms: [propext, Classical.choice, Quot.sound]
    'RRFusion.C6_mul_A6' depends on axioms: [propext, Classical.choice, Quot.sound]
    'RRFusion.A6_mul_C6' depends on axioms: [propext, Classical.choice, Quot.sound]
    'RRFusion.Q_fuse' depends on axioms: [propext, Classical.choice, Quot.sound]
    'RRFusion.fuse_symm' depends on axioms: [propext, Classical.choice, Quot.sound]
    'RRFusion.E8_chain' depends on axioms: [propext, Classical.choice, Quot.sound]
    'RRFusion.E8_chain_linear' depends on axioms: [propext, Classical.choice, Quot.sound]
    'RRFusion.E6_chain' depends on axioms: [propext, Classical.choice, Quot.sound]
    'RRFusion.E6_chain_linear' depends on axioms: [propext, Classical.choice, Quot.sound]
    'RRFusion.nodeE8_bijective' depends on axioms: [propext, Classical.choice, Quot.sound]
    'RRFusion.nodeE6_bijective' depends on axioms: [propext, Classical.choice, Quot.sound]
    'RRFusion.durfee_gauss' depends on axioms: [propext, Classical.choice, Quot.sound]
    'RRFusion.durfee_inv' depends on axioms: [propext, Classical.choice, Quot.sound]
    'RRFusion.durfee_field' depends on axioms: [propext, Classical.choice, Quot.sound]
    'RRFusion.durfee_ratFunc' depends on axioms: [propext, Classical.choice, Quot.sound]

These are the three standard axioms of Lean/Mathlib; there is no `sorryAx` and no `Lean.ofReduceBool`
/ `Lean.trustCompiler` (which `native_decide` would introduce).

## What is NOT formalized

The project certifies the finite matrix computations and the two algebraic identities behind the proof.
It does **not** formalize the theorems of the document themselves. In particular:

* **Lemma 2 for infinite sums.** Only its algebraic core is proved: the exponent identity (1.2)
  (`Q_fuse`) and the inner finite sum (Lemma 1, `durfee_*`). The Nahm sums `F_{A,B}` as formal power
  series (in a fractional power of `q`), admissibility, the summability of the families involved, the
  finiteness of the fibres of `(m, c) ↦ m + c(e_a + e_b)`, and the regrouping of the summable family
  that yields `F_{A°,B°} = F_{A,B}` are not formalized. `Q_fuse` and `durfee_inv` are also not combined
  in Lean into a single term-by-term statement.
* **The Andrews–Gordon identities** (and the other known identities quoted in the document: Rogers–
  Ramanujan, the rank 3 identities used for E6, …) are taken from the literature; none is formalized.
* **The final power-series identities**: Theorems 1 and 2 and Corollaries 1 and 2 of the document
  (the Nahm sums of `C(E8)⁻¹` and `C(E6)⁻¹` as sums/products) are not stated or proved in Lean. Positive
  definiteness of the matrices and the modularity statements are likewise outside the formalization.

What Lean does guarantee: the matrices printed in the document are the inverses of the E8 and E6 Cartan
matrices; iterating the fusion operation of Lemma 2 along the stated chains transforms `2 min(i,j)`
(rank 4) into `C(E8)⁻¹` and `E6start` into `C(E6)⁻¹`, with the stated linear terms and a bijective
relabelling of the nodes; each fusion step satisfies the quadratic-form identity (1.2); and the
Durfee rectangle identity (1.1) holds.
