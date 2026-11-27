# lean-krtw: the first Kanade–Russell–Tsuchioka–Warnaar Andrews–Gordon-type identity, in Lean 4

Source of the mathematics: `../rr-identities/proofs/KRTW-general.md` (referee reports
`../rr-identities/REFEREE-KRTW-general.md`, `../rr-identities/REFEREE-2-KRTW.md`).
The identity (S. Kanade, M. C. Russell, S. Tsuchioka, S. O. Warnaar, *Remarks on the conjectures of Capparelli,
Meurman, Primc and Primc*, arXiv:2404.03851v2, displayed after Conjecture 4.9; conjectural there for k ≥ 3):

    Σ_{r_1≥…≥r_k≥0, 0≤s_1≤…≤s_k} ∏_{i=1}^k q^{(r_i+s_i)²+s_i²+s_i} / ((q;q)_{r_i−r_{i+1}} (q²;q²)_{s_i−s_{i−1}})
        = (q, q^{k+1}, q^{k+2}; q^{k+2})_∞ / ((q;q²)_∞ (q;q)_∞)            (r_{k+1} = s_0 = 0).

## Status (2026-11-26)

**The identity is machine-checked for every k ≥ 0** (`KRTW.krtw_identity`), in the product form printed by KRTW.
No `sorry`, no `axiom`, no `native_decide`; `#print axioms` gives `[propext, Classical.choice, Quot.sound]` for
every statement (`logs/print_axioms.txt`).

```lean
theorem krtw_identity (k : ℕ) :
    krtwSum (1 : ℚ⟦X⟧) k
      = qPochInf 1 (k + 2) * qPochInf (k + 1) (k + 2) * qPochInf (k + 2) (k + 2)
        * (qPochInf 1 2 * qPochInf 1 1)⁻¹
```

* **Stage 1 (the new mathematics)**: the 2k-fold KRTW sum equals the k-fold "Melzer-type" sum, with the
  refinement by `x` (Theorem 2, first equality, of the write-up): `krtw_sum_eq_melzer_sum_refined`,
  `krtw_sum_eq_melzer_sum`. Sections 1, 2.
* **Stage 2 (the known part, re-proved)**: the k-fold sum is the product, by Bailey's lemma with the unit
  Bailey pair, the limit of the Bailey chain, and the Jacobi triple product (proved here for the instance
  needed; Mathlib has no triple product): `krtw_identity_theta`, `theta_eq_triple_product`, `melzer_identity`,
  `krtw_identity`. Section 3.

## 1. What is proved (Stage 1)

The statements are in `KRTWStatement.lean` (imports only Mathlib; definitions + statements, proofs `sorry` on
purpose) and are proved verbatim in `KRTWSolution.lean`.

```lean
-- for an arbitrary field K and an ARBITRARY power series x : K⟦X⟧   (q = X)
theorem krtw_sum_eq_melzer_sum_refined {K : Type*} [Field K] (x : K⟦X⟧) (k : ℕ) :
    krtwSum x k = melzerSum x k

-- x = 1, K = ℚ
theorem krtw_sum_eq_melzer_sum (k : ℕ) :
    krtwSum (1 : ℚ⟦X⟧) k = melzerSum (1 : ℚ⟦X⟧) k

-- the finite sums in the definitions lose nothing
theorem coeff_krtwTerm_eq_zero {K : Type*} [Field K] (x : K⟦X⟧) {k : ℕ} (r s : Fin k → ℕ) (N : ℕ)
    (h : r ∉ box k N ∨ s ∉ box k N) : coeff N (krtwTerm x r s) = 0
theorem coeff_melzerTerm_eq_zero {K : Type*} [Field K] (x : K⟦X⟧) {k : ℕ} (M : Fin k → ℕ) (N : ℕ)
    (h : M ∉ box k N) : coeff N (melzerTerm x M) = 0
```

In words: for every k ≥ 0, as formal power series in q,

    Σ_{r_1≥…≥r_k≥0, 0≤s_1≤…≤s_k} x^{Σ(r_i+2s_i)} ∏_i q^{(r_i+s_i)²+s_i²+s_i} / ((q;q)_{r_i−r_{i+1}} (q²;q²)_{s_i−s_{i−1}})
        = Σ_{M_1≥…≥M_k≥0} x^{M_1+…+M_k} ∏_i q^{M_i(M_i+1)/2} / (q;q)_{M_i−M_{i+1}} .

`x` is an arbitrary element of `K⟦q⟧`; with `K = ℚ(x)` and the constant series `x` this is the identity of
`ℤ[x][[q]]` of the write-up (the inclusion `ℤ[x][[q]] ⊂ ℚ(x)[[q]]` is injective); `x = 1` is the unrefined identity.

### The definitions (all of them; `KRTW/Defs.lean`, inlined verbatim in `KRTWStatement.lean`)

```lean
variable {K : Type*} [Field K]
noncomputable def qPoch (n : ℕ) : K⟦X⟧ := ∏ i ∈ range n, (1 - X ^ (i + 1))            -- (q;q)_n
noncomputable def q2Poch (n : ℕ) : K⟦X⟧ := ∏ i ∈ range n, (1 - X ^ (2 * (i + 1)))     -- (q²;q²)_n
def nxt {k : ℕ} (r : Fin k → ℕ) (i : Fin k) : ℕ :=                                     -- r_{i+1}, r_{k+1} := 0
  if h : (i : ℕ) + 1 < k then r ⟨i + 1, h⟩ else 0
def prv {k : ℕ} (s : Fin k → ℕ) (i : Fin k) : ℕ :=                                     -- s_{i-1}, s_0 := 0
  if h : (i : ℕ) = 0 then 0 else s ⟨i - 1, by omega⟩
def box (k N : ℕ) : Finset (Fin k → ℕ) := Fintype.piFinset fun _ => range (N + 1)      -- tuples in {0..N}^k
noncomputable def krtwTerm (x : K⟦X⟧) {k : ℕ} (r s : Fin k → ℕ) : K⟦X⟧ :=
  ∏ i : Fin k, x ^ (r i + 2 * s i) * X ^ ((r i + s i) ^ 2 + s i ^ 2 + s i)
    * (qPoch (r i - nxt r i))⁻¹ * (q2Poch (s i - prv s i))⁻¹
noncomputable def krtwSum (x : K⟦X⟧) (k : ℕ) : K⟦X⟧ :=
  PowerSeries.mk fun N => ∑ r ∈ box k N, ∑ s ∈ box k N,
    if (∀ i, nxt r i ≤ r i) ∧ (∀ i, prv s i ≤ s i) then coeff N (krtwTerm x r s) else 0
noncomputable def melzerTerm (x : K⟦X⟧) {k : ℕ} (M : Fin k → ℕ) : K⟦X⟧ :=
  ∏ i : Fin k, x ^ (M i) * X ^ (M i * (M i + 1) / 2) * (qPoch (M i - nxt M i))⁻¹
noncomputable def melzerSum (x : K⟦X⟧) (k : ℕ) : K⟦X⟧ :=
  PowerSeries.mk fun N => ∑ M ∈ box k N,
    if ∀ i, nxt M i ≤ M i then coeff N (melzerTerm x M) else 0
```

Reading guide. A tuple `(r_1,…,r_k)` is a function `r : Fin k → ℕ` (`r 0 = r_1`, …). `⁻¹` is the inverse in
`K⟦X⟧` (Mathlib's `PowerSeries.inv`; the series inverted have constant term 1). Each side is defined **by its
coefficients as a finite sum**: the coefficient of `q^N` is the sum over the index tuples with all entries
`≤ N` (explicit bound) that satisfy the chain conditions `r_{i+1} ≤ r_i`, `s_{i−1} ≤ s_i` (resp.
`M_{i+1} ≤ M_i`) of the coefficient of `q^N` of the summand. The bound is harmless: a summand with an entry `> N`
is divisible by `q^{N+1}` (`coeff_krtwTerm_eq_zero`, `coeff_melzerTerm_eq_zero`; because
`(r_i+s_i)²+s_i²+s_i ≥ max(r_i,s_i)` and `M_i(M_i+1)/2 ≥ M_i`). Natural-number subtraction is truncated, but
it is only evaluated under the chain conditions, where it is the true difference.
KRTW print the sum over all `r_i, s_i ≥ 0` with the convention `1/(q)_n = 0` for `n < 0`; this is the same sum.

### Sanity test of the definitions

`scripts/Sanity.lean` (run by `lake env lean scripts/Sanity.lean`) is a computable mirror of the definitions
(it calls the very `KRTW.nxt`, `KRTW.prv` of the statement file and copies the formulas of `krtwTerm`,
`krtwSum`, `melzerTerm`, `melzerSum`, `qPochInf`), with `#guard`s: for k = 1, 2, 3 the first 31 coefficients of
the 2k-fold sum, of the k-fold sum and of the product `(q,q^{k+1},q^{k+2};q^{k+2})_∞/((q;q²)_∞(q;q)_∞)` equal
the integers computed independently by `scripts/sanity_independent.py` (plain Python, depth-first enumeration
of the printed sums), e.g. for k = 3: `1, 1, 2, 4, 5, 8, 12, 17, 24, 34, 46, 62, 84, 111, 146, 192, …, 5632`;
the theta form (`thetaSum`) is checked in the same way; two negative controls are detected. (The mirror is linked to the noncomputable definitions by inspection, not
by a proof.)

## 2. The proof, and how it differs from the write-up

All sums live in `K⟦X⟧` with the coefficientwise (q-adic) notion of convergence (Mathlib's `tsum` for the
product topology of the discrete topology on `K`; `KRTW/Summation.lean`); `x` is an arbitrary power series, so
that the substitutions `x ↦ xq^c` of the write-up are just evaluations of a function at another point.

| file | content |
|---|---|
| `KRTW/Basic.lean` | triangular numbers; `(q)_n`, `1/(q)_n` for an arbitrary base `q` with zero constant term; Gaussian binomials as power series, q-Pascal, `[R;r] ≡ 1/(q)_r mod q^{R−r+1}` |
| `KRTW/Summation.lean` | summability from divisibility by powers of `X`; truncation and congruence lemmas; the triangular rearrangement `Σ_b Σ_m F(b, m+σ(b)) = Σ_M Σ_{σ(b)≤M} F(b,M)` |
| `KRTW/Lemma4.lean` | Lemma 4 in the two forms that are used: (∗∗) `star2` and (∗∞) `starInf`; Euler's series `Einf a x = (xq^{a+1};q)_∞` with its functional equation |
| `KRTW/Theorem1.lean` | the chain sums defined one level at a time (`Hf`, `Bf`, `Tf`), their shift identities, and Theorem 1 `thm1 : Tf (j+1) R 0 x = Bf j R x` by induction on `j` |
| `KRTW/Limit.lean` | `R → ∞`: `Kf_eq_Mel : Kf k x = Mel k x` |
| `KRTW/Defs.lean`, `KRTW/Tuples.lean`, `KRTW/Bridge.lean`, `KRTW/Trunc.lean` | the statement-level definitions (all 11 of them are in `Defs.lean`); the recursive sums are the explicit sums over tuples (induction on `k`, peeling off the first resp. last entry, modulo `X^{N+1}`); the truncation lemmas |

Theorem 1, the limit `R → ∞` and the regrouping (3.1) are formalized as in the write-up. **Lemma 4 is proved
by a different route** (found during the formalization; shorter to check by machine): instead of the
q-binomial theorem and Lemmas 2 and 3, three two-term recurrences, each verified termwise:

* with `c(R,n) = Σ_j (−1)^j q^{Rj}/((q)_j(q)_{n−j})`: `c(R+1,n+1) = c(R,n+1) + q^R c(R,n)` and
  `(1−q^{n+1}) c(R,n+1) = c(R+1,n) − q^R c(R,n)`, hence `(1−q^{n+2}) c(0,n+2) = c(0,n)`, i.e. `c(0,2b) = 1/(q²;q²)_b`,
  `c(0,2b+1) = 0` (this is Lemma 2 at `R = 0`);
* both sides of (∗∗) and of (∗∞) satisfy `F(R+1;x) = F(R;x) + x q^{R+1} F(R;xq)` (left sides: q-Pascal; right
  sides: `(xq^{a+1})_∞ = (1−xq^{a+1})(xq^{a+2})_∞`), which reduces both identities to `R = 0`;
* at `R = 0`: with `Z(R,M;x) = Σ_{n≤M} x^n q^{n(n+1)/2} c(R,n)/(q)_{M−n}` (the form (2.3) of the write-up) the right
  side `U` of (∗∗) and `(xq^{R+M+1})_∞ Z` satisfy the same recurrence
  `(1−q^{M+1}) F(R,M+1) = F(R,M) + x q^{M+1} F(R+1,M)`, so `U(R,M;x) = (xq^{R+M+1})_∞ Z(R,M;x)` by induction on `M`,
  and `Z(0,M;x)` is the left side by the first item; for (∗∞) one triangular rearrangement and the first item.

So Lemmas 2 (for `R > 0`) and 3 and the finite form (∗) of the write-up are **not** formalized as stated (they are
not needed); what is checked is that their consequences (∗∗), (∗∞) hold. **No gap or error was found in the
write-up** in the part that was formalized (Theorem 1, the limit, the exponent bookkeeping (0.1), the summability
claims).

### What is NOT machine-checked

* Corollaries 4, 5 and §8 of the write-up (Hall–Littlewood form, the parameter `y`, identities (ii)–(iv), (F_k)).
* The `x`-refined product-side formula (5.1) (second equality of Theorem 2); only its `x = 1` case.
* The statement that the Lean definitions transcribe the printed identity: this is a matter of reading the 11
  definitions (above and in §3) against KRTW's display; the sanity test below supports it numerically.
* Novelty/priority claims.

## 3. Stage 2: the product formula

```lean
-- (q^a; q^d)_∞ and the theta series, by coefficients          (KRTW/Defs.lean, in KRTWStatement.lean)
noncomputable def qPochInf (a d : ℕ) : K⟦X⟧ :=
  PowerSeries.mk fun N => coeff N (∏ i ∈ range (N + 1), (1 - X ^ (a + d * i) : K⟦X⟧))
noncomputable def thetaSum (k : ℕ) : K⟦X⟧ :=
  PowerSeries.mk fun N => ∑ n ∈ range (N + 1),
    coeff N ((-1 : K⟦X⟧) ^ n * X ^ (((k + 2) * n ^ 2 + k * n) / 2) * (1 - X ^ (2 * n + 1)))

theorem krtw_identity (k : ℕ) :                     -- the printed identity
    krtwSum (1 : ℚ⟦X⟧) k
      = qPochInf 1 (k + 2) * qPochInf (k + 1) (k + 2) * qPochInf (k + 2) (k + 2)
        * (qPochInf 1 2 * qPochInf 1 1)⁻¹
theorem krtw_identity_theta (k : ℕ) :               -- output of the Bailey chain
    krtwSum (1 : ℚ⟦X⟧) k * (qPochInf 1 2 * qPochInf 1 1) = thetaSum k
theorem theta_eq_triple_product (k : ℕ) :           -- the instance of the Jacobi triple product
    (thetaSum k : ℚ⟦X⟧) = qPochInf 1 (k + 2) * qPochInf (k + 1) (k + 2) * qPochInf (k + 2) (k + 2)
theorem melzer_identity (k : ℕ) :                   -- the known evaluation of the k-fold sum
    melzerSum (1 : ℚ⟦X⟧) k
      = qPochInf 1 (k + 2) * qPochInf (k + 1) (k + 2) * qPochInf (k + 2) (k + 2)
        * (qPochInf 1 2 * qPochInf 1 1)⁻¹
theorem coeff_qPochInf {K : Type*} [Field K] (a d : ℕ) (ha : 1 ≤ a) (hd : 1 ≤ d) (N M : ℕ) (h : N ≤ M) :
    coeff N (qPochInf a d : K⟦X⟧) = coeff N (∏ i ∈ range (M + 1), (1 - X ^ (a + d * i) : K⟦X⟧))
```

Reading guide. `qPochInf a d` is `(q^a;q^d)_∞ = ∏_{i≥0}(1 − q^{a+di})` for `a, d ≥ 1`: its coefficient of `q^N`
is that of the finite product over `i ≤ N`, and `coeff_qPochInf` says that every longer finite product has the
same coefficient of `q^N`. So the right side of `krtw_identity` is
`(q;q^{k+2})_∞ (q^{k+1};q^{k+2})_∞ (q^{k+2};q^{k+2})_∞ · ((q;q²)_∞ (q;q)_∞)⁻¹`, where `⁻¹` is the inverse in
`ℚ⟦X⟧` of a series with constant term 1. `thetaSum k` is `Σ_{n∈ℤ} (−1)^n q^{((k+2)n²+kn)/2}` with the terms
`n ≥ 0` and `−n−1` paired: `Σ_{n≥0} (−1)^n q^{((k+2)n²+kn)/2}(1 − q^{2n+1})` (the exponent at `−n−1` is the
exponent at `n` plus `2n+1`); this pairing is done by hand in the definition, the genuinely bilateral form is
not formalized. All four statements are also proved over an arbitrary field in place of `ℚ`
(`krtw_product_field`, `krtw_theta_field` in `KRTWSolution.lean`; `jtp_theta`).

| file | content |
|---|---|
| `KRTW/Bailey.lean` | `gp y n = (y;q)_n` for arbitrary `y`; q-Chu–Vandermonde `cv` (induction on `N`, all `b`); the unit Bailey pair `ubp` (telescoping); Bailey's lemma (S2) `bailey_step` and the chain `bailey_chain`, relative to `a = x²q` for an arbitrary `x`, **cleared of the denominators** `(aq;q)_m`, `(−xq;q)_m`: `(x²q²;q)_{2n} G_i(n) = (−xq;q)_n Σ_{r≤n} x^{(i+1)r} q^{(i+1)r(r+1)/2} α_r (x²q^{n+r+2};q)_{n−r}/(q)_{n−r}` with `G_i(n) = Σ_{n≥M_1≥…≥M_i≥0} ∏ x^{M_l}q^{M_l(M_l+1)/2}/((q)_{n−M_1}⋯(q)_{M_i})` |
| `KRTW/JTP.lean` | homogeneous q-binomial theorem `qbinomial`; finite triple product `finite_tp`; its limit; `jtp_theta` |
| `KRTW/Stage2.lean` | `(q)_n G_k(n) ≡ Mel_k(x) mod q^{n+1}` (`poch_Gc_cong`, via the tuple form); at `x = 1` the limit `limit_one : (q)_n Mel_k(1) ≡ (−q;q)_n Σ_{r≤n}(−1)^r q^{((k+2)r²+kr)/2}(1−q^{2r+1}) mod q^{n+1}`; `(−q;q)_n (q;q²)_n = (q^{n+1};q)_n`; `mel_theta`, `krtw_theta'`, `krtw_product'` |

Differences from §5 of the write-up: (a) the Bailey pair relation is multiplied by `(aq;q)_{2n}` and the
prefactor `1/(−xq;q)_n` of `β^{(i)}_n` is kept on the other side, so that no inverse of a Pochhammer symbol
with a generic argument is needed; Bailey's lemma is derived from q-Chu–Vandermonde (which is proved), not quoted;
(b) the chain is formalized for arbitrary `x`, but **the limit `n → ∞` only at `x = 1`**: the `x`-refined closed
form (5.1) of the write-up (second equality of Theorem 2) is NOT formalized; (c) the Jacobi triple product is
proved only in the instance `(q, q^{k+1}, q^{k+2}; q^{k+2})_∞` (for all `k`), from the q-binomial theorem.
No gap or error was found in §5 either.

## 4. Pins and how to check

Lean `leanprover/lean4:v4.33.0-rc2`, Mathlib `51e6992efd06126df61a496bebf8f49482a4e129` (the pins of
`../rr-identities/lean`).

    curl -sSfL https://raw.githubusercontent.com/leanprover/elan/master/elan-init.sh | sh -s -- -y --default-toolchain none
    export PATH="$HOME/.elan/bin:$PATH"
    cd lean/lean-krtw
    lake exe cache get          # Mathlib binaries; never build Mathlib from source
    lake build                  # the warnings "declaration uses sorry" come from KRTWStatement.lean only (deliberate)
    lake env lean scripts/PrintAxioms.lean
    ./scripts/compare_statements.sh
    lake env lean scripts/Sanity.lean        # about 20 s
    python3 scripts/sanity_independent.py

The outputs of these commands are saved in `logs/` (`print_axioms.txt`, `compare_statements.txt`, `sanity.txt`,
`sanity_independent.txt`). The whole library is about 2700 lines and builds in a few minutes.

`scripts/compare_statements.sh` checks that `KRTWStatement.lean` imports only Mathlib, that its definition block
is byte for byte the block of `KRTW/Defs.lean`, that each statement (docstring included) occurs byte for byte in
`KRTWSolution.lean`, that the solution side contains no `sorry` / `axiom` / `native_decide`, and that `#check`
/ `#print` give the same output from both modules. It is a weak substitute for `leanprover/comparator`, which is
not installed here; `comparator.json` is the configuration for that tool (9 theorems, 11 definitions, the three
permitted axioms). `KRTWStatement.lean` is generated by `scripts/gen_statement.py`.
