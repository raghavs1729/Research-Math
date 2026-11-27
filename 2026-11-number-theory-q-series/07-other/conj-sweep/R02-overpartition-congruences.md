# R02. Four conjectural Ramanujan-type congruences for overpartitions (mod 7, 17, 19, 23): PROVED

**Source.** XuanLing Wei, *New Ramanujan-type Congruences for Overpartitions Modulo 11 and 13*,
arXiv:2603.08510 (9 March 2026), Section "Closing remarks". `p̄(n)` is the number of
overpartitions of `n`, `Σ p̄(n) qⁿ = ∏ (1+qⁿ)/(1−qⁿ)`.

**Statements as printed.**

* (closing remarks, first display) "Numerical computations suggest that the overpartition function
  `p̄(n)` satisfies the following conjectural Ramanujan-type congruences:
  `p̄(17 × 11² n) ≡ 0 (mod 17)`, where `n ≡ 3 (mod 8)` and `(n/11) = −1`. Equivalently,
  `p̄(17 × 11²(88n+m)) ≡ 0 (mod 17)` for all `n ≥ 0` and `m ∈ {19,35,43,51,83}`."
* (second display) "`p̄(23 × 13² n) ≡ 0 (mod 23)`, where `n ≡ 5 (mod 8)` and `(n/13) = 1`.
  Equivalently, `p̄(23 × 13²(104n+l)) ≡ 0 (mod 23)`, for all `n ≥ 0` and `l ∈ {29,53,61,69,77,101}`."
  The author writes that these two "are expected to be provable by the methods developed in this
  paper", reduces them to checking `n ≤ 1,226,649,601` resp. `n ≤ 3,968,520,193` coefficients, and
  states: "This verification requires a substantial amount of computation and has not been completed."
* **Conjecture 1 (printed).** "`p̄(2⁴ n) ≡ 0 (mod 7)`, where `n ≡ 3 (mod 8)` and `(n/7) = 1`.
  Equivalently, `p̄(2⁴(56n+k)) ≡ 0 (mod 7)` for all `n ≥ 0` and `k ∈ {11,43,51}`."
* **Conjecture 2 (printed).** "`p̄(17²·4n) ≡ 0 (mod 19)`, where `n ≡ 3 (mod 8)`, `(n/17) = −1`, and
  `(n/19) = 1`. It involves 72 residue classes modulo 2584".
  Both conjectures are introduced by "conjectural Ramanujan-type congruences for overpartitions
  which cannot be proved by the methods in this paper".

## Results

**Theorem A.** `p̄(17·11²·n) ≡ 0 (mod 17)` for all `n ≡ 3 (mod 8)` with `(n/11) = −1`.

**Theorem B.** `p̄(23·13²·n) ≡ 0 (mod 23)` for all `n ≡ 5 (mod 8)` with `(n/13) = +1`.

**Theorem C** (= Conjecture 1). `p̄(16n) ≡ 0 (mod 7)` for all `n ≡ 3 (mod 8)` with `(n/7) = +1`.

**Theorem D** (= Conjecture 2). `p̄(4·17²·n) ≡ 0 (mod 19)` for all `n ≡ 3 (mod 8)` with
`(n/17) = −1` and `(n/19) = +1`.

Each proof is a reduction, by standard theorems on modular forms, to a finite computation in exact
arithmetic, which was carried out (certificates in `cert/`). Nothing below is "numerical evidence":
the finite checks cover the full range required by the cited bounds.

## Notation and the facts used

`φ(q) = Σ_{n∈ℤ} q^{n²}`, `F(q) = Σ_{n odd ≥ 1} σ(n) qⁿ = η(4τ)⁸/η(2τ)⁴`, `q = e^{2πiτ}`.

* (F1) `Σ p̄(n)(−q)ⁿ = 1/φ(q)`. (Because `Σ p̄(n)qⁿ = (−q;q)_∞/(q;q)_∞ = 1/φ(−q)`, Gauss.)
* (F2) For every integer `K ≥ 0`, `M_{K/2}(Γ̃₀(4))` (Shimura's space of weight `K/2`; for even `K`
  this is `M_{K/2}(Γ₀(4), χ_{−4}^{K/2})`) has the basis `B_j = F^j φ^{K−4j}`, `0 ≤ j ≤ ⌊K/4⌋`
  (Cohen; Koblitz, *Introduction to Elliptic Curves and Modular Forms*, Ch. IV, Prop. 4). Since
  `B_j = q^j + O(q^{j+1})` with integer coefficients, the forms with integer coefficients are exactly
  `Λ_{K/2} := ⊕ ℤ B_j`, an element of `Λ_{K/2}` (or of `Λ_{K/2}⊗𝔽_ℓ`) is determined by its first
  `⌊K/4⌋+1` coefficients, and `Λ_{K/2}⊗𝔽_ℓ → 𝔽_ℓ[[q]]` is injective.
* (F3) Hecke operators of half-integral weight (Shimura, Ann. of Math. 97 (1973), Thm. 1.7; Koblitz
  IV, Prop. 13): for an odd prime `p`, `λ ≥ 1` and `f = Σ a(n)qⁿ ∈ M_{λ+1/2}(Γ̃₀(4))`,
  `T(p²)f = Σ b(n) qⁿ ∈ M_{λ+1/2}(Γ̃₀(4))` with
  `b(n) = a(p²n) + (−1/p)^λ (n/p) p^{λ−1} a(n) + p^{2λ−1} a(n/p²)`.
  Hence `T(p²)` preserves `Λ_{λ+1/2}` and acts on `Λ_{λ+1/2}⊗𝔽_ℓ` by the same formula.
* (F4) `U_4 : Σ a(n)qⁿ ↦ Σ a(4n)qⁿ` maps `M_{K/2}(Γ̃₀(4))` to itself, hence preserves `Λ_{K/2}`.
  (This is the operator `U_4` of Kohnen's theory on level 4. The usual citation, Shimura 1973,
  Prop. 1.5 / Ono, *Web of Modularity*, Prop. 3.7(2), is stated for `U(d)` with `d | N` on level `4N`
  and does not literally cover level 4, as the referee pointed out (`AUDIT-R02.md`, G1); so here is
  a direct proof, due to the referee and re-checked.)

  *Proof.* `U_4f(τ) = ¼ Σ_{j mod 4} f((τ+j)/4)`. Let `γ = (a b; c d) ∈ Γ₀(4)` and `σ_j = (1 j; 0 4)`.
  Since `4 | c` and `a` is odd, `j ↦ j′ := a^{−1}(b+jd) mod 4` is a bijection of `ℤ/4`, and
  `σ_j γ = γ′ σ_{j′}` with `γ′ = (a+jc, (b+jd−j′(a+jc))/4; 4c, d−cj′) ∈ Γ₀(16)`. With `τ′ = σ_{j′}τ`:
  `4cτ′ + (d−cj′) = cτ+d`; `d′ := d−cj′ ≡ d (mod 4)`, so `ε_{d′} = ε_d`; and
  `(4c/d′) = (c/d′) = (c/d)` because for `4 | c`, `c ≠ 0`, the symbol `d ↦ (c/d)` (Shimura's
  convention) on odd `d` is periodic modulo `|c|` (for `c = 0`, `d = d′ = ±1`). Hence
  `j(γ′,τ′) = j(γ,τ)` and `U_4f(γτ) = ¼Σ_j f(γ′τ′_j) = j(γ,τ)^K·¼Σ_{j′} f((τ+j′)/4) = j(γ,τ)^K U_4f(τ)`.
  Each `f((τ+j)/4)` is holomorphic at all cusps. ∎
  (Numerical confirmation by the referee: `audit/f4test.gp`, all basis elements for `K ≤ 48`.)
* (F5) Integer-weight Hecke operator: for a prime `ℓ ∤ 4` and `f = Σ a(n)qⁿ ∈ M_w(Γ₀(4), χ)`,
  `T_ℓ f = Σ (a(ℓn) + χ(ℓ)ℓ^{w−1}a(n/ℓ)) qⁿ ∈ M_w(Γ₀(4), χ)`; so `U_ℓ f ≡ T_ℓ f (mod ℓ)` if `w ≥ 2`.
* (F6) Sturm's bound (Sturm, LNM 1240 (1987); Stein, *Modular Forms: a Computational Approach*,
  Thm. 9.18): if `G ∈ M_w(Γ₁(N))` has integer coefficients `a(n)` with `a(n) ≡ 0 (mod ℓ)` for all
  `n ≤ w·[SL₂(ℤ):Γ₁(N)]/12`, then `G ≡ 0 (mod ℓ)`. For `N = 256` the index is `256²·(3/4) = 49152`.
* (F7) **Sieving lemma.** Let `g ∈ Λ_{λ+1/2}` and `r mod 8`. If `g(n) ≡ 0 (mod ℓ)` for all
  `n ≡ r (mod 8)` with `n ≤ (λ+1)·49152/12`, then `g(n) ≡ 0 (mod ℓ)` for all `n ≡ r (mod 8)`.

  *Proof.* `G(τ) := g(τ)φ(8τ)` has integer weight `λ+1`. For `γ = (a b; c d) ∈ Γ₀(32)`,
  `g(γτ) = j(γ,τ)^{2λ+1} g(τ)` with `j(γ,τ) = (c/d) ε_d^{−1} (cτ+d)^{1/2}`, and
  `φ(8γτ) = j(γ′, 8τ) φ(8τ)` with `γ′ = (a 8b; c/8 d) ∈ Γ₀(4)`, `j(γ′,8τ) = (2/d)(c/d)ε_d^{−1}(cτ+d)^{1/2}`.
  So `G(γτ) = (2/d)(−1/d)^{λ+1}(cτ+d)^{λ+1} G(τ)`; in particular `G ∈ M_{λ+1}(Γ₁(32))`.
  Put `α_j = (1 j/8; 0 1)`. Then `S(τ) := Σ_{n≡r (8)} [qⁿ]G · qⁿ = (1/8) Σ_{j=0}^{7} e^{−2πi rj/8} G|α_j`,
  and for `γ ∈ Γ₁(256)` the matrix
  `α_j γ α_j^{−1} = (a + jc/8, b + j(d−a)/8 − j²c/64; c, d − jc/8)` is integral and lies in `Γ₁(32)`.
  Hence `S ∈ M_{λ+1}(Γ₁(256))`, with integer coefficients. Since `φ(q⁸)` is a series in `q⁸`,
  `S = (Σ_{n≡r (8)} g(n)qⁿ)·φ(q⁸)`. The hypothesis gives `[qⁿ]S ≡ 0` for `n ≤ (λ+1)·49152/12`, so
  `S ≡ 0 (mod ℓ)` by (F6), and `φ(q⁸) = 1 + 2q⁸ + …` is invertible in `𝔽_ℓ[[q]]`. ∎
* (F8) **Θ-lemma.** Let `ℓ ≥ 5` be prime and `Θ = q d/dq`. Then `Θ` maps `Λ_k⊗𝔽_ℓ` into
  `Λ_{k+ℓ+1}⊗𝔽_ℓ` (as subsets of `𝔽_ℓ[[q]]`), and `Λ_k⊗𝔽_ℓ ⊂ Λ_{k+ℓ−1}⊗𝔽_ℓ`.

  *Proof.* From `η`-product formulas, `ΘF/F = (4E₂(q⁴) − E₂(q²))/3` and
  `Θφ/φ = (10E₂(q²) − 2E₂(q) − 8E₂(q⁴))/24` (`φ = η(2τ)⁵/(η(τ)²η(4τ)²)`, `Θη(dτ)/η(dτ) = dE₂(q^d)/24`).
  For `ℓ ≥ 5`, `E₂ ≡ E_{ℓ+1}` and `E_{ℓ−1} ≡ 1 (mod ℓ)` (Kummer, von Staudt–Clausen; Swinnerton-Dyer),
  and `E_{ℓ+1}(q^d) ∈ M_{ℓ+1}(Γ₀(4))` for `d | 4` has `ℓ`-integral coefficients, so it lies in
  `Λ_{ℓ+1}⊗ℤ_{(ℓ)}` by (F2). Hence `Θ(F^jφ^m) ≡ F^jφ^m·H_{j,m}` with `H_{j,m} ∈ Λ_{ℓ+1}⊗𝔽_ℓ`. The
  second claim is multiplication by `E_{ℓ−1} ∈ Λ_{ℓ−1}⊗ℤ_{(ℓ)}`. ∎

## Proof of Theorems A and B

Let `(ℓ, p, r, s) = (17, 11, 3, −1)` for A and `(23, 13, 5, +1)` for B; `λ = (ℓ−3)/2`.

1. By (F1), `1/φ(q) = φ(q)^{ℓ−1}/φ(q)^ℓ ≡ φ(q)^{ℓ−1}/φ(q^ℓ) (mod ℓ)`, so
   `Σ p̄(ℓn)(−1)ⁿqⁿ = U_ℓ(1/φ) ≡ U_ℓ(φ^{ℓ−1})/φ (mod ℓ)` (`ℓ` odd).
   `φ^{ℓ−1} ∈ M_{(ℓ−1)/2}(Γ₀(4), χ_{−4}^{(ℓ−1)/2})`, so by (F5) and (F2)
   `U_ℓ(φ^{ℓ−1}) ≡ T_ℓ(φ^{ℓ−1}) = Σ_j α_j F^j φ^{ℓ−1−4j}` with integers `α_j` determined by the first
   `⌊(ℓ−1)/4⌋+1` coefficients. Computed (exactly):
   `α mod 17 = (1, 13, 13, 0, 0)`; `α mod 23 = (1, 9, 5, 14, 17, 20)`.
   For `ℓ = 17` the coefficient of `F⁴` (the only basis element not divisible by `φ`) vanishes mod 17.
   Hence, with `c₀ := Σ_j γ_j F^j φ^{ℓ−2−4j} ∈ Λ_{λ+1/2}`,
   `γ = (1,13,13,0)` resp. `(1,9,5,14,17,20)`:  `p̄(ℓn)(−1)ⁿ ≡ c₀(n) (mod ℓ)` for all `n`.
   (This agrees with the two displayed formulas for `c`, `d` in the source.)
2. By (F3), `b := T(p²)c₀ ∈ Λ_{λ+1/2}`; its coordinates `β` follow from `b(0),…,b(⌊(ℓ−2)/4⌋)`:
   `β mod 17 = (8, 14, 14, 0)` (`λ = 7`, `(−1/11)⁷ = −1`);
   `β mod 23 = (3, 20, 8, 10, 12, 6)` (`λ = 10`, `(−1/13)^{10} = +1`).
3. For `p ∤ n` with `(n/p) = s`: `c₀(p²n) = b(n) − s(−1/p)^λ p^{λ−1} c₀(n)`. Here
   `s(−1/p)^λ = +1` in both cases, `11⁶ ≡ 8 (mod 17)`, `13⁹ ≡ 3 (mod 23)`. So `c₀(p²n) ≡ h(n)` with
   `h := Σ η_j F^jφ^{ℓ−2−4j}`, `η = β − p^{λ−1}γ`: `η = (0,12,12,0)` mod 17, `η = (0,16,16,14,7,15)` mod 23.
4. By the exact computation, `h(n) ≡ 0 (mod ℓ)` for every `n ≡ r (mod 8)` with
   `n ≤ (λ+1)·49152/12` (= 32768 for A, 45056 for B). By (F7), `h(n) ≡ 0` for all `n ≡ r (mod 8)`.
5. So for `n ≡ r (mod 8)`, `(n/p) = s`: `p̄(ℓp²n) ≡ ±c₀(p²n) ≡ ±h(n) ≡ 0 (mod ℓ)`. ∎

*Remark.* For A one can avoid step 4's computation: `h ≡ 12·Fφ⁷(φ⁴+F) ≡ 12·Fφ(q)⁷φ(−q)⁴ (mod 17)`
because `φ(q)⁴ − φ(−q)⁴ = 16F`, and the coefficients of `Fφ(q)⁷φ(−q)⁴` at `n ≡ 3 (mod 8)` vanish
exactly (checked to `n < 400`; not used).

Certificate: `cert/overpartition_hecke_cert.txt` (program `code/overp_hecke_proof.py`).

## Proof of Theorem D

`ℓ = 19`, `p = 17`. Let `W(q) := Σ_{(N/19)=1} p̄(N)(−q)^N`.

1. `1/φ(q) ≡ φ(q)^{18}/φ(q^{19}) (mod 19)`, and `φ(q^{19})` is a series in `q^{19}`; taking the terms
   whose exponent is a nonzero square mod 19 gives `W(q)·φ(q^{19}) ≡ A₊(q) (mod 19)`, where
   `A = φ^{18} = Σ r₁₈(n)qⁿ` and `A₊ := Σ_{(n/19)=1} r₁₈(n)qⁿ`.
2. `A₊ ≡ ½(Θ^{18}A + Θ⁹A) (mod 19)` since `n^{18} ≡ [19∤n]`, `n⁹ ≡ (n/19)`. `A ∈ Λ₉`; by (F8),
   `Θ^{18}A ∈ Λ_{369}⊗𝔽₁₉` and `Θ⁹A ∈ Λ_{189}⊗𝔽₁₉ ⊂ Λ_{369}⊗𝔽₁₉` (`369 − 189 = 10·18`). So
   `A₊ mod 19 ∈ Λ_{369}⊗𝔽₁₉`, with coordinates `u_0,…,u_{184}` in the basis `F^jφ^{738−4j}` determined by
   its first 185 coefficients. Computed: `u_{180} = … = u_{184} = 0`. Therefore
   `A₊ ≡ φ^{19}·w₀` with `w₀ := Σ_{j<180} u_j F^jφ^{719−4j} ∈ Λ_{719/2}⊗𝔽₁₉`, and since
   `φ(q^{19}) ≡ φ^{19}` is invertible, **`W ≡ w₀ (mod 19)` as power series.**
3. `b := T(17²)w₀ ∈ Λ_{719/2}⊗𝔽₁₉` by (F3) (`λ = 359`, `(−1/17)^{359} = 1`, `17^{358} ≡ 5`); its
   coordinates follow from `b(n)`, `n < 180`, i.e. from `w₀(m) = W(m)`, `m ≤ 289·179`.
4. For `(n/17) = −1`: `w₀(289n) = b(n) + 17^{358}w₀(n) = h(n)`, `h := b + 5w₀ ∈ Λ_{719/2}⊗𝔽₁₉`.
5. `g := U₄h ∈ Λ_{719/2}⊗𝔽₁₉` by (F4); coordinates from `h(4n)`, `n < 180`.
6. Computed: `g(n) ≡ 0 (mod 19)` for all `n ≡ 3 (mod 8)`, `n ≤ 360·49152/12 = 1 474 560`. By (F7)
   (applied to an integral lift of `g`), `g(n) ≡ 0` for all `n ≡ 3 (mod 8)`.
7. Let `n ≡ 3 (mod 8)`, `(n/17) = −1`, `(n/19) = 1`. Then `N := 4·289·n` has `(N/19) = 1`, so
   `±p̄(N) = W(N) ≡ w₀(289·4n) = h(4n) = g(n) ≡ 0 (mod 19)` (step 4 with `4n`, `(4n/17) = −1`). ∎

Certificate: `cert/overpartition_mod19_cert.txt` (program `code/overp19_proof.py`, 138 s).

## Proof of Theorem C

`Σ p̄(n)qⁿ = f₂/f₁²` (`f_k = (q^k;q^k)_∞`) is an eta-quotient, and the three progressions
`896n + 16k`, `k ∈ {11,43,51}`, form exactly the set `P_{m,r}(t)` of Radu's algorithm for `m = 896`,
`t = 176`. Radu's Lemma 4.5 (Ramanujan J. 20 (2009); in the form of Radu–Sellers, Int. J. Number
Theory 7 (2011), Lemma 2.4; statement reproduced in `code/radu.py`) applies with `N = 56`
(`(m, M, N, t, r) ∈ Δ*` verified) and an explicit integer vector `a` on the divisors of 56 for which
`p_{m,r}(γ) + p*_a(γ) ≥ 0` at all cusps (verified exactly), giving `ν = 88009/28`, `⌊ν⌋ = 3143`.
The congruence `p̄(896n + t′) ≡ 0 (mod 7)` was verified for the three `t′` and all `0 ≤ n ≤ 3143`. ∎

Certificate: `cert/radu_overpartition_mod7.txt` (program `code/radu_overp.py`; engine sanity tests
in `cert/radu_sanity.txt`).

## Independent numerical test (not part of the proofs)
`code/overp_direct.py`: direct computation of `p̄ mod ℓ` up to `8·10⁶`: the four congruences hold for
all 220, 119, 26 785 and 196 admissible `n` in range.

## Status, audit, novelty
Independent adversarial referee: `AUDIT-R02.md` (verdicts: A, B, C proof complete; D complete
modulo the citation for (F4), repaired above by a direct proof; every computed number reproduced by
independently written PARI/GP and Python programs over the full Sturm ranges). Remark R1 of the
audit: the sieved form in (F7) already lives on `Γ₀(64) ∩ Γ₁(8)`, so the ranges checked here are
larger than necessary.

**Novelty: NOT NEW as statements.** `NOVELTY-1.md`, item 1: Sumit Kumar Jha, *On Four Conjectural
Congruence Families for the Overpartition Function*, Research Square preprint
DOI 10.21203/rs.3.rs-11336379/v1 (posted 26 Oct 2026, unrefereed, under review at The Ramanujan
Journal, with an AI-assistance note) claims all four congruences, and more general Hecke relations
modulo 17 and 23 (via Cohen–Eisenstein series). That preprint was not read here beyond its first
pages by the search agent, and it was found only after the proofs above were complete; the proofs
above are independent of it and come with exact certificates and an independent audit. The older
literature (Ryan–Sirolli–Villegas-Morales–Zheng, Ramanujan J. 2024, and earlier) does not contain
these congruences.
