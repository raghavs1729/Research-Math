# AUDIT-R02. Adversarial referee report on `R02-overpartition-congruences.md`

Scope: the four claimed theorems A, B, C, D (overpartition congruences mod 17, 23, 7, 19 conjectured in
arXiv:2603.08510). I read the note and the session's programs, checked every "fact used" and every numbered step by
hand, and recomputed every finite computation with programs written from scratch (PARI/GP, and Python with
python-flint used only as a polynomial-arithmetic library), using different formulas and algorithms from the
session's code wherever that was possible. Programs: `audit/`; their complete output: `audit/OUTPUT.txt`.

## Verdicts

| Theorem | Verdict |
|---|---|
| A (mod 17) | **PROOF COMPLETE.** No error found. All numbers reproduced. |
| B (mod 23) | **PROOF COMPLETE.** No error found. All numbers reproduced. |
| C (mod 7)  | **PROOF COMPLETE** (as an application of the cited Radu / Radu–Sellers lemma, whose statement `radu.py` implements faithfully). All numbers reproduced. |
| D (mod 19) | **COMPLETE MODULO ONE LISTED GAP** (G1: the *citation* for fact (F4)). The fact itself is true; a self-contained proof is given below, and it was also confirmed numerically. With the citation replaced by that proof, D is complete. All numbers reproduced over the full Sturm range. |

I found **no mathematical error** and **no discrepancy in any computed number**. The statements proved are exactly
the statements printed in the source.

## Gaps and remarks (complete list)

* **G1 (Theorem D, step 5; fact F4) — citation not established, statement true.** (F4) says `U_4` maps
  `M_{K/2}(Γ̃₀(4))` to itself and cites "Shimura, Prop. 1.5 with m = 4". I could not consult Shimura's paper here.
  In the form in which that proposition is normally quoted (I cite it from memory: Ono, *The Web of Modularity*, Prop. 3.7(2):
  `f ∈ M_{λ+1/2}(Γ₀(4N), ψ)`, `d | N` ⟹ `f|U(d) ∈ M_{λ+1/2}(Γ₀(4N), (4d/·)ψ)`) the hypothesis is `d | N` with the
  level written as `4N`; for level 4 (`N = 1`) and `d = 4` that hypothesis is **not** satisfied. So the citation, as
  given, does not visibly cover the case used, and the note should not rest on it. The statement is nevertheless
  true (it is the operator `U_4` of Kohnen's theory on level 4). Direct proof:

  *Proof of (F4).* `U_4 f(τ) = ¼ Σ_{j mod 4} f((τ+j)/4)`. Let `γ = (a b; c d) ∈ Γ₀(4)`. For each `j` choose `j′`
  with `j′(a+jc) ≡ b+jd (mod 4)`; since `4 | c` and `a` is odd this is `j′ ≡ a⁻¹b + j (mod 4)` (as `ad ≡ 1`,
  `a² ≡ 1 mod 4`), a bijection of `ℤ/4`. Then
  `(1 j; 0 4)·γ = γ′·(1 j′; 0 4)` with `γ′ = (a+jc, (b+jd−j′(a+jc))/4; 4c, d−cj′) ∈ Γ₀(16)`.
  With `τ′ = (τ+j′)/4`: `4cτ′ + (d−cj′) = cτ+d`, `d′ := d−cj′ ≡ d (mod 4)` so `ε_{d′} = ε_d`, and
  `(4c/d′) = (c/d′) = (c/d)`, because for `4 | c`, `c ≠ 0`, the map `d ↦ (c/d)` on odd integers (Shimura's
  convention, `(c/−1) = sgn c`) is periodic modulo `|c|` (its conductor is that of `ℚ(√c)`, which divides `c` when
  `4 | c`; checked by brute force for `|c| ≤ 400`, `|d| ≤ 801` in `audit/misc.gp`), and for `c = 0` there is
  nothing to prove. Hence `j(γ′, τ′) = j(γ, τ)` and
  `U_4f(γτ) = ¼ Σ_j f(γ′τ′_j) = j(γ,τ)^K · ¼ Σ_{j′} f((τ+j′)/4) = j(γ,τ)^K U_4f(τ)`. Holomorphy at the cusps is clear
  since each `f((τ+j)/4)` is a holomorphic modular form on a congruence subgroup. ∎

  Numerical confirmation: `audit/f4test.gp` (exact rational arithmetic): for all `K ≤ 48` and all 324 basis
  elements `B_j = F^jφ^{K−4j}`, `U_4B_j` lies in the span of the basis with *integer* coordinates, verified on
  `dim + 40` coefficients; negative control `U_2`: in the span for all 168 basis elements of even `K` (integer
  weight) and for **none** of the 156 of odd `K`. In the actual proof: `(U_4h)/φ^{719}` is a polynomial of degree
  `≤ 179` in `t = F/φ⁴` up to `O(t^{700})` (`audit/d.gp`), and `g(n) = h(4n)` was confirmed directly from
  overpartition numbers for all `n ≤ 100000` (`audit/direct_cd.py`).

* **R1 (not a gap; the bound is conservative).** In (F7) the sieved form `S` is in fact on `Γ₀(64) ∩ Γ₁(8) ⊇ Γ₁(64)`
  (the conjugation only needs `64 | c` and `a ≡ d ≡ 1 mod 8`); the note uses `Γ₁(256)` and the index in `SL₂(ℤ)`
  (49152, not the projective index 24576). Both choices only make the checked range larger than necessary.
* **R2 (program, not proof).** In `code/overp_hecke_proof.py` the "extra consistency" assertion after step 2 is
  vacuous: `D2 = p²·dim2 + 5` gives `nchk = dim2`, so it re-checks only the coefficients that were used to solve
  for `β`. This does not affect the proof (membership of `T(p²)c₀` in `Λ` is (F3)). My program checks the expansion
  of `T(p²)c₀` on 25 further coefficients and of `T_ℓ(φ^{ℓ−1})` on 40 further coefficients (all agree).
* **R3 (cosmetic).** In the proof of C the tuple is written `(m, M, N, t, r)`; the lemma writes
  `(m, M, N, (r_δ), t)`. The display for `P_{m,r}(t)` in arXiv:2302.11830 omits "`∈ 𝕊_{24m}`" (it only says
  "∃ [s]_{24m}"), although `𝕊_{24m}` (squares in `ℤ*_{24m}`) is introduced in the sentence before; `radu.py` uses
  the squares, which is the definition of Radu–Sellers (and the only reading for which `(s−1)/24` is an integer).
* **R4 (references I could not open).** Shimura (Thm 1.7), Koblitz (IV, Prop. 4 and Prop. 13), Stein (Thm 9.18),
  Radu–Sellers "Lemma 2.4": statements checked from memory / from the verbatim copy in arXiv:2302.11830, numbering
  not checked against the originals. (F3) was additionally tested numerically (`T(9)` on all basis elements of odd
  `K ≤ 47`; `T(11²)`, `T(13²)`, `T(17²)` on the forms of the proofs, each on more coefficients than the dimension).

## 1. Mathematical audit, item by item

**(F1)** `(−q;q)_∞/(q;q)_∞ = 1/φ(−q)` (Gauss), so `Σ p̄(n)(−q)ⁿ = 1/φ(q)`. Correct.

**(F2)** Koblitz IV Prop. 4: `M_{K/2}(Γ̃₀(4))` = polynomials in `Θ` (weight ½) and `F` (weight 2) of pure weight
`K/2`. `B_j = q^j + …` with integer coefficients, hence triangularity, the description of the integral forms, the
determination by the first `⌊K/4⌋+1` coefficients and injectivity mod `ℓ` are all correct. `F = η(4τ)⁸/η(2τ)⁴`,
`φ = η(2τ)⁵/(η(τ)²η(4τ)²)` verified against `Σ_{n odd}σ(n)qⁿ`, `Σ q^{n²}` (n < 3000).

**(a) (F3).** Shimura Thm 1.7 / Koblitz IV Prop. 13 with trivial character, level 4, `p` odd:
`b(n) = a(p²n) + (−1/p)^λ (n/p) p^{λ−1} a(n) + p^{2λ−1} a(n/p²)` — sign `((−1)^λ n / p) = (−1/p)^λ (n/p)`, exponents
`λ−1` and `k−2 = 2λ−1`: correct. Level 4 is allowed (`4 | N`, `p ∤ N`). `λ ≥ 1` in all uses (7, 10, 359), so the
formula has integer coefficients, `T(p²)` maps the lattice `Λ_{λ+1/2}` to itself (integral coefficients + (F2)),
hence induces an `𝔽_ℓ`-linear map of `Λ⊗𝔽_ℓ` given by the same coefficient formula (lift, apply, reduce;
well defined because `ℓΛ` is mapped into `ℓΛ`). Legitimate. Signs used: `(−1/11)⁷ = −1`, `(−1/13)^{10} = +1`,
`(−1/17)^{359} = +1`: correct.

**(b) (F7).** Checked line by line.
* `g(γτ) = j(γ,τ)^{2λ+1}g(τ)` on `Γ₀(4)`; for `γ ∈ Γ₀(32)`, `8γτ = γ′(8τ)`, `γ′ = (a, 8b; c/8, d) ∈ Γ₀(4)`,
  `j(γ′,8τ) = ((c/8)/d) ε_d⁻¹ (cτ+d)^{1/2}` and `((c/8)/d) = (c/d)(2/d)` (also for `d < 0`, `c = 0`). Correct.
* `G(γτ) = (2/d) [ (c/d)² ε_d⁻² (cτ+d) ]^{λ+1} G = (2/d)(−1/d)^{λ+1}(cτ+d)^{λ+1}G` (`ε_d² = (−1/d)`). Correct; the
  character is trivial on `Γ₁(8)`, so `G ∈ M_{λ+1}(Γ₀(32)∩Γ₁(8)) ⊂ M_{λ+1}(Γ₁(32))`.
* `α_jγα_j⁻¹ = (a+jc/8, b+j(d−a)/8−j²c/64; c, d−jc/8)`: recomputed, correct; integral and in `Γ₁(32)` for
  `γ ∈ Γ₁(256)`. So `S ∈ M_{λ+1}(Γ₁(256))`, holomorphic at cusps, integer coefficients, and `S = (sieve of g)·φ(q⁸)`.
* `[SL₂(ℤ):Γ₁(256)] = 256²(1−¼) = 49152`. Correct. Sturm (Stein Thm 9.18, any congruence subgroup of index `m`:
  `ord_ℓ > km/12 ⟹ f ≡ 0`) is applied to `S`: integer weight `λ+1`, integer coefficients, hypothesis
  `[qⁿ]S ≡ 0` for `n ≤ (λ+1)·4096`, which follows from the hypothesis on `g` because `[qⁿ]S` is a ℤ-combination of
  `g(n−8m²)`. Right object, right bound: `8·4096 = 32768`, `11·4096 = 45056`, `360·4096 = 1474560`.
* In D the lemma is applied to an integral lift of `g ∈ Λ⊗𝔽₁₉`; fine.

**(c) (F8) and D step 2.** `ΘF/F = (4E₂(q⁴)−E₂(q²))/3`, `Θφ/φ = (10E₂(q²)−2E₂(q)−8E₂(q⁴))/24`: rederived and
checked as power series (`audit/misc.gp`). `E₂ ≡ E_{ℓ+1}`, `E_{ℓ−1} ≡ 1 (mod ℓ)`, both `ℓ`-integral (checked for
`ℓ = 19`). `E_{ℓ±1}(q^d)`, `d | 4`, are forms of *even* weight and trivial character on `Γ₀(4)`, i.e. elements of
`Λ_{ℓ±1}⊗ℤ_{(ℓ)}` (`χ_{−4}^{even} = 1`). Since `Λ_k·Λ_{k′} ⊂ Λ_{k+k′}` for all (half-)integers (monomials in `φ`, `F`),
multiplying keeps the right space, including the character for odd integer weight: `Λ₉ = M₉(Γ₀(4),χ_{−4})`,
`Λ₉·Λ₂₀ ⊂ Λ₂₉ = M₂₉(Γ₀(4),χ_{−4})`, etc. Bookkeeping: `9 + 18·20 = 369`, `9 + 9·20 = 189`,
`369 − 189 = 180 = 10·18` (multiply by `E₁₈^{10}`): correct. `A₊ ≡ ½(Θ^{18}A + Θ⁹A)`: correct (`n^{18} ≡ [19∤n]`,
`n⁹ ≡ (n/19)`, constant term killed). `dim Λ₃₆₉ = 185`; the last five basis elements `F^jφ^{738−4j}`,
`j = 180..184`, are exactly those with `φ`-exponent `< 19` (18, 14, 10, 6, 2), so `u₁₈₀ = … = u₁₈₄ = 0` gives
`A₊ ≡ φ^{19}w₀`, `w₀ ∈ Λ_{719/2}⊗𝔽₁₉` (`dim = 180`). Numerically: the `t`-expansions (`t = F/φ⁴`) of
`Θ^{18}A/φ^{738}`, `Θ⁹A/φ^{378}` and `A₊/φ^{738}` are polynomials of degree `≤ 184`, `≤ 94`, `≤ 179` up to `O(t^{1200})`.

**(d) (F4).** True; citation not established — see G1.

**(e) Step 1 and the final chains.**
* A, B step 1: `1/φ(q) ≡ φ(q)^{ℓ−1}/φ(q^ℓ)`; `U_ℓ(f(q)g(q^ℓ)) = U_ℓ(f)·g`; `U_ℓ(Σp̄(n)(−q)ⁿ) = Σp̄(ℓn)(−1)ⁿqⁿ` (`ℓ` odd).
  `φ^{ℓ−1} ∈ M_{(ℓ−1)/2}(Γ₀(4),χ_{−4}^{(ℓ−1)/2})`, weight 8 resp. 11 `≥ 2`, so `U_ℓ ≡ T_ℓ`. For `ℓ = 17` the only
  basis element of `Λ₈` not divisible by `φ` is `F⁴` and `α₄ ≡ 0`; for `ℓ = 23` every element of `Λ₁₁` is divisible
  by `φ`. Correct. `γ` agrees with the formulas for `c`, `d` printed in the source.
* A, B steps 3–5: `s(−1/p)^λ = +1` in both cases, `11⁶ ≡ 8 (17)`, `13⁹ ≡ 3 (23)`, `η = β − p^{λ−1}γ` recomputed by
  hand and by machine; `(n/p) = s ≠ 0` kills the third Hecke term; `(−1)^{p²n} = (−1)ⁿ` only affects a sign. Correct.
* D step 1: sieving by a set of residues mod 19 commutes with multiplication by the `q^{19}`-series `φ(q^{19})`;
  `φ(q^{19}) ≡ φ^{19}` is a unit. Correct.
* D steps 3–7: `λ = 359`, `17^{358} ≡ 5`, `(n/17) = −1 ⟹ w₀(289n) = b(n) + 5w₀(n) = h(n)`; `g(n) = h(4n)` by
  definition of `U₄`; `(4n/17) = (n/17)` (4 is a square, `17 ∤ 4`); `(4·289·n/19) = (n/19)` (`4·289 = 34²`,
  `19 ∤ 34`); `N = 1156n` is even so `W(N) = +p̄(N)`; the condition `n ≡ 3 (mod 8)` enters only in step 6, on `g`,
  with the same `n`. Correct.

**(f) Statements.** Compared with the TeX source of arXiv:2603.08510 (downloaded; Section "Closing remarks"):
the two displays, the two "Equivalently" lists, Conjecture 1 and Conjecture 2 are quoted correctly, and Theorems
A–D are exactly these statements. The residue lists were recomputed: `{19,35,43,51,83}` mod 88,
`{29,53,61,69,77,101}` mod 104, `{11,43,51}` mod 56, 72 classes mod 2584.

## 2. Independent recomputation

All programs are in `audit/`; none imports or copies the session's code. `φ`, `F`, `p̄` are built from Euler's
pentagonal series / eta-quotients (`p̄ = f₂/f₁²`), not from `Σq^{n²}`, divisor sums or `1/φ(−q)`.

**Theorems A, B — `audit/ab.gp` (PARI/GP, exact integers, 4 s).** Coordinates by solving a linear system and
verifying the expansion on 25–40 further coefficients.

| | A (ℓ=17, p=11) | B (ℓ=23, p=13) |
|---|---|---|
| α mod ℓ | (1,13,13,0,0) | (1,9,5,14,17,20) |
| β mod ℓ | (8,14,14,0) | (3,20,8,10,12,6) |
| η | (0,12,12,0) | (0,16,16,14,7,15) |
| `h(n) ≡ 0`, `n ≡ r (8)`, `n ≤` bound | all 4096 values, `n ≤ 32768` | all 5632 values, `n ≤ 45056` |

All identical to the certificate (also the comparison counts of nonzero `h(n)` by residue class). `c₀(n) =
(−1)ⁿp̄(ℓn)` checked for `n < 3000`.

*Second, basis-free check* (`audit/direct_ab.py`, `audit/direct_b_full.py`): from overpartition numbers alone,
`hd(n) := c₀(p²n) + ((−1/p)^λ(n/p) − s(−1/p)^λ)p^{λ−1}c₀(n) + p^{2λ−1}c₀(n/p²)`, `c₀(n) = (−1)ⁿp̄(ℓn)`. For A:
`hd(n) = h(n)` for **all** `n ≤ 32768` (0 mismatches), so `hd` vanishes on all 4096 `n ≡ 3 (8)`; this needs `p̄ mod 17`
up to `6.74·10⁷`. For B: `hd(n) = h(n)` for **all** `n ≤ 45056` (0 mismatches), so `hd` vanishes on all 5632
`n ≡ 5 (8)` (see the note under the table below for how `p̄(23k) mod 23` was obtained).

**Theorem D — `audit/d.gp` (PARI/GP, 3.5 min).** Coordinates are computed as coefficients of the expansion of
`f/φ^K` in powers of `t = F/φ⁴` (repeated division by `t`), always to more terms than the dimension; `W` is taken
directly from overpartition numbers.
* `u₁₈₀..u₁₈₄ = 0`; moreover all `t`-coefficients of `A₊/φ^{738}` for `180 ≤ j < 1200` vanish.
* `W(n) = w₀(n)` (basis expansion) for all `n < 63590`.
* `b = T(17²)W`: `t`-coefficients `180 ≤ j < 220` vanish; `h` agrees with `(T(17²)+5)W` for `n < 220`;
  `(U₄h)/φ^{719}`: `t`-coefficients `180 ≤ j < 700` vanish.
* The four coordinate vectors `u (= w₀)`, `b`, `h`, `g` (180 entries each) are **identical** to those printed in
  `cert/overpartition_mod19_cert.txt` (`audit/compare_d_coords.py`).
* **Full range:** `g(n) ≡ 0 (mod 19)` for all 184320 values `n ≡ 3 (mod 8)`, `n ≤ 1474560`. Nonzero counts for
  `n < 4000` by class mod 8: `[211,220,224,0,217,223,223,223]`, as in the certificate.
* Basis-free check (`audit/direct_cd.py D 100000`, `p̄ mod 19` up to `1.156·10⁸`):
  `g(n) = W(1156n) + 5((4n/17)+1)W(4n) + 17^{717}W(4n/289)` for **all** `n ≤ 100000` (0 mismatches).

**Direct tests of the congruences themselves** (own `p̄`, python-flint series arithmetic):

| Theorem | range of `p̄` | admissible `n` tested | failures |
|---|---|---|---|
| A | `n < 6.74·10⁷` | 1861 (all admissible `n ≤ 32768`) | 0 |
| B | `p̄(23k)`, `k ≤ 7614464` (arguments up to `1.75·10⁸`) | 2598 (all admissible `n ≤ 45056`) | 0 |
| C | `n < 2²⁶` | 224695 (all admissible `n ≤ 4194299`) | 0 |
| D | `n ≤ 1.156·10⁸` | 2790 (all admissible `n ≤ 99987`) | 0 |

(For B a plain series inversion to `1.75·10⁸` terms did not fit into the 15 GB of this machine. The plain
inversion was run to `1.166·10⁸` (`n ≤ 30000`: 0 mismatches, 1729 admissible `n`, 0 failures); the full range was
then done with `direct_b_full.py`, which gets `p̄(23k) mod 23` from the elementary identity
`Σ_k p̄(23k)q^k ≡ U₂₃(G²²)/G (mod 23)`, `G = f₁²/f₂`, and was cross-checked against the plain inversion for
`k ≤ 130434`.)

Negative controls behave as expected (C with `(n/7) = −1`: 175 nonzero values for `n < 4000`; D with
`(n/17) = +1`: 523 nonzero values for `n < 20000`). The counts 220, 119, 26785, 196 quoted in the note for the
range `8·10⁶` are the correct numbers of admissible `n`.

## 3. Theorem C: `code/radu.py` against the Radu–Sellers lemma

Source of the statement: TeX of arXiv:2302.11830, subsection "An algorithmic approach by Radu and Sellers"
(downloaded into a fresh directory, read only). Line-by-line comparison:

| lemma | `radu.py` | |
|---|---|---|
| `κ = gcd(m²−1, 24)` | `kappa` | ✓ |
| `P_{m,r}(t) = {ts + (s−1)/24·Σδr_δ mod m : [s]_{24m} ∈ 𝕊_{24m}}` | `P_set` (squares of units mod `24m`) | ✓ (see R3) |
| (a) `p | m ⟹ p | N` | `any(N % p …)` | ✓ |
| (b) `r_δ ≠ 0 ⟹ δ | mN` | ✓ | ✓ |
| (c) `κN Σ r_δ mN/δ ≡ 0 (24)` | exact fraction, `% 24` | ✓ |
| (d) `κN Σ r_δ ≡ 0 (8)` | ✓ | ✓ |
| (e) `24m / gcd(−24κt − κΣδr_δ, 24m)` divides `N` | `gcd(k*(-24*t - S), 24*m)` | ✓ |
| (f) `2 | m`: (`4 | κN` and `8 | sN`) or (`2 | s` and `8 | (1−j)N`), `Πδ^{|r_δ|} = 2^s j` | ✓ | ✓ |
| `p_{m,r}(γ) = min_{0≤λ<m} (1/24) Σ r_δ gcd²(δa+δκλc, mc)/(δm)` | `p_mr` | ✓ |
| `p*_a(γ) = (1/24) Σ_{δ|N} a_δ gcd²(δ,c)/δ` | `p_star` | ✓ |
| `ν = (1/24){(Σr_δ+Σa_δ)[Γ:Γ₀(N)] − Σδa_δ} − (1/24m)Σδr_δ − t_min/m` | `nu_bound` | ✓ |
| check `c_r(mn+t′) ≡ 0 (u)`, all `t′ ∈ P`, `0 ≤ n ≤ ⌊ν⌋` | loop in `certificate`, series to `m(⌊ν⌋+1)` | ✓ |

*Representatives.* The double coset `Γ₀(N)γΓ_∞` is determined by the cusp `γ∞ = a/c` of `Γ₀(N)`; every cusp is
`Γ₀(N)`-equivalent to some `a/c` with `c | N`, `gcd(a,c) = 1`, `1 ≤ a ≤ c`, and `p_{m,r}`, `p*` depend on `γ` only
through its first column. So the 56 pairs used by `radu.py` contain a complete set of representatives (8 cusps for
`N = 56`), and checking the inequality on a superset is harmless. Valid. (The ILP that *finds* `a` uses floating
point, but the inequality is re-verified exactly with `Fraction`s, which is all the lemma needs.)

*Recomputation (`audit/radu_c.gp`, PARI/GP, written from the TeX statement).*
* `(m, M, N, (r_δ), t) = (896, 2, 56, (−2, 1), 176)`: conditions (a)–(f) all hold (`κ = 3`; (e): the quotient is
  exactly 56; (f): `s = 1`, `j = 1`, first alternative).
* `P_{m,r}(t) = {176, 688, 816} = {16k : k = 11, 43, 51}`.
* With `a = (a₁,a₂,a₄,a₇,a₈,a₁₄,a₂₈,a₅₆) = (1829, −915, 7143, −256, −14284, 128, −49993, 99986)`:
  `p_{m,r} + p*_a` at the 8 cusps `1/c`: `3/224, 0, 1/8, 3/224, 3/56, 0, 1/8, 1371237/8` for
  `c = 1, 2, 4, 7, 8, 14, 28, 56`; minimum 0 over the 8 cusps, over the session's 56 pairs, and over **all** 4387
  first columns `(a,c)` with `c ≤ 120` (which also tests the double-coset invariance).
* `[Γ:Γ₀(56)] = 96`, `ν = 88009/28`, `⌊ν⌋ = 3143`.
* `p̄(896n + t′) ≡ 0 (mod 7)` for `t′ ∈ {176, 688, 816}`, `0 ≤ n ≤ 3143` (9432 values; `p̄` from PARI's `eta`): all
  zero. Negative control `896n + 48`: 171 nonzero among `n ≤ 200`.

Everything agrees with `cert/radu_overpartition_mod7.txt`. The three progressions `896n + 16k` are exactly
`16(56n + k)`, `k ∈ {11,43,51}`, i.e. Conjecture 1 as printed.

## 4. Programs (all in `audit/`)

| file | what it does | run time |
|---|---|---|
| `ab.gp` | A, B: α, β, η; `h` to the Sturm bound; writes `h_17.txt`, `h_23.txt` | 4 s |
| `direct_ab.py` | A (and B on a partial range): `h` and the congruence straight from `p̄` | 2 min (A) |
| `direct_b_full.py` | B: the same over the full Sturm range, low-memory (23-dissection of `G¹⁶·G⁶`) | 17 min |
| `d.gp` | D: `u`, `w₀`, `b`, `h`, `g`, tests of F3/F4/F8, `g` to the Sturm bound; writes `d_coords.txt`, `d_g.txt` | 3.5 min |
| `compare_d_coords.py` | compares `d_coords.txt` with the session's certificate | — |
| `direct_cd.py` | C, D: congruences from `p̄`; D: `g` against its definition from `p̄` | 2–3 min each |
| `f4test.gp` | (F4), (F3) on all basis elements, `K ≤ 48`, exact arithmetic; `U₂` control | 4 min |
| `misc.gp` | Kronecker periodicity (for G1), Remark after Thm A, constants, residue lists, (F8) identities | 1 s |
| `radu_c.gp` | Theorem C: Δ*, `P`, cusp inequalities, ν, the finite check | 11 s |
| `OUTPUT.txt` | concatenated output of all runs (the `*.log` files are git-ignored) | |

(The Remark after Theorem A — exact vanishing of the coefficients of `Fφ(q)⁷φ(−q)⁴` at `n ≡ 3 (mod 8)` — was
confirmed for `n < 2000`; it is not used in the proof.)
