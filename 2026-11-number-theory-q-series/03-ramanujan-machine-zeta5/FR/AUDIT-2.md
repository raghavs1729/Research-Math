# AUDIT-2: hostile referee report on THEOREM.md §5–§7 (Theorems C, D, E, Lemma 2, §6, §7)

**Verdict: complete after listed fixes.** I found no mathematical error in Theorem C, Lemma 2 or the main lines of
Theorems D and E, but there are three major problems:
* The written proof of Theorem D (1) has a gap when `ω = 0`. The statement is still true and the fix is one line.
* Theorem E misidentifies its condition (a) with the families *as printed*, because the printed (Eq. 6) requires
  `k ≠ 0`. This also produces a counterexample to the printed "only if" half of 1.4 that THEOREM.md misses, and it
  is non-terminating.
* The document does not say that `P_n` is the paper's own `GCD2_{n+1}`. As a result the remaining gap is exactly
  Conjecture 3 of the paper, which THEOREM.md says it does not treat.

Sections 1–4 (Lemma 1, the proofs of Theorems A and B) were re-read only as far as Theorem E uses them. The
consequences claimed there for `P_n` do follow (finding 5 is a small presentational gap).

## Findings

### 1. (major) Theorem D (1): the proof breaks when `ω = 0`
The proof sets `Q := {p ∉ S ∪ S_0, p ∤ 2Bω : p | some P_n}` and says that "the finitely many other primes"
contribute `O(N)`. If `ω = 0`, then `p ∤ 2Bω` fails for **every** prime, so the "other primes" are all primes and
the argument proves nothing. Case (1) claims "every `A, k`", so it includes `k = A(B+C)/(2B)`, which is the
(Eq. 7)-shaped `a_n` for irreducible `b`. Examples: `PCF[2n+1, n^2+1]`, `PCF[2n+1, n^2−2]`, and the Euler-shaped
`PCF[2n+1, −(n^2+1)]`, `PCF[4n+2, −(4n^2+3)]`.

The statement is still true. Case (1) never uses Theorem C, only the fact that `b` has a root mod `p`, and that does
not need `p ∤ ω`. Numerically, `(ln P_n − ½ ln n!)/n` stays in `[−0.47, −0.10]` for `n = 100, …, 1600` in all four
examples (`audit2/a2_irreducible_omega0.py`), so the `n!^{1/2}` bound looks sharp there.

**Fix:** in case (1), define `Q` without the condition `p ∤ ω`: `Q := {p ∉ S ∪ S_0 ∪ {2}, p ∤ B : p | some P_n}`
`⊂ {p : (Δ_b/p) = 1}`. Keep `p ∤ 2Bω` only in cases (2) and (3).

### 2. (major) The printed (Eq. 6) excludes `k = 0`
This makes Theorem E (a)'s "i.e." false, and it gives a non-terminating counterexample to the printed "only if"
half of 1.4.

The paper prints the families with "for any `k ∈ Q∖{0}`, and `m ∈ Q∖{0, ±√|B|}`" (fr.txt l. 715–718; THEOREM.md
l. 46 quotes this correctly). So `a_n = An` (`k = 0`) with `D` a nonzero square is **not** in the printed (Eq. 6).
It is not in (Eq. 7) either unless `ω = 0`. Yet by Theorem B (i) ("any `k`") such fractions have FR.

* `PCF[n, 2(n+1)(n+2)]`: `A = 1`, `B = 2`, `m = 1`, `D = 9`, `ω = −8`, `b_n > 0`. (Eq. 7) would need
  `a_n = k'(2n+4)`. Exactly, `(ln GCD_n − ln n!)/n = −1.23, −1.42, −1.47, −1.49` at `n = 100, 200, 400, 799`
  (bounded, as Theorem B says). A control outside the families decreases steadily
  (`PCF[n, 3(n+1)(n+2)]`: −2.23, −2.72, −3.07, −3.46).
* `PCF[7n, 8n^2 − 2]`: this is the paper's *own* (Eq. 6) example "for `m = 1` and any `k ∈ Z`" (fr.txt l. 762)
  at `k = 0`. It has FR (−1.30, −1.40, −1.39, −1.48).

These are balanced-degree, non-terminating fractions with `b` split over `Q`. They contradict the printed "only if"
half of 1.4. The paper's own examples show `k ∈ Q∖{0}` is a misprint for (Eq. 6), but it is exactly as fair to call
this a counterexample "to the printed statement" as the examples of §6.

The algebraic identification "(Eq. 6) ⇔ `D` is a nonzero rational square" is correct when `A ≠ 0`, `k` is arbitrary
and `m ≠ ±√|B|`:
* `A = B/m − m` ⇔ `m = (−A ± √D)/2`.
* `m = ±√|B|` forces `A = 0` or `D = 0`.

So the mathematics of Theorem E is right. Its sentence "(a) `a_n` belongs to family (Eq. 6) or (Eq. 7) of
Conjecture 1.4, i.e. …" is not. The search's class `F6` (search/fr_search12.py l. 20) also includes `k = 0`. The
search therefore met such counterexamples and filed them as "in a family".

**Fix:**
* State (a) as "`D` is a nonzero rational square or `ω = 0`", adding "(the families of 1.4 with `k = 0` allowed
  in (Eq. 6), as in the paper's examples)".
* Add `PCF[n, 2(n+1)(n+2)]` / `PCF[7n, 8n^2−2]` to §6 as a fourth, **non-terminating** discrepancy.
* Correct the sentence "These do not touch the intended content" (l. 328).
* Correct the headline claims in the summary (l. 13–19), "the 'only if' half of Conjecture 1.4", and the §7 bullet
  "Proved: outside the families …".

### 3. (major, framing and attribution) `P_n` is the paper's `GCD2`, and the open gap is Conjecture 3
The paper defines `GCD2_n := GCD[p_n, q_n, p_{n−1}, q_{n−1}]` (fr.txt l. 1136), so `P_n = GCD2_{n+1}`. It notes
that `GCD2_n | GCD2_{n+1}` (l. 1144) and states **Conjecture 3**: "For any PCF: `GCD2_n ≐ GCD_n`" (l. 1155).

THEOREM.md says "Conjectures 2 and 3: not treated" (l. 27, l. 92) and presents the persistent gcd as its own
device. In fact:
* Theorems D and E are theorems about the paper's `GCD2`.
* Theorem E together with Conjecture 3 restricted to degree (1,2) would give the full "only if" half of 1.4 (with
  `k = 0` allowed). The open "transient factor" problem of §7 is exactly Conjecture 3 in this family.

Moreover, §6.1's `PCF[n+1, n(n−1)]` appears to contradict Conjecture 3 *as printed* ("for any PCF"):
* `GCD_n ≥ (n+1)!`, while `(ln P_n − ln n!)/n = −2.09, −2.39, −2.80` at `n = 50, 100, 199`.
* So `(1/n) ln(GCD_n/GCD2_{n+1})` grows (2.6, 2.9, 3.3) instead of tending to 0.

This is numerical only; a proof would need an upper bound on `gcd(t_n, t_{n+1})`. The same holds for §6.2.

**Fix:**
* Identify `P_n = GCD2_{n+1}` and cite l. 1136–1155.
* Say that `P_n | P_{n+1}`. This also simplifies Theorem E (a)⇒(b): `P_n ≥ c^n n!` for one index propagates.
* Rephrase §7's open problem as "Conjecture 3 for degree (1,2), non-terminating".
* Replace "Conjectures 2 and 3: not treated" accordingly.
* Optionally, report the §6.1 numerics against Conjecture 3 as printed.

### 4. (minor) Theorem C: two of the four hypotheses are redundant
Theorem C's proof is correct. I checked Steps 1–7 in full:
* Minimality: `n = 0` is impossible since `V_0 = (a_0, 1)`, so in fact `n ≥ 1`. For `n' − 1 ≥ 0` the pair
  `(n'−1, n')` precedes `(n, n+1)`.
* `V_{n'−1} ≢ 0`.
* Step 4, including `n + 1 = n' + p ⇒ x' = x`.
* (R) at `j < n'`, `j = n'`, `j = n'+1` (via `t_{n'} = a_{n'}`), `n'+1 < j ≤ n`, `j = n+1` (in both cases
  `n' < n` and `n' = n`), `j = n+2` and `j > n+2`.
* The operator identity, coefficientwise.
* The logarithmic derivative in characteristic `p`. Multiplicities divisible by `p` give zero coefficients, which
  is harmless. A pole of `N/Δ` at `1/λ_i` must be a root of `Γ`. If `Γ(1/λ_i) ≠ 0`, then `α_i = 0 ∈ F_p`. The
  constant term needs no condition. The decomposition is unique.
* The double-root case: `p ∤ B` forces `p ∤ A` when `p | D`, so `λ ≠ 0`, and `N(1/λ) = ω/(2Bλ)`.
* `α_1 − α_2 = ω/(B(λ_1 − λ_2))`, verified symbolically in `audit2/a2_algebra.py`.

However, "`p` odd" and "`p ∤ B`" follow from the other two hypotheses:
* If `p | B`, then `D ≡ A^2`, which is not a nonzero square only if `p | A`, and then `p | ω`.
* If `p = 2` and `D` is not a nonzero square mod 2, then `D` is even, so `A` is even and `ω` is even.

So the exact hypotheses are "`p ∤ ω` and `D` is not a nonzero square mod `p`". Both are needed: persistence occurs
when either one is dropped (see the tests). This matches the brute force: the classes "only `p | B` fails" and
"only `p = 2` fails" are empty.

**Fix:** state this (optional). Also change "so `n ≥ 0`" to "`n ≥ 1`" (cosmetic).

### 5. (minor) Theorem E, Euler subcase: integrality of `H_j`
The subcase `D = 0 = ω` uses Theorem A's formula `GCD_n = |Π_1(n+1)|/den(σ_n)`. Theorem A proves that formula only
for integer-valued `H_j = κ(n − x_j)`. With integer coefficients `κ = A/2` is an integer, but `κ(n − x_1) ∈ Z[n]`
is not argued.

**Fix:** add "after inflation by a suitable `c` (which multiplies `P_n` by a factor between `c^n` and `c^{n+2}`)".
Then `Π_1(n+1) | Π_1(n+2)` and the lcm argument go through as written.

The Theorem B subcase is fine:
* `M_n = ∏_{j≤n+1}(je − a)`, so `M_n | M_{n+1}`.
* `M_n | e^{n+1} D_n GCD_n` and `M_n | M_{n+1} | e^{n+2} D_{n+1} GCD_{n+1}`, hence
  `M_n | e^{n+2} lcm(D_n, D_{n+1}) P_n`.
* `lcm ≤ C_2^{2n+3}`.
* `P_n ≠ 0` because `Φ_n ≠ 0`.

Coverage of (a) is complete:
* `D` a nonzero square → Theorem B (i).
* `ω = 0`, `D ≠ 0` → Theorem B (ii).
* `D = 0 = ω` → Theorem A with `d = 1`, `m = −κ`.

"`x_1, x_2` not positive integers" suffices everywhere: Theorems A and B need only one root, Theorem D needs both.

### 6. (minor) §7 overstatement: "which is how it arises in every known case"
The terminating examples of §6 (`PCF[n+1, n(n−1)]`, `PCF[n^2+1, (n−1)^2(n^2−2)]`) have FR with `P_n` far from
factorial (finding 3).

**Fix:** write "every known non-terminating case".

### 7. (minor/cosmetic) Theorem D and the summary say "outside the families" with no condition on `A`
For `A = 0` (constant `a_n`), `D = 4B` can be a nonzero square. Such fractions are outside both printed families
(`m = ±√B` is excluded), yet `P_n` is factorial: `PCF[1, (n+1)^2]` gives `(ln P_n − ln n!)/n → −0.975`
(`audit2/a2_examples.py`). §7 says "degree (1,2)", so it is not wrong, but summary l. 15–17 is not explicit.

**Fix:** say "for `A ≠ 0`" in the summary (Theorem D itself is fine as stated).

### 8. (cosmetic) §6 and the paper's standing assumptions
* §2.1 of the paper says "Unless stated otherwise, we focus on 'balanced-degree' PCFs" (fr.txt l. 364). The
  remark "for a constant `a_n` the definition of FR is vacuous" is therefore outside the conjectures' scope; say so.
* The paper assumes `b_n ≠ 0` "to guarantee that the PCF is infinite", but only in Appendix B (l. 1551, proofs of
  Theorems 1 and 2), not for §2.3. Calling §6 counterexamples "to the printed statements" is fair, but this
  assumption should be mentioned.
* §6.3: `PCF[n−1, −(n−1)(n−2)]` already has `p_n = q_n = 0` for `n ≥ 1`, not only for `n ≥ 2`.
* The family membership in §6 checks out:
  * 6.1: `D = 5` and `a ≠ 2kn`.
  * 6.3a: `m ∈ {1, −2}`, `≠ ±√2`, `k = −3 ≠ 0`.
  * 6.3b: `k = 1/2 ≠ 0`.
  * 6.2: `b` has the irreducible factor `n^2 − 2`, and the degrees `(2, 4)` are balanced.

### 9. (cosmetic) Smaller points
* Theorem D cites "Montgomery–Vaughan §4.3" for Mertens' theorem in progressions. Verify the section number, or
  derive it directly: `(4d/·)` is a real non-principal character mod `4|d|` for every non-square integer `d`
  (negative and non-squarefree `d` included), so its kernel is half of the reduced classes, and
  `Σ_{p≤y} χ(p) log p/p = O(1)` because `L(1, χ) ≠ 0`. The claim itself is correct.
* Lemma 2 (c): `i < −C/(2B)` is the sharp form of the bound, but `i < −C/B` is also true.
* The "vanishing `p`-curvature" parenthetical is a heuristic aside. Label it as such.
* §7: "v_p ≈ 2n/p" should read `2n/(p−1)`.

## Checks with no problem found
* **Lemma 2:**
  * (b): Hensel gives at most 2 roots mod `p^t` for `p ∤ 2BΔ_b`; the count is `≤ 2(N/p^t + 1)` and is cut off at
    `p^t ≤ c_2 N^2`.
  * The case `Δ_b = 0` reduces to `B_0(ei − a)^2`.
  * (a): `4Bb = Y^2 − Δ_b`; for `t > v`, `v_p(Y^2) = v`; there are at most 4 square roots of a unit mod `p^s`
    (`p = 2` included); and `i` lies in at most 4 classes mod `p^{t−c}`.
  * (c): `b_i − b_j = (i−j)(B(i+j)+C)`; in the exceptional case `p` divides one of finitely many nonzero `b_i`.
* **Theorem D:**
  * `P_n^2 | GCD_n GCD_{n+1} | Φ_{n+1}`.
  * Steps 1–2 of Theorem C use nothing about `a_n`, and the two zeros of `b` they give have indices `< N`.
  * In case (2), `Q ⊂ {(D/p) = 1}`, and `Q = ∅` in case (3).
  * The exponent: `N Σ_{p≤c_1N, χ(p)=1} log p/p = ½ N log N + O(N)`, and `(n+1)^{(n+1)/2} ≤ C^n n!^{1/2}`.
* **§7 search:** spot check of the "slowest" fraction `(8, 10, −7, 1, 4)`: `s = −2.04, −2.58, −2.84, −3.08` at
  `n = 80, …, 640` reproduces exactly (`audit2/a2_spot_search.py`).
* **"Transient":** correctly described. `p | GCD_n` alone means `t_n ≡ 0` with `t_{n+1} ≢ 0`, and
  `v_p(GCD_n) ≤ v_p(∏ b_i)`. "Heuristically `n^{O(1)}`" is labelled as a heuristic, which is acceptable.

## Test summary (scripts in `audit2/`, my own code)
* **`a2_algebra.py` (sympy):** the operator identity of Step 6 (coefficientwise, generic degree-6 `Γ`); the
  formula `α_1 − α_2 = ω/(B(λ_1−λ_2))`; and `N(1/λ) = ω/(2Bλ)` when `D = 0`. All true.
* **`a2_thmC_bruteforce.py` (seeds 1, 2):**
  * Setup: 2 × 108,000 pairs (fraction, prime), primes `{2, 3, 5, 7, 11, 13}` plus 6 random primes `< 400`,
    `n ≤ 5000`, computed mod `p`.
  * Kinds of fraction: generic; `b` irreducible; `b` with a double root; `D = 0`; `b` vanishing at a positive
    integer; `ω = 0`.
  * **No violation:** 0 persistent among 46,891 + 46,756 pairs satisfying the hypotheses (about 8,000 per kind,
    about 14,700 with `D = 0`).
  * Dropping one hypothesis:
    * Only `D` a nonzero square mod `p`: persistent in 19,810/26,177 and 20,077/26,370 pairs (e.g.
      `(A,k,B,C,E) = (−5, 6, −7, −7, −12)`, `p = 13`, `n = 10`).
    * Only `p | ω`: persistent in 8,769/13,970 and 8,721/13,868 pairs (e.g. `(7, −12, 4, −2, 4)`, `p = 5`,
      `n = 8`).
    * Only `p | B` or only `p = 2`: never occurs (finding 4).
* **`a2_examples.py` (exact integers):**
  * All §6 examples verified: `p_n = q_n`, `GCD_n ≥ (n+1)!` (respectively `≥ n!^2`), and vanishing from `n = 2`
    (respectively `n = 1`).
  * The new `k = 0` examples of finding 2.
  * The `A = 0` example of finding 7.
  * `P_n` for §6.1 (finding 3).
* **`a2_irreducible_omega0.py`:** finding 1.
* **`a2_spot_search.py`:** the §7 spot check.

## Prior art (searched 2026-11-26)
No proof of the "only if" halves, and no result on the gcd of consecutive convergent pairs for these families, was
found.
* The **paper itself** defines `GCD2` and states Conjecture 3 (finding 3). This is the closest prior formulation
  and must be cited: https://arxiv.org/abs/2111.04468 ; journal version https://link.springer.com/article/10.1007/s40598-024-00250-z
* Elimelech, David, …, Kaminer, *Algorithm-assisted discovery of an intrinsic order among mathematical constants*:
  factorial reduction is treated as an empirical property and used as a search filter; no proof.
  https://arxiv.org/abs/2308.11829
* O. David, *On Euler polynomial continued fractions*: https://arxiv.org/abs/2308.02567 (v1 Example 35 as already
  cited). O. David, conservative matrix fields: https://arxiv.org/abs/2303.09318 . Neither treats non-Euler
  `a_n` or the "only if" direction.
* MathWorld "Factorial Reduction" (descriptive only): https://mathworld.wolfram.com/FactorialReduction.html
* Nothing relevant from 2024–2026 turned up. arXiv:2511.10851 is unrelated (polynomial factorization).
