# Referee report on FR/THEOREM.md (independent, adversarial)

Date: 2026-11-26. Sources:
- arXiv:2111.04468, the 57-page PDF fetched from arXiv's e-print endpoint (journal ref. Arnold Math. J. 10,
  529–566 (2024)).
- My scripts in `/tmp/claude-0/-home-user-weighted-prime-races/7b9224b4-6348-502d-a364-b106f526c180/scratchpad/u2ref/`
  (`fr_test.py`, `fr_degenerate.py`, `fr_degen_check.py`, `fr_thmB.py`, `fr_thmB2.py`).

## Verdict

**COMPLETE MODULO LISTED GAPS.**

- Theorem A, as stated in FR/THEOREM.md, is **correctly proved**. Its proof has no gap. Its exact gcd formula holds
  in every test I ran.
- Theorem A is not literally Conjecture 1.3. It proves the conjecture on the range `B ≠ 0` and
  `x_1, …, x_d ∉ {1,2,3,…}`, with any `m ≠ 0`, which is a superset of the conjecture's `m`.
- On the excluded range the printed conjecture is **false**. I give an explicit counterexample below (finding 1).
  So the exclusion is necessary, and FR/THEOREM.md proves the correct form of the conjecture. The document should
  say this explicitly instead of calling the excluded range a "degenerate case we do not include". It should also
  stop labelling Theorem A simply "Conjecture 1.3".
- The case `B = 0` is excluded silently (finding 2; trivial).
- The side claims about Conjectures 1.1, 1.2 and 1.4 are overstated in the header (findings 3–4).

## The conjecture as printed (arXiv:2111.04468, §2.3)

> "Conjecture 1.3. For any 𝑑 ∈ ℕ, and 𝐵, 𝑥1, … , 𝑥2𝑑 ∈ ℚ, the following PCF[𝑎𝑛, 𝑏𝑛] has FR:
> 𝑏𝑛 = 𝐵 ⋅ ∏_{𝑖=1}^{2𝑑}(𝑛 − 𝑥𝑖), 𝑎𝑛 = (𝐵/𝑚) ⋅ ∏_{𝑖=1}^{𝑑}(𝑛 + 1 − 𝑥𝑖) − 𝑚 ⋅ ∏_{𝑖=𝑑+1}^{2𝑑}(𝑛 − 𝑥𝑖)
> for any 𝑚 ∈ ℚ\{0, ±√|𝐵|} (since then 𝐵 = −𝐴²/4, or 𝐴 = 0)."

FR is defined in §2.1 as `limsup (1/n) ln(GCD_n / n!^{deg a_n}) > −∞`. The paper says the rational formulation
is equivalent to the integer one "by multiplying a_n and b_n by a constant (see inflation process in Appendix C)".

FR/THEOREM.md §0 quotes the conjecture correctly, with one omission: it drops the parenthetical "(since then
B = −A²/4, or A = 0)". Its definition of FR and its convention for rational parameters match the paper.
Inflation by `c` sends `(p_n, q_n)` to `(c^{n+1} p_n, c^n q_n)`; I checked this, so the definition does not
depend on `c`. The family and the roles of `x_1..x_d` (shifted) versus `x_{d+1}..x_{2d}` are the same as in the
paper.

## Findings

### Fatal
None against Theorem A as stated.

### Gaps (in the claim "Theorem A = Conjecture 1.3")

1. **The excluded range `x_i ∈ Z_{≥1}` for some `i ≤ d`: the conjecture is false there; FR/THEOREM.md does not
   say so.**
   - Counterexample: `d = 2`, `B = −3/4`, `m = 1` (allowed: `m ∉ {0, ±√(3/4)}`), `(x_1,x_2,x_3,x_4) = (−2, 1, −2, 2)`.
     Then `a_n = 4, 0, −15/2, −37/2, −33, …` and `b_n = −6, 0, 0, −75/2, …` for `n = 0, 1, 2, …`.
   - This gives `p_n = q_n = 0` for every `n ≥ 1`. So `GCD_n = 0` for all `n ≥ 1`, the "PCF" is undefined (0/0),
     and FR fails (`ln 0 = −∞`).
   - The integer inflation (`B = −12`, `m = 4`; a_n and b_n multiplied by 4 and 16) behaves identically. Both
     were checked in exact arithmetic (`fr_degen_check.py`).
   - A search over small integer roots found 16 such collapses for `d = 2` and none for `d = 1`, where the
     analogous collapse forces `m^2 = −B`, which the conjecture excludes.
   - Remark 2 of FR/THEOREM.md does mention the possibility that "q_n = 0 for all large n". But it treats the
     excluded range only by an informal sketch, and it presents the restriction as a convenience rather than as
     necessary.
   - *What closes it:* state explicitly that Conjecture 1.3 as printed fails on part of the excluded range (give
     the example), and present Theorem A as "Conjecture 1.3 in its correct form (`x_1..x_d ∉ Z_{≥1}`), plus a
     counterexample to the literal statement".
   - Optional: turn Remark 2 into a proved statement for the excluded range. After the last positive-integer
     root `N` of `H_1`, `(p_n, q_n) = (p_{N−1}, q_{N−1})·X_n`. Here `X_n` is the numerator-type sequence of the
     shifted Euler fraction, to which Theorem A's formula applies. FR then holds exactly when `X_n` is not
     eventually 0.
2. **`B = 0` is excluded silently.** The conjecture allows `B ∈ Q` and then only `m ≠ 0`. With `B = 0` one has
   `b_n ≡ 0`, `q_n = ∏_{k=1}^n a_k` and `GCD_n = |q_n|`. FR is trivial unless some `a_k = 0`, which happens iff
   some `x_i` (`i > d`) is a positive integer; then the convergents are 0/0 again. Add one sentence.

### Minor (side claims, Theorem B, presentation)

3. **Header and Corollary: "the 'if' halves of Conjectures 1.1, 1.2".**
   - The Corollary covers only `b` of even degree `2d` with at least `d` roots, counted with multiplicity, that
     are not positive integers.
   - Conjectures 1.1 and 1.2 as printed cover every `b` that splits over `Q`. The paper's §2.1 "unless stated
     otherwise, balanced-degree" arguably restricts to even degree, but the root condition is an extra restriction.
   - The Corollary's own text is accurate. The file title (line 1) is not: change it to "partial 'if' halves".
4. **Theorem B versus the "if" half of Conjecture 1.4.**
   - Theorem B assumes `A^2 + 4B ≠ 0` and "not both `x_i` positive integers". The conjecture's family (Eq. 7),
     `a_n = k(2n − x_1 − x_2 + 1)` with `k ≠ 0`, includes `B = −k^2` (double root `λ_1 = λ_2`). That case is not
     covered.
   - Numerically it does seem to have FR. With `k = 3/2`, `x = (1/3, −2)`, `B = −9/4`, I get
     `(1/n) log(GCD_n/n!) = 0.37, 0.27, 0.22, 0.19` at `n = 50, 100, 200, 400`.
   - So this is a genuine uncovered case, not a counterexample. Either add it (the generating-function method
     adapts, with `W = (1−λz)^{−α}·exp(…)`-type factors) or state the restriction in the header.
5. **Theorem B, "GCD_n ≥ c^n n! for all n".**
   - This fails at any `n` with `p_n = q_n = 0`. That can only happen at isolated `n` once some `b_k = 0`, i.e.
     when `x_2 ∈ Z_{≥1}`.
   - *Fix:* say "for all n with (p_n, q_n) ≠ (0,0)". FR, being a limsup, is unaffected.
6. **Theorem A, "GCD_n ≥ c^n n!^d … for all n".** This is correct: `p_n` and `q_n` cannot both vanish because
   `Π_1(n+1) ≠ 0`. The same is true when `deg a_n < d` (`m^2 = B`), since then `n!^d ≥ n!^{deg a}`. So the claim
   "no condition `m ≠ ±√|B|` is needed" is right.
7. **Novelty.**
   - Theorem A is elementary. The identity `q_n = Π_1(n+1)σ_n` is Euler's classical conversion between such
     continued fractions and series, and is used throughout O. David's "Euler family" papers.
   - A gcd lower bound of the form `n!^d/lcm[n]^d | gcd(p_n, q_n)` for the special case `h_1 = h_2 = n^d` already
     appears in O. David, arXiv:2308.02567**v1** ("Introduction to polynomial continued fractions", the
     Euler-family example before Part I/Appendix).
   - The doc's statement that 2308.02567 "does not discuss FR" holds for v2 but **not for v1**, which defines
     factorial reduction and proves this special case. Correct the citation.

## Line-by-line on the proof of Theorem A (§1–§2)

- **Lemma 1:** correct.
  - The `p`-adic count `#{j<k : p^t | b+jf} ≤ ⌊k/p^t⌋ + 1` holds, and the count is 0 for `p^t > |b|+kf`, since
    `b + jf ≠ 0` because `β ∉ Z_{≤0}`.
  - The numerator bound `≥ v_p(k!)` for `p ∤ e` is fine (if `(α)_k = 0` the claim is vacuous).
  - The cases `p | f` and `p | e` are fine, and so is `v_p(k!) ≤ k`.
  - `lcm(1..M) ≤ 3^M` is Hanson's bound; `4^M` would also do.
  - The constant `C = e·rad(e)·3^{2(|b|+f)}` suffices, and the extra factor for `(β)_{k+1}` is accounted for.
- **Decomposition:** correct. `a_n = H_1(n+1) + H_2(n)` and `b_n = −H_1(n)H_2(n)`. Inflation `H_j → cH_j`
  corresponds to `(B, m) → (c^2B, cm)`.
- **`U_n = Π_1(n+1)` is a solution:** correct.
- **Reduction to first order:** correct. `H_1(n+1)(s_n − s_{n−1}) = H_2(n)(s_{n−1} − s_{n−2})`; dividing by
  `Π_1(n) ≠ 0` is legitimate.
- **Initial values for `p` and `q`:** correct. The formulas hold for all `n ≥ −1` by induction.
- **The gcd step:** correct. `gcd(p_n, q_n) = gcd(N, q_n)` uses `H_2(0) ∈ Z`. Then `D | N` and the gcd is
  `|N|/D`; the case `σ_n = 0` (`D = 1`) is also fine.
- **Size of `D`:** correct. Pair the `d` numerator Pochhammers with the `d` denominator Pochhammers (Lemma 1); the
  powers of `−m` and `m/B` contribute `den^{n+1}`.
- **Lower bound for `|Π_1(n+1)|`:** correct.
- **Undoing inflation:** this is not even needed, because FR for rational data is *defined* through inflation.

## Theorem B (checked less exhaustively)

- The recurrence for `c_j` is correct.
- The operator identity, including the boundary terms `r_0, r_1`, is correct.
- `ΔθW/W = (k + Ax_1)z + B(1+x_1−x_2)z^2` with `λ_1 + λ_2 = A`, `λ_1λ_2 = −B`: correct.
- `C = WV` and `V_j = e_j/(j − x_1)`: correct.
- Case (i): correct, since rational `λ_i` give rational `α_i`.
- Case (ii): correct. `α_1 = α_2` iff `k = (A/2)(1 − x_1 − x_2)`, which is family (Eq. 7) with `A = 2k`.
- The denominator bounds (Lemma 1 with `β = 1`, and the lcm of `|je − a|` bounded by `lcm(1..ne+|a|)`) are
  correct, and so is the final divisibility `M_n | e^{n+1} D_n GCD_n`.
- Remaining issues: findings 4–5.
- Numerically, `(1/n)log(GCD_n/n!)` for an Eq. 6 instance stays bounded (`0.26, 0.23, 0.18, 0.19` at
  `n = 400, …, 3200`). A control outside both families decreases (`−1.24, −1.68, −2.18, −2.55` at
  `n = 50, …, 400`).

## What I recomputed

- **Theorem A, exact integer arithmetic (`fr_test.py`):**
  - 40 random parameter sets: `d ∈ {1,2,3,4}`, rational `B`, `m`, `x_i` with denominators up to 4, and
    `n ≤ 300` (`d ≤ 2`) or `n ≤ 150` (`d = 3, 4`).
  - Four further special cases (`m^2 = B`, `m^2 = −B`, the paper's `PCF[2n^2+2n+1, −n^4]`, `d = 1`, `B = m = 1`)
    up to `n = 200`.
  - In every case `q_n = Π_1(n+1)σ_n`, `p_n = Π_1(n+1)(1 + H_2(0)σ_n)` and
    `GCD_n = |Π_1(n+1)|/den(σ_n)` hold exactly at every `n`.
  - `(1/n) log(GCD_n/n!^d)` stays bounded.
- **The counterexample in finding 1**, exact.
- **Theorem B** spot checks, as above.

## Prior-art search

- Searches: two web searches ("factorial reduction" with Ben David/Kaminer; "factorial reduction" with Euler
  family / conservative matrix field, 2024–2025), Semantic Scholar's citation list for arXiv:2111.04468, and full
  text of the citing papers 2502.17533, 2412.16818, 2211.01058 and 2210.15669.
- **No proof of Conjecture 1.3 was found.** The citing papers (2210.15669, 2211.01058, 2212.09470, 2301.01624,
  2308.11829, 2412.12361, 2412.16818, 2502.17533) cite FR only as a definition or as a heuristic.
- Closest prior art: O. David, arXiv:2308.02567v1, which proves the special case noted in finding 7.
- Not checked: the authors' `verify/` scripts (I wrote my own tests instead), and the journal version's text
  (I used the arXiv PDF, which carries the journal reference).
