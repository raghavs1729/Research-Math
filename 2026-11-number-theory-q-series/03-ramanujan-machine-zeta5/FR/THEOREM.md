# Factorial reduction: Conjecture 1.3 of arXiv:2111.04468 on its non-degenerate range (and a counterexample outside it); partial "if" halves of Conjectures 1.1, 1.2, 1.4

**Exactly what is proved (summary).**
* *Theorem A:* Conjecture 1.3 for all `d ≥ 1`, `B, m ∈ Q∖{0}`, `x_i ∈ Q` with `x_1, ..., x_d ∉ {1,2,3,...}` — with the exact
  value of `GCD_n` for every `n`.
* *Counterexample:* without that restriction the conjecture **as printed is false**: for `d = 2`, `B = −3/4`, `m = 1`,
  `(x_1,x_2,x_3,x_4) = (−2, 1, −2, 2)` one has `p_n = q_n = 0` for all `n ≥ 1`, so `GCD_n = 0` (§4, Remark 2; found by the
  parent session's referee, re-verified here exactly).
* *Corollary:* Conjecture 1.2 and the "if" half of 1.1 **only for `b` of even degree `2d` having at least `d` roots
  (with multiplicity) that are not positive integers**.
* *Theorem B:* the "if" half of Conjecture 1.4 when `A^2 + 4B ≠ 0` and not both roots are positive integers; the
  remaining case `A^2 + 4B = 0` of family (Eq. 7) (`B = −k^2`) is an Euler family and is covered by Theorem A.
* *Theorems C, D, E (added 2026-11-26; "only if" direction, degree (1,2); refereed, `AUDIT-2.md`: "complete after
  listed fixes", fixes incorporated).* Let `a_n = An + k` (`A ≠ 0`), `b_n = Bn^2 + Cn + E`, `ω := 2kB − A(B+C)`, and let
  `GCD2_n := gcd(p_n, q_n, p_{n−1}, q_{n−1})` — **the paper's own quantity** (its §2.6, where Conjecture 3 states
  `GCD2_n ≐ GCD_n`). *Theorem C:* a prime `p ∤ ω` modulo which `A^2 + 4B` is not a nonzero square never divides two
  consecutive `GCD_n`. *Theorem D:* for non-terminating fractions, `GCD2_n ≤ C^n (n!)^{1/2}` if `b` is an irreducible
  quadratic (every `a_n`), or if `ω ≠ 0` and `A^2+4B` is not a square (`≤ C^n` if `ω ≠ 0`, `A^2 + 4B = 0`). *Theorem E:* for
  non-terminating fractions with `b` split over `Q`: **`GCD2_n` has factorial size iff `A^2+4B` is a nonzero rational
  square or `ω = 0`**, i.e. iff `a_n` is in family (Eq. 6) or (Eq. 7) of Conjecture 1.4 *with `k = 0` allowed in
  (Eq. 6)*. This is Conjecture 1.4 with `GCD_n` replaced by `GCD2_n`; consequently **Conjecture 3 of the paper,
  restricted to non-terminating fractions of degree (1,2), implies the "only if" halves of 1.4 and (for `deg b = 2`)
  of 1.1** (in the corrected form of the next item).
* *Discrepancies with the printed statements (§6):* (i) non-terminating: the printed (Eq. 6) requires `k ≠ 0`, but
  `PCF[n, 2(n+1)(n+2)]` and the paper's own example `PCF[7n + k, 8n^2 − 2]` at `k = 0` have FR by Theorem B and are in
  neither printed family (evidently a misprint: `k = 0` must be allowed in (Eq. 6)); (ii) terminating fractions have
  FR trivially — `PCF[n+1, n(n−1)]` (in neither family) and `PCF[n^2+1, (n−1)^2(n^2−2)]` (`b` not split over `Q`)
  contradict the printed "only if" halves of 1.4 and 1.1; (iii) `PCF[n−3, 2(n−1)(n−3)]` (family (Eq. 6)) has
  `p_n = q_n = 0` for `n ≥ 2`, contradicting the printed "if" half of 1.4.
* **Still open:** the "only if" halves of 1.1 and 1.4 for the paper's definition (`limsup` of a single `GCD_n`), for
  non-terminating fractions — equivalently, by Theorem E, Conjecture 3 in degree (1,2): what is missing is a bound
  for common factors of `(p_n, q_n)` that do not persist to `n+1` ("transient" factors, §7). A systematic exact
  search (1.39 million fractions of degree (1,2), §7) found no non-terminating fraction with FR outside the
  (corrected) families. Conjecture 2: not treated. Conjecture 3: only as just described.

## 0. Source and statements

**Source (read 2026-11-25).** N. Ben David, G. Nimri, U. Mendlovic, Y. Manor, C. De la Cruz Mengual, I. Kaminer,
*On the Connection Between Irrationality Measures and Polynomial Continued Fractions*, arXiv:2111.04468v2
(Arnold Math. J. 2024), §2.1–§2.3. Convention: `PCF[a_n,b_n] = a_0 + b_1/(a_1 + b_2/(a_2 + ...))`, convergents
`p_n/q_n` (`u_n = a_n u_{n−1} + b_n u_{n−2}`, `(p_{−1},p_0) = (1,a_0)`, `(q_{−1},q_0) = (0,1)`), `GCD_n := gcd(p_n,q_n)`.
*Factorial reduction (FR)* (§2.1): `limsup_{n→∞} (1/n) ln(GCD_n / n!^{deg a_n}) > −∞`.

> "**Conjecture 1.1.** Given a polynomial b_n, there exists an a_n for which PCF[a_n, b_n] has FR if and only if b_n
> is decomposable to linear factors over Q."
> "**Conjecture 1.2.** For each b_n that is decomposable to linear factors over Q, there exists at least one infinite
> family of a_n s for which every PCF[a_n, b_n] has FR. Such a family is presented in Conjecture 1.3."
> "**Conjecture 1.3.** For any d ∈ N, and B, x_1, ..., x_{2d} ∈ Q, the following PCF[a_n, b_n] has FR:
> b_n = B·∏_{i=1}^{2d}(n − x_i), a_n = (B/m)·∏_{i=1}^{d}(n + 1 − x_i) − m·∏_{i=d+1}^{2d}(n − x_i)
> for any m ∈ Q∖{0, ±√|B|}."
> "**Conjecture 1.4.** For every b_n of the form b_n = B(n − x_1)(n − x_2) with B, x_1, x_2 ∈ Q, PCF[a_n, b_n] has FR
> if and only if a_n belongs to one (or more) of the families a_n^{(k)} = (B/m − m)n + k (Eq. 6),
> a_n^{(k)} = k(2n − x_1 − x_2 + 1) (Eq. 7), for any k ∈ Q∖{0}, and m ∈ Q∖{0, ±√|B|}."

The paper works with integer polynomials and says that the rational formulation is equivalent "by multiplying a_n
and b_n by a constant (see inflation process in Appendix C)". Accordingly, for `a, b ∈ Q[n]` we say that
`PCF[a,b]` has FR if `PCF[c·a, c^2·b]` has FR for one (equivalently every) integer `c ≥ 1` with `ca, c^2 b`
integer-valued; this is independent of `c`, because inflation by `c` multiplies `(p_n, q_n)` by `(c^{n+1}, c^n)`.

**Openness.** Listed as conjectures in the journal version; our survey (`../STATUS-MAP-RM.md` §2.4) found no proof.
The family of 1.3 is the "Euler family" of O. David, arXiv:2308.02567. **Prior partial result:** the first arXiv
version of that paper (v1, Aug 2023; the passage is not in v2) introduces the term "factorial reduction" informally
and proves in its Example 35 the special case `h_1 = h_2 = x^d` (the `ζ(d)` fractions `PCF[n^d + (n+1)^d, −n^{2d}]`):
`(n!/lcm[n])^d` divides `gcd(p_n, q_n)` — by an eigenvector argument. Theorem A extends this to the whole family of
Conjecture 1.3 (arbitrary rational roots, `B`, `m`), with the exact gcd; in David's case it reads
`GCD_n = (n+1)!^d/den(Σ_{k≤n} (k+1)^{−d})`, consistent with his divisibility. A second search on 2026-11-26 (web, arXiv; journal version Arnold Math. J. 2024, doi 10.1007/s40598-024-00250-z,
which still states them as conjectures) found no proof either. "Open" means "not found after these searches". The results below are elementary.

> **Theorem A (Conjecture 1.3 on the range `x_1..x_d ∉ Z_{≥1}`).** Let `d ≥ 1`, `B, m ∈ Q∖{0}`, `x_1, ..., x_{2d} ∈ Q` with
> `x_1, ..., x_d ∉ {1, 2, 3, ...}`, and `a_n, b_n` as in Conjecture 1.3. Put `H_1(n) = (B/m)∏_{i≤d}(n − x_i)`,
> `H_2(n) = −m∏_{i>d}(n − x_i)`, `Π_j(n) = ∏_{i=1}^n H_j(i)`, `σ_n = Σ_{k=0}^n Π_2(k)/Π_1(k+1)`. Then
> `q_n = Π_1(n+1) σ_n`, `p_n = Π_1(n+1)(1 + H_2(0) σ_n)`; if the `H_j` are integer-valued,
> `GCD_n = |Π_1(n+1)| / den(σ_n)`; and there is `C > 0` with `den(σ_n) ≤ C^n`, so that
> `GCD_n ≥ c^n n!^d` for some `c > 0` and **all** `n`: the fraction has FR. (No condition `m ≠ ±√|B|` is needed.)
> Here `p_n` and `q_n` never vanish simultaneously (`p_n − H_2(0) q_n = Π_1(n+1) ≠ 0`).
> *Excluded cases.* `B = 0`: then `b ≡ 0`, the fraction is the constant `a_0` and `q_n = ∏_{i≤n} a_i`; this trivial
> case is not treated. Some `x_i` (`i ≤ d`) a positive integer: the fraction terminates (`b_{x_i} = 0`), and here the
> conjecture as printed can fail — see the counterexample in Remark 2.

> **Corollary (Conjecture 1.2 and the "if" half of 1.1, restricted to: `b` of even degree — the paper's
> "balanced-degree" setting — with at least `d` roots that are not positive integers).** If `b_n = B∏_{i=1}^{2d}(n − x_i)` with rational `B ≠ 0`, `x_i`, and at least `d` of the `x_i`
> (with multiplicity) are not positive integers, then the `a_n` of Theorem A, `m ∈ Q∖{0}`, form an infinite family
> with FR.

> **Theorem B (the "if" half of Conjecture 1.4).** Let `A, k, B, x_1, x_2 ∈ Q`, `A ≠ 0`, `B ≠ 0`, `A^2 + 4B ≠ 0`,
> not both `x_i` positive integers, `a_n = An + k`, `b_n = B(n − x_1)(n − x_2)`. If
> (i) `A^2 + 4B` is the square of a rational number (equivalently `A = B/m − m` with `m ∈ Q∖{0}`: family (Eq. 6),
> any `k`), or (ii) `k = (A/2)(1 − x_1 − x_2)` (family (Eq. 7)), then `PCF[a_n, b_n]` has FR; in fact
> `GCD_n ≥ c^n n!` for every `n` with `(p_n, q_n) ≠ (0,0)`. If neither `x_i` is a positive integer this is every `n`
> (`p_n q_{n−1} − p_{n−1} q_n = ±∏ b_i ≠ 0`); if one of them is a positive integer `N`, then `(p_n,q_n) ≠ (0,0)` for at
> least every second `n` (for `n ≥ N` the pair is `t_n · (p_{N−1}, q_{N−1})` with `t_{N−1} = 1` and
> `t_n = a_n t_{n−1} + b_n t_{n−2}`, `b_n ≠ 0` for `n > N`, so `t_n, t_{n+1}` are not both 0), which suffices for the limsup.
> *The case `A^2 + 4B = 0` of family (Eq. 7)* (`A = 2κ`, `B = −κ^2`, `a_n = κ(2n − x_1 − x_2 + 1)`) is not covered by
> the proof below but is the Euler family with `d = 1`, `m = −κ` (`a_n = κ(n+1−x_1) + κ(n−x_2)`), so FR holds by
> Theorem A when `x_1` or `x_2` is not a positive integer. The case `A^2 + 4B = 0` of family (Eq. 6) is excluded
> in the conjecture itself (`m ≠ ±√|B|`).

The "only if" direction is treated in §5–§7 (Theorems C, D, E): proved for the paper's two-step gcd
`GCD2_{n+1} = gcd(p_n,q_n,p_{n+1},q_{n+1})`, open for the `limsup` definition with `GCD_n` (it would follow from the
paper's Conjecture 3 in degree (1,2)); Conjecture 2 is not treated.

**Verification.** `verify/fr_euler_family.py`: for 60 random parameter sets (`d ≤ 3`) the formulas for `p_n, q_n`
and `GCD_n = |Π_1(n+1)|/den(σ_n)` hold exactly for all `n ≤ 50`. `verify/fr_degree1_families.py`: the closed form
of Theorem B's proof equals the recurrence exactly (`j ≤ 120`) for 25 + 25 random parameter sets, and
`(1/n) log(GCD_n/n!)` stays bounded (≥ −1.91 at n = 30, 60, 120), while in a control outside both families it
decreases steadily (−1.03, −1.60, −2.05, −2.47 at n = 30, 60, 120, 240). `verify/fr_only_if.py`: the algebra of
Theorem C (sympy); Theorem C by brute force (11064 pairs (fraction, prime `p < 400`), `n ≤ 3000`: no violation;
for primes modulo which `A^2+4B` is a nonzero square, two consecutive gcds divisible by `p` do occur, in 6621 of
11769 pairs); the examples of §6 exactly. `search/fr_search12.py`, `search/fr_stage2.py`: the search of §7.

---

## 1. A lemma on denominators

**Lemma 1.** Let `α ∈ Q`, `β ∈ Q∖Z_{≤0}`. There is an integer `C ≥ 1` such that for every `n ≥ 0` the numbers
`(α)_k/(β)_k` and `(α)_k/(β)_{k+1}` (`0 ≤ k ≤ n`) have a common denominator `≤ C^{n+1}`.

*Proof.* Write `α = a/e`, `β = b/f` in lowest terms, `e, f ≥ 1`. Then
`(α)_k/(β)_k = (f/e)^k ∏_{j<k}(a + je)/∏_{j<k}(b + jf)`. Let `p` be a prime. If `p ∤ f`, then for each `t ≥ 1`
exactly one residue class of `j mod p^t` has `p^t | b + jf`, so `#{j < k : p^t | b + jf} ≤ ⌊k/p^t⌋ + 1`, and it is
`0` when `p^t > |b| + kf`; hence `v_p(∏_{j<k}(b+jf)) ≤ v_p(k!) + ⌊log_p(|b| + kf)⌋`. If `p | f` then
`v_p(∏(b+jf)) = 0`. Likewise `v_p(∏_{j<k}(a+je)) ≥ v_p(k!)` if `p ∤ e` (and `≥ 0` always). Therefore

    v_p((α)_k/(β)_k) ≥ −k v_p(e) − [p | e] v_p(k!) − [p ∤ f] ⌊log_p(|b| + kf)⌋ ≥ −k (v_p(e) + [p | e]) − ⌊log_p(|b| + kf)⌋.

So the denominator of `(α)_k/(β)_k` divides `(e·rad(e))^k · lcm(1, 2, ..., |b| + kf)`, and
`lcm(1..M) ≤ 3^M`. For `(α)_k/(β)_{k+1} = ((α)_k/(β)_k) · f/(b + kf)` multiply by one more factor
`lcm(1..|b|+nf)`. Take `C = e·rad(e)·3^{2(|b|+f)}`. ∎

Products of `d` such ratios and powers of fixed rationals again have common denominators `≤ C^{n+1}`.

## 2. Proof of Theorem A

`a_n = H_1(n+1) + H_2(n)` and `b_n = −H_1(n) H_2(n)`. By inflation we may assume `H_1, H_2` integer-valued (multiply
`a` by `c`, `b` by `c^2`, `H_j` by `c`). By hypothesis `Π_1(n) ≠ 0` for all `n`.

*Solutions.* `U_n := Π_1(n+1)` satisfies `U_n = a_n U_{n−1} + b_n U_{n−2}`, since
`a_n Π_1(n) − H_1(n)H_2(n) Π_1(n−1) = (H_1(n+1) + H_2(n) − H_2(n)) Π_1(n)`. Writing a solution as `u_n = Π_1(n+1) s_n`
gives `H_1(n+1)(s_n − s_{n−1}) = H_2(n)(s_{n−1} − s_{n−2})`, i.e.
`s_n − s_{n−1} = (s_0 − s_{−1}) · H_1(1) Π_2(n)/Π_1(n+1)`. With the initial values of `q` (`s_{−1} = 0`,
`s_0 = 1/H_1(1)`) and of `p` (`s_{−1} = 1`, `s_0 = a_0/H_1(1) = 1 + H_2(0)/H_1(1)`):

    q_n = Π_1(n+1) σ_n,       p_n = Π_1(n+1) (1 + H_2(0) σ_n),      σ_n = Σ_{k=0}^n Π_2(k)/Π_1(k+1).

*The gcd.* `gcd(p_n, q_n) = gcd(Π_1(n+1), q_n)` (`H_2(0) ∈ Z`). Write `σ_n = P/D` in lowest terms and
`N = Π_1(n+1) ∈ Z∖{0}`. As `q_n = NP/D ∈ Z` and `gcd(P,D) = 1`, `D | N`; then
`gcd(N, (N/D)P) = (N/D) gcd(D, P) = |N|/D`.

*Size.* `Π_2(k)/Π_1(k+1) = (−m)^k (m/B)^{k+1} ∏_{i>d}(1−x_i)_k / ∏_{i≤d}(1−x_i)_{k+1}`, a product of `d` ratios as
in Lemma 1 (`1 − x_i ∉ Z_{≤0}` for `i ≤ d`) and of powers of two fixed rationals; so `D = den(σ_n) ≤ C^{n+1}`. And
`|Π_1(n+1)| = |B/m|^{n+1} ∏_{i≤d} |(1−x_i)_{n+1}| ≥ c_1^{n+1} n!^d`, because `|(1−x)_{n+1}| = ∏_{j=1}^{n+1}|j − x|`
with all factors `≥ δ > 0` and `|j − x| ≥ j/2` for `j ≥ 2|x|`. Hence `GCD_n ≥ (c_1/C)^{n+1} n!^d`. Undoing the
inflation changes `GCD_n` by a factor between `c^{−(n+1)}` and `1`. ∎

*Proof of the Corollary.* Order the roots so that `x_1, ..., x_d ∉ Z_{≥1}`; distinct `m` give distinct `a_n`
(the leading coefficients `B/m − m` take infinitely many values). ∎

## 3. Proof of Theorem B

By symmetry assume `x_1 ∉ Z_{≥1}`. Inflation maps `(A, k, B) ↦ (cA, ck, c^2B)` and preserves (i) and (ii), so
assume `a_n, b_n` integer-valued. For a solution `u` of the recurrence put `c_j := u_{j−1}/(1−x_1)_j` (`j ≥ 0`).
Dividing the recurrence by `(1−x_1)_{j−1}`:

    (j − x_1) c_j = (A(j−1) + k) c_{j−1} + B (j − 1 − x_2) c_{j−2}       (j ≥ 2).

For `C(z) = Σ c_j z^j` (formal) and `θ = z d/dz`, `Δ(z) := 1 − Az − Bz^2 = (1−λ_1 z)(1−λ_2 z)` (`λ_1 ≠ λ_2` the roots
of `λ^2 = Aλ + B`), this reads

    Δ · (θ − x_1) C = ((k + A x_1) z + B(1 + x_1 − x_2) z^2) C + r_0 + r_1 z,    r_0 = −x_1 c_0,  r_1 = (1−x_1)c_1 − k c_0.

Let `W := (1−λ_1 z)^{−α_1}(1−λ_2 z)^{−α_2}` with `α_1 + α_2 = 1 + x_1 − x_2`, `α_1λ_1 + α_2λ_2 = k + A x_1`; then
`Δ θW = ((k + A x_1) z + B(1+x_1−x_2) z^2) W`, and `C = W V` with `(θ − x_1) V = (r_0 + r_1 z)/(Δ W) =: Σ_j e_j z^j`,
i.e.

    c_n = Σ_{i+j=n} g_i · e_j/(j − x_1),      W = Σ g_i z^i,     Σ e_j z^j = (r_0 + r_1 z)(1−λ_1 z)^{α_1−1}(1−λ_2 z)^{α_2−1}

(if `x_1 = 0` then `r_0 = 0 = e_0` and the term `j = 0` is replaced by `c_0`). (Checked exactly against the
recurrence: `verify/fr_degree1_families.py`.)

*Case (i):* `λ_1, λ_2 ∈ Q`, hence `α_1, α_2 ∈ Q`, and `(1−λz)^{−α} = Σ (α)_i λ^i z^i/i!`.
*Case (ii):* the two linear equations give `α_1 = α_2 = α := (1 + x_1 − x_2)/2 ∈ Q` exactly when
`α(λ_1+λ_2) = k + A x_1`, i.e. `k = (A/2)(1 − x_1 − x_2)`; then `W = Δ^{−α} = Σ_i ((α)_i/i!)(Az + Bz^2)^i` and
`1/(ΔW) = Δ^{α−1}`.
In both cases Lemma 1 (with `β = 1`) shows that `g_0..g_n` and `e_0..e_n` have common denominators `≤ C^{n+1}`; the
numerators of `j − x_1` (`j ≤ n`) have lcm `≤ C^{n+1}`; and `c_0, c_1` are fixed rationals for `u = p` and for
`u = q`. So there are integers `D_n ≤ C_2^{n+1}` with `D_n c^{(p)}_j, D_n c^{(q)}_j ∈ Z` for `j ≤ n+1`.

Write `x_1 = a/e` in lowest terms; `(1−x_1)_{n+1} = M_n/e^{n+1}` with `M_n = ∏_{j=1}^{n+1}(je − a) ∈ Z∖{0}`,
`|M_n| ≥ c_1^{n+1} n!`. From `p_n = M_n c^{(p)}_{n+1}/e^{n+1}`, `q_n = M_n c^{(q)}_{n+1}/e^{n+1}` we get that `M_n` divides
`e^{n+1} D_n p_n` and `e^{n+1} D_n q_n`, hence `e^{n+1} D_n GCD_n`. Therefore, whenever `(p_n, q_n) ≠ (0,0)`,
`GCD_n ≥ |M_n|/(e^{n+1} D_n) ≥ c^{n+1} n!`. ∎

## 4. Remarks

1. In Theorem A the conclusion holds for every `n` (not only as a limsup) and gives the exact value of `GCD_n`.
   For the paper's example `PCF[2n^2+2n+1, −n^4]` (`B = −1`, `m = −1`, `x_i = 0`): `GCD_n = (n+1)!^2/den(Σ_{k≤n} 1/(k+1)^2)`.
2. **The excluded range, and a counterexample to Conjecture 1.3 as printed.** If `x_i ∈ Z_{≥1}` for some `i ≤ d`,
   then `b_{x_i} = 0` and the fraction is finite. In this range the printed statement ("for any ... `x_1, ..., x_{2d} ∈ Q`")
   is **false**: take `d = 2`, `B = −3/4`, `m = 1`, `(x_1, x_2, x_3, x_4) = (−2, 1, −2, 2)`, i.e.
   `b_n = −(3/4)(n+2)^2(n−1)(n−2)`, `a_n = −(3/4)(n+3)n − (n+2)(n−2)`. Then `a_0 = 4`, `a_1 = 0`, `b_1 = b_2 = 0`, so
   `p_1 = q_1 = 0`, and then `p_n = a_n p_{n−1} + b_n p_{n−2} = 0`, `q_n = 0` for all `n ≥ 1` (`b_2 = 0` kills the term with
   `p_0, q_0`). Hence `GCD_n = 0` for all `n ≥ 1` and there is no factorial reduction; the same happens for the
   integer inflation `B = −12`, `m = 4` (`verify/fr_counterexample.py`, exact). (This example is due to the parent
   session's referee, `REFEREE-FR.md` on the parent branch.) When the pair `(p_n, q_n)` does not vanish
   identically for large `n`, the argument of §2 applied to the shifted fraction still gives FR in the limsup
   sense; we do not formulate this as a theorem.
3. For `deg a = 1`, `deg b = 2` and `λ_1 ≠ λ_2`: if neither (i) nor (ii) holds, then `λ_1, λ_2` are conjugate quadratic
   irrationals and `α_1 ≠ α_2` are conjugate, irrational. The coefficients `(α)_i/i!` then have unbounded
   denominators at every prime inert in `Q(√(A^2+4B))`, which is the expected mechanism for the failure of FR
   (the "only if" half of Conjecture 1.4). Turning this into a proof requires a lower bound for the denominators
   of the convolution sums for all large `n`, which we did not obtain.

## 5. The "only if" direction in degree (1,2): Theorems C, D, E

Throughout: `a_n = An + k`, `b_n = Bn^2 + Cn + E` with integers `A, k, B, C, E`, `B ≠ 0` (the rational case reduces
to this by inflation, which multiplies `(p_n, q_n)` by `(c^{n+1}, c^n)`),

    D := A^2 + 4B,        Δ_b := C^2 − 4BE,        ω := 2kB − A(B + C).

If `b_n = B(n−x_1)(n−x_2)`, then `(B+C)/B = 1 − x_1 − x_2`, so `ω = 0` iff `k = (A/2)(1 − x_1 − x_2)` (family (Eq. 7)),
and, for `A ≠ 0`, `A = B/m − m` with `m ∈ Q∖{0, ±√|B|}` iff `D` is a nonzero rational square (`m = (−A ± √D)/2`;
`m = ±√|B|` forces `A = 0` or `D = 0`) — family (Eq. 6) **with every `k ∈ Q` allowed**; the printed family excludes
`k = 0`, see §6. `Φ_n := p_n q_{n−1} − p_{n−1} q_n = (−1)^{n−1} ∏_{i≤n} b_i`.

**The paper's `GCD2`.** §2.6 of the paper defines `GCD2_n := GCD[p_n, q_n, p_{n−1}, q_{n−1}]`, notes `GCD2_n | GCD2_{n+1}`
and `GCD2_n | GCD_n`, and states "Conjecture 3. For any PCF: GCD2_n ≐ GCD_n" (`≐`: same exponential order). Below,
`P_n := gcd(p_n, q_n, p_{n+1}, q_{n+1}) = GCD2_{n+1}`; so `P_n | P_{n+1}`.

> **Theorem C (no persistence).** Let `p` be a prime with `p ∤ ω`, and suppose that `D` is not a nonzero
> square modulo `p` (i.e. `p | D` or `(D/p) = −1`). Then there is no `n ≥ 0` with `p | GCD_n` and `p | GCD_{n+1}`;
> i.e. `p ∤ GCD2_n` for all `n`.
> (These hypotheses imply `p` odd and `p ∤ B`: if `p | B` then `D ≡ A^2`, so `p | A` and `p | ω`; if `p = 2` then `D`,
> hence `A`, hence `ω`, is even. Both hypotheses are needed: the referee's tests find persistence in most cases
> when either is dropped.)

*Proof.* Work in `F = F_p`; congruences are mod `p`. Suppose the contrary and take `n` minimal. `V_m := (p_m, q_m)`.
`V_{−1} = (1,0)` and `V_0 = (a_0, 1)`, so `n ≥ 1`, and `V_n ≡ V_{n+1} ≡ 0`.
*Step 1.* `V_{n+1} = a_{n+1} V_n + b_{n+1} V_{n−1}` gives `b_{n+1} V_{n−1} ≡ 0`; by minimality `V_{n−1} ≢ 0`; so `b_{n+1} ≡ 0`.
*Step 2.* Some `i ∈ [1, n]` has `b_i ≡ 0`: otherwise `Φ_n ≢ 0`, contradicting `V_n ≡ 0`. Let `n'` be the largest such `i`.
*Step 3.* Since `b_{n'} ≡ 0`, `V_{n'} ≡ a_{n'} V_{n'−1}`, hence `V_m ≡ t_m V_{n'−1}` for all `m ≥ n'−1`, where
`t_{n'−1} = 1`, `t_{n'} = a_{n'}`, `t_m = a_m t_{m−1} + b_m t_{m−2}` (`m > n'`). `V_{n'−1} ≢ 0` (otherwise `V_{n'−1} ≡ V_{n'} ≡ 0`
with `n'−1 < n`). Hence `t_n = 0`.
*Step 4.* As `p ∤ B`, `b ≡ B(m − x)(m − x')` in `F[m]` with `x ≡ n'` (`b` has the root `n'`, hence splits). The indices
`m` with `b_m ≡ 0` are those `≡ x` or `≡ x'`; `n' < n+1` are consecutive such indices, so `n + 1 ≤ n' + p` and
`n + 1 ≡ x'` (if `n+1 = n'+p`, no index of the class `x'` lies strictly between, so `x' = x`).
*Step 5.* Put `γ_j := t_{j−1} / ∏_{i=n'+1}^{j} (i − x)` for `n' ≤ j ≤ n` (the factors are nonzero since `n' < i < n'+p`)
and `γ_j := 0` otherwise; `γ_{n'} = 1`. Then for **all** `j ∈ Z`

    (j − x) γ_j = (A(j−1) + k) γ_{j−1} + B (j − 1 − x') γ_{j−2}.                                          (R)

Indeed: `j < n'` trivial; `j = n'`: both sides 0 (`n' ≡ x`); `n' < j ≤ n`: divide the recurrence of `t` at `m = j−1`
by `∏_{i=n'+1}^{j−1}(i−x)` (for `j = n'+1` use `t_{n'} = a_{n'}`); `j = n+1`: the right side is `t_n/∏_{i=n'+1}^{n}(i−x) = 0`;
`j = n+2`: the right side is `B(n+1−x')γ_n = 0` by Step 4; `j > n+2` trivial.
*Step 6.* `Γ(z) := Σ γ_j z^j ∈ F[z]` is a nonzero polynomial, and (R) says, with `θ = z d/dz`,

    Δ(z) · θΓ = N(z) · Γ,        Δ := 1 − Az − Bz^2,      N := x + kz + B(1 − x')z^2.

*Step 7.* Over an algebraic closure, `Γ = z^{n'} ∏_ζ (1 − z/ζ)^{m_ζ}`, so
`N/Δ = θΓ/Γ = n' + Σ_ζ m_ζ − Σ_ζ m_ζ/(1 − z/ζ)`: a rational function with at most simple poles whose coefficients
at `1/(1 − z/ζ)` lie in the prime field `F_p`.
If `p | D`: `Δ = (1 − λz)^2`, `λ = A/2 ≠ 0` (`λ^2 = −B`), so `(1 − λz) | N`, i.e. `0 = N(1/λ) = ω/(2Bλ)` (using
`x + x' = −C/B`): contradiction with `p ∤ ω`.
If `(D/p) = −1`: `Δ = (1−λ_1z)(1−λ_2z)` with `λ_1 ≠ λ_2` in `F_{p^2}`, `(λ_1 − λ_2)^2 = D`, and
`N/Δ = c + α_1/(1−λ_1z) + α_2/(1−λ_2z)` with `α_1 − α_2 = ω/(B(λ_1 − λ_2))` (`verify/fr_only_if.py` (a)). By the above
`α_1, α_2 ∈ F_p`, so `λ_1 − λ_2 = ω/(B(α_1−α_2)) ∈ F_p` (note `α_1 ≠ α_2` as `ω ≢ 0`), i.e. `D` is a square mod `p`:
contradiction. ∎

(Heuristic aside: a prime can divide two consecutive gcds only if the rank-one differential equation
`Δ θW = N W` — the one solved by `W = z^{x}(1−λ_1z)^{−α_1}(1−λ_2z)^{−α_2}` in the proof of Theorem B — has a polynomial
solution mod `p`, i.e. vanishing `p`-curvature.)

**Lemma 2 (valuations of `∏ b_i`).** Assume `b_i ≠ 0` for all `i ≥ 1`. (a) For every prime `p` there is `C_p` with
`v_p(∏_{i≤N} b_i) ≤ C_p N` for all `N ≥ 1`. (b) There are a finite set `S` of primes and `c_2` such that for `p ∉ S`:
`v_p(∏_{i≤N} b_i) ≤ 2N/(p−1) + 2 log_p(c_2 N^2)`. (c) There are a finite set `S_0` and `c_1` such that every prime `p`
dividing `b_i` and `b_j` for two indices `1 ≤ i < j ≤ N` satisfies `p ≤ c_1 N` or `p ∈ S_0`.

*Proof.* (b) If `Δ_b ≠ 0` and `p ∤ 2BΔ_b`, the congruence `b_i ≡ 0 (mod p^t)` has at most 2 solutions mod `p^t`
(simple roots mod `p`, Hensel), so `#{i ≤ N : p^t | b_i} ≤ 2(N/p^t + 1)`, and it is 0 when `p^t > max_{i≤N}|b_i|`,
which is `≤ c_2 N^2`; sum over `t ≥ 1`. If `Δ_b = 0`, `b_i = B_0 (ei − a)^2` with `gcd(a,e) = 1`, and for `p ∤ e·num(B_0)den(B_0)`
the same count holds with `ei − a`. (a) For the finitely many remaining primes: if `Δ_b = 0` or `b` is split, use
`b_i = B_0(e_1 i − a_1)(e_2 i − a_2)` (`p | e_j` gives `v_p(e_j i − a_j) = 0`; otherwise as in (b); `B_0` contributes
`N v_p(B_0)`). In general, `4B b_i = Y_i^2 − Δ_b` with `Y_i = 2Bi + C`; let `v := v_p(Δ_b)`. If `p^t | Y^2 − Δ_b` and `t > v`,
then `v_p(Y^2) = v`, `Y = p^{v/2} Y_1`, and `Y_1^2 ≡ Δ_b/p^v (mod p^{t−v})` has at most 4 solutions mod `p^{t−v}`; so `i`
lies in at most 4 classes modulo `p^{t−c}`, `c := v/2 + v_p(2B)` (when `t ≥ c`). Hence
`#{i ≤ N : p^t | 4Bb_i} ≤ 4(N p^{c−t} + 1)` for `t > v` and `≤ N` for `t ≤ v`; sum over `t ≤ log_p(4|B| c_2 N^2)`.
(c) `p | b_i − b_j = (i − j)(B(i+j) + C)`. If `B(i+j) + C ≠ 0`, then `p ≤ max(N, 2|B|N + |C|)`. If `B(i+j) + C = 0`, then
`i < −C/B`, and `p` divides one of the finitely many nonzero integers `b_i`, `1 ≤ i < −C/B`. ∎

> **Theorem D (upper bounds for `GCD2`).** Assume `b_n ≠ 0` for all `n ≥ 1`, and let
> `P_n := gcd(p_n, q_n, p_{n+1}, q_{n+1}) = gcd(GCD_n, GCD_{n+1}) = GCD2_{n+1}`. There is a constant `C` such that for all `n`:
> (1) if `b` is irreducible over `Q` (`Δ_b` not a square): `P_n ≤ C^n (n!)^{1/2}` — for every `A, k` (the proof does not
>     use the shape of `a_n` at all);
> (2) if `ω ≠ 0` and `D` is not a square in `Q`: `P_n ≤ C^n (n!)^{1/2}`;
> (3) if `ω ≠ 0` and `D = 0`: `P_n ≤ C^n`.

*Proof.* `P_n^2 | GCD_n GCD_{n+1} | Φ_{n+1} = ± ∏_{i ≤ n+1} b_i`, so `v_p(P_n) ≤ (1/2) v_p(∏_{i≤N} b_i)`, `N := n+1`.
If `p | P_n`, then `p | GCD_n, GCD_{n+1}`, and Steps 1–2 of the proof of Theorem C (which use nothing about `a_n`)
applied to the minimal such index show that `b` has two zeros mod `p` at indices `≤ N`; so `p ≤ c_1 N` or `p ∈ S_0`
(Lemma 2 (c)), and `b` has a root mod `p`. Let `Q` be the set of primes `p ∉ S ∪ S_0 ∪ {2}`, `p ∤ B`, that can divide
some `P_n` and, in cases (2) and (3) only, do not divide `ω` (in case (1) `ω` may be 0 and plays no role). The
excluded primes are finitely many in each case. By the last remark, in case (1) `Q ⊂ {p : (Δ_b/p) = 1}`; by
Theorem C, in case (2) `Q ⊂ {p : (D/p) = 1}` and in case (3) `Q = ∅`. By Lemma 2,

    log P_n ≤ Σ_{p ∈ Q, p ≤ c_1 N} ( N/(p−1) + log_p(c_2 N^2) ) log p + O(N) = N Σ_{p ∈ Q, p ≤ c_1 N} (log p)/p + O(N)

(the finitely many other primes contribute `O(N)` by Lemma 2 (a); `Σ_{p ≤ c_1 N} log(c_2N^2) = O(N)`). For a
non-square integer `d` (negative or non-squarefree allowed), `χ := (4d/·)` is a real non-principal Dirichlet
character modulo `4|d|` with `χ(p) = (d/p)` for odd `p ∤ d`; by Mertens' theorem `Σ_{p≤y} (log p)/p = log y + O(1)` and
`Σ_{p≤y} χ(p)(log p)/p = O(1)` (a standard consequence of `L(1,χ) ≠ 0`; e.g. Davenport, *Multiplicative Number
Theory*, ch. 7, or Mertens' theorem for arithmetic progressions), so
`Σ_{p ≤ y, (d/p)=1} (log p)/p = (1/2) log y + O(1)`. Hence `log P_n ≤ (1/2) N log N + O(N)` in cases (1), (2), and
`log P_n = O(N)` in case (3); and `N^{N/2} ≤ C^n (n!)^{1/2}`. ∎

> **Theorem E (Conjecture 1.4 for `GCD2`).** Let `A ≠ 0`, `b_n = B(n − x_1)(n − x_2)` with rational `x_1, x_2`
> that are not positive integers, `B ≠ 0`. Then the following are equivalent:
> (a) `D` is a nonzero rational square or `ω = 0`; equivalently, `a_n = An + k` belongs to family (Eq. 6) or (Eq. 7)
>     of Conjecture 1.4, where in (Eq. 6) `k = 0` is allowed (as in the paper's examples, but not in the printed
>     conjecture, see §6);
> (b) there is `c > 0` with `P_n ≥ c^n n!` for all `n`;
> (c) `limsup_n (1/n) ln(P_n/n!) > −∞`.
> If they fail, `P_n ≤ C^n (n!)^{1/2}` for all `n`. Moreover, if `b` is an irreducible quadratic, (b) and (c) fail for
> every `a_n = An + k`.

*Proof.* (a) ⇒ (b): in the proof of Theorem B, `M_n | e^{n+1} D_n GCD_n` with `M_n | M_{n+1}` and `D_n ≤ C_2^{n+1}`, so
`M_n` divides `e^{n+2} lcm(D_n, D_{n+1}) · gcd(GCD_n, GCD_{n+1})`, whence `P_n ≥ |M_n|/(e^{n+2} C_2^{2n+3}) ≥ c^n n!`
(here `(p_n,q_n) ≠ (0,0)` for all `n` because `Φ_n ≠ 0`). The remaining case of (a), `D = 0 = ω`, is the Euler family
(Theorem A with `d = 1`, `m = −κ`, `κ = A/2`); after inflation by a suitable integer `c` (which makes `H_1, H_2`
integer-valued and multiplies `P_n` by a factor between `c^n` and `c^{n+2}`), `GCD_n = |Π_1(n+1)|/den(σ_n)`, `Π_1(n+1) | Π_1(n+2)`, so
`P_n ≥ |Π_1(n+1)|/lcm(den σ_n, den σ_{n+1}) ≥ c^n n!`. (b) ⇒ (c) is trivial. (c) ⇒ (a): if (a) fails, then `ω ≠ 0` and
`D` is either 0 or not a rational square, and Theorem D (2), (3) gives `P_n ≤ C^n (n!)^{1/2}`. The last statement is
Theorem D (1). ∎

## 6. Discrepancies between the printed statements and what is true

The paper says (§2.1) "Unless stated otherwise, we focus on 'balanced-degree' PCFs" (`deg b = 2 deg a`), and assumes
`b_n ≠ 0` ("to guarantee that the PCF is infinite") only in Appendix B, for Theorems 1 and 2; the conjectures of
§2.3 are printed "for every `b_n`". All examples below are balanced-degree and are checked exactly
(`verify/fr_only_if.py` (c); the referee's `audit2/a2_examples.py`).

0. **Non-terminating: `k = 0` in (Eq. 6).** The conjecture prints the families "for any `k ∈ Q∖{0}`, and
   `m ∈ Q∖{0, ±√|B|}`". But `PCF[n, 2(n+1)(n+2)]` (`A = 1`, `B = 2 = m(A+m)` with `m = 1`, `D = 9`, `k = 0`, `b_n > 0`) has FR
   by Theorem B (i), and is not in the printed (Eq. 6) (`k = 0`) nor in (Eq. 7) (`a_n ≠ k'(2n+4)`); likewise
   `PCF[7n, 8n^2 − 2]`, the `k = 0` member of the paper's own example "for `m = 1` and any `k ∈ Z`". So the printed
   "only if" half of 1.4 is false as stated; the exclusion `k ≠ 0` is evidently meant for (Eq. 7) only (where
   `k = 0` gives `a_n ≡ 0`). All statements in this file use the corrected family. (Found by the referee.)
1. **Terminating: "only if" of Conjecture 1.4.** `PCF[n+1, n(n−1)]` (`b_1 = 0`, value 1): `p_n = q_n = t_n` with
   `t_0 = 1, t_1 = 2`, `t_n = (n+1)t_{n−1} + n(n−1)t_{n−2} ≥ (n+1)!`, so `GCD_n ≥ (n+1)!`: FR. But `A^2 + 4B = 5` is not a
   square (not (Eq. 6)) and `a_n = n+1 ≠ 2kn` (not (Eq. 7)). More generally, if `b_N = 0` for a positive integer `N`,
   then `(p_n, q_n) = t_n (p_{N−1}, q_{N−1})` for `n ≥ N−1` with an integer sequence `t_n` satisfying the same recurrence,
   and `GCD_n ≥ |t_n|`, which typically has factorial size whatever `a_n` is (search of §7: all 83548 terminating
   fractions outside the families that were examined at `n ≤ 640` have `s` bounded).
   *Numerical remark (not proved):* for this fraction `(ln GCD2_{n+1} − ln n!)/n = −2.09, −2.39, −2.80` at
   `n = 50, 100, 199` (referee's computation), so it also appears to contradict Conjecture 3 as printed ("for any
   PCF").
2. **Terminating: "only if" of Conjecture 1.1.** `PCF[n^2+1, (n−1)^2(n^2−2)]`: again `p_n = q_n = t_n ≥ ∏_{i≤n}(i^2+1) ≥ n!^2`,
   so FR with `deg a = 2`, although `b` is not a product of linear factors over `Q`.
3. **Terminating: "if" of Conjecture 1.4.** `PCF[n−3, 2(n−1)(n−3)]` is in family (Eq. 6) (`B = 2`, `m = 1`, `k = −3`), and
   `p_n = q_n = 0` for all `n ≥ 2`; `PCF[n−1, −(n−1)(n−2)]` is in family (Eq. 7) (`k = 1/2`), and `p_n = q_n = 0` for all
   `n ≥ 1`. So `GCD_n = 0`. (Both roots of `b` are positive integers; this is the case excluded in Theorem B.)

Items 1–3 concern terminating fractions and do not touch the intended content of the conjectures; item 0 is a
misprint-level correction of the family (Eq. 6). (Outside the balanced-degree setting the definition degenerates:
for constant `a_n`, `n!^{deg a} = 1` and FR is vacuous; e.g. `PCF[1, (n+1)^2]` has `A = 0`, `D = 4` and `GCD2` of
factorial size. Theorems C–E assume `A ≠ 0` where it matters, as stated.)

## 7. What remains open, and the search

For the paper's definition, FR means `GCD_n ≥ c^n n!` for infinitely many `n`. Theorem D bounds only
`GCD2_{n+1} = gcd(GCD_n, GCD_{n+1})`, the part of `GCD_n` that persists. A prime `p` as in Theorem C can still divide a
single `GCD_n` ("transient" factor: in the notation of the proof, `t_n ≡ 0` with `t_{n±1} ≢ 0`), and `v_p(GCD_n)` can then
be as large as `v_p(∏_{i≤n} b_i) ≈ 2n/(p−1)` as far as we can prove; heuristically such factors contribute `n^{O(1)}` to
`GCD_n`. So, for non-terminating fractions of degree (1,2) with `A ≠ 0`:

* **Proved:** outside the (corrected) families, `GCD2_{n+1} = gcd(GCD_n, GCD_{n+1}) ≤ C^n (n!)^{1/2}`; in particular
  `GCD_n ≥ c^n n!` cannot come from a common factor shared by two consecutive convergent pairs, which is how it
  arises in every known non-terminating case (Theorems A, B; Apéry-type examples).
* **Reduction:** the paper's Conjecture 3 (`GCD2_n ≐ GCD_n`) for these fractions, together with Theorem E, gives
  the "only if" half of 1.4 (corrected family) and, with Theorem D (1), of 1.1 for `deg b = 2`, `deg a = 1`.
* **Not proved:** `limsup (1/n) ln(GCD_n/n!) = −∞` outside the families — i.e. Conjecture 3 in this range.
* **Search (exact integer arithmetic; `search/`).** All `a_n = An + k`, `b_n = Bn^2 + Cn + E` with `1 ≤ A ≤ 12`,
  `|k| ≤ 12`, `1 ≤ |B| ≤ 8`, `|C|, |E| ≤ 8` (1.39 million fractions). With
  `s(n) := max_{n−3 ≤ m ≤ n} (ln GCD_m − ln m!)/m` (a window of 4 consecutive `m`, to catch transient factors),
  evaluated at `n = 40, 80`: about 1.2 million fractions have `s(80) ≤ s(40) − 0.40` and are discarded as
  non-FR; the 186,913 others are re-examined at `n = 160, 320, 640`. Result: every non-terminating fraction in a
  (corrected) family — `k = 0` included in (Eq. 6) — has `s` bounded (16254 cases); every remaining non-terminating fraction outside the families — 23975
  with `b` split over `Q`, 55088 with `b` irreducible — has `s(n)` decreasing steadily (by 0.25–0.7 per doubling
  of `n`; the slowest, `(A,k,B,C,E) = (8,10,−7,1,4)`: −2.04, −2.58, −2.84, −3.08). No counterexample to the "only if" halves was found in the non-terminating range. The observed size is
  `GCD_n ≈ C^n (n!)^θ` with `θ ≈ 1/2` for split `b` outside the families, in agreement with Theorem D being sharp
  (`verify/fr_only_if.py` (d): `(log P_n − (1/2) log n!)/n` is stable for `PCF[3n+1, n^2+n]`).
