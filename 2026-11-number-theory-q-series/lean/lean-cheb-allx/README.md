# lean-cheb-allx: Chebyshev's 1853 assertion at every x > 0, and for every convex weight, under GRH for L(s, χ₋₄) (Lean 4)

A machine-checked proof, in Lean 4 with Mathlib, the public library `zeta23`
(`anthropics/formal-math`), and `monotonicity/lean-allx`, `monotonicity/lean-chebyshev` of this
repository, of:

> **Theorem (`grh_chebyshev_sum_neg`).** Assume that every zero of the Dirichlet L-function
> L(s, χ₋₄) in the strip 0 < Re s < 1 has real part 1/2. Then for **every** real x > 0
>
> C(x) = Σ_p χ₋₄(p) e^{−p/x} < 0,
>
> the sum over all primes p, with χ₋₄(p) = (−1)^{(p−1)/2} for odd p and χ₋₄(2) = 0.

Chebyshev asserted in 1853 that C(x) → −∞; Hardy–Littlewood and Landau showed that this is
equivalent to the Riemann hypothesis for L(s, χ₋₄), and the implication "GRH ⟹ C(x) → −∞" is
formalized in `monotonicity/lean-chebyshev` (`grh_chebyshev_bias`). The statement here is for every
single x. The paper proof is `monotonicity/chebyshev1853/CHEBYSHEV1853.md`, Theorem 1 (first
proof); the Lean proof uses the same idea (second-order partial summation against the one-signed
function Θ₁) in a discrete form, which turns out to need fewer hypotheses; see "How the Lean proof
differs from the paper proof".

More generally (`grh_chi4_convex_weight_sum_neg`, `grh_chi4_convexOn_weight_sum_neg`): under the
same hypothesis, **Σ_p χ₋₄(p) w(p) < 0 for every weight w that is convex on [3, ∞)** (as a sequence
on the integers, or as a real function), with Σ_p w(p) convergent and w(3) > 0. This contains
Fujii's weights e^{−(p/x)^α}, 0 < α ≤ 1 (`grh_fujii_sum_neg`), the weights p^{−σ}, σ > 1
(`grh_chi4_prime_dirichlet_sum_neg`), and the same statements with an extra factor log p
(`grh_chebyshev_log_sum_neg`, `grh_chi4_convex_log_weight_sum_neg`).

Without any hypothesis (`chebyshev_sum_neg_of_le`): **C(x) < 0 for every 0 < x ≤ 250000**.

The Riemann hypothesis for L(s, χ₋₄) is an **explicit hypothesis** of the conditional theorems
(`GRHFor chi4`), never an axiom. There is **no other hypothesis**: no numerical certificate, no
data on zeros, no linear independence, no hypothesis on ζ.

Status: `lake build` succeeds; no `sorry` outside the challenge file; all 23 statements depend
only on `propext`, `Classical.choice`, `Quot.sound`; no `native_decide`; no `axiom`. Nobody has
reviewed this development other than the Lean kernel and the agent that wrote it; the statements
below are short and are what a reader should check.

## The exact Lean statements

From `ChebAllxChallenge.lean` (which imports only Mathlib; 4 definitions, 23 statements; docstrings
omitted here):

```lean
def GRHFor {q : ℕ} [NeZero q] (χ : DirichletCharacter ℂ q) : Prop :=
  ∀ s : ℂ, 0 < s.re → s.re < 1 → χ.LFunction s = 0 → s.re = 1 / 2

def chi4 : DirichletCharacter ℂ 4 := ZMod.χ₄.ringHomComp (Int.castRingHom ℂ)

def chebyshevSum (x : ℝ) : ℝ :=
  ∑' p : Nat.Primes, ((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ) * Real.exp (-((p : ℕ) : ℝ) / x)

def chebyshevLogSum (x : ℝ) : ℝ :=
  ∑' p : Nat.Primes, ((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ) * Real.log (p : ℕ)
    * Real.exp (-((p : ℕ) : ℝ) / x)

open Finset MeasureTheory

-- (1) the headline and its logarithmic companion
theorem grh_chebyshev_sum_neg (hGRH : GRHFor chi4) (x : ℝ) (hx : 0 < x) :
    chebyshevSum x < 0

theorem grh_chebyshev_log_sum_neg (hGRH : GRHFor chi4) (x : ℝ) (hx : 0 < x) :
    chebyshevLogSum x < 0

-- (1') a finite range with no hypothesis at all
theorem chebyshev_sum_neg_of_le (x : ℝ) (hx : 0 < x) (hN : x ≤ 250000) :
    chebyshevSum x < 0

-- (2) every convex weight
theorem grh_chi4_convex_weight_sum_neg (hGRH : GRHFor chi4) (w : ℕ → ℝ)
    (hconv : ∀ n, 3 ≤ n → 2 * w (n + 1) ≤ w n + w (n + 2))
    (hsum : Summable fun p : Nat.Primes => w p) (hw3 : 0 < w 3) :
    ∑' p : Nat.Primes, ((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ) * w p < 0

theorem grh_chi4_convexOn_weight_sum_neg (hGRH : GRHFor chi4) (w : ℝ → ℝ)
    (hconv : ConvexOn ℝ (Set.Ici 3) w)
    (hsum : Summable fun p : Nat.Primes => w (p : ℕ)) (hw3 : 0 < w 3) :
    ∑' p : Nat.Primes, ((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ) * w (p : ℕ) < 0

theorem grh_chi4_convex_log_weight_sum_neg (hGRH : GRHFor chi4) (f : ℕ → ℝ)
    (hconv : ∀ n, 3 ≤ n → 2 * f (n + 1) ≤ f n + f (n + 2))
    (hsum : Summable fun p : Nat.Primes => Real.log (p : ℕ) * f p) (hf3 : 0 < f 3) :
    ∑' p : Nat.Primes, ((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ) * Real.log (p : ℕ) * f p < 0

theorem grh_fujii_sum_neg (hGRH : GRHFor chi4) (x α : ℝ) (hx : 0 < x) (hα0 : 0 < α)
    (hα1 : α ≤ 1) :
    ∑' p : Nat.Primes, ((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ)
      * Real.exp (-(((p : ℕ) : ℝ) / x) ^ α) < 0

theorem fujii_summable (x α : ℝ) (hx : 0 < x) (hα0 : 0 < α) :
    Summable (fun p : Nat.Primes => ((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ)
      * Real.exp (-(((p : ℕ) : ℝ) / x) ^ α))

theorem grh_chi4_prime_dirichlet_sum_neg (hGRH : GRHFor chi4) (σ : ℝ) (hσ : 1 < σ) :
    ∑' p : Nat.Primes, ((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ) * ((p : ℕ) : ℝ) ^ (-σ) < 0

-- (3) the identities behind them (no hypothesis on zeros)
theorem chi4_log_sum_eq_theta1_series (f : ℕ → ℝ)
    (hconv : ∀ n, 3 ≤ n → 2 * f (n + 1) ≤ f n + f (n + 2))
    (hsum : Summable fun p : Nat.Primes => Real.log (p : ℕ) * f p) :
    HasSum (fun k : ℕ =>
        (∑ p ∈ (Icc 0 (k + 4)).filter Nat.Prime,
          -((ZMod.χ₄ (p : ZMod 4) : ℤ) : ℝ) * Real.log p * (((k + 4 : ℕ) : ℝ) - p))
          * (f (k + 3) - 2 * f (k + 4) + f (k + 5)))
      (-∑' p : Nat.Primes, ((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ) * Real.log (p : ℕ) * f p)

theorem chebyshevSum_eq_neg_integral (x : ℝ) (hx : 0 < x) :
    chebyshevSum x = -∫ t in Set.Ioi (2 : ℝ),
      (∑ p ∈ (Icc 0 ⌊t⌋₊).filter Nat.Prime,
        -((ZMod.χ₄ (p : ZMod 4) : ℤ) : ℝ) * Real.log p * (t - p))
        * (Real.exp (-t / x) * (1 / (x ^ 2 * Real.log t) + 2 / (x * t * Real.log t ^ 2)
          + 1 / (t ^ 2 * Real.log t ^ 2) + 2 / (t ^ 2 * Real.log t ^ 3)))

theorem chebyshevSum_integrand_integrableOn (x : ℝ) (hx : 0 < x) :
    IntegrableOn (fun t : ℝ =>
      (∑ p ∈ (Icc 0 ⌊t⌋₊).filter Nat.Prime,
        -((ZMod.χ₄ (p : ZMod 4) : ℤ) : ℝ) * Real.log p * (t - p))
        * (Real.exp (-t / x) * (1 / (x ^ 2 * Real.log t) + 2 / (x * t * Real.log t ^ 2)
          + 1 / (t ^ 2 * Real.log t ^ 2) + 2 / (t ^ 2 * Real.log t ^ 3)))) (Set.Ioi 2)

theorem chebyshev_kernel_hasDerivAt (x t : ℝ) (hx : x ≠ 0) (ht : 1 < t) :
    HasDerivAt (fun t : ℝ => Real.exp (-t / x) / Real.log t)
        (-Real.exp (-t / x) * (1 / (x * Real.log t) + 1 / (t * Real.log t ^ 2))) t ∧
      HasDerivAt
        (fun t : ℝ => -Real.exp (-t / x) * (1 / (x * Real.log t) + 1 / (t * Real.log t ^ 2)))
        (Real.exp (-t / x) * (1 / (x ^ 2 * Real.log t) + 2 / (x * t * Real.log t ^ 2)
          + 1 / (t ^ 2 * Real.log t ^ 2) + 2 / (t ^ 2 * Real.log t ^ 3))) t

theorem chebyshev_kernel_second_deriv_pos (x t : ℝ) (hx : 0 < x) (ht : 1 < t) :
    0 < Real.exp (-t / x) * (1 / (x ^ 2 * Real.log t) + 2 / (x * t * Real.log t ^ 2)
      + 1 / (t ^ 2 * Real.log t ^ 2) + 2 / (t ^ 2 * Real.log t ^ 3))

-- (4) from the sum without logarithms to the sums with logarithms
theorem chebyshevSum_eq_integral_log_sums (x : ℝ) (hx : 0 < x) :
    chebyshevSum x = ∫ u in Set.Ioi (0 : ℝ), ∑' p : Nat.Primes,
      ((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ) * Real.log (p : ℕ)
        * (((p : ℕ) : ℝ) ^ (-u) * Real.exp (-((p : ℕ) : ℝ) / x))

theorem grh_chi4_log_sum_rpow_exp_neg (hGRH : GRHFor chi4) (x u : ℝ) (hx : 0 < x) (hu : 0 ≤ u) :
    ∑' p : Nat.Primes, ((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ) * Real.log (p : ℕ)
      * (((p : ℕ) : ℝ) ^ (-u) * Real.exp (-((p : ℕ) : ℝ) / x)) < 0

-- (5) the two general summation-by-parts identities over the infinite range (pure analysis)
theorem second_order_summation_by_parts_infinite (a f : ℕ → ℝ)
    (hconv : ∀ k, 2 * f (k + 1) ≤ f k + f (k + 2))
    (h0 : Filter.Tendsto f Filter.atTop (nhds 0)) (hs : Summable fun n => |a n| * f n) :
    HasSum (fun k : ℕ => (∑ n ∈ range (k + 1), a n * ((k : ℝ) + 1 - n))
        * (f k - 2 * f (k + 1) + f (k + 2))) (∑' n, a n * f n)

theorem second_order_partial_summation_infinite (a : ℝ) (ha : 0 ≤ a) (c : ℕ → ℝ)
    (hc : ∀ n : ℕ, (n : ℝ) < a → c n = 0) (F F' F'' : ℝ → ℝ)
    (hF : ∀ t, a ≤ t → HasDerivAt F (F' t) t) (hF' : ∀ t, a ≤ t → HasDerivAt F' (F'' t) t)
    (hF'' : ContinuousOn F'' (Set.Ici a)) (hF''0 : ∀ t, a ≤ t → 0 ≤ F'' t)
    (hlim : Filter.Tendsto F Filter.atTop (nhds 0)) (hs : Summable fun n : ℕ => |c n| * F n) :
    IntegrableOn (fun t => (∑ n ∈ Icc 0 ⌊t⌋₊, c n * (t - n)) * F'' t) (Set.Ioi a) ∧
      ∑' n, c n * F n = ∫ t in Set.Ioi a, (∑ n ∈ Icc 0 ⌊t⌋₊, c n * (t - n)) * F'' t

-- (6) conventions
theorem chebyshevSum_summable (x : ℝ) (hx : 0 < x) :
    Summable (fun p : Nat.Primes =>
      ((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ) * Real.exp (-((p : ℕ) : ℝ) / x))

theorem chebyshevLogSum_summable (x : ℝ) (hx : 0 < x) :
    Summable (fun p : Nat.Primes => ((ZMod.χ₄ ((p : ℕ) : ZMod 4) : ℤ) : ℝ) * Real.log (p : ℕ)
      * Real.exp (-((p : ℕ) : ℝ) / x))

theorem chi4_natCast_two : ((ZMod.χ₄ ((2 : ℕ) : ZMod 4) : ℤ) : ℝ) = 0

theorem chi4_natCast_of_odd (n : ℕ) (hn : n % 2 = 1) :
    ((ZMod.χ₄ (n : ZMod 4) : ℤ) : ℝ) = (-1) ^ ((n - 1) / 2)

theorem grhFor_of_dirichlet_zeros_on_line
    (h : ∀ (N : ℕ) [NeZero N] (χ : DirichletCharacter ℂ N) (s : ℂ),
      0 < s.re → s.re < 1 → DirichletCharacter.LFunction χ s = 0 → s.re = 1 / 2)
    {q : ℕ} [NeZero q] (χ : DirichletCharacter ℂ q) : GRHFor χ
```

Reading notes (what a reader must check to trust that these say what the prose says):

- `chi4` is Mathlib's quadratic character `ZMod.χ₄` (values 0, 1, 0, −1 on 0, 1, 2, 3 mod 4) as a
  `ℂ`-valued Dirichlet character; `χ.LFunction` is Mathlib's `DirichletCharacter.LFunction` (the
  analytic continuation of the L-series). `GRHFor`, `chi4`, `chebyshevSum`, `chebyshevLogSum` are
  character for character those of `monotonicity/lean-chebyshev`.
- `∑' p : Nat.Primes, …` is Mathlib's unconditional sum over the primes; it is `0` for a
  non-summable family, so each statement about a `∑'` should be read together with the
  corresponding summability statement: `chebyshevSum_summable`, `chebyshevLogSum_summable`,
  `fujii_summable`, and the hypotheses `hsum` of the general theorems (for `p^{−σ}` summability is
  Mathlib's `Real.summable_nat_rpow`). `Summable` means absolute convergence here.
- In the general theorems only `w(n)` for `n ≥ 3` matters (`χ₋₄(2) = 0`); the hypotheses "w ≥ 0",
  "w nonincreasing", "w(t) → 0" of the paper statement are not assumed because they follow from
  convexity and `Summable` (a convex sequence that is summable along the primes is nonincreasing,
  nonnegative and tends to 0). If `w(3) = 0` the sum is `0`, so `0 < w 3` cannot be dropped.
- `(Icc 0 ⌊t⌋₊).filter Nat.Prime` is the set of primes p ≤ t; the inner sum in
  `chebyshevSum_eq_neg_integral` and `chi4_log_sum_eq_theta1_series` is
  Θ₁(t) = Σ_{p ≤ t} (−χ₋₄(p)) log p · (t − p), the function of `chebyshev_race_theta1_pos` in
  `monotonicity/lean-allx`. `∫ t in Set.Ioi 2, …` is the Lebesgue integral over (2, ∞).
- `(p : ℝ) ^ (-σ)`, `(p / x) ^ α` are real powers (`Real.rpow`).
- `chebyshev_sum_neg_of_le` has no hypothesis on zeros: it rests on the kernel-checked finite
  computation of lean-allx (Θ₁(t) ≥ 0 for 2 ≤ t ≤ 10⁷, by `decide +kernel`) and on real analysis.
- The hypothesis `GRHFor chi4` has exactly the shape that a theorem "all zeros of all Dirichlet
  L-functions in the open critical strip are on the line" would discharge
  (`grhFor_of_dirichlet_zeros_on_line`). Nothing here depends on, or says anything about, any
  announced proof of such a theorem.

## Structure of the proof

Write a(n) = [n prime] (−χ₋₄(n)) log n and Θ₁(m) = Σ_{n ≤ m} a(n)(m − n).

1. **Second-order summation by parts over the infinite range** (`ChebAllxLib/Discrete.lean`;
   Mathlib only). For a real sequence f put Δ²f(k) = f(k) − 2f(k+1) + f(k+2). If Δ²f ≥ 0 and
   f → 0, then f(n) = Σ_{j ≥ 0} (j+1) Δ²f(n+j) for every n (`hasSum_mul_secondDiff_shift`; the
   boundary term K(f(K) − f(K+1)) tends to 0 because it is ≤ 2f(⌊K/2⌋)). Hence for every sequence
   a with Σ|a(n)| f(n) < ∞ the double family a(n)(j+1)Δ²f(n+j) is absolutely summable with row sums
   a(n)f(n); summing it along the antidiagonals n + j = k gives
   Σ_n a(n)f(n) = Σ_k A₁(k+1) Δ²f(k), A₁(m) = Σ_{n<m} a(n)(m − n)
   (`hasSum_secondSum_mul_secondDiff`). No boundary terms, no growth condition on a.
2. **Positivity** (`tsum_mul_pos_of_secondSum_pos`): if all A₁(k+1) > 0 and f(0) > 0 then some
   Δ²f(k) > 0 and the sum is positive.
3. **The input from lean-allx** (`ChebAllxLib/Main.lean`): with the coefficients shifted by 3,
   A₁(k+1) = Θ₁(k+4) (`secondSum_coeff`), which is positive for every k ≥ 0 by
   `Allx.theta1_wlog_pos` (= `chebyshev_race_theta1_pos`, the theorem "Θ₁(t) > 0 for every real
   t > 3 under GRH for L(s, χ₋₄)" of `monotonicity/lean-allx`: a kernel-checked computation up to
   10⁷ and an explicit-formula argument beyond). Only the integers t = 4, 5, 6, … are used.
4. **Convex sequences along the primes** (`tendsto_zero_of_secondDiff_nonneg`,
   `exists_prime_small`): if f is convex on n ≥ 3 and Σ_p log p · f(p) converges, then f(j+3) → 0,
   so f is nonincreasing and nonnegative there. This gives `tsum_log_mul_neg`
   (Σ_p χ₋₄(p) log p · f(p) < 0).
5. **Removing the logarithm** (`tsum_mul_neg`): apply 4 to f(n) = w(n)/log n. The sequence
   1/log(j+3) is convex because log n + log(n+2) ≤ 2 log(n+1) and 1/a + 1/b ≥ 4/(a+b)
   (`inv_log_convex`), and the product of two nonnegative, nonincreasing, convex sequences is
   convex (`secondDiff_mul_nonneg`).
6. **Instances**: e^{−t/x} (`convexOn_exp_neg_div`), e^{−(t/x)^α} for 0 < α ≤ 1 and t^{−σ}
   (`ChebAllxLib/Weights.lean`: e^{−g} is convex when g is concave).
7. **The integral form** (`ChebAllxLib/Integral.lean`, `Weights.lean`), which is the identity of
   the paper proof and is not needed for the sign: for F twice differentiable on [a, ∞) with
   F'' ≥ 0 continuous and F → 0, F(p) = ∫_p^∞ (t − p) F''(t) dt (`integral_Ioi_sub_mul_deriv2`;
   Taylor's formula on [p, T], and (T − p)F'(T) → 0 because −F'(T)(T − p) ≤ 2F((T+p)/2)), then
   Tonelli for the kernels c(n)(t − n)₊F''(t) (`tsum_mul_eq_integral`). For F(t) = e^{−t/x}/log t
   the derivatives are computed in closed form (`hasDerivAt_chebF`, `hasDerivAt_chebF'`,
   `chebF''_pos`), which gives `chebyshevSum_eq_neg_integral`.
8. **From C to the log-weighted sums** (`chebyshevSum_eq_integral_log_sums`): 1/log p = ∫₀^∞ p^{−u} du
   and one exchange of sum and integral.
9. **The unconditional finite range** (`ChebAllxLib/Finite.lean`). Let f(n) = e^{−n/x}/log n,
   M = 10⁷, δ = f(M) − f(M+1), and split f = h + r with h(n) = f(n) − f(M) − δ(M − n) for n ≤ M
   and h(n) = 0 for n > M (f minus its secant through M and M+1, cut off at M). Then h is convex
   on n ≥ 3 (`headPart_convex`), so the series identity of step 3 applies to h, and only
   Θ₁(4), …, Θ₁(M) occur in it; these are ≥ 0 by lean-allx's kernel check
   (`Allx.theta1_wlog_nonneg_of_le`), and the first term is Θ₁(4)·(f(3) − 2f(4) + f(5)) with
   Θ₁(4) = log 3 (`log_three_mul_le_neg_tsum_head`). The remainder is bounded trivially:
   Σ_p log p · r(p) ≤ ((M+1)² + x) e^{−M/x} (`tsum_log_mul_tail_le`). For x ≤ 250000 the first
   quantity is ≥ 0.03 e^{−4/x} (Mathlib's 9-digit bounds for log 2, log 3, log 5) and the second
   is ≤ 1.1·10¹⁴ e^{−39} e^{−4/x} < 0.002 e^{−4/x} (`finite_numeric`).

## How the Lean proof differs from the paper proof

`CHEBYSHEV1853.md` (Theorem 1, first proof) writes −C(x) = ∫_{3}^∞ Θ₁(t) F''(t) dt with
F = e^{−t/x}/log t, F'' > 0, and uses Θ₁(t) > 0 for every real t > 3; the boundary terms of the
partial summation are removed with t·F(t) → 0, t²·F'(t) → 0 and the Chebyshev bound θ(t) ≪ t; the
general statement is for C² weights w ≥ 0, nonincreasing, convex, with Σ_p w(p) < ∞ and
t·w(t) → 0.

The Lean proof of the sign theorems is the discrete analogue (second differences on the integers
instead of F''). Consequences:

- no differentiability of the weight: the theorem holds for every convex weight;
- no boundary terms: the exchange is Tonelli's theorem for a nonnegative kernel, so neither
  t·w(t) → 0 nor any bound for θ(t) is used; "w ≥ 0" and "w nonincreasing" are consequences of
  convexity and summability and are not assumed;
- only Θ₁ at the integers ≥ 4 is used;
- the convexity of w/log is a two-line inequality instead of a computation with derivatives.

The integral identity itself, with the closed form of F'', is formalized as well
(`chebyshevSum_eq_neg_integral`, `chebyshev_kernel_hasDerivAt`, and the general
`second_order_partial_summation_infinite`), again in the Tonelli form (so its hypotheses are
"F'' ≥ 0, F → 0, Σ|c(n)|F(n) < ∞" rather than the decay conditions of the note's Lemma 1, which is
stated there for F not necessarily convex). Together with `chebyshev_race_theta1_pos` and
`chebyshev_kernel_second_deriv_pos` it gives the paper's first proof verbatim; the Lean file does
not spell out that second derivation of the sign, since the discrete one already proves it.

## What is not formalized

- **The converse / the equivalence with GRH** (Theorem 2 of the note: "C(x) ≤ 0 for all large x"
  ⟹ GRH for L(s, χ₋₄), by Landau's theorem). Not attempted, because none of its ingredients is
  available in Mathlib or in zeta23, and each is a separate formalization: (i) Landau's theorem
  for a Mellin integral ∫₁^∞ f(x) x^{−s−1} dx with f of one sign, with analytic continuation only
  to a half-plane (Mathlib has the corresponding statement only for L-series that extend to an
  *entire* function, `LSeries.positive_of_differentiable_of_eqOn`); (ii) the Mellin transform
  ∫₀^∞ C(x) x^{−s−1} dx = Γ(s) Σ_p χ₋₄(p) p^{−s}; (iii) the decomposition
  Σ_p χ(p) p^{−s} = log L(s,χ) − ½ log(ζ(2s)(1 − 2^{−2s})) − H(s) with H analytic in Re s > 1/3,
  and the passage from an analytic logarithm to "no zeros"; (iv) L(σ, χ₋₄) ≠ 0 for real
  1/2 < σ < 1 (which the argument needs as an unconditional input); (v) the reflection of zeros
  by the functional equation. So **only the implication GRH ⟹ (C(x) < 0 for all x > 0) is
  machine-checked here**; "GRHFor chi4 ↔ ∀ x > 0, chebyshevSum x < 0" is not.
- **The quantitative statements** of the note: Theorem 3 (two-sided bounds for the log-weighted
  sum with explicit constants), Theorem 4 and Corollary 4a (−C(x) ≥ 0.88 √x / log x for x ≥ 70),
  Proposition 8. They need explicit prime-counting bounds (Broadbent–Kadiri–Lumley–Ng–Wilk) and
  certified numerics of Γ at the zeros that are not available in Lean. Only signs are proved here;
  the asymptotic statement C(x) → −∞ is `grh_chebyshev_bias` of `monotonicity/lean-chebyshev`.
- **Lemma 5 (Volterra's kernel)** of the note, −C(x) = ∫₀^∞ W((λ + 1/x)^{−1}) μ(λ) dλ. A simpler
  relation of the same kind is formalized instead (`chebyshevSum_eq_integral_log_sums`: the weights
  p^{−u} e^{−p/x}, u ≥ 0, in place of e^{−(λ+1/x)p}).
- Of **Theorem 5** (no hypothesis, finite range) the statement for C is formalized, on the larger
  range 0 < x ≤ 250000 and by a different method (`chebyshev_sum_neg_of_le`; the note certifies
  0 < x ≤ 69 905 by interval arithmetic); the statement for the log-weighted sum W on a finite
  range is not. The range 250000 is what the crude remainder bound gives with Θ₁ known up to 10⁷;
  it is not optimized.
- **Theorem 6** (the character mod 3; lean-allx and lean-chebyshev only treat χ₋₄).
- Fujii's weights for α > 1 (not convex; the note only explores them numerically).
- The official Comparator run (`leanprover/comparator` with `lean4export`, `nanoda`, `landrun`)
  was not performed; `comparator.json` is provided, and `scripts/compare_statements.sh` checks
  that the 23 statements and the 4 definitions elaborate identically in the challenge and solution
  modules.
- No human has read the proofs. The trusted base is: the Lean kernel (including, through
  lean-allx, its GMP-backed evaluation of natural-number arithmetic used by `decide +kernel`),
  Mathlib, and the definitions and statements in `ChebAllxChallenge.lean`.

## Files

| file | lines | content |
|---|---|---|
| `ChebAllxChallenge.lean` | 320 | **trusted**: definitions and the 23 statements with `sorry`; imports only Mathlib (generated) |
| `ChebAllxSolution.lean` | 280 | the same statements, proved by one-line delegation |
| `ChebAllxLib/Discrete.lean` | 287 | discrete second-order summation by parts over the infinite range; convex sequences (Mathlib only) |
| `ChebAllxLib/Integral.lean` | 406 | the integral form; the kernel e^{−t/x}/log t and its derivatives (Mathlib only) |
| `ChebAllxLib/Main.lean` | 389 | link with Θ₁ (lean-allx) and `chebyshevSum` (lean-chebyshev); the theorems for convex weights; the headline |
| `ChebAllxLib/Weights.lean` | 262 | the integral identity for C(x); Fujii's weights, powers; from C to the log-weighted sums |
| `ChebAllxLib/Finite.lean` | 477 | the unconditional range 0 < x ≤ 250000 |
| `scripts/gen_challenge.py`, `challenge_header.txt` | | generator of the challenge file, of `PrintAxioms.lean` and of `comparator.json` |
| `scripts/PrintAxioms.lean`, `compare_statements.sh`, `comparator.json` | | audit |
| `logs/` | | build output, axiom audit, statement comparison |

## Pins and build

| | |
|---|---|
| Lean | `leanprover/lean4:v4.33.0-rc2` (`lean-toolchain`) |
| Mathlib | `51e6992efd06126df61a496bebf8f49482a4e129` |
| zeta23 | `https://github.com/anthropics/formal-math`, commit `11b6c7e3c927df2658f513b6b68ea7f8929c62c5`, `subDir = "zeta23"` (Apache-2.0) |
| lean-allx | `../lean-allx` in this repository (path dependency; `AllxLib`) |
| lean-chebyshev | `../lean-chebyshev` in this repository (path dependency; `ChebLib`) |

```
cd lean/lean-cheb-allx
lake exe cache get                               # Mathlib oleans; never build Mathlib from source
lake build ChebAllxChallenge ChebAllxSolution    # zeta23 modules, ChebLib, AllxLib, then this project
lake env lean scripts/PrintAxioms.lean           # axiom audit, expected output in logs/print-axioms.txt
./scripts/compare_statements.sh                  # challenge vs solution statements
python3 scripts/gen_challenge.py --check         # the challenge file is the generated one
```

Do not run `lake clean` (it removes the Mathlib cache). The first build compiles the needed zeta23
modules, ChebLib and AllxLib (including lean-allx's 1000 kernel checks, about 620 MB per chunk
file); see the lean-allx README for timings. This project's own modules take about 45 seconds
(`logs/build.txt`).

Output of the axiom audit (`logs/print-axioms.txt`), identical for all 23 statements:

```
'grh_chebyshev_sum_neg' depends on axioms: [propext, Classical.choice, Quot.sound]
'grh_chebyshev_log_sum_neg' depends on axioms: [propext, Classical.choice, Quot.sound]
'chebyshev_sum_neg_of_le' depends on axioms: [propext, Classical.choice, Quot.sound]
'grh_chi4_convex_weight_sum_neg' depends on axioms: [propext, Classical.choice, Quot.sound]
...
```

## Provenance

Written by an AI agent (Claude Code) in one session on 2026-11-24; all files were compiled
together and the axiom audit was run on the final build. zeta23 is Apache-2.0.
