# Chebyshev's 1853 assertion, pointwise: the Abel sum is negative for every x > 0 (GRH for χ_{−4} only)

Working notes, 2026-11-24 (work order T8 of `../OPEN-PROBLEMS-TRIAGE.md`). Everything in
`monotonicity/chebyshev1853/` except `REFEREE-REPORT.md` and `referee_checks/` was produced for this note; nothing
outside it was modified. Revised the same day after the referee report (§11, "Review history").

Status in one paragraph. Chebyshev asserted in 1853 that `C(x) = Σ_{p>2} (−1)^{(p−1)/2} e^{−p/x} → −∞`.
Hardy–Littlewood (1917) and Landau (1918) proved that this is *equivalent* to GRH for `L(s,χ_{−4})`. This
note proves, under GRH for that one L-function, the pointwise form: **`C(x) < 0` for every real `x > 0`**,
with the explicit bound **`−C(x) ≥ 0.88 √x / log x` for every `x ≥ 70`** (no constant above `0.885836` is
admissible for all large x: under GRH alone `liminf −C(x) log x/√x` lies in `[0.8858349, 0.8858360]`, Proposition 8),
and conversely "`C(x) ≤ 0` for
all large x" implies GRH(χ_{−4}). Without any hypothesis, `C(x) < 0` is certified for `0 < x ≤ 69 905`.
The sign statement is a three-line corollary of the (in-house refereed) Theorem 1 of `../explicit/EXPLICIT.md`
(`Θ_1 > 0`) and of the convexity of `e^{−t/x}/log t`; a second proof, independent of the analytic part of
that theorem, goes through the explicit formula with the weights `Γ(ρ)`, whose zero sum is `< 4·10^{-4}√x`.
The equivalence with GRH is classical; what is new here is only the all-x statement, the explicit
constant, and the two devices that make the log-free sum accessible without loss (§7). The same is done
for `χ_{−3}`. An in-house adversarial referee found no fatal error and certified the theorems; its four
repairable and nine minor items are all applied (§11). For Fujii's weights `e^{−(p/x)^α}` the all-x statement is proved for `0 < α ≤ 1`, and the
truth for larger α is explored numerically (§6).

---------------------------------------------------------------------------

## 0. Notation

As in `../explicit/EXPLICIT.md` (cited below as **[E]**): `q ∈ {4, 3}`, `χ = χ_{−q}` the nontrivial
character mod q (real, primitive, odd), `s(p) = −χ(p)`, `ℓ` the prime dividing q, `p_0` the least prime
not dividing q (`p_0 = 3` for q = 4, `p_0 = 2` for q = 3).

    C(x)    = Σ_{p∤q} χ(p) e^{−p/x}                  (x > 0;  q = 4: Chebyshev's sum Σ_{p>2} (−1)^{(p−1)/2} e^{−p/x})
    W(y)    = Σ_p s(p) log p · e^{−p/y} = −Σ_p χ(p) log p · e^{−p/y}
    N_C(x)  = −C(x) log x / √x                       (x > 1)
    ψ_e(y,χ)= Σ_n χ(n) Λ(n) e^{−n/y}
    Θ(t)    = Σ_{p≤t} s(p) log p,    Θ_1(t) = ∫_0^t Θ(u) du = Σ_{p≤t} s(p) log p (t − p)
    ϑ_o(u)  = Σ_{p≤u, p∤q} log p
    S_Γ     = Σ_ρ |Γ(ρ)|          (all nontrivial zeros of L(s,χ), with multiplicity)
    μ(λ)    = ∫_0^∞ λ^{u−1}/Γ(u) du   (λ > 0)

Platt–Trudgian and Fujii write the sum with `e^{−xp}`, `x → 0`; replace x by 1/x. "GRH(χ)" means: all
zeros of `L(s,χ)` in `0 < Re s < 1` have real part 1/2. No hypothesis on ζ, no linear independence.

"Certified" below means: decided by Arb ball arithmetic (python-flint) or exact integers, with no
floating-point decision. "Floating point" means IEEE double / x87 long double, with or without an a-priori
rounding budget, as stated each time (§5). Since the revision, the Taylor-grid certificates called
"certified" compute the remainder, the prime tails and every comparison in Arb as well (§5).

---------------------------------------------------------------------------

## 1. Results

Constants (χ_{−4} / χ_{−3}); all are certified enclosures or cited, see §5 and the table in §8:

    S_Γ  ≤ 3.9199328·10^-4 / 1.6531314·10^-5        (under GRH; the first pair of zeros gives 3.9147·10^-4 / 1.6425·10^-5)
    √π/2 = 0.886226925…,   √π/2 − S_Γ ≥ 0.88583493 / 0.88621039
    c_1  = 2.071193 / 2.081942       (y − c_1√y ≤ ϑ_o(y) < y for 0 < y ≤ 2^22;  [E] Lemma 4(a), re-certified)
    c_3  = 2.153716634 / 2.074662259 (|Θ(u)| < c_3 √u for 0 < u ≤ 2^22; new, certified)
    ε_a  = 0.2/log²(2^22) = 8.6007·10^-4,  c_θ = 1 + ε_a      (BKLNW Cor. 11.2, row k = 2)
    U_1  = 2^22 = 4 194 304          (not a prime; the last prime below it is 4 194 301, and the certificates include the endpoint)

**Theorem 1 (pointwise form of Chebyshev's assertion; convex weights).** Assume GRH(χ_{−4}). Then for
every real `x > 0`

    C(x) = Σ_{p>2} (−1)^{(p−1)/2} e^{−p/x} < 0     and     W(x) = −Σ_{p>2} (−1)^{(p−1)/2} log p · e^{−p/x} > 0.

More generally, let `w : [3,∞) → [0,∞)` be C², nonincreasing and convex, not identically 0 on `(3,∞)`,
with `Σ_p w(p) < ∞` and `t·w(t) → 0`. Then `Σ_{p>2} (−1)^{(p−1)/2} w(p) < 0`, and in fact

    −Σ_p χ(p) w(p) = ∫_3^∞ Θ_1(t) (w(t)/log t)'' dt  ≥  ∫_3^∞ Θ_1(t) [ w''(t)/log t + w(t) (1/log t)'' ] dt  >  0.

In particular `Σ_{p>2} (−1)^{(p−1)/2} e^{−(p/x)^α} < 0` for every `x > 0` and every `0 < α ≤ 1` (Fujii's
weights), and `Σ_{p>2} (−1)^{(p−1)/2} p^{−σ} < 0` for every `σ > 1`. (For `χ_{−3}`, Theorem 6, the same holds
with w defined on `[2,∞)`, the integrals starting at 2 and the sums over `p ≠ 3`.)

**Theorem 2 (equivalence).** The following are equivalent:

    (a)  GRH(χ_{−4});
    (b)  C(x) < 0 for every x > 0;
    (b') C(x) ≤ 0 for all sufficiently large x;
    (c)  C(x) → −∞ as x → ∞                                   (Chebyshev's assertion);
    (d)  W(y) > 0 for every y > 0;
    (d') W(y) ≥ 0 for all sufficiently large y.

(a)⇔(c) is the theorem of Hardy–Littlewood and Landau; according to Suzuki (arXiv:2411.07436v3, §1.2, read)
they proved the same for the log-weighted sum (`W(y) → +∞`). New: (b) and (d) (the all-x forms). The
implications (b')⟹(a), (d')⟹(a) are Landau's argument (§2.7), written out. Close relatives of (b'), i.e.
"one sign for all `x ≥ x_0` ⇔ GRH(χ_4)" for other weights: Suzuki's Theorem 3(2) (weight
`Λ(n)χ_4(n) n^{−1/2} log(x/n)`, `n ≤ x`), and Knapowski–Turán for the weight `exp(−log²(p/x))` (as quoted by
Suzuki; not read). Neither contains a statement for all `x > 0`.

**Theorem 3 (the log-weighted sum, two-sided, explicit).** Assume GRH(χ), `χ ∈ {χ_{−4}, χ_{−3}}`. For
every `y ≥ 4`

    κ(y) √y  ≤  W(y)  ≤  κ^up(y) √y,

where

    κ(y) √y    = (√π/2 − S_Γ) √y − E(y) − c_1 Γ(5/4) y^{1/4} + max{0, (1−ε_a) Γ(5/4) y^{1/4} − c_1 Γ(9/8) y^{1/8}}
                 − Σ_{k odd, k≥3, p_0^k<y} B_k(y) − 1,
    κ^up(y) √y = (√π/2 + S_Γ) √y + E(y) + b + 0.2555 y^{-2} + Σ_{k odd, k≥3, p_0^k<y} B_k(y)
                 + Σ_{k even, k≥4, p_0^k<y} c_θ Γ(1+1/k) y^{1/k} + 2,
    B_k(y)     = min{ c_θ Γ(1+1/k) y^{1/k},  c_3 Γ(1+1/(2k)) y^{1/(2k)} + c_θ U_1 e^{−U_1^k/y} (1 + y/(k U_1^k)) },
    E(y)       = √y Σ_i ε_i [ F(U_i/√y) − F(U_{i+1}/√y) ],     F(v) = v e^{−v²} + (√π/2) erfc(v),

`b = (L'/L)(0,χ) = 0.78319 / 0.94820`, and `(U_i, ε_i)` is the ladder `U_1 = 2^22, ε_1 = ε_a`;
`U_2 = 89 967 803`, and `U_i ∈ {10^10, 10^12, 10^15, 10^20, 10^30, 10^50, 10^100, 10^200, 10^400}` with
`ε_i = 1/log³ U_i` for `i ≥ 2` (BKLNW Cor. 11.2, rows k = 2 and k = 3; `F(U_{last+1}) := 0`). The ladder may be
continued by any further points `U > 10^400` with `ε = 1/log³ U` (the inequalities remain true). `E(y)/√y` is
`< 10^-10` for `y ≤ 10^12` and `≤ ε_a √π/2 = 7.622·10^-4` for all y; for the finite ladder as printed its limit
is `√π/(2 log³ 10^400) = 1.13·10^-9`, and it tends to 0 when the ladder is continued (§2.6).
Consequences: `κ(y) > 0` on the grid `10^{j/200}` from `y = 412.1` (χ_{−4}) resp. `749.9` (χ_{−3}) on (a grid
statement, not used; positivity is proved for `y ≥ 69 905` in §2.9, which is all that is needed); and

    √π/2 − S_Γ  ≤  liminf W(y)/√y  ≤  limsup W(y)/√y  ≤  √π/2 + S_Γ          (0.88583 … 0.88662 for χ_{−4}).

**Theorem 4 (the sum without logarithms, explicit).** Assume GRH(χ). Let `κ_*(y) = max{0, inf_{y'≥y} κ(y')}`.
For `x > 1`, integers `1 = V_0 < V_1 < … < V_n` with `x/(1+V_n) ≥ 4`, and `A_n = log(x/V_n) > 0`,

    N_C(x) ≥ κ_*(x/2) (1 − c_a/log x) + Σ_{i=1}^n κ_*(x/(1+V_i)) [A(V_i) − A(V_{i−1})] log x / log²(x/V_{i−1})
             − √π √(1+V_n) (1 + A_n) log x / (A_n² x),

with `c_a = 2 log((1+√2)/2) = 0.37645`, `A(V) = −log((√(1+V)+1)/(√(1+V)−1))` (so `A(1) = −2 log(1+√2)`,
`A(∞) = 0`). For n → ∞ and x → ∞ the right side is `κ_*·(1 + (2 log 2 + o(1))/log x)`; the true expansion
is `(√π/2)(1 + (γ + 2 log 2)/log x + …)` up to the oscillating term of size `S_Γ`.

**Corollary 4a (uniform constant).** Assume GRH(χ_{−4}). Then

    −C(x) = −Σ_{p>2} (−1)^{(p−1)/2} e^{−p/x}  ≥  0.88 · √x / log x      for every x ≥ 70,

and `liminf_{x→∞} N_C(x) ≥ √π/2 − S_Γ ≥ 0.885834`. (`N_C(68.5) < 0.88`: the crossing is at `x ≈ 68.86`.
The same computation gives the constant `0.883` for `x ≥ 4·10^6`; by Proposition 8 no constant above `0.885836`
is admissible for all large x.) In particular `C(x) → −∞`.

**Proposition 8 (asymptotic shape).** Assume GRH(χ), and let `Z(y) = Σ_ρ Γ(ρ) y^{iγ}` (real, `|Z| ≤ S_Γ`). Then

    W(y)/√y = √π/2 + Z(y) + o(1)   (y → ∞),          N_C(x) = √π/2 + Z(x) + o(1)   (x → ∞).

Consequently, under GRH(χ_{−4}) alone (`ρ_1` the lowest zero, `2|Γ(ρ_1)| ≥ 3.9147116·10^-4`):

    0.8858349 ≤ √π/2 − S_Γ ≤ liminf N_C ≤ √π/2 + S_Γ − 4|Γ(ρ_1)| ≤ 0.8858360,
    0.8866178 ≤ √π/2 − S_Γ + 4|Γ(ρ_1)| ≤ limsup N_C ≤ √π/2 + S_Γ ≤ 0.8866190

(χ_{−3}: `liminf N_C ∈ [0.88621039, 0.88621061]`); the same holds for `W(y)/√y`. In particular `N_C` has no
limit, and with any constant above `0.885836` the inequality of Corollary 4a fails at some arbitrarily large x. If in addition the
positive ordinates are linearly independent over Q (LI), then `liminf N_C = √π/2 − S_Γ` and
`limsup N_C = √π/2 + S_Γ` exactly. (Proof in §2.11; added in the revision, replacing a sketch.)

**Theorem 5 (no hypothesis; finite range).** `C(x) < 0` and `W(x) > 0` for every `0 < x ≤ 69 905`
(= 2^22/60), for both characters; certified.

**Theorem 6 (χ_{−3}).** Theorems 1–5 hold for `χ_{−3}`, i.e. for
`C_3(x) = Σ_{p≠3} χ_{−3}(p) e^{−p/x}` (classes 1 against 2 mod 3; `3` replaced by `p_0 = 2` in the
integrals; the weights w of Theorem 1 are defined, nonincreasing and convex on `[2,∞)`), under GRH(χ_{−3}): `C_3(x) < 0` for every `x > 0`, the equivalence of Theorem 2, the bounds of
Theorems 3 and 4 with the χ_{−3} constants, and

    −C_3(x) ≥ 0.88 · √x / log x      for every x ≥ 14            (crossing at x ≈ 13.43).

**Corollary 7 (Chebyshev's assertion, given the announced proof of GRH).** A Lean 4 formalisation of a
proof of GRH for Dirichlet L-functions was announced on 2026-11-20 (see `../GRH-PROOF-DIGEST.md`, which
records that the statement was matched against Mathlib and the proof accepted by proof checkers, and that
no human-readable or refereed proof exists as of that digest; I have not examined the announcement
myself). The digest also records that the final default-kernel stage of `lake comparator` was still running
and that a PhilArchive paper alleging that the formal statement is vacuous exists and was not read.
*If* that proof is correct, then GRH(χ_{−4}) holds, hence by Hardy–Littlewood–Landau Chebyshev's
assertion of 1853 is true, and by Theorem 1 and Corollary 4a it is true in the pointwise form: `C(x) < 0` for
every `x > 0`, and `−C(x) ≥ 0.88 √x/log x` for `x ≥ 70`. Nothing here is claimed unconditionally: **GRH(χ)
remains a hypothesis of every theorem of this note**, and the first implication is the classical one —
the only thing we did not find written down is the remark itself (§7).

Numerical truth (floating point, §3; not part of any theorem): `C(x) < 0` at all 473 grid points with
`0.06 ≤ x ≤ 10^9` (grid `10^{j/50}` up to `10^8`, `10^{j/10}` up to `10^9`), and, with the rounding budget of §5,
on every real x in `[0.103, 1.02·10^8]`; C is strictly decreasing along
the grid; `N_C` rises from 0.88 at x ≈ 68.9 to its maximum `1.1087617` at `x ≈ 6356` (grid maximum `1.108761` at
`x = 6309.6`) and then decreases
slowly (`1.0129` at `10^8`, `0.9956` at `10^9`), consistent with `(√π/2)(1 + 1.9635/log x)`;
`W(y)/√y` rises to `0.88467` at `10^9`, with the oscillation of the first zero (amplitude `3.9·10^-4`) visible.

---------------------------------------------------------------------------

## 2. Proofs

### 2.1 Lemma 1 (kernel identity; no hypothesis)

Let `F ∈ C²[p_0, ∞)` with `t F(t) → 0` and `t² F'(t) → 0`. Then

    lim_{T→∞} Σ_{p≤T} s(p) log p · F(p) = ∫_{p_0}^∞ Θ_1(t) F''(t) dt

whenever one side converges.

*Proof.* [E] Lemma 2 on `[p_0, T]` gives `Σ_{p≤T} s(p) log p F(p) = F(T)Θ(T) − F'(T)Θ_1(T) + ∫_{p_0}^T Θ_1 F''`.
Unconditionally `|Θ(T)| ≤ ϑ_o(T) ≤ c_θ T` and `|Θ_1(T)| ≤ c_θ T²/2` (§2.5), so the boundary terms tend to 0. ∎

Two instances. (i) `F(t) = e^{−t/y}`: `W(y) = y^{-2} ∫_{p_0}^∞ Θ_1(t) e^{−t/y} dt`.
(ii) `F(t) = e^{−t/x}/log t`:

    −C(x) = ∫_{p_0}^∞ Θ_1(t) F''(t) dt,
    F''(t) = e^{−t/x} [ 1/(x² log t) + 2/(x t log² t) + 1/(t² log² t) + 2/(t² log³ t) ]  >  0     (t > 1).

(Both checked numerically to 15 digits at x = 1, 7.5, 40: `identities.py`, item 4.)

### 2.2 Proof of Theorem 1

By [E] Theorem 1 (GRH(χ_{−4})), `Θ_1(t) > 0` for every `t > 3`. Let w be as in the theorem and
`F = w/log t`. Then

    F'' = w''/log t + 2 w' · (1/log t)' + w · (1/log t)''

and each of the three terms is `≥ 0` on `(1,∞)`: `w'' ≥ 0`; `w' ≤ 0` and `(1/log t)' < 0`; `w ≥ 0` and
`(1/log t)'' = 1/(t² log² t) + 2/(t² log³ t) > 0`. The hypotheses of Lemma 1 hold: `tF → 0`; and by
convexity `|w'(T)| ≤ 2 w(T/2)/T`, so `T² F'(T) → 0`. The left side of Lemma 1 converges absolutely. Hence
`−Σ χ(p) w(p) = ∫ Θ_1 F'' ≥ ∫ Θ_1 (w''/log t + w (1/log t)'')`, and this is `> 0` because w is positive
on an interval `(3, 3+δ)`. For `w(t) = e^{−(t/x)^α}`: `w'' = e^{−(t/x)^α} α t^{α−2} x^{−α} [α (t/x)^α − (α−1)] ≥ 0`
when `α ≤ 1`. The statement for W is Lemma 1(i). ∎

*Remark.* This is the same mechanism as [E] Theorem 5 (the cone of weights): `Θ_1` is one-signed, so
every nonnegative, nonincreasing, convex C² weight w with `Σ_p w(p) < ∞` and `t·w(t) → 0` inherits the sign,
with no loss of a logarithm. The two decay hypotheses are needed: `Σ_p w(p) < ∞` is not implied by the others
(`w = 1/(t log log t)`), and without them the sum need not converge.
For `w = t^{−σ}` with `1/2 < σ ≤ 1` the same argument works under GRH if one uses `Θ(T) ≪ √T log² T` and
`Θ_1(T) ≪ T^{3/2}` for the boundary terms; I have not written this out.

### 2.3 Lemma 3 (explicit formula with the weights Γ(ρ))

Let χ be primitive and odd, `b = (L'/L)(0,χ)`, `b_0` as in [E] Lemma 3. For every `y > 0`

    ψ_e(y,χ) = − Σ_ρ Γ(ρ) y^ρ − b (1 + 1/y) e^{−1/y} + [ (1 + b_0) e^{−1/y} + E_1(1/y) ] / y + R_1(y) − R_4(y),

    R_1(y) = y^{-2} Σ_ρ (ρ(ρ+1))^{-1} ∫_0^1 t^{ρ+1} e^{−t/y} dt,       |R_1(y)| ≤ (2/5) σ_0² y^{-2}   under GRH(χ),
    R_4(y) = y^{-2} ∫_1^∞ Σ_{n≥1} t^{−2n}/(2n(2n+1)) · e^{−t/y} dt,      0 ≤ R_4(y) ≤ (log 2 − 1/2) y^{-2},

`E_1` the exponential integral, `σ_0² = Σ_ρ |ρ|^{-2}` ([E] Lemma 3(b): 0.155568 / 0.113230). All sums over ρ
converge absolutely. Under GRH(χ), `|Γ(1/2 + iγ)| = √(π / cosh πγ)`, so `|Σ_ρ Γ(ρ) y^ρ| ≤ S_Γ √y`.

*Proof.* By partial summation twice, `ψ_e(y,χ) = y^{-2} ∫_1^∞ ψ_1(t,χ) e^{−t/y} dt` (`ψ_1(1,χ) = 0`).
Insert the explicit formula [E] Lemma 3(a), valid for `t ≥ 1`,
`ψ_1(t,χ) = −Σ_ρ t^{ρ+1}/(ρ(ρ+1)) − b t + log t + 1 + b_0 − Σ_{n≥1} t^{−2n}/(2n(2n+1))`, and integrate term by
term (the zero sum converges absolutely and uniformly on compacta, with `Σ |ρ(ρ+1)|^{-1} ∫ t^{Re ρ+1} e^{−t/y} < ∞`):
`y^{-2} ∫_0^∞ t^{ρ+1} e^{−t/y} dt = Γ(ρ+2) y^ρ` and `Γ(ρ+2)/(ρ(ρ+1)) = Γ(ρ)`; the piece `∫_0^1` is `R_1`,
bounded using `|t^{ρ+1}| = t^{3/2}`, `∫_0^1 t^{3/2} = 2/5` and `|ρ+1| > |ρ|`;
`y^{-2} ∫_1^∞ t e^{−t/y} dt = (1 + 1/y) e^{−1/y}`; `y^{-2} ∫_1^∞ e^{−t/y} dt = e^{−1/y}/y`;
`y^{-2} ∫_1^∞ log t · e^{−t/y} dt = E_1(1/y)/y`; and
`∫_1^∞ Σ_n t^{−2n}/(2n(2n+1)) dt = Σ_n 1/((2n−1)2n(2n+1)) = log 2 − 1/2`. ∎

This is the Mellin-transform formula `ψ_e = (1/2πi) ∫ (−L'/L)(s) Γ(s) y^s ds` with the contour moved to the
left; deriving it from [E] Lemma 3(a) avoids a new contour estimate. Numerical check (`identities.py`,
item 1, 511 zeros): direct sum minus formula is `−3.3·10^-9, −1.6·10^-9, −3.6·10^-10, −2.7·10^-11, −1.1·10^-12,
−1.2·10^-14` at `y = 1, 2, 5, 20, 100, 1000` (the residual `≈ 3·10^-9/y²` is the truncation of `R_1` at 511
zeros). A wrong `b` would show as an error of order 1, a wrong `b_0` as an error of order 1/y.

Used below, for `y ≥ 4` (with `(1+u)e^{−u} ≥ 1 − u²/2`, `0 < E_1(u) < e^{−u} log(1 + 1/u)`, `1 + b_0 > 0`):

    −ψ_e(y,χ) ≥ −S_Γ √y + b − r^+(y),    r^+(y) = b/(2y²) + (1 + b_0 + log(1+y))/y + (2/5)σ_0²/y²,
    −ψ_e(y,χ) ≤  S_Γ √y + b + 0.2555/y².

`r^+` is decreasing and `b − r^+(4) = 0.1175 / 0.2274 > 0`, so **`−ψ_e(y,χ) ≥ −S_Γ √y` for `y ≥ 4`**.

### 2.4 Lemma 4 (the zero sum S_Γ; GRH(χ))

`S_Γ ≤ 3.9199328·10^-4` (χ_{−4}) and `≤ 1.6531314·10^-5` (χ_{−3}).

*Proof* (`cert_abel.py`, part B; Arb, 160 bits). (1) For each of the 511 (537) ordinates `γ_j` listed in
`prototype/explore/chi{4,3}_zeros.txt`, the Hardy Z-function changes sign on `(γ_j − 10^-20, γ_j + 10^-20)`
(certified), and these intervals are disjoint; so there are zeros `1/2 ± iγ_j^*` with `γ_j^*` in these
intervals. (2) Under GRH, `Σ_ρ |ρ|^{-2} = σ_0²` exactly ([E] Lemma 3(b); closed form
`log(q/π) − γ + 2(L'/L)(1,χ)`). Hence the zeros *not* accounted for in (1) — unlisted zeros and extra
multiplicities — satisfy `Σ' |ρ|^{-2} = d := σ_0² − Σ_j 2/(1/4 + γ_j^{*2}) = 0.0034836` (χ_{−3}: `0.0030943`),
certified `> 0`. They come in conjugate pairs (`L(1/2,χ) ≠ 0`), so each has `2/(1/4+γ²) ≤ d`, i.e.
`|γ| ≥ g := √(2/d − 1/4) = 23.955` (χ_{−3}: `25.418`). (3) `t ↦ (1/4 + t²) √(π/cosh πt)` is decreasing for
`t ≥ 2`, so `Σ' |Γ(ρ)| ≤ d (1/4 + g²) √(π/cosh πg) < 2.3·10^-16`. (4) `2 Σ_j √(π/cosh πγ_j^*) ≤ 3.91993271295·10^-4`
(χ_{−3}: `1.65313139265·10^-5`). ∎

So completeness of the zero list is not assumed: it is replaced by the identity for `σ_0²`. (That the
listed zeros are all the zeros up to height `23.955` (χ_{−3}: `25.418`), and are simple there, follows from (2),
given GRH. Independently, Platt–Trudgian report a Turing-method verification of all zeros of `L(s,χ_4)` up to
height 1127; that is not used here.)

### 2.5 Inputs on primes (no hypothesis)

(P1) `u − c_1 √u ≤ ϑ_o(u) < u` for `0 < u ≤ U_1 = 2^22` ([E] Lemma 4(a); re-certified in `cert_abel.py`, A2,
at every prime, and A0 at the endpoint: `2^22` is not prime, the last prime is `4 194 301`, and
`(u − ϑ_o(u))/√u` increases on the last segment up to `(U_1 − ϑ_o(U_1))/√U_1 = 1.01442 / 1.01462 < c_1`).
(P2) `|θ(u) − u| < 0.2 u/log² u` for `u ≥ 3 594 641` and `< u/log³ u` for `u ≥ 89 967 803`: Broadbent–Kadiri–
Lumley–Ng–Wilk, arXiv:2002.11068v2, Corollary 11.2 (eq. (3.39)), rows `(k, M_k, X_0) = (2, 0.2, 3 594 641)`
and `(3, 1, 89 967 803)`; **re-read for this note** in the arXiv PDF (the table and the proof paragraph for
each row; their input is the RH verification height `2 445 999 556 030`, which their Appendix A.1.1 says is
confirmed by Platt–Trudgian's `3 000 175 332 800`). (P3) `|Θ(u)| < c_3 √u` for `0 < u ≤ U_1`
(`cert_abel.py`, A3; the sup is `2.15371663398` at `u = 126 227` for χ_{−4}, `2.07466225820` at
`u = 1 419 839` for χ_{−3}; `|Θ(u)|/√u` decreases between primes, so the last segment needs nothing). Consequences, for all `u > 0`, with `ε̄(u) = 0` for `u ≤ U_1` and
`ε̄(u) = ε_i` on `[U_i, U_{i+1})`:

    ϑ_o(u) ≤ (1 + ε̄(u)) u ≤ c_θ u,      ϑ_o(u) ≥ u − c_1 √u − ε̄(u) u,      |Θ(u)| ≤ min(c_3 √u [u ≤ U_1], c_θ u).

(For `u > U_1`: `ϑ_o = θ − log ℓ ≥ u − ε̄u − log ℓ` and `log ℓ ≤ c_1 √u`.)

For a nonincreasing C¹ function g with `u g(u) → 0`: `Σ_{p∤q} log p · g(p) = ∫_0^∞ ϑ_o(u)(−g'(u)) du` and
`Σ_p s(p) log p · g(p) = ∫_0^∞ Θ(u)(−g'(u)) du`. With `g_k(u) = e^{−u^k/y}`:
`∫_0^∞ g_k = Γ(1+1/k) y^{1/k}`, `∫_0^∞ √u (−g_k') = Γ(1+1/(2k)) y^{1/(2k)}`.

### 2.6 Proof of Theorem 3

Since `χ(p^k) = χ(p)` for odd k and `= 1` (p ∤ q) for even k,

    W(y) = −ψ_e(y,χ) + Q_2(y) − Σ_{k≥3 odd} W_k(y) + Σ_{k≥4 even} Q_k(y),
    Q_k(y) = Σ_{p∤q} log p · e^{−p^k/y} ≥ 0,     W_k(y) = Σ_p s(p) log p · e^{−p^k/y}.

*Zeros:* §2.3. *Squares:* by §2.5,

    Q_2(y) ≥ ∫_0^∞ (u − c_1√u)(−g_2') du − ∫_{U_1}^∞ ε̄(u) u (−g_2'(u)) du = (√π/2) √y − c_1 Γ(5/4) y^{1/4} − E(y),
    Q_2(y) ≤ (√π/2) √y + E(y),

because `∫_a^b u(−g_2') du = a g_2(a) − b g_2(b) + ∫_a^b g_2 = √y [F(a/√y) − F(b/√y)]`. F decreases from `√π/2`
to 0 (`F' = −2v² e^{−v²}`), so `E(y) ≤ ε_a F(U_1/√y) √y`. *Fourth powers:* `Q_4(y) ≥ max{0, (1−ε_a) Γ(5/4) y^{1/4} − c_1 Γ(9/8) y^{1/8}}`
(from `ϑ_o(u) ≥ (1−ε_a)u − c_1√u`); the other even powers are dropped in the lower bound, and for the
upper bound `Q_k(y) ≤ c_θ Γ(1+1/k) y^{1/k}`. *Odd powers with `p_0^k < y`:*
`|W_k(y)| ≤ ∫ |Θ|(−g_k')`, which is at most `c_θ Γ(1+1/k) y^{1/k}` (trivially) and at most
`c_3 Γ(1+1/(2k)) y^{1/(2k)} + c_θ [U_1 g_k(U_1) + ∫_{U_1}^∞ g_k]` (by (P3)), where
`∫_{U_1}^∞ e^{−u^k/y} du ≤ ∫_{U_1}^∞ (u/U_1)^{k−1} e^{−u^k/y} du = y e^{−U_1^k/y}/(k U_1^{k−1})`: this is `B_k(y)`.
*Powers with `p_0^k ≥ y`:* let `k_1 ≥ 3` be the least such k of the given parity. Since `e^{−v} ≤ 1/(e v)`,

    Σ_{k≥k_1, same parity} Σ_p log p · e^{−p^k/y} ≤ (y/e) Σ_p log p · p^{−k_1}/(1 − p^{−2})
        ≤ (1/e) Σ_{p∤q} log p · (p_0/p)³/(1 − p^{−2}) ≤ 0.695 (q = 4), 0.411 (q = 3)  ≤ 1

(`bounds_abel.py` prints the constant). Adding up gives the two inequalities, for `y ≥ 4`.

*Positivity threshold and limits.* `bounds_abel.py`, item 1: `κ(y) > 0` on the grid `10^{j/200}` from
`y = 412.1` (χ_{−4}), `749.9` (χ_{−3}) on; that κ stays positive for all larger y follows from the
minorant `κ^♭` of §2.9 (positive at 69 905 and beyond; below 69 905 positivity of W is Theorem 5 anyway).
*Limits.* For any `U ≥ 89 967 803` the ladder may be cut at U with the single value `1/log³ U` beyond, which
gives `E(y)/√y ≤ ε_a [F(U_1/√y) − F(U/√y)] + (√π/2)/log³ U`; the bracket tends to 0 as `y → ∞`, so
`limsup E(y)/√y ≤ (√π/2)/log³ U` for every U, i.e. the E-term disappears in the limit. All other correction
terms are `O(y^{-1/6} log y)` relative to `√y`. This gives the liminf/limsup. ∎

### 2.7 Proof of Theorem 2 (Landau's argument)

(a)⟹(b),(d): Theorem 1. (a)⟹(c): Corollary 4a (or Theorem 3 and §2.8). (b)⟹(b'), (c)⟹(b'), (d)⟹(d')
are trivial. It remains to prove (b')⟹(a) and (d')⟹(a).

*Step 1 (Mellin transforms).* For `Re s > 1`, by absolute convergence (`|C(x)| ≤ Σ_p e^{−p/x} ≤ x`),

    ∫_0^∞ C(x) x^{−s−1} dx = Σ_p χ(p) ∫_0^∞ e^{−p/x} x^{−s−1} dx = Γ(s) P(s),    P(s) = Σ_p χ(p) p^{−s},
    ∫_0^∞ W(y) y^{−s−1} dy = Γ(s) Z(s),                                           Z(s) = −Σ_p χ(p) log p · p^{−s}.

As `x → 0`, `|C(x)| ≤ Σ_{n≥2} e^{−n/x}` and `|W(x)| ≤ Σ_{n≥2} log n · e^{−n/x}` are `O(e^{−2/x})`, so `∫_0^1 C(x) x^{−s−1} dx` and the
analogue for W are entire functions of s. Hence `F_C(s) := ∫_1^∞ C(x) x^{−s−1} dx = Γ(s)P(s) − (entire)` and
`F_W(s) := ∫_1^∞ W(y) y^{−s−1} dy = Γ(s)Z(s) − (entire)`, for `Re s > 1`.

*Step 2 (the functions P and Z).* From the Euler product, for `Re s > 1`,

    P(s) = log L(s,χ) − (1/2) log[ζ(2s)(1 − ℓ^{−2s})] − H(s),      Z(s) = (L'/L)(s,χ) + Σ_{p∤q} log p · p^{−2s} + H_1(s),

with `H(s) = Σ_p Σ_{k≥3} c_k(p) p^{−ks}` and `H_1(s) = Σ_p Σ_{k≥3} χ(p)^k log p · p^{−ks}` analytic in
`Re s > 1/3` (bounded coefficients `|c_k(p)| ≤ 1`), and
`Σ_{p∤q} log p · p^{−2s} = −(ζ'/ζ)(2s) − Σ_p log p/(p^{2s}(p^{2s}−1)) − log ℓ · ℓ^{−2s}`, which is analytic in
`Re s > 1/2` except for the poles at the points where `ζ(2s) = 0`; there are none with `Re s > 1/2`.

*Step 3 (the real axis).* For real `σ > 1/2`: `Γ(σ)` is analytic; `ζ(2σ)(1 − ℓ^{−2σ}) > 0`; and
`L(σ,χ) > 0`, because `L(σ,χ_{−4}) = Σ_{n≥0} ((4n+1)^{−σ} − (4n+3)^{−σ})` and
`L(σ,χ_{−3}) = Σ_{n≥0} ((3n+1)^{−σ} − (3n+2)^{−σ})` are series of positive terms for every `σ > 0`. Hence
`Γ(s)P(s)` and `Γ(s)Z(s)` continue analytically to a neighbourhood of every real point `σ > 1/2`. (This is
where the absence of real zeros in `(1/2, 1)` is used; exactly as in [E] Theorem 6 it cannot be dropped for
a general real character.)

*Step 4 (Landau).* Suppose `C(x) ≤ 0` for `x ≥ x_0` (the case `C ≥ 0` is identical). By Landau's theorem
(the abscissa of convergence of `∫_1^∞ f(x) x^{−s−1} dx`, for f of constant sign near infinity, is a
singularity of the function it represents; Montgomery–Vaughan, *Multiplicative Number Theory I*, Lemma
15.1) and Step 3, the abscissa of convergence of `F_C` is `≤ 1/2`. So `F_C`, and with it
`P(s) = (F_C(s) + entire)/Γ(s)` (`1/Γ` is entire), is analytic in `Re s > 1/2`. Then
`log L(s,χ) = P(s) + (1/2) log[ζ(2s)(1−ℓ^{−2s})] + H(s)` has an analytic continuation to `Re s > 1/2`
(`ζ(2s) ≠ 0` there), and its exponential is `L(s,χ)`, which therefore has no zero with `Re s > 1/2`; by the
functional equation there is none with `0 < Re s < 1/2`. For W: the same argument makes `Z(s)`, hence
`(L'/L)(s,χ) = Z(s) − Σ log p · p^{−2s} − H_1(s)`, analytic in `Re s > 1/2`, so L has no zero there. ∎

*Remark (unconditional).* It follows that `C(x) ≥ 0` for all large x is impossible: it would imply GRH
(same proof), which implies `C(x) < 0`. So unconditionally `C(x) < 0` for arbitrarily large x.

### 2.8 Lemma 5 (removing the logarithm: Volterra's kernel) and proof of Theorem 4

**Lemma 5.** (i) `μ(λ) = ∫_0^∞ λ^{u−1}/Γ(u) du` is finite and positive for `λ > 0`, and for every `t > 1`

    1/log t = ∫_0^∞ e^{−λt} μ(λ) dλ.

(ii) For `0 < λ < 1`: `μ(λ) ≥ (1 − λ + λ log λ)/(λ log² λ) > 0`.
(iii) (no hypothesis) For every `x > 0`:   `−C(x) = ∫_0^∞ W( (λ + 1/x)^{-1} ) μ(λ) dλ`.
More generally, for every `α > 0`, with `C_α(x) = Σ χ(p) e^{−(p/x)^α}` and `W_α(y) = Σ s(p) log p · e^{−(p/y)^α}`:
`−C_α(x) = α ∫_0^∞ W_α( (λ + x^{−α})^{−1/α} ) μ(λ) dλ`.

*Proof.* (i) By Tonelli, `∫_0^∞ e^{−λt} μ(λ) dλ = ∫_0^∞ (1/Γ(u)) ∫_0^∞ λ^{u−1} e^{−λt} dλ du = ∫_0^∞ t^{−u} du = 1/log t`.
(ii) `1/Γ(u) = u/Γ(1+u) ≥ u` for `0 ≤ u ≤ 1`, so `μ(λ) ≥ ∫_0^1 u λ^{u−1} du`, which is the stated expression;
`1 − λ + λ log λ` vanishes at 1 and has derivative `log λ < 0`. (iii) Multiply `1/log p = ∫ e^{−λp} μ(λ) dλ`
by `s(p) log p · e^{−p/x}` and sum over p; the interchange is Fubini's theorem (the terms change sign, but
`Σ_p log p · e^{−p/x} ∫ e^{−λp} μ(λ) dλ = Σ_p e^{−p/x} < ∞`). For general α use
`1/log p = α/log(p^α) = α ∫ e^{−λ p^α} μ(λ) dλ`. ∎

`μ = ν'` where `ν(λ) = ∫_0^∞ λ^u/Γ(u+1) du` is Volterra's function, whose Laplace transform `1/(s log s)` is
classical (Erdélyi et al., *Higher Transcendental Functions* III, §18.3 — quoted from memory; the proof above
is self-contained). Numerical checks (`identities.py`, items 2, 3): `∫ e^{−λt} μ = 0.91014, 0.43419` against
`1/log t = 0.91024, 0.43429` at t = 3, 10 (the quadrature neglects `λ < e^{−10^4}`, worth `10^-4`); the
mixture identity at x = 2 and 10 to 6 digits; `μ(λ)·λ log² λ = 1.051, 1.081, 1.148` at `λ = 10^-8, 10^-4, 10^-2`
(the lower bound (ii) has 1 there; the true expansion is `1 + 2γ/log(1/λ) + …`).

So `−C(x)` is an average, with a *positive* measure, of the values `W(y)`, `0 < y < x`: positivity of W on
`(0, x]` gives `C(x) < 0`, and a lower bound for W transfers to C without any further hypothesis.

*Proof of Theorem 4.* Put `λ = v/x`, so `y_λ = x/(1+v)`. By Theorem 1 (or Theorem 5 and Theorem 3),
`W ≥ 0` everywhere, so the range `v > V_n` may be dropped, and for `v ≤ V_n < x`:
`W(y_λ) ≥ κ_*(x/(1+v)) √x (1+v)^{−1/2}`. With (ii), `μ(v/x)/x ≥ m(v/x)/(v log²(x/v))`,
`m(λ) = 1 − λ(1 + log(1/λ)) ∈ (0,1)`. Hence

    −C(x) ≥ √x ∫_0^{V_n} κ_*(x/(1+v)) (1+v)^{−1/2} m(v/x) dv/(v log²(x/v)).

On `(0,1]`: `κ_* ≥ κ_*(x/2)` and
`∫_0^1 (1+v)^{−1/2} dv/(v log²(x/v)) ≥ ∫_0^1 dv/(v log²(x/v)) − (1/log² x) ∫_0^1 (1 − (1+v)^{−1/2}) dv/v = 1/log x − c_a/log² x`.
On `(V_{i−1}, V_i]`: `κ_* ≥ κ_*(x/(1+V_i))`, `log(x/v) ≤ log(x/V_{i−1})`, and
`∫ dv/(v√(1+v)) = A(v)`. The factor m costs at most
`(√π/2) x^{-1} ∫_0^{V_n} (1+v)^{−1/2} (1 + log(x/v))/log²(x/v) dv ≤ (√π/2) x^{-1} · 2√(1+V_n) (1+A_n)/A_n²`
(`κ_* ≤ √π/2`; `(1+a)/a²` is decreasing). Multiply by `log x/√x`. ∎

### 2.9 Proof of Corollary 4a

*A minorant with controlled monotonicity.* Let `Y_E = 10^12`, `Y_I = 10^16` and, for `y ≥ 500`,

    κ^♭(y) = √π/2 − S_Γ − Ē(y) − (c_1 − 1 + ε_a) Γ(5/4) y^{−1/4} − c_1 Γ(9/8) y^{−3/8} − [ Σ_{k odd ≥3, p_0^k<y} b_k(y) + J ] / √y,

`Ē = 8.5·10^-11` for `y ≤ Y_E` and `ε_a √π/2` for `y > Y_E`; for `y ≤ Y_I`: `b_k = c_3 Γ(1+1/(2k)) y^{1/(2k)}`,
`J = c_3 √p_0/(1 − 1/p_0) + 1 + 10^{-3000}`; for `y > Y_I`: `b_k = c_θ Γ(1+1/k) y^{1/k}`, `J = c_θ p_0/(1−1/p_0) + 1`.
Then **`κ(y') ≥ κ^♭(y)` whenever `y' ≥ y ≥ 500` lie in the same regime** (`[500, Y_E]`, `(Y_E, Y_I]`, `(Y_I, ∞)`):
`Ē ≥ E/√y` on each regime; the maximum in κ is attained by its second argument from `y = 464` (χ_{−4}) resp.
`484` (χ_{−3}) on; every term
with a fixed k is decreasing in y; and an odd power k that is absent at y (`p_0^k ≥ y`) but present at `y'`
contributes at `y'` at most `c_3 (p_0^k)^{1/(2k) − 1/2} = c_3 √p_0 · p_0^{−k/2}` (resp. `c_θ p_0 · p_0^{−k/2}`), and
the sum of these over the odd k with `p_0^k ≥ y` is at most `c_3 √p_0 y^{−1/2}/(1 − 1/p_0)` (resp.
`c_θ p_0 y^{−1/2}/(1−1/p_0)`). Therefore `κ_*(y) ≥ κ^⋆(y)`, where
the regime-II formula at `Y_E` and the regime-III formula at `Y_I` give lower bounds for κ on `(Y_E, Y_I]` and on
`(Y_I, ∞)`: `K_E = 0.883736` and `K_I = 0.883030` (χ_{−3}: `0.884090`, `0.883404`), and
`κ^⋆(y) = min(κ^♭(y), K_E, K_I)` for `y ≤ Y_E`.

*Range `x ≥ 10^14`.* Take `V = 1, 2, 4, …, 64` in Theorem 4. All arguments of `κ_*` are `≥ x/65 > 10^12`, so
`κ_* ≥ min(K_E, K_I) = 0.883030`; `log x/log²(x/V) ≥ 1/log x`; `A(64) − A(1) − c_a = 1.137 > 0`; the last term is `< 2·10^-13`.
Hence `N_C(x) ≥ 0.883030 − 2·10^-13 > 0.88`.

*Range `X_A ≤ x ≤ 10^14`* (grid argument; `bounds_abel.py`, item 3). For `x ∈ [a, ρa]`, `ρ = 10^{1/100}`,
Theorem 4 with `V_i = 2^i` (as long as `a/(1+V_i) ≥ 500`) gives, since `κ_*(y') ≥ κ_*(y) ≥ κ^⋆(y)` for `y' ≥ y`,
`log x/log²(x/V)` is decreasing in x, and `log x/x ≤ log(ρa)/a`:

    N_C(x) ≥ β(a) := κ^⋆(a/2)(1 − c_a/log a) + Σ_i κ^⋆(a/(1+V_i)) [A(V_i) − A(V_{i−1})] log a/log²(ρa/V_{i−1})
                     − √π √(1+V_n)(1 + log(a/V_n)) log(ρa)/(log²(a/V_n) · a).

With `X_A = 10^{6.6} = 3 981 072`: `β(a) ≥ 0.89604` for every `a = 10^{j/100}`, `X_A ≤ a ≤ 10^14` (the minimum is at
the left end; margin 0.016). (`β(a) ≥ 0.88` already holds from `a = 2 137 963` on, with margin `1.4·10^-4` there.
30-digit mpmath evaluation of closed forms, not interval arithmetic.)

*Range `70 ≤ x ≤ X_A`* (certified computation, §5): `N_C(x) ≥ 0.88` for every real x in `[69.79, 2001]` and in
`[1995, 5.01·10^6]`, by Taylor expansion around grid points in ball arithmetic (`certify_arb_q4.out`,
`certify_arb_q4_5e6.out`). The two certified ranges overlap, and `X_A = 3.98·10^6 < 5.01·10^6`. ∎

For χ_{−3}: the same with `X_A = 10^{6.6}`: `β(a) ≥ 0.88607` on `[X_A, 10^14]` (`β ≥ 0.88` from `a = 3 311 312` on),
`N_C ≥ 0.883404` for `x ≥ 10^14`; certified ranges `[13.53, 2001]` and `[1995, 5.01·10^6]`
(`certify_arb_q3.out`, `certify_arb_q3_5e6.out`).

### 2.10 Lemma 7 and proof of Theorem 5 (no hypothesis)

Let `T = 2^22`. `Θ_1(p) > 0` at every race prime `p_0 < p ≤ T` (certified: `cert_abel.py` A1 =
[E] `cert_small.py` B; for `p ≤ 10^4` also by exact integers, [E]) and at the endpoint T itself, which is not
a prime (A0: `Θ_1(2^22) = 5 892 140 560.76` for χ_{−4}, `6 049 850 364.72` for χ_{−3}), hence `Θ_1 > 0` on `(p_0, T]`
(piecewise linear with nodes at the primes; on the last segment `[4 194 301, T]` it is linear, positive at both
ends). For `0 < y ≤ T/60`, by Lemma 1(i),

    y² W(y) = ∫_{p_0}^T Θ_1 e^{−t/y} dt + ∫_T^∞ Θ_1 e^{−t/y} dt ≥ A e^{−a/y} − (c_θ/2) y e^{−T/y} (T² + 2Ty + 2y²),

where `A e^{−a/y}` bounds the first integral from below by its part over `[5,7]` (q = 4: `Θ_1 ≥ Θ_1(7) =
4 log 3 − 2 log 5`, `A = 2Θ_1(7) = 2.3511`, `a = 7`) resp. `[3,5]` (q = 3: `Θ_1 ≥ log 2`, `A = 2 log 2`,
`a = 5`), and `|Θ_1(t)| ≤ c_θ t²/2`. The quotient `(c_θ/2) y (T² + 2Ty + 2y²) e^{−(T−a)/y}/A` is increasing in y
and equals `2.4·10^-9` (`4.0·10^-9`) at `y = T/60` (certified, A4). So `W(y) > 0` for `0 < y ≤ 69 905`, and
by Lemma 5(iii), `C(x) < 0` for `0 < x ≤ 69 905`. ∎

*Second proof of Theorem 1 for the exponential weight (independent of the analytic part of [E]
Theorem 1).* `W(y) > 0` for `y ≤ 69 905` (Theorem 5) and for `y ≥ 69 905` (Theorem 3 and §2.9, GRH:
`κ(y) ≥ min(κ^♭(69 905), K_E, K_I) ≥ 0.63`); so `W > 0` on `(0,∞)`, and `C < 0` on `(0,∞)` by Lemma 5(iii). This proof uses from [E] only the explicit formula
(Lemma 3(a)), `σ_0²`, and the certified prime data.

*Elementary range.* For `x ≤ 2` (q = 4): `C(x) ≤ −e^{−3/x} + Σ_{n≥5} e^{−n/x} = −e^{−3/x}(1 − u²/(1−u))`,
`u = e^{−1/x} < 0.607`, and `u² + u < 1`.

### 2.11 Proof of Proposition 8 (added in the revision)

Throughout GRH(χ); `γ_0 = 0.5772…` is Euler's constant (the letter γ is used for ordinates). `Z(y) = Σ_ρ Γ(ρ) y^{iγ}`
converges absolutely, is real (the zeros are symmetric and `Γ(ρ̄) = conj Γ(ρ)`), `|Z| ≤ S_Γ`, and
`|Z(y') − Z(y)| ≤ S' |log(y'/y)|` with `S' = Σ_ρ |γ| |Γ(ρ)| < ∞`.

*(i) W.* By §2.6, `W = −ψ_e + Q_2 − Σ_{k≥3 odd} W_k + Σ_{k≥4 even} Q_k`. By Lemma 3, `−ψ_e(y) = √y Z(y) + O(1)`
for `y ≥ 1` (every term other than the zero sum is bounded). By §2.6 (Squares and Limits),
`|Q_2(y) − (√π/2)√y| ≤ c_1 Γ(5/4) y^{1/4} + E(y)` and `E(y) = o(√y)` when the ladder is continued. The powers
`k ≥ 3` contribute at most `Σ_{k≥3, p_0^k<y} c_θ Γ(1+1/k) y^{1/k} + 2 ≤ c_θ y^{1/3} (1 + log y/log p_0) + 2 = o(√y)`.
Hence `W(y) = √y (√π/2 + Z(y) + ε(y))` with `ε(y) → 0`. Also `M := sup_{y>0} |W(y)|/√y < ∞`: for `y ≤ 4`,
`|W(y)| ≤ Σ_{n≥2} log n · e^{−n/y}` is bounded by a constant times `e^{−2/y}`.

*Three upper bounds for μ.* (μ1) For `0 < λ ≤ 1/4`, with `a = log(1/λ)`: `μ(λ) ≤ 1/(λ (a − γ_0)²) + 3.1`.
(μ2) For `1/4 ≤ λ ≤ 1`: `μ(λ) ≤ 8`. (μ3) For `λ ≥ 1`: `μ(λ) ≤ 0.75 + 1.13 λ e^λ`.
Proof: `log Γ(1+u) = −γ_0 u + Σ_{k≥2} (−1)^k ζ(k) u^k/k ≥ −γ_0 u` for `0 ≤ u ≤ 1` (alternating series with
decreasing terms), so `1/Γ(u) = u/Γ(1+u) ≤ u e^{γ_0 u}` on `[0,1]`; and `Γ(u) ≥ 0.8856 (n−1)!` on `[n, n+1]`,
`n ≥ 1` (since `Γ ≥ 0.8856` on `[1,2]`). For `λ ≤ 1`: `∫_0^1 λ^{u−1}/Γ(u) du ≤ λ^{-1} ∫_0^1 u e^{−(a−γ_0)u} du`, which
is `≤ 1/(λ(a−γ_0)²)` if `a > γ_0` and `≤ 4 ∫_0^1 u e^{γ_0 u} du ≤ 3` if `λ ≥ 1/4`; and
`∫_1^∞ λ^{u−1}/Γ(u) du ≤ Σ_{n≥1} 1/(0.8856 (n−1)!) ≤ 3.1`. For `λ ≥ 1`: `∫_0^1 ≤ ∫_0^1 u e^{γ_0 u} du ≤ 0.75` and
`∫_n^{n+1} λ^{u−1}/Γ(u) du ≤ λ^n/(0.8856 (n−1)!)`. Together with Lemma 5(ii): `λ log²(1/λ) μ(λ) → 1` as `λ → 0`.

*(ii) N_C.* Let `x ≥ 16`, `0 < η ≤ 1`, and split `−C(x) = ∫_0^∞ W(y_λ) μ(λ) dλ`, `y_λ = (λ + 1/x)^{-1}`, at
`λ = η/x, 1/4, 1`; put `v = λx`, so `y_λ = x/(1+v)`.

* `λ ≥ 1`: `y_λ ≤ 1/λ ≤ 1`, so `|W(y_λ)| ≤ C_0 e^{−2λ}` with `C_0 = Σ_{n≥2} log n · e^{−(n−2)}`; with (μ3) the integral is `O(1)`.
* `1/4 ≤ λ ≤ 1`: `|W(y_λ)| ≤ M √4` and (μ2): `O(1)`.
* `η/x < λ ≤ 1/4`: `|W(y_λ)| ≤ M √x (1+v)^{−1/2}`, and by (μ1) `μ(v/x)/x ≤ 1/(v (log(x/v) − γ_0)²) + 3.1/x`. The
  second term gives `O(1)`. For `v ≤ √x`: `log(x/v) − γ_0 ≥ (1/4) log x`, and `∫_η^∞ dv/(v √(1+v)) = −A(η) ≤ log(6/η)`:
  contribution `≤ 16 M √x log(6/η)/log² x`. For `√x ≤ v ≤ x/4`: `(1+v)^{−1/2} ≤ x^{−1/4}` and
  `∫ dv/(v (log(x/v) − γ_0)²) ≤ 1/(log 4 − γ_0) ≤ 1.24`: contribution `O(x^{1/4})`.
* `0 < λ ≤ η/x`: here `y_λ ∈ [x/2, x]`, `|Z(y_λ) − Z(x)| ≤ S' log(1+v) ≤ S' η`, `|(1+v)^{−1/2} − 1| ≤ η/2`, so
  `W(y_λ) = √x [ √π/2 + Z(x) + θ (c η + 2 ε^*(x)) ]`, `|θ| ≤ 1`, `c = S' + (√π/2 + S_Γ + 1)/2`, `ε^*(x) = sup_{y ≥ x/2} |ε(y)| → 0`.
  The mass `m_η(x) = ∫_0^{η/x} μ(λ) dλ` satisfies, by Lemma 5(ii) and (μ1),
  `1/log(x/η) − 2η/x ≤ m_η(x) ≤ 1/(log(x/η) − γ_0) + 3.1 η/x`, hence `log x · m_η(x) → 1` as `x → ∞` (η fixed).

Multiplying by `log x/√x`: `N_C(x) = (√π/2 + Z(x)) log x · m_η(x) + O(c η + ε^*(x)) + O_η(1/log x)`. Let `x → ∞`:
`limsup |N_C(x) − √π/2 − Z(x)| ≤ 2 c η`. As η is arbitrary, `N_C(x) = √π/2 + Z(x) + o(1)`.

*(iii) Bounds for liminf and limsup under GRH alone.* `Z ≥ −S_Γ` gives the lower bound for the liminf. Let
`ρ_1 = 1/2 + iγ_1` be the lowest zero (certified to exist, Lemma 4), of multiplicity `m_1 ≥ 1`. There are
arbitrarily large x with `Γ(ρ_1) x^{iγ_1} = −|Γ(ρ_1)|`; at such x,
`Z(x) ≤ −2 m_1 |Γ(ρ_1)| + (S_Γ − 2 m_1 |Γ(ρ_1)|) ≤ S_Γ − 4 |Γ(ρ_1)|`. With `S_Γ ≤ 3.9199328·10^-4` and
`2|Γ(ρ_1)| ≥ 3.9147116·10^-4` (`cert_abel.out`, B3): `liminf Z ≤ −3.90949·10^-4`. The limsup is symmetric.
(χ_{−3}: `S_Γ ≤ 1.6531314·10^-5`, `2|Γ(ρ_1)| ≥ 1.6425034·10^-5`.)

*(iv) Under LI.* If the positive ordinates `γ_1 < γ_2 < …` are linearly independent over Q (so the zeros are
simple), Kronecker's theorem (for linearly independent reals `γ_1,…,γ_J`, any targets and any `ε > 0`, there
are arbitrarily large t with `‖(γ_j t − φ_j)/2π‖ < ε` for all `j ≤ J`; Hardy–Wright, Theorem 444, quoted from
memory) gives arbitrarily large `x = e^t` at which each of the first J terms `2 Re Γ(ρ_j) x^{iγ_j}` is within
`2|Γ(ρ_j)|(1 − cos 2πε)` of `−2|Γ(ρ_j)|`; the remaining terms are at most `2 Σ_{j>J} |Γ(ρ_j)|` in absolute value.
Letting `ε → 0`, `J → ∞`: `liminf Z = −S_Γ`, and likewise `limsup Z = S_Γ`. ∎

*Consistency with the data.* (i) is what the table of §4 shows between `10^8` and `10^9`. For (ii) the
convergence is only logarithmic: `(N_C(x) − √π/2) log x = 2.39, 2.33, 2.27` at `x = 10^7, 10^8, 10^9`, against the
limit `(√π/2)(γ_0 + 2 log 2) = 1.74` of the `1/log x` term.

---------------------------------------------------------------------------

## 3. Numerics (the truth)

Program `abel_scan.c` (segmented sieve; `gcc -O2 abel_scan.c -lm`): for each grid point `x_j = 10^{j/M}` it
sums `e^{−p/x_j}` over all race primes `p ≤ 50 x_j` (the neglected tail is `< (x+1) e^{−50} < 2·10^-22 (x+1)`:
`Σ_{n>P} e^{−n/x} ≤ e^{−P/x}/(1 − e^{−1/x}) ≤ (x+1) e^{−P/x}`), separately for the two classes, together with
`W`, the prime-square sum, and the moments `Σ (p/x)^k e^{−p/x}`, `k ≤ 8`. Runs: `q = 4`, `10^-2 ≤ x ≤ 10^8`,
M = 50 (primes to `5·10^9`, 4.6 min); `q = 4`, `10^8 ≤ x ≤ 10^9`, M = 10 (primes to `5·10^10`, 8.7 min);
`q = 3`, `10^-2 ≤ x ≤ 10^8`, M = 50 (4.1 min).

χ_{−4} (floating point):

    x        C(x)               N_C(x)     W(x)/√x    (prime squares)/√x
    1        −0.043975409       —          0.045661   0.000137
    10       −0.727045820       0.529392   0.327200   0.187611
    10^2     −2.029394816       0.934571   0.558237   0.489184
    10^3     −4.850496994       1.059554   0.691396   0.657423
    10^4     −12.014849566      1.106609   0.803107   0.760480
    10^5     −29.953353835      1.090514   0.849241   0.815012
    10^6     −76.908394976      1.062529   0.870169   0.848382
    10^7     −202.979782396     1.034586   0.878160   0.865439
    10^8     −549.890650713     1.012936   0.882741   0.874914
    10^9     −1519.254802874    0.995609   0.884672   0.880102

* `C(x) < 0` at every grid point with `x ≥ 0.0603` (462 points up to `10^8`, 11 from `10^8` to `10^9`; no
  positive value). At the 39 grid points below 0.06 the program's cut `p ≤ 50x` excludes p = 3 and prints 0;
  there `C(x) = −e^{−3/x}(1 + o(1))`, see the elementary range in §2.10. As `x → 0+`, `C(x) → 0−`:
  the supremum of C over `(0,∞)` is 0 and is not attained; on `[1,∞)` the grid maximum is `C(1) = −0.04398`.
* C is strictly decreasing along the grid, and `x C'(x) = Σ χ(p)(p/x)e^{−p/x} < 0` at all grid points
  with `x ≥ 0.07` (observation; monotonicity of C is not proved here). `W > 0` at all grid points.
* `N_C`: below 0.88 for `x < 68.86`; maximum over the grid `1.108761` at `x = 6309.6` (the true maximum, on a
  finer mesh, is `1.1087617` at `x ≈ 6356`); minimum over
  `[10^4, 10^9]` at the right end, `0.995609`. The expected leading shape is
  `(√π/2)(1 + (γ + 2 log 2)/log x)`, which is `0.9702` at `10^9` (observed 0.9956; the difference is of the
  order of the next term, `1/log² x`).
* `W(y)/√y` increases towards `√π/2 = 0.886227` up to `y ≈ 7·10^7`; from there on the oscillation of the
  first zero is visible (table in §4).

χ_{−3} (floating point): `C_3(x) < 0` at all grid points with `x ≥ 0.05`; `N_C` crosses 0.88 at `x ≈ 13.43`,
has its grid maximum `1.122095` at `x = 2089.3` (true maximum `1.122098` at `x ≈ 2129`), and equals `1.013309` at `10^8`; `C_3(10^8) = −550.093087594`.

---------------------------------------------------------------------------

## 4. Cross-check of the analytic bounds against the scan

`bounds_abel.py`, item 2 (`bounds_q4.out`, `bounds_q3.out`). `κ`, `κ^up`: Theorem 3. "bound for N_C":
Theorem 4 with `V_i = 2^i` and the rigorous minorant `κ^⋆`.

χ_{−4}:

    x = y      κ(y)       W/√y (scan)  κ^up(y)    κ^⋆(y)     | bound for N_C   N_C (scan)   ratio
    10^2       −0.2227    0.558237     1.8667     —          | —               0.934571     —
    10^3        0.2197    0.691396     1.5443     0.0291     | —               1.059554     —
    10^4        0.5379    0.803107     1.2572     0.4768     | 0.3683          1.106609     0.333
    10^5        0.7128    0.849241     1.0778     0.6933     | 0.6931          1.090514     0.636
    10^6        0.8017    0.870169     0.9810     0.7954     | 0.8503          1.062529     0.800
    10^7        0.8450    0.878160     0.9324     0.8430     | 0.9167          1.034586     0.886
    10^8        0.8659    0.882741     0.9088     0.8653     | 0.9395          1.012936     0.927
    10^9        0.8760    0.884672     0.8974     0.8758     | 0.9446          0.995609     0.949
    10^12       0.8845    —            —          0.8830     |
    10^50       0.88583   —            —          0.88507    |

χ_{−3}:

    10^4        0.4743    0.784857     1.3657     0.3866     | 0.2736          1.113664     0.246
    10^6        0.7936    0.867815     0.9958     0.7840     | 0.8321          1.063046     0.783
    10^8        0.8652    0.882810     0.9105     0.8642     | 0.9370          1.013309     0.925

`κ ≤ W/√y ≤ κ^up` holds at all 380 (χ_{−4}) resp. 370 (χ_{−3}) grid points with `y ≥ 4`, and the bound for `N_C`
lies below the scan value at every grid point `x ≥ 2000`. The bounds are on the correct side and not
absurdly far: for W the lower bound is 67% of the truth at `10^4`, 92% at `10^6`, 98% at `10^8`, 99% at `10^9`; for
`N_C` it is 33% at `10^4`, 80% at `10^6`, 93% at `10^8`, 95% at `10^9`.

*Component by component* (`crosscheck.py`, `crosscheck_q4.out`, `crosscheck_q3.out`; double precision). So that
a wrong constant in one component cannot hide behind slack in another, each inequality of §2.6 was tested
separately at `y = 10, 10^2, …, 10^6` against direct sums: squares (two-sided; also at all 381 scan points),
fourth powers, `|W_k| ≤ B_k` for k = 3, 5, 7, the high powers (`≤ 0.23` odd, `≤ 0.01` even, against the bound 1),
and the zero term `−S_Γ√y ≤ −ψ_e ≤ S_Γ√y + b + 0.2555/y²`: all hold. The zero bound is nearly attained:
at `y = 10^6` (χ_{−4}), `−ψ_e − b = +0.391229` against `S_Γ √y = 0.391993`. Likewise Theorem 4 was tested with
the *true* W (scan values) in place of κ, which isolates Lemma 5 and the constants `c_a`, `A(V)`: the right
side is below `N_C` at all 236 grid points in `[2000, 10^8]`, at 69% (x = 3162), 79% (`10^4`), 92% (`10^6`),
96% (`10^8`) of it.

*The zero term is visible in the data.* By Lemma 3 the oscillating part of `W(y)/√y` is
`Σ_ρ Γ(ρ) y^{iγ} ≈ 2 Re Γ(ρ_1) y^{iγ_1}` (`γ_1 = 6.0209…`, amplitude `3.91·10^-4`, period 0.453 in `log_10 y`).
Scan values for χ_{−4} between `10^8` and `10^9`:

    y            W/√y (scan)   2 Re Γ(ρ_1) y^{iγ_1}   difference
    1.000·10^8   0.882741      −0.000337              0.883078
    1.259·10^8   0.883060      −0.000258              0.883317
    1.585·10^8   0.883787      +0.000243              0.883544
    1.995·10^8   0.884107      +0.000347              0.883760
    2.512·10^8   0.883852      −0.000115              0.883967
    3.162·10^8   0.883776      −0.000389              0.884165
    3.981·10^8   0.884324      −0.000027              0.884351
    5.012·10^8   0.884908      +0.000379              0.884529
    6.310·10^8   0.884858      +0.000166              0.884692
    7.943·10^8   0.884517      −0.000318              0.884835
    1.000·10^9   0.884672      −0.000283              0.884955

The raw values are not monotone; after subtracting the first-zero term they increase smoothly. This
confirms sign and size of the zero sum in Lemma 3 at large y, independently of the small-y check of §2.3.

Where the loss comes from (χ_{−4}, `y = 10^8`, relative to √y): zeros `0.0004` (worst-case bound for an
oscillating term of true size `≤ 0.0004`); prime squares `c_1 Γ(5/4) y^{−1/4} = 0.0188` minus the fourth
powers' `0.0090 − 0.0020` (the constant `c_1 = 2.07` is a worst case; the mean of `(u − ϑ_o(u))/√u` is
about 1); odd powers `0.0043 + 0.0013 + …`. For `N_C` the further loss is in `μ ≥ …` (Lemma 5(ii) captures
`2 log 2` of the `γ + 2 log 2` in the `1/log x` term) and in using `κ_*` at the lower end of each piece.

---------------------------------------------------------------------------

## 5. What is certified and what is floating point

| item | how | where |
|---|---|---|
| `Θ_1(p) > 0`, `p ≤ 2^22` | Arb, 256 bits (and exact integers for `p ≤ 10^4` in [E]) | `cert_abel.py` A1 |
| endpoint `u = 2^22` (not a prime): `Θ_1(2^22) > 0`, `(2^22 − ϑ_o(2^22))/2^11 ≤ c_1` | Arb, 256 bits (added in the revision, REP-1) | A0 |
| `c_1`, `ϑ_o(y) < y`, `y ≤ 2^22` | Arb | A2 |
| `c_3` (`|Θ(u)| < c_3√u`, `u ≤ 2^22`) | Arb | A3 |
| Theorem 5 (`x ≤ 69 905`, no hypothesis) | Arb + Lemma 7 | A4 |
| zeros: sign changes of Z, `d`, `g`, `S_Γ` | Arb, 160 bits | B1–B3 |
| `N_C > 0.88` on `[69.79, 2001]` (χ_{−4}), `[13.53, 2001]` (χ_{−3}) | Arb, 80 bits, grid `10^{j/400}`, Taylor order 8; moments, prime tails, remainder R and the comparisons all in Arb (revision, REP-2) | `certify_grid.py arb`, `certify_arb_q*.out` |
| `N_C > 0.88` on `[1995.27, 5.01·10^6]` (both characters) | the same, grid `10^{j/25}`, Taylor order 5 | `certify_arb_q*_5e6.out` |
| closed-form functions `κ`, `κ^♭`, `β(a)`; thresholds 413 / 750, `X_A` | mpmath, 30 digits, **not** interval arithmetic; smallest margins: `β − 0.88 ≥ 0.016` (χ_{−4}), `0.006` (χ_{−3}) on `[10^{6.6}, 10^14]`, `0.003` beyond | `bounds_abel.py` |
| scan of `C`, `W`, moments to `10^8`/`10^9` | long double accumulation of double `exp`; **floating point** | `abel_scan.c` |
| `C < 0` on every x in `[0.103, 1.02·10^8]`, `N_C ≥ 0.88` on `[74.2, 1.02·10^8]` (χ_{−3}: `[0.068, 1.02·10^8]`, `[14.2, 1.02·10^8]`) | from the scan, with an a-priori rounding budget (below); **not** interval arithmetic | `certify_grid.py scan`, `certify_scan_q*.out` |
| identities (explicit formula, μ, mixture, kernel) | mpmath quadrature, numerical checks only | `identities.py` |
| Fujii weights, §6 | double precision, exploration only | `fujii_scan.c` |

*The Taylor certificates.* `C(x_j e^h) = Σ_p χ(p) exp(−z e^{−h})`, `z = p/x_j`. Its k-th derivative in h is
`Σ_p χ(p) P_k(w) e^{−w}`, `w = z e^{−h}`, with `P_0 = 1`, `P_{k+1}(w) = w (P_k(w) − P_k'(w))` (signed Stirling
numbers of the second kind). So the derivatives at `h = 0` are integer combinations of the moments
`M_m^± = Σ_{χ(p)=±1} z^m e^{−z}`, and for `|h| ≤ δ`

    −C(x_j e^h) ≥ −D_0 − Σ_{k=1}^K |D_k| δ^k/k! − R,     R = (δ^{K+1}/(K+1)!) Σ_m |c_{K+1,m}| e^{2mδ} [X m! + (m/e)^m],  X = x_j e^δ

(the bracket bounds `Σ_{n≥1} v^m e^{−v}`, `v = n/X`, for the unimodal function `v^m e^{−v}`; the factor
`e^{2mδ}` passes from `sup_{|h|≤δ} w^m e^{−w}` to `v^m e^{−v}`). On a cell with left end `≥ e²`,
`log x/√x ≥` its value at the right end. In `arb` mode everything that enters a decision is an Arb ball:
the moments (primes up to an integer `P_j ≈ 60 x_j`); the tails, bounded for each m by
`Σ_{n>P_j} (n/x)^m e^{−n/x} ≤ c^m e^{−c} + x Γ(m+1, c)`, `c = P_j/x_j` (the function `v^m e^{−v}` decreases for
`v ≥ m`, and `c > m` is checked), which are added to `|D_k|`; the remainder R; the weight `log X/√X`; the test
"left end `> e²`"; and the comparison with 0.88. Floating point is used only to choose `P_j` (any choice is
covered by the tail bound) and to print rounded decimals; printed interval ends are rounded inwards. (Before
the revision R and the tails were IEEE doubles inflated by `10^-6`; the referee's item REP-2. The certified
ranges did not change.) R is at most `9·10^-20 |C|` for the fine grid (order 8) and `3.2·10^-3 |C|` for the coarse one (order 5, at
the top of the range), the tail at most `3·10^-22 |C|`. The loss of the method is about δ
relative (`δ = log 10/(2M)`; observed 0.2% for M = 400, 2% for M = 50, 4% for M = 25), small against the
margins except near the crossing, hence the fine grid there. Validation of the
machinery: the order-8 Taylor polynomial at `x_j` predicts the scanned values at the neighbouring grid points
`x_{j±1}` (`h = ±2δ`) to within `4·10^-12` relative, at all grid points `x_j ≥ 1` of `out_q4_1e8.txt`.

*Rounding budget of the scan mode.* For each moment, `|computed − true| ≤ ETA · (M_m^+ + M_m^−) + tail_m`,
`ETA = n·2^-64 + 2^-45 + 2^-47` with `n ≤ 2.35·10^8` the number of race primes (`ETA = 1.28·10^-11`):
each summand `z^m e^{−z}` is formed from the double `z = p·(1/x_j)` (relative error `≤ 2^-52`, hence
`≤ 50·2^-52` in `e^{−z}`), the libm `exp` (relative error assumed `≤ 2^-46`, i.e. 64 ulp; glibc documents
`< 1` ulp), and at most 8 long double multiplications; the long double accumulation of n nonnegative
terms loses at most `n·2^-64` relative. This is a worst-case a-priori bound under the stated assumption on
`exp`, not a machine-checked enclosure. Resulting relative error budget of `−C(x_j)`, including truncation
and Taylor remainder: `≤ 1.4·10^-7` for all `x_j ≥ 1`. Agreement of the two code paths: `−C(100) =
2.029394815655` (Arb) against `2.0293948157` (scan); `−C(10^3) = 4.850496993866` against `4.8504969939`;
`−C(1 737 800.83) = 96.824879497983` (Arb).

None of the theorems uses a floating-point sign decision. Theorem 1 uses [E] Theorem 1 (exact integers +
analysis). Theorem 5 is certified (Arb, including the endpoint `2^22`). Corollary 4a uses the certified ranges
up to `5.01·10^6` (Arb throughout) and, beyond `X_A`, 30-digit evaluations of closed forms with the margins
stated; the referee re-implemented those closed forms independently (`referee_checks/beta.py`) and obtained
the same `K_E`, `K_I` and minima of β.

---------------------------------------------------------------------------

## 6. Fujii's exponent α

`C_α(x) = Σ_p χ_{−4}(p) e^{−(p/x)^α}`, `W_α(y) = Σ_p s(p) log p · e^{−(p/y)^α}`. Fujii (1988) proved that
"`C_α(x) → −∞`" is equivalent to GRH(χ_{−4}) for `0 < α < 4.19`; Platt–Trudgian (arXiv:1802.00089, Theorem 1)
for `0 < α < 20.40442` (according to Suzuki, §1.2, Fujii conjectured that it holds for every `α > 0`), the
criterion being their inequality (3):
`Σ_ρ |Γ((1/2 + iγ)/α)| < (1/2) Γ(1/(2α))` (in our normalisation: the zero term
`(1/α) Σ_ρ |Γ(ρ/α)| √y` against the bias `Γ(1 + 1/(2α)) √y` of the prime squares), and they note that their
method cannot go beyond 20.40443.

**(F1) Proved (GRH(χ_{−4})): `C_α(x) < 0` for every `x > 0` and every `0 < α ≤ 1`** (Theorem 1; the weight
`e^{−(t/x)^α}` is convex on `[3,∞)` for every `x > 0` if and only if `α ≤ 1`. For `α > 1` it is still convex on
`[3,∞)` when x is small — for `α = 2` when `x ≤ 3√2` — but for large x the kernel `w/log t` fails to be convex
near `t ≍ x`, so the method stops at `α = 1`). The same for χ_{−3}. (For `α ≤ 1` one could also transfer the
bound of Corollary 4a, since `e^{−t^α}` is then a positive mixture of exponentials; not done.)

**(F2) Reduction for every α (no hypothesis): `−C_α(x) = α ∫_0^∞ W_α((λ + x^{−α})^{−1/α}) μ(λ) dλ`**
(Lemma 5(iii)). So `C_α < 0` on `(0, X]` follows from `W_α > 0` on `(0, X]`, and lower bounds transfer.

**(F3) The criterion, computed** (`fujii.py`, `fujii_crit_q4.out`; mpmath, 511 zeros, not certified). With
`β(α) = Γ(1 + 1/(2α))` (prime squares: `Σ_p log p · e^{−(p²/y)^α} ~ β(α) √y`) and
`S(α) = (1/α) Σ_ρ |Γ(ρ/α)|` (zeros), under GRH `liminf W_α(y)/√y ≥ β(α) − S(α)` (the proof of Theorem 3 with
`Γ(s)` replaced by `Γ(s/α)/α`; only the asymptotic statement is claimed here):

    α        β(α)        S(α)          S/β
    1        0.886227    0.000392      0.00044
    2        0.906402    0.017428      0.0192
    3        0.927719    0.063345      0.0683
    4        0.941743    0.124104      0.1318
    4.19     0.943840    0.136393      0.1445        (Fujii's range ends here)
    6        0.958286    0.255830      0.2670
    8        0.967580    0.382842      0.3957
    12       0.977604    0.608730      0.6227
    16       0.982901    0.801943      0.8159
    20       0.986174    0.970517      0.9841
    20.40442 0.98643651  0.98643638    1 − 1.4·10^-7
    20.40443 0.98643652  0.98643677    1 + 2.5·10^-7
    24       0.988396    1.120325      1.1335
    64       0.995550    2.111006      2.1204
    256      0.998876    4.001647      4.0062

The root of `S = β` is `α_* = 20.4044235` (listed zeros; the largest listed zero contributes `2·10^-22`).
This reproduces Platt–Trudgian's `20.40442 / 20.40443`; at `α = 20.40442` I get `Σ_ρ |Γ(ρ/α)| = 20.1276621`
against their bound `< 20.1276643`, and `Γ(1/(2α))/2 = 20.1276649` as they state. `S/β` is increasing at all
sampled α, so the set where the criterion holds is the interval `(0, α_*)` as far as sampled; monotonicity of
`S/β` is not proved here. Heuristically (GRH and LI; the analogue of Proposition 8 for `W_α`, `C_α`, not written out) the criterion is
sharp: for `α > α_*` one expects `C_α(x) > 0` for some arbitrarily large x, which would contradict the
conjecture attributed to Fujii above. This is an expectation, not a result of this note.

**(F4) An unconditional obstruction for the route through `W_α`.** `W_α(y) > 0` for all `y > 0` is *false* for
every sampled `α ≥ 6.77` (grid computation; e.g. `W_8(6.07) = −0.122`, `W_16(6.28) = −0.463`): for large α the
sum approaches `Θ(y)`, and `Θ(5) = log 3 − log 5 < 0`. The least such α is `α_W ≈ 6.76` (bisection on a grid). So (F2) can give all-x
statements for `C_α` at most for `α < 6.76`; beyond that one has to work with `C_α` itself.

**(F5) The truth for `x ≤ 10^7`** (`fujii_scan.c`, `fujii_q4_1e7.out`; double precision, grid `10^{j/M}` with
`M = max(50, 20α)` points per decade; exploration, nothing certified). "min N" is the minimum over the grid
points `x ≥ 100` of `−C_α(x) log x/√x`.

    α       C_α > 0 on the grid?   min N (at x)             min W_α/√x on x ≥ 100
    2       no                     0.9485 (100)             0.5352
    3       no                     0.9669 (103.9)           0.5416
    4       no                     0.9666 (115.5)           0.5355
    6       no                     0.8927 (413.7)           0.4591
    8       no                     0.8001 (427.8)           0.3699
    12      no                     0.6648 (643 181)         0.2410
    16      no                     0.5718 (644 726)         0.1609
    20      no                     0.4979 (641 948)         0.1141
    24      no                     0.4363 (640 103)         0.0885
    32      no                     0.3512 (640 103)         0.0552
    64      no                     0.2101 (629 823)         −0.0662
    128     no                     0.0819 (629 257)         −0.1389
    256     no                     0.0460 (630 674)         −0.2286
    512     no                     0.0292 (628 974)         −0.3262
    1024    no                     0.0041 (623 693)         −0.4005
    2048    YES, 75 grid points    −0.0276 (623 658)        −0.4181
    4096    YES, 200 grid points   −0.0470 (623 675)        −0.4187

(For `α ≥ 96` some grid points in the tie intervals of `π(x;4,3) − π(x;4,1)`, such as `[5,7)`, have `|C_α| ≤ 10^-9`,
below double resolution; there `C_α < 0` by Abel summation, see below. They are not counted.)
Independent recomputation in 30-digit arithmetic (window of the sum treated exactly):
`C_2048(623 658) = +1.63350307888` and `C_4096(26 862.745) = +0.00889782500`; the latter is the first sign
change of Chebyshev's race at 26 861, the former the region `616 841 … 633 798` where the class 1 mod 4
leads. So **the all-x statement is false for α = 2048 and 4096, without any hypothesis**, and on `x ≤ 10^7` it
holds at every grid point for all sampled `α ≤ 1024`.

*Why.* By Abel summation `C_α(x) = −∫_0^∞ D_0(t) (−∂_t e^{−(t/x)^α}) dt`, `D_0(t) = π(t;4,3) − π(t;4,1)`: a
positive average of `−D_0` over a window of relative width about `1/α` around x. Hence `C_α(x) < 0` whenever
`D_0 ≥ 0` on the window (in particular, for every α, for all x up to nearly 26 861), and `C_α(x) > 0`
needs a window on which class 1 leads on average: relative width `7·10^-5` at 26 861, about `3·10^-3` near
623 000.

*Summary for the question "for which α does the all-x statement hold".* Proved under GRH: `0 < α ≤ 1`.
False unconditionally: `α = 2048, 4096` (presumably from some α between 1024 and 2048 on). Expected false under GRH and
LI for every `α > α_* = 20.4044`, but the first counterexample is then beyond `10^7` for `α ≤ 1024`. For
`1 < α < α_*`: no counterexample up to `10^7`, the minimum of N is comfortably positive, and I expect the
statement to be true; a proof for a fixed `α < 6.76` would follow the plan of this note (explicit version of
(F3) for `y ≥ Y_0(α)`, a Taylor-grid certificate for `W_α` on `(0, Y_0(α)]`, then (F2)); for `6.76 < α < α_*` the
certificate has to be made for `C_α` directly. None of this was carried out, and a uniform statement for a
*range* of `α > 1` would need certificates in two variables.

---------------------------------------------------------------------------

## 7. Relation to prior work; what is and is not new

Sources and how much of each was read:

* **D. Platt, T. Trudgian, "Fujii's development on Chebyshev's conjecture", arXiv:1802.00089v1.** Read in
  full (5 pages). History as stated there: Chebyshev's 1853 conjecture `Σ_{p>2} (−1)^{(p−1)/2} e^{−xp} → −∞`
  as `x → 0`; Hardy–Littlewood [Acta Math. 41 (1917), 119–196] and Landau [Math. Z. 1 (1918), 1–24, 213–219]
  showed it equivalent to GRH for `χ_4`; Fujii [Proc. Japan Acad. Ser. A 64 (1988), 260–263, Thm 1] for
  `e^{−(xp)^α}`, `0 < α < 4.19`, using only `γ_1 > 6` and `Σ_{γ>0} γ^{-2} < 1/5`; their Theorem 1 for
  `0 < α < 20.40442`, with 1000 zeros from lcalc checked by Arb and Turing's method. Their decomposition
  (from Fujii §3) is `S = S_1 + S_2 + S_3 + S_4` with `S_1 = −(1/2)Γ(1/(2α)) x^{−1/2} + o(x^{−1/2})` and
  `|S_3| ≤ x^{−1/2} Σ_ρ |Γ(1/(2α) + iγ/α)|` — an asymptotic statement, with `o(·)` terms and without the
  factor `1/log`: it is the log-weighted sum, in our notation `W_α`, that has this clean form. No statement
  for all x, and no explicit version of the `o(·)` terms, is made there.
* **Hardy–Littlewood, Landau 1918, Fujii 1988.** Not read in the original; known through Platt–Trudgian
  and Suzuki. (Platt–Trudgian date Hardy–Littlewood's paper, Acta Math. 41, 119–196, to 1917; Suzuki's
  reference list says 1916.) In particular I do not know how the originals pass from the log-weighted to the
  unweighted sum, nor whether a pointwise remark occurs in them.
* **M. Suzuki, "On variants of Chebyshev's conjecture", arXiv:2411.07436v3.** Read for the revision: abstract
  and §1.2 (pp. 4–6), i.e. Theorem 3, the historical paragraph and Theorem 4. Used for: Hardy–Littlewood and
  Landau proved the equivalence with GRH(χ_4) for both the plain sum (his (16)) and the log-weighted sum (his
  (17)); Knapowski–Turán [Number Theory and Analysis (papers in honor of E. Landau), Plenum 1969, 157–171; not
  read] for the weight `exp(−log²(p/x))`; Fujii's conjecture that his statement holds for every `α > 0`. His
  Theorem 3 — GRH(χ_4) ⇔ `Σ_{n≤x} Λ(n)χ_4(n) n^{−1/2} log(x/n) ≤ 0` for all `x ≥ x_0` — and his Theorem 4 are
  the closest published relatives of Theorem 2 (b'), (d') that I know: sign for all large x, other weights.
  None of these is a statement for every `x > 0`. The rest of the paper was not read.
* **BKLNW, arXiv:2002.11068v2.** Corollary 11.2 and its proof paragraph, Appendix A.1.1: read for this note.
* **[E] = `../explicit/EXPLICIT.md`** and its referee report: Lemma 2, Lemma 3(a),(b), Lemma 4, Theorem 1 are
  used as stated there (the referee's verdict: proved as written, given the exact check of `Θ_1` for
  `p ≤ 10^4`, which is in `cert_small.py`).
* `../OPEN-PROBLEMS-TRIAGE.md` §3.1 reports that the pointwise all-x form was not found in Suzuki
  (arXiv:2411.07436), Platt–Trudgian, Ford–Konyagin; the referee checked Suzuki's text independently and found
  no statement for all `x > 0`. I did not extend that search.

Classical (not ours): the equivalence (a)⇔(c) of Theorem 2 and its log-weighted form (`W → +∞`); Landau's theorem and the Mellin-transform
mechanism; the explicit formula with `Γ(ρ)` and the observation that the zero sum is tiny
(`|Γ(1/2 + iγ_1)| ≈ 2·10^-4`); the criterion `Σ|Γ(ρ/α)| < (1/2)Γ(1/(2α))`; second-order partial summation
(Ingham; Johnston's Lemma 2.7); Volterra's function and its Laplace transform.

Apparently new (not found in the sources above; a limited search): (1) `C(x) < 0` for **every** `x > 0`
under GRH(χ_{−4}), and its equivalence with GRH (Theorem 2 (b), (d)); the same for every convex decreasing
weight, in particular Fujii's weights with `α ≤ 1` (Theorem 1), and for `χ_{−3}`; (2) the explicit
two-sided bound for W and the explicit lower bound `0.88 √x/log x` for `x ≥ 70`, with a constant within
0.7% of the best possible; (3) the unconditional certified range `x ≤ 69 905`; (4) as devices: the remark
that the kernel `e^{−t/x}/log t` is convex, so that [E] Theorem 1 applies verbatim, and the positive
mixture `−C(x) = ∫ W((λ+1/x)^{-1}) μ(λ) dλ`, which removes the logarithm exactly. An expert would rightly
call (1) a corollary of Hardy–Littlewood–Landau-type analysis plus a finite verification; the content is
that the finite verification can be done and that the margin is enormous.

---------------------------------------------------------------------------

## 8. Constants and where each comes from

| constant | value (χ_{−4} / χ_{−3}) | source | used in |
|---|---|---|---|
| `√π/2 = ∫_0^∞ e^{−u²} du` | 0.886226925452758 | closed form | bias of the prime squares |
| `S_Γ` | ≤ 3.9199328·10^-4 / 1.6531314·10^-5 | Lemma 4 (Arb; GRH) | Theorems 3, 4 |
| `σ_0²` | 0.155567979923586 / 0.113229969857472 | [E] Lemma 3(b); recomputed as an Arb ball | Lemma 3 (`R_1`), Lemma 4 |
| `b`, `b_0` | 0.783188785 / 0.948198827; −0.060266914 / 0.143831645 | [E]; b recomputed as an Arb ball; formula check to `10^-9` | §2.3 (only `b − r^+(4) > 0`, `1 + b_0 > 0`) |
| `c_1` | 2.071193 / 2.081942 | Arb, sup as `y → 1423^-` | squares, fourth powers |
| `c_3` | 2.153716634 / 2.074662259 | Arb, sup at 126 227 / 1 419 839 | odd powers `k ≥ 3` |
| `ε_a`, ladder `ε_i` | `0.2/log²(2^22)`; `1/log³ U_i` | BKLNW Cor. 11.2, rows k = 2, 3 (read) | `E(y)`, `c_θ` |
| `Γ(5/4), Γ(9/8), Γ(1+1/k)` | 0.906402477, 0.941742700, … | closed form (Arb values printed) | prime powers |
| 0.2555 | `(2/5)σ_0² + log 2 − 1/2 ≤ 0.2554` | Lemma 3 | `κ^up` |
| 1 (twice in `κ^up`) | `(1/e) Σ log p (p_0/p)³/(1−p^{−2}) ≤ 0.695 / 0.411` | §2.6 | powers with `p_0^k ≥ y` |
| `c_a`, `A(V)` | `2 log((1+√2)/2)`, `−2 arcoth √(1+V)` | elementary integrals | Theorem 4 |
| 69 905 | `2^22/60` | Lemma 7 | Theorem 5 |
| 70 / 14 | first integers after the crossings 68.86 / 13.43; certified from 69.79 / 13.53 | Arb | Corollary 4a, Theorem 6 |

---------------------------------------------------------------------------

## 9. What I am not fully sure of

1. **Dependence on [E].** Theorem 1 (first proof) is [E] Theorem 1 plus a convexity remark; [E] has an
   in-house referee report, no external review. The second proof (§2.10) avoids the analytic part of [E]
   Theorem 1 but still uses [E] Lemma 3(a) (explicit formula for `ψ_1`; its contour estimate quotes Davenport
   from memory) and Lemma 3(b) (`Σ|ρ|^{-2} = σ_0²`, Davenport Ch. 12, quoted from memory there). Here the
   resulting formula of Lemma 3 is checked numerically to `3·10^-9/y²` for both characters, which tests `b`,
   `b_0` and the shape of the formula, but is not a proof.
2. **Lemma 4 (S_Γ).** The completeness step uses the identity for `σ_0²` (GRH). If one distrusts it, the
   winding-number count of `prototype/explore/certify_arb.py` (not rerun here) is an independent route; and
   even the crude information "no zero with `|γ| < 6`" gives `S_Γ ≤ 1.15·10^-3`, which would change
   `√π/2 − S_Γ` to `0.8851` and leave Corollary 4a intact (`K_E`, `K_I` would drop by `7.5·10^-4`, to `0.8823`).
3. **Cited inequalities.** BKLNW Corollary 11.2: the two rows were re-read in arXiv v2 (table and the proof
   paragraph for each row); the published version and the tables it relies on were not checked.
   `E_1(u) < e^{−u} log(1 + 1/u)` (Abramowitz–Stegun 5.1.20) and Landau's theorem (Montgomery–Vaughan Lemma
   15.1) are quoted from memory; the first is checked numerically (`selfcheck.py`), and is used only to see
   that `b − r^+(y) ≥ 0`. The reference for Volterra's function is from memory (the proof of Lemma 5 is
   self-contained).
4. **Closed forms evaluated without interval arithmetic** (`bounds_abel.py`, mpmath 30 digits): κ, `κ^♭`,
   `β(a)`, the thresholds 413/750. Margins: `β − 0.88 ≥ 0.016` (χ_{−4}), `0.006` (χ_{−3}) on the grid,
   `0.003` for `x ≥ 10^14`. The grid argument rests on the monotonicity statements of §2.9, which are argued
   by hand and tested on a grid (`min_{y'≥y} κ(y') ≥ κ^♭(y)` at all grid points in `[500, 10^6]`), not
   machine-checked. The referee re-derived them and re-implemented `κ^♭`, `K_E`, `K_I`, β independently, with
   the same values; they are still not interval arithmetic.
5. **Taylor certificates.** The Arb-mode results depend on python-flint/Arb and on the hand-derived
   remainder and tail bounds of §5 (re-derived by the referee); since the revision no floating-point number
   enters a decision there. The scan-mode results (to `10^8`) depend on the assumption about libm
   `exp` and are **not** interval arithmetic; no theorem uses them (they overlap the certified ranges and
   extend the numerical evidence).
6. **Hand computations.** The elementary integrals (`c_a`, `A(V)`, `F`, the Γ-integrals, `log 2 − 1/2`, the
   constant `≤ 1` for the high prime powers) were each checked numerically (`selfcheck.py`: all ok;
   `bounds_abel.py` item 1), and each component inequality of Theorem 3 and the reduction of Theorem 4 were
   tested against direct sums (`crosscheck.py`: all ok), but the chain of inequalities in §2.6, §2.8, §2.9 has not been read by a second
   person. The direction of every rounding (S_Γ, c_1, c_3 rounded up) was checked by me only.
7. **Statements not proved here:** monotonicity of `C(x)` (observed); the remark after Theorem 1 for
   `p^{−σ}`, `1/2 < σ ≤ 1`; everything in §6 except (F1), (F2) and the two counterexamples (in particular the
   expected failure for `α > α_*` under LI, and monotonicity of `S/β`).
   **Proposition 8 and its proof (§2.11) were written after the referee report and have not been refereed.**
   The proof replaces the earlier sketch; its ingredients are Theorem 3's decomposition, three elementary
   upper bounds for μ (spot-checked numerically in `selfcheck.py`), and Kronecker's theorem
   for the LI part (Hardy–Wright Theorem 444, quoted from memory). The numerical bounds
   `liminf N_C ∈ [0.8858349, 0.8858360]` need only GRH, the certified `S_Γ` and the certified first zero.
8. **Classical attributions.** Hardy–Littlewood (1917) and Landau (1918) are cited through Platt–Trudgian;
   the originals and Fujii's paper were not read. Whether a pointwise statement or the mixture device occurs
   in them or elsewhere is unknown to me; the novelty claims of §7 are limited accordingly.
9. **Corollary 7.** I rely on `../GRH-PROOF-DIGEST.md` for the status of the announced proof and did not
   examine it; the corollary is a conditional statement and should be worded as such wherever it is used.
10. **χ_{−3}.** Same proofs and scripts with `q = 3`; the Arb certificates and the bounds were run, the
    identities check only for Lemma 3 (`identities_q3.out`).
11. **Theorem 3, limit of `E(y)/√y`.** The first version said that `E(y)/√y → 0` for the ladder as printed;
    with a finite ladder the limit is `1.13·10^-9`, not 0. The statement and the proof of the liminf/limsup
    now use the continued ladder (§2.6, "Limits"). Found while writing §2.11, not by the referee.

---------------------------------------------------------------------------

## 10. Files

    abel_scan.c              scan of C, W, prime squares, moments (C; abel_scan is the gcc -O2 binary)
    out_q4_1e8.txt           q = 4, 10^-2..10^8, 50 points per decade   (time_q4_1e8.txt)
    out_q4_1e9.txt           q = 4, 10^8..10^9, 10 points per decade    (time_q4_1e9.txt)
    out_q3_1e8.txt           q = 3, 10^-2..10^8, 50 points per decade   (time_q3_1e8.txt)
    out_q4_1e5.txt           first test run (superseded by out_q4_1e8.txt)
    cert_abel.py / .out      Arb: prime data to 2^22, Lemma 7, zeros and S_Gamma, constants
    bounds_abel.py           kappa, kappa^up, kappa-flat, Theorem 4, grid argument; bounds_q4.out, bounds_q3.out
    certify_grid.py          Taylor certificates (scan mode: rounding budget; arb mode: balls)
    certify_scan_q4.out, certify_scan_q3.out
    certify_arb_q4.out, certify_arb_q4_5e6.out, certify_arb_q3.out, certify_arb_q3_5e6.out
    identities.py            numerical checks of Lemma 1, Lemma 3, Lemma 5; identities_q4.out, identities_q3.out
    fujii_scan.c, fujii.py   Fujii weights: scan (double) and the criterion; fujii_q4_1e7.out, fujii_crit_q4.out
    selfcheck.py / .out      numerical checks of the elementary integrals and inequalities used by hand
    crosscheck.py            component-wise tests of Theorems 3 and 4 against direct sums; crosscheck_q4.out, crosscheck_q3.out
    REFEREE-REPORT.md        the referee's report (not written by me)
    referee_checks/          the referee's scripts nc.py, ef.py, beta.py, fujii.py (not written by me; untouched)

Reproduce: `gcc -O2 -o abel_scan abel_scan.c -lm && ./abel_scan 4 -2 8 50 > out_q4_1e8.txt && ./abel_scan 4 8 9 10 > out_q4_1e9.txt && ./abel_scan 3 -2 8 50 > out_q3_1e8.txt && python3 cert_abel.py > cert_abel.out && python3 bounds_abel.py 4 > bounds_q4.out && python3 bounds_abel.py 3 > bounds_q3.out && python3 certify_grid.py scan out_q4_1e8.txt > certify_scan_q4.out && python3 certify_grid.py arb 4 60 2000 400 > certify_arb_q4.out && TAYLOR_K=5 python3 certify_grid.py arb 4 2000 5000000 25 > certify_arb_q4_5e6.out && python3 identities.py 4 > identities_q4.out` ; for q = 3: `python3 certify_grid.py scan out_q3_1e8.txt`, `python3 certify_grid.py arb 3 8 2000 400`, `TAYLOR_K=5 python3 certify_grid.py arb 3 2000 5000000 25`, `python3 identities.py 3 part1`; Fujii: `gcc -O2 -o fujii_scan fujii_scan.c -lm && ./fujii_scan 4 0 7 20 2 3 4 6 8 12 16 20 24 32 48 64 96 128 192 256 512 1024 2048 4096 > fujii_q4_1e7.out && python3 fujii.py 4 > fujii_crit_q4.out`; `python3 selfcheck.py`; `python3 crosscheck.py 4`; `python3 crosscheck.py 3`. Total CPU time about 45 minutes on one core.

---------------------------------------------------------------------------

## 11. Review history

Referee report: `REFEREE-REPORT.md` (2026-11-24, adversarial, in-house). Verdict: no fatal error; Lemma 1,
Theorems 1–6, Corollary 4a and its χ_{−3} analogue, Theorem 5 and the Fujii counterexamples certified as
proved (given [E] Theorem 1 and [E] Lemma 3(a)); four repairable items and nine minor ones. Resolution:

| item | referee's point | resolution |
|---|---|---|
| REP-1 | the certificates A1, A2 test only at primes, and the last prime is `4 194 301 < 2^22`; the segment `[4 194 301, 2^22]` is not covered for `Θ_1 > 0` and for `(u − ϑ_o(u))/√u ≤ c_1` | `cert_abel.py` has a new part A0 at the endpoint `T = 2^22`, in 256-bit Arb: `Θ_1(T) = 5 892 140 560.76018` (χ_{−4}), `6 049 850 364.71940` (χ_{−3}), both certified `> 0`; `(T − ϑ_o(T))/√T = 1.01442484837129` resp. `1.01462282938111`, certified `≤ c_1`; `ϑ_o(T) < T`. These agree with the referee's double-precision values (`5.89·10^9`, `6.05·10^9`, `1.0144`, `1.0146`). `cert_abel.out` regenerated (ALL CERTIFIED: True). (P1), Lemma 7 and the §5 table say so. A3 needs no endpoint (`|Θ|/√u` decreases between primes). |
| REP-2 | the "certified" Taylor certificates evaluated the remainder R, the prime tails and the tail argument in IEEE double inflated by `10^-6` | Done the preferred way: `certify_grid.py arb` now computes the tails (`c^m e^{−c} + x Γ(m+1,c)`, `c = P_j/x_j` a ball), the remainder R, the weight and every comparison (`low > 0`, left end `> e²`, `N_C > 0.88`) in Arb; floats only choose the integer cut-off `P_j` and print. All four certificates were rerun: the certified ranges are unchanged (χ_{−4}: `[69.783059, 2001.013414]` and `[1995.262315, 5 011 872.336]`; χ_{−3}: `[13.528512, 2001.013414]` and `[1995.262315, 5 011 872.336]`; 5 s and 5.3 min, resp. 7 s and 5.1 min). §0, §5 updated. The scan mode is unchanged and remains labelled floating point. |
| REP-3 | the LI statements appeared as facts in the status paragraph and in Corollary 4a, while §9.7 called the argument a sketch | Proved: new Proposition 8 with a full proof (§2.11): under GRH, `N_C(x) = √π/2 + Z(x) + o(1)`; hence under GRH *alone* `liminf N_C ∈ [0.8858349, 0.8858360]` (so no constant above `0.885836` is admissible — LI is not needed for that), and under LI `liminf N_C = √π/2 − S_Γ` by Kronecker's theorem. The status paragraph and Corollary 4a now quote Proposition 8. This proof was written after the report and is itself unrefereed (§9.7). |
| REP-4 | weight-class wording | Remark after Theorem 1 now states the decay hypotheses `Σ_p w(p) < ∞`, `t·w(t) → 0` and why they are needed; (F1) reads "convex on `[3,∞)` for every `x > 0` if and only if `α ≤ 1`", with the referee's example `α = 2`, `x ≤ 3√2`; Theorem 1 and Theorem 6 state the χ_{−3} domain `[2,∞)`; "height 24" → `23.955`; "Tonelli" → "Fubini" in Lemma 5(iii); the floating-point parentheticals were removed from Theorem 5. |
| M1 | floating-point claims inside the statement of Theorem 5 | moved to the "Numerical truth" paragraph after the theorems, marked "not part of any theorem". |
| M2 | "Tonelli" | "Fubini's theorem", with the absolute convergence stated. |
| M3 | `ε_a √π/2 = 7.63·10^-4` | corrected to `7.622·10^-4`. |
| M4 | "height 24"; Platt–Trudgian's Turing verification | `23.955` (χ_{−3}: `25.418`); Platt–Trudgian's verification up to height 1127 is mentioned as independent and not used. |
| M5 | Corollary 7 cited Theorems 1 and 4 for the bound | now "Theorem 1 and Corollary 4a". |
| M6 | cite Suzuki's Theorem 3(2) and Knapowski–Turán | done in the remark after Theorem 2 and in §7; I read Suzuki's §1.2 for this (Theorem 3, the historical paragraph, Theorem 4); Knapowski–Turán is cited through Suzuki, not read. Two further facts from that paragraph were added: the log-weighted limit form is also due to Hardy–Littlewood and Landau, and Fujii conjectured his statement for all `α > 0`. |
| M7 | the digest's caveats (comparator stage still running; unread vacuity allegation) | one sentence added to Corollary 7 (checked against `../GRH-PROOF-DIGEST.md`, lines on `lake comparator` and the PhilArchive paper). |
| M8 | the maximum of `N_C` quoted was the grid maximum | true maximum `1.1087617` at `x ≈ 6356` (my recomputation on a mesh of step 1; the referee has 6360) now given next to the grid maximum, also for χ_{−3} (`1.122098` at `x ≈ 2129`). |
| M9 | "κ > 0 for y ≥ 413 / 750" is a grid statement | reworded in Theorem 3 as a grid statement that is not used. |

Also changed, not requested by the referee: Theorem 3's sentence "`E(y)/√y` tends to 0" was wrong for the finite
ladder as printed (limit `1.13·10^-9`); the statement now allows the ladder to be continued and §2.6 proves the
limit (§9.11). `certify_grid.py arb` now reports `>` instead of `≥` 0.88 (strict ball comparison).

Not done, and why: (1) the closed-form evaluations of `bounds_abel.py` are still 30-digit mpmath, not interval
arithmetic (margins `≥ 0.003`; the referee considers this sufficient and re-implemented them); (2) the
scan-mode certificate (to `10^8`) is still floating point with an a-priori budget — it would need directed
rounding in the C program — and no theorem uses it; (3) Hardy–Littlewood, Landau, Fujii and Knapowski–Turán
were not read in the original; (4) the analogue of Proposition 8 for Fujii's weights (which would turn the
"expected failure for `α > α_*` under LI" into a proof) was not written; (5) Proposition 8 has not been seen by
a referee.

