# The Cesàro mean of Chebyshev's race is positive for every x > 3 (GRH only)

Working notes, 2026-11-24; revised the same day after the referee report
`REFEREE-REPORT.md` (see §11, "Review history"). Everything in this directory
except that report was produced for this note; nothing outside
`monotonicity/explicit/` was modified.

Status in one paragraph. Under GRH for `L(s,χ_4)` alone, the asymptotic
"for all sufficiently large x" statement of Theorem A (pair m = 0, m' = 1) is
turned into a statement for **every** `x > 3`, with an explicit lower bound.
The analytic part of the proof works for `t ≥ 10^4` (it would work from
`t ≥ 4892`); the finite part needed by the proof (`t ≤ 10^4`) is done in exact
integer arithmetic, and a floating-point sieve covers `3 < x ≤ 10^10`, so the two
ranges overlap by six orders of magnitude. The same is done for `χ_3` (every
`x > 2`). The key step is Ingham's second-order partial summation (Lemma 2, the
identity Johnston uses as his Lemma 2.7) together with the observation that the
integrated function `Θ_1` is one-signed, so that the positive kernel loses
nothing — in particular no `1/log x`. The same observation gives the whole cone
of weight pairs `m ≤ 0, m < m' ≤ 1` (Theorem 5). What is needed for a general
pair `(m, m')` is stated in §8.

---------------------------------------------------------------------------

## 0. Notation

`q ∈ {4, 3}`, `χ = χ_q` the nontrivial (real, primitive, odd) character mod q,
`s(p) = −χ(p)`:

    q = 4:  s(p) = +1 (p ≡ 3 mod 4), −1 (p ≡ 1 mod 4), 0 (p = 2)
    q = 3:  s(p) = +1 (p ≡ 2 mod 3), −1 (p ≡ 1 mod 3), 0 (p = 3)

`ℓ` = the prime dividing q (2 resp. 3); `p_0` = the least prime not dividing q
(3 resp. 2).

    D_m(x)  = Σ_{p≤x} s(p) p^m                      (D_0 = π(x;4,3) − π(x;4,1) for q = 4)
    R(x)    = Σ_{p≤x} s(p) (1 − p/x) = D_0(x) − D_1(x)/x = (1/x) ∫_0^x D_0(t) dt
    N(x)    = R(x) log x / √x
    Θ(t)    = Σ_{p≤t} s(p) log p = −θ(t,χ)
    Θ_1(t)  = ∫_0^t Θ(u) du = Σ_{p≤t} s(p) log p · (t − p)
    G(t)    = Θ_1(t) / t^{3/2}
    θ(y)    = Σ_{p≤y} log p,     ϑ_o(y) = Σ_{p≤y, p∤q} log p = θ(y) − (log ℓ)·1_{y≥ℓ}
    ψ(t,χ)  = Σ_{n≤t} χ(n)Λ(n),  ψ_1(t,χ) = ∫_0^t ψ(u,χ) du = Σ_{n≤t} χ(n)Λ(n)(t − n)

`R` is continuous (the term of p vanishes at x = p) and equals `A − B/x` between
consecutive primes; `Θ_1` is continuous and piecewise linear. "GRH(χ)" means:
all zeros of `L(s,χ)` in `0 < Re s < 1` have real part 1/2. No hypothesis on
`ζ` is used anywhere, and no linear-independence hypothesis.

---------------------------------------------------------------------------

## 1. Results

Constants (all digits shown are correct; see §6 for their provenance):

                         χ_4                       χ_3
    σ_0²(χ)              0.155567979923586…        0.113229969857472…
    c_0 = 2/3 − σ_0²     0.511098686743…           0.553436696809…
    c_1                  2.071193                  2.081942
    Y_1                  10^10                     10^10
    ε_1 = 1/log³ Y_1     8.191·10^-5               8.191·10^-5
    c_θ = 1 + ε_1        1.0000820                 1.0000820

`σ_0²(χ) = log(q/π) − γ + 2 (L'/L)(1,χ)` is a closed form; under GRH it equals
`Σ_ρ |ρ|^{-2}` (Lemma 3(b)). For χ_4 it is `γ + 6 log 2 + 5 log π − 8 log Γ(1/4)`.

Define, for `t ≥ p_0`,

    E_odd(t) = (3/4) t^{-1/6} + (5/6) t^{-3/10} + t^{-5/14} · max(0, (log t / log p_0 − 5)/2),

    L_χ(t)   = c_0 − (4/5) c_1 t^{-1/4} − c_θ E_odd(t)
                   − 1_{t > Y_1²} · [ (2/3) ε_1 + (log ℓ) t^{-1/2} ].

`L_χ` is nondecreasing on `[4000, Y_1²]` and on `(Y_1², ∞)`;
`L_{χ_4}(10^4) = 0.0681`, `L_{χ_3}(10^4) = 0.0182`, `L_χ(Y_1²+) = 0.5107` resp.
`0.5530`, and `L_χ(t) → c_0 − (2/3)ε_1 = 0.51104` resp. `0.55338`.

**Theorem 1 (integrated θ-race).** Assume GRH(χ_4). Then

    Θ_1(t) ≥ L_{χ_4}(t) · t^{3/2}   for every t ≥ 3,     and     Θ_1(t) > 0   for every t > 3.

**Theorem 2 (Cesàro mean of Chebyshev's race).** Assume GRH(χ_4). Then for every `x > 3`

    R(x) = (1/x) ∫_0^x (π(t;4,3) − π(t;4,1)) dt  ≥  Θ_1(x) / (x log x)  >  0,

and for every `x ≥ 10^4` and every `t_1 ∈ [10^4, x]`

    R(x) ≥ (√x / log x) · [ L(x) + (L_*(t_1, x) / log x) · (8/3 − 2 (t_1/x)^{1/2} − (2/3)(t_1/x)^{3/2}) ],

where `L = L_{χ_4}` and `L_*(t_1,x) = inf_{t_1≤t≤x} L(t)` (`= L(t_1)` if `x ≤ 10^20`,
`= min(L(t_1), L(10^20+))` otherwise). In particular `N(x) ≥ L(x)`.

**Corollary 3.** Assume GRH(χ_4). Then `R(x) > √x / (2 log x)` for every `x ≥ 23`.
(Under GRH and LI in the usual sense — the positive ordinates are linearly
independent over Q, in particular the zeros are simple —
`liminf N = liminf G = 2/3 − Σ_{γ>0} 2/|ρ(ρ+1)| ≈ 0.5129`, so for all large x the
constant 1/2 cannot be replaced by anything larger than that. Here
`liminf N = liminf G` because `N − G = (log x / x^{3/2}) ∫ Θ_1 h''` is `≥ 0` and
`O(1/log x)`. The theorem's own limiting constant is `c_0 − (2/3)ε_1 = 0.51104`.)

**Corollary 4 (inheritance of leads, m = 0 versus m' = 1).** Assume GRH(χ_4).
For every `x > 3`:

    D_1(x) ≥ 0  ⟹  D_0(x) > 0,          D_0(x) ≤ 0  ⟹  D_1(x) < 0,

and more precisely `D_0(x) ≥ D_1(x)/x + √x/(2 log x)` for `x ≥ 23`. So `x_1 = 3^+`:
every lead (or tie) of the class 3 mod 4 in the race weighted by p is a strict
lead in Chebyshev's classical race, and every lead (or tie) of the class 1 mod 4
in the classical race is a strict lead of class 1 in the p-weighted race — at
every x > 3, not just for large x.

**Theorem 5 (the cone `m ≤ 0, m < m' ≤ 1`; m' may be negative).** Assume GRH(χ_4). For all real
`m ≤ 0` and `m < m' ≤ 1`, with `κ = m' − m`, and every `x > 3`:

    D_m(x) − x^{-κ} D_{m'}(x)  ≥  κ x^{m−1} Θ_1(x) / log x  ≥  κ x^{m+1/2} L(x) / log x,

the first lower bound being `> 0`. Hence `D_{m'}(x) ≥ 0 ⟹ D_m(x) > 0` for every `x > 3`.

**Theorem 6 (equivalence).** The following are equivalent:
(a) GRH(χ_4); (b) `∫_0^x (π(t;4,3) − π(t;4,1)) dt > 0` for every `x > 3`;
(c) the same integral is `≥ 0` for all sufficiently large x.

**Theorem 7 (χ_3).** Theorems 1, 2, 5, 6 and Corollary 4 hold verbatim for
`χ_3` (classes 2 versus 1 mod 3, the race starting at p = 2) with `x > 3`,
`t > 3`, `t ≥ 3` replaced by `x > 2`, `t > 2`, `t ≥ 2`, with `L = L_{χ_3}`, and
Corollary 3 holds for every `x ≥ 6` (it fails at `x = 5`, where `N = 0.4319`, and
holds from `x ≈ 5.343` on).

Numerical truth (sieve, §5): for q = 4, `R(x) > 0` on `3 < x ≤ 10^10` (exact
integer test) with `min_{x≥5} R = 2/7` at x = 7; `inf N` over `[10^4, 10^10]` is
`0.682694` at the prime `x = 1 070 765 639`; `inf N` over `[11, 10^10]` is
`N(11) = 0.394360`; `inf G` over `[2^13, 10^10]` is `0.5196`.

---------------------------------------------------------------------------

## 2. Proofs

### Lemma 2 (summation identity; no hypothesis)

Let `x > p_0` and `F ∈ C²[p_0, x]`. Then

    Σ_{p≤x} s(p) log p · F(p) = F(x) Θ(x) − F'(x) Θ_1(x) + ∫_{p_0}^x Θ_1(t) F''(t) dt.

*Proof.* The left side is `∫_{p_0^-}^x F dΘ = F(x)Θ(x) − ∫_{p_0}^x Θ F' dt`. Since
`Θ_1' = Θ` a.e., `Θ_1` is absolutely continuous and `Θ_1(p_0) = 0`, a second
integration by parts gives `∫ Θ F' = Θ_1(x)F'(x) − ∫ Θ_1 F''`. ∎

Take `F(t) = h(t) := (x − t)/log t`. Then `h(x) = 0`, `h'(x) = −1/log x`, and

    h''(t) = 2/(t log² t) + (x − t) [ 1/(t² log² t) + 2/(t² log³ t) ]
           = (x + t)/(t² log² t) + 2(x − t)/(t² log³ t)   > 0     (1 < t ≤ x).

Since `Σ_{p≤x} s(p) log p · h(p) = Σ s(p)(x − p) = x R(x)`:

    x R(x) = Θ_1(x)/log x + ∫_{p_0}^x Θ_1(t) h''(t) dt.                                  (2.1)

(Checked numerically to 22 digits at x = 200.5 and x = 1000.25, `verify_small.py`.)

**Consequence.** If `Θ_1 ≥ 0` on `[p_0, x]` then `R(x) ≥ Θ_1(x)/(x log x)`. So the
positivity of `R` — a sum with weights `1 − p/x` and no `log p` — follows from the
positivity of the *one* function `Θ_1`, with no loss.

*Attribution.* Lemma 2 is not new. It is the second-order partial summation of
Ingham (*The distribution of prime numbers*, pp. 103–104, as cited by Johnston),
and Johnston's Lemma 2.7 is the same identity, with the same positive kernel
`1/(t² log² t) + 2/(t² log³ t)`, acting on the integrated error
`R_1(x) = ∫_2^x (ψ(t) − t) dt`. In his setting `R_1` is not one-signed, so the
kernel term has to be estimated ((3.1)–(3.6) of his paper, applied for
`x > 10^19`, with Büthe's verification used below). The only new observation
here is that the function integrated against the kernel, `Θ_1`, *is* one-signed,
so the kernel term can simply be dropped; this is why the analytic argument can
start at `10^4`.

### Lemma 3 (zeros)

Let χ be a primitive odd character mod q (so `L(0,χ) ≠ 0` and the trivial zeros
are simple, at `s = −1, −3, −5, …`). Put `b = (L'/L)(0,χ)` and let `b_0` be the
constant term of the Laurent expansion `(L'/L)(s,χ) = 1/(s+1) + b_0 + O(s+1)`.

(a) (Explicit formula for ψ_1.) For `x ≥ 1`,

    ψ_1(x,χ) = − Σ_ρ x^{ρ+1}/(ρ(ρ+1)) − b x + log x + 1 + b_0 − Σ_{n≥1} x^{-2n}/(2n(2n+1)),

the sum over all nontrivial zeros (with multiplicity) being absolutely convergent.

(b) Under GRH(χ), `Σ_ρ 1/|ρ(ρ+1)| < Σ_ρ 1/|ρ|² = log(q/π) − γ + 2 (L'/L)(1,χ) =: σ_0²`
(χ real).

*Proof of (a).* Let `x > 1` and `c = 2`. Then
`ψ_1(x,χ) = (1/2πi) ∫_{(c)} (−L'/L)(s,χ) x^{s+1}/(s(s+1)) ds`
(Mellin inversion of `(x−n)_+`, term by term; the integral converges absolutely).
Integrate over the rectangle with vertices `c ± iT`, `−U ± iT`, where U is an
even integer (so the left side passes midway between the trivial zeros
`−U ± 1`) and `T = T_k → ∞` is chosen with `|T_k − γ| ≫ 1/log T_k` for every
ordinate γ (possible because there are `≪ log(qT)` ordinates in `[T, T+1]`).
Two standard bounds (Davenport, *Multiplicative Number Theory*, 3rd ed.,
Ch. 19, the paragraphs between eqs. (2) and (7), which treat the harder kernel
`x^s/s`; the same as for ζ in Ch. 17) are used:

    (i)   (L'/L)(σ ± iT_k, χ) ≪ log²(qT_k)            uniformly for −1 ≤ σ ≤ 2;
    (ii)  (L'/L)(s, χ) ≪ log(q(|s| + 2))               for Re s ≤ −1, |s + (2n+1)| ≥ 1/2 for all n ≥ 0

((ii) follows from the functional equation). *Horizontal sides.* On `Im s = ±T_k`
we have `|s(s+1)| ≥ T_k²`, so by (i) and (ii) the two integrals are

    ≪ T_k^{-2} ∫_{-∞}^{2} [ log²(qT_k) + log(|σ| + 2) ] x^{σ+1} dσ ≪_x log²(qT_k) / T_k²,

uniformly in U (here `x > 1` makes `∫_{-∞}^2 log(|σ|+2) x^{σ+1} dσ` finite); this
tends to 0 as `k → ∞`. *Left side.* On `Re s = −U`, by (ii),

    | ∫ | ≪ x^{1−U} ∫_{-∞}^{∞} log(q(U + |t| + 2)) / (U² + t²) dt ≪ x^{1−U} log(qU) / U → 0   (U → ∞).

Letting first `U → ∞` and then `k → ∞`, `ψ_1(x,χ)` equals the sum of all residues
to the left of `Re s = c`; the sum over ρ needs no grouping, since
`Σ_ρ |x^{ρ+1}/(ρ(ρ+1))| ≤ x² Σ_ρ |ρ|^{-2} < ∞`. The residues are: at `s = ρ`,
`−x^{ρ+1}/(ρ(ρ+1))`; at `s = 0`, `−b x`; at `s = −(2n+1)`, `n ≥ 1`,
`−x^{-2n}/(2n(2n+1))`; at the double pole `s = −1`, with
`−L'/L = −1/(s+1) − b_0 + …` and `x^{s+1}/s = −1 − (s+1)(log x + 1) + …`, the
residue is `log x + 1 + b_0`. Both sides of the formula are continuous on
`[1, 2]` (the two series converge uniformly there), so it also holds at `x = 1`. ∎

Three independent checks of this formula (`constants.py`, `constants.out`):
(i) its x-derivative is Davenport's formula for `ψ(x,χ)`, Ch. 19 eq. (3) with
`a = 1`: `ψ(x,χ) = −Σ x^ρ/ρ − b(χ) + Σ_{m≥1} x^{1−2m}/(2m−1)`;
(ii) evaluating at `x = 1` (`ψ_1(1,χ) = 0`) the constant must equal
`Σ_ρ 1/(ρ(ρ+1)) + b + 1 − log 2`, and the Hadamard product gives
`Σ_ρ 1/(ρ(ρ+1)) = b_0 − b + log 2`, so the constant is `1 + b_0` — the same
value; numerically `Σ_ρ 1/(ρ(ρ+1))` from the certified zero lists is
`−0.150312` (χ_4), `−0.111224` (χ_3) against `b_0 − b + log 2 = −0.150309`,
`−0.111220`;
(iii) direct comparison with the prime-power sum, using 511 / 537 zeros, at
small x, where the truncation of the zero sum is small enough for the test to
say something ("bound" is `x^{3/2} Σ_{γ>T} 2/γ²` from the Riemann–von Mangoldt
density; a wrong constant term would show up as a difference of order 1):

    x       ψ_1 direct − formula (χ_4)   bound     ψ_1 direct − formula (χ_3)   bound
    1       +3.48·10^-3                  0.0035    +3.09·10^-3                  0.0031
    3       +1.64·10^-3                  0.018     −6.5·10^-7                   0.016
    5.5     −1.2·10^-6                   0.045     +1.6·10^-5                   0.040
    10      +7.5·10^-5                   0.11      +1.7·10^-4                   0.098
    20      −1.6·10^-3                   0.31      −8.9·10^-4                   0.28
    50      −3.2·10^-4                   1.2       +3.1·10^-3                   1.1
    200     −8.3·10^-2                   9.9       +3.6·10^-3                   8.8

(At `x = 1` the difference is the non-oscillating tail `Σ_{γ>T} 1/(ρ(ρ+1))`. The
first version of this note quoted comparisons at `x = 10^3, 10^4`; there the
allowed error is `10^2`–`10^3` and the comparison is vacuous, so they were
removed.)

Values: χ_4: `b = 0.78318878541`, `b_0 = −0.06026691418`; χ_3: `b = 0.94819882667`,
`b_0 = 0.14383164456`, from the functional equation
(`b = −log(q/π) + γ + log 2 − (L'/L)(1,χ)`,
`b_0 = −log(q/π) + γ/2 − ψ(3/2)/2 − (L'/L)(2,χ)`, χ real) and confirmed by
direct evaluation of `L'/L` near 0 and −1.

*Proof of (b).* `|ρ + 1| > |ρ|` since `Re ρ > 0`. By Davenport Ch. 12, eqs.
(17)–(18), `(L'/L)(s,χ) = −(1/2)log(q/π) − (1/2)(Γ'/Γ)((s+1)/2) + B(χ) + Σ_ρ (1/(s−ρ) + 1/ρ)`
and `Re B(χ) = −Σ_ρ Re(1/ρ)`. At `s = 1` (real parts; χ real so the left side is
real): `Σ_ρ Re 1/(1−ρ) = (L'/L)(1,χ) + (1/2)log(q/π) − γ/2`. Under GRH,
`Re 1/(1−ρ) = (1/2)/|ρ|²`. ∎

(Check: 511 zeros of χ_4 plus the Riemann–von Mangoldt tail give 0.1555712 against
0.1555680; `RESEARCH-CONTEXT.md` §5 has the same certified value
0.15556797992358592886. For χ_3: 0.1132341 against 0.1132300.)

*Remark.* The true value of `S = Σ_ρ 1/|ρ(ρ+1)|` is about `0.15380` (χ_4) and
`0.11256` (χ_3); the proof uses the closed-form majorant `σ_0²`, losing 0.0018.

### Lemma 4 (inputs on θ)

(a) (Computation.) For `0 < y ≤ Y_1 = 10^10`:  `y − c_1 √y ≤ ϑ_o(y) < y`, with
`c_1 = 2.071193` (q = 4), `2.081942` (q = 3).

Basis of (a). The supremum of `(y − ϑ_o(y))/√y` is attained as `y → 1423^-`
(`2.0711926180` resp. `2.0819411967`). For `y ≤ 2^22` both inequalities are
certified in ball arithmetic (`cert_small.py`, part C). For `2^22 < y ≤ 10^10`
they rest on the `long double` sieve (§5): there the supremum of
`(y − ϑ_o(y))/√y` over `y ≥ 1427` is `1.98197` resp. `1.99271` (as `y → 1427^-`),
so the lower inequality holds with margin `≥ 0.089√y ≥ 180`, and
`ϑ_o(y) − y ≤ −0.26√y`; the worst-case accumulated rounding error of θ is `0.21`.
This part of (a) is therefore safe but not certified in interval arithmetic;
the Remark after this lemma (variant C) shows how to do without it.

Remarks on the literature. Subtracting the `log ℓ/√1423` contributed by the
missing prime gives `2.05282` for θ itself, which I recall as the classical
Rosser–Schoenfeld constant (no equation number checked). J. Büthe, "An analytic
method for bounding ψ(x)", Math. Comp. 87 (2018), Theorem 2 (read in
arXiv:1511.02032): (1.7) `x − ϑ(x) > 0.05√x` for `1 ≤ x ≤ 10^19` is an
independent source for the upper inequality `ϑ_o(y) < y`. Its (1.6),
`x − ϑ(x) ≤ 1.95√x` "for `1423 ≤ x ≤ 10^19`", is **not** used: as printed in the
arXiv text it fails on the short interval `(1426.48, 1427)`, where
`(x − θ(x))/√x` rises to `1.96363` (my computation: `θ(1423) = 1352.82275`);
on `[1427, 3·10^6]`, the range I checked for θ itself, it holds (next largest
value `1.93722` as `x → 19373^-`), so the threshold is presumably a misprint for
1427.
In any case these cited bounds plus the shift by `log ℓ` do not give `c_1` for
small y (at `y = 4` the shift is `log 2/√4 = 0.35`): below 1427, (a) rests on the
computation alone.

(b) (Cited.) For `y ≥ Y_1`: `|θ(y) − y| < y/log³ y ≤ ε_1 y`, `ε_1 = 1/log³(10^10)`.
Source: S. Broadbent, H. Kadiri, A. Lumley, N. Ng, K. Wilk, "Sharper bounds for
the Chebyshev function θ(x)", Math. Comp. 90 (2021); read in arXiv:2002.11068v2.
Corollary 11.2, eq. (3.39): `|ϑ(x) − x| < M_k x/log^k x` for `x ≥ X_0`, with the
row `k = 3, M_3 = 1, X_0 = 89 967 803`. (This is the statement of Dusart's
Theorem 4.2 / arXiv Theorem 5.2, for which BKLNW say they "provide valid
proofs".) On the input: their Appendix A.1.1 says the computations use the
partial verification of RH up to `H = 2 445 999 556 030` "announced by Gourdon
and Wedeniwski", and that Platt and Trudgian, "The Riemann hypothesis is true
up to 3·10^12" (Bull. London Math. Soc.; their reference [38]), verify
`H = 3 000 175 332 800` with "the more rigorous techniques" and "thus confirm the
earlier calculations". Since `3 000 175 332 800 > 2 445 999 556 030`, the only
unproved-at-the-time input of BKLNW is a consequence of Platt–Trudgian's
rigorous verification. (The same remark applies to Dusart's arXiv:1002.0442,
Theorem 5.2, which the first version of this note cited and which rests on
Gourdon's computation; BKLNW is now the cited source because it is published,
states its input explicitly, and re-proves the row we need.)

(b') (Fallback, classical.) `x − x/log x < θ(x) < x + x/(2 log x)` for `x ≥ 41`:
this is the form in which Johnston's Lemma 2.2, (2.2), quotes Rosser–Schoenfeld,
"Approximate formulas for some functions of prime numbers", Illinois J. Math. 6
(1962), (3.15)–(3.16); the original was not consulted. It is used only for
`y ≥ 10^10`. With (b') instead of (b): `ε_low = 1/log Y_1 = 0.0434`,
`ε_up = 0.0217`, `c_θ = 1.0217`, and `(2/3)ε_1` is replaced by `(2/3)ε_low = 0.0290`.

From (a) and (b): `ϑ_o(y) ≤ c_θ y` for all `y > 0`, and
`ϑ_o(y) ≥ (1 − ε_1) y − log ℓ` for `y > Y_1`.

*Remark (variant C: certified inputs only).* The uncertified part of (a) can be
avoided altogether. Take `Y_1^C = 2^22 = 4 194 304`, where (a) is certified in
ball arithmetic, and for `y ≥ 2^22` use another row of the same Corollary 11.2
of BKLNW, `k = 2, M_2 = 0.2, X_0 = 3 594 641`:
`|θ(y) − y| < 0.2 y/log² y ≤ ε_1^C y` with `ε_1^C = 0.2/log²(2^22) = 8.601·10^-4`,
`c_θ^C = 1.000860`. Let `L^C` be the function L of §1 with `(Y_1, ε_1, c_θ)`
replaced by `(Y_1^C, ε_1^C, c_θ^C)` (so the indicator becomes `1_{t > 2^44}`,
`2^44 = 1.759·10^13`). The proofs below go through word for word, so **every
statement of §1 holds with `L^C` in place of L** (and `2^44` in place of `10^20`
in Theorem 2). Numbers (`bounds.out`, variant C, χ_4 / χ_3):
`L^C(10^4) = 0.0679 / 0.0179`; positivity thresholds `4904 / 8548`;
`L^C(2^44 +) = 0.5048 / 0.5469`; limit `0.51053 / 0.55286`; in the analytic part of
Corollary 3 (`x ≥ 10^10`) the margins are `0.0287 / 0.0726` on the grid and
`0.0105 / 0.0529` in the tail. `L^C ≤ L`, and
`L − L^C` is below `3·10^-4` for `10^4 ≤ t ≤ 2^44` and below `8·10^-4` for all
`t ≥ 10^4`. In
this variant Theorems 1, 2, 5, 6, 7 and Corollaries 3, 4 use only: the exact and
ball-certified finite computations of `cert_small.py` and `cert_gaps.py`, BKLNW
Corollary 11.2, and the analysis. The main text keeps `Y_1 = 10^10` only because the constants
are slightly better.

### Proof of Theorem 1 (and of its χ_3 analogue)

*Step 1: decomposition.* Since `ψ(u,χ) = Σ_{k≥1} θ(u^{1/k}, χ^k)` and `χ^k` is χ
for odd k and the principal character for even k,

    Θ(u) = −ψ(u,χ) + ϑ_o(u^{1/2}) + Σ_{k≥3 odd} θ(u^{1/k},χ) + Σ_{k≥4 even} ϑ_o(u^{1/k}).

Integrate over `[0,t]`, `t ≥ p_0`. The last sum is `≥ 0` and is dropped.

*Step 2: zeros.* By Lemma 3, under GRH (`|t^{ρ+1}| = t^{3/2}`):

    −ψ_1(t,χ) ≥ −σ_0² t^{3/2} + (b t − log t − 1 − b_0).

The bracket is increasing for `t > 1/b` and positive at `t = p_0` (0.311 for
χ_4 at t = 3; 0.059 for χ_3 at t = 2), hence `−ψ_1(t,χ) ≥ −σ_0² t^{3/2}` for `t ≥ p_0`.

*Step 3: prime squares (the bias).* `B_2(t) := ∫_0^t ϑ_o(u^{1/2}) du = ∫_0^{√t} ϑ_o(y) 2y dy`.
If `√t ≤ Y_1`, Lemma 4(a) gives `B_2(t) ≥ ∫_0^{√t} (2y² − 2c_1 y^{3/2}) dy = (2/3)t^{3/2} − (4/5)c_1 t^{5/4}`.
If `√t > Y_1`, split at `Y_1` and use Lemma 4(b) on `[Y_1, √t]`:

    B_2(t) ≥ (2/3)Y_1³ − (4/5)c_1 Y_1^{5/2} + (2/3)(1 − ε_1)(t^{3/2} − Y_1³) − (log ℓ)(t − Y_1²)
           ≥ (2/3)(1 − ε_1) t^{3/2} − (4/5) c_1 t^{5/4} − (log ℓ) t.

*Step 4: odd prime powers k ≥ 3.* `|θ(y,χ)| ≤ ϑ_o(y) ≤ c_θ y`, and `θ(y,χ) = 0` for
`y < p_0`. So `|∫_0^t θ(u^{1/k},χ) du| ≤ c_θ (k/(k+1)) t^{1+1/k}`, and the term vanishes
unless `k ≤ log t/log p_0`. For k = 3 and k = 5 keep the constants 3/4 and 5/6;
for odd `k ≥ 7` use `(k/(k+1)) t^{1+1/k} ≤ t^{8/7}`; there are at most
`max(0, (log t/log p_0 − 5)/2)` such k. Hence the total is at most
`c_θ E_odd(t) t^{3/2}`.

*Step 5.* Adding up: `Θ_1(t) ≥ L_χ(t) t^{3/2}` for all `t ≥ p_0`.

*Step 6: monotonicity and positivity.* `t^{-1/4}`, `t^{-1/6}`, `t^{-3/10}` decrease,
and `t^{-5/14}(log t/log p_0 − 5)` decreases once `log t > 5 log p_0 + 14/5`, i.e.
`t > 3996` (q = 4), `t > 527` (q = 3). So `L_χ` is nondecreasing on `[4000, Y_1²]`
and on `(Y_1², ∞)`, and `inf_{t ≥ 10^4} L_χ(t) = min(L_χ(10^4), L_χ(Y_1²+)) = 0.0681`
(χ_4), `0.0182` (χ_3). Hence `Θ_1(t) > 0` for `t ≥ 10^4`. (The exact thresholds of
the analytic bound are `t ≥ 4892` for χ_4 and `t ≥ 8526` for χ_3.)

*Step 7: small t.* `Θ_1` is piecewise linear with nodes at the primes p ∤ q,
`Θ_1(p_0) = 0`, and its first slope is `log p_0 > 0`; so `Θ_1 > 0` on `(p_0, 10^4]`
as soon as `Θ_1(p) > 0` at every race prime `p_0 < p ≤ 10^4`. This is checked
**exactly** (`cert_small.py`, part A): `Θ_1(p) > 0` is the inequality between
integers

    Π_{r ≤ p, s(r) = +1} r^{p−r}  >  Π_{r ≤ p, s(r) = −1} r^{p−r},

verified with integers of up to `3.6·10^7` bits at all 1227 race primes
`p_0 < p ≤ 10^4`, for both characters: no failure. The smallest values are
`Θ_1(7) = 4 log 3 − 2 log 5 = 1.17557` (χ_4; `G(7) = 0.0635`) and
`Θ_1(5) = 3 log 2 = 2.07944` (χ_3; `G(5) = 0.1860`). (Beyond what the proof needs:
`Θ_1(p) > 0` is certified in 256-bit ball arithmetic for `p ≤ 2^22`, part B of the
same script, and found in `long double` for `p ≤ 10^10`, §5.) ∎

### Proof of Theorem 2

By Theorem 1 and (2.1), for `x > p_0`: `x R(x) ≥ Θ_1(x)/log x > 0`, i.e. `N(x) ≥ G(x) ≥ L(x)`.
For the refinement let `10^4 ≤ t_1 ≤ x`. On `[t_1, x]`, `Θ_1(t) ≥ L_*(t_1,x) t^{3/2}`
with `L_* > 0`, and `h''(t) ≥ (x + t)/(t² log² x)`. Therefore

    ∫_{p_0}^x Θ_1 h'' ≥ (L_*/log² x) ∫_{t_1}^x t^{-1/2}(x + t) dt
                      = (L_*/log² x) [ (8/3) x^{3/2} − 2 x t_1^{1/2} − (2/3) t_1^{3/2} ].  ∎

(The asymptotic expansion of the mean of N is `(2/3)(1 + (8/3)/log x + …)`,
which is why the refinement matters numerically: at `x = 10^10` it is worth 0.048.)

### Proof of Corollary 3

*x ≥ 10^10 (analytic).* Put `φ(u) = 8/3 − 2√u − (2/3)u^{3/2}` and `ρ = 10^{1/50}`.
For `x ∈ [a, ρa]` with `a ≥ 10^6`, Theorem 2 with `t_1 = a/100` gives

    N(x) ≥ β(a) := inf_{[a, ρa]} L  +  ( inf_{t ≥ a/100} L ) · φ(1/100) / log(ρ a),

because `t_1/x ≤ 1/100` and `log x ≤ log(ρa)`. Both infima are explicit, since L is
nondecreasing on `[4000, Y_1²]` and on `(Y_1², ∞)`: the infimum over an interval is
the value at its left end, or the smaller of that and `L(Y_1²+)` if the interval
contains `Y_1²`. `bounds.py` evaluates β on the grid `a = 10^10 ρ^j`, `a < 10^30`:

    β(a) − 1/2 ≥ 0.0293 on the whole grid   (0.0346 at a = 10^10; the minimum is at the top end),

and for `x ≥ 10^30`, `N(x) ≥ inf_{t ≥ 10^30} L(t) = 0.51104`, i.e. margin `0.0110`.
Hence `N(x) > 1/2` for all `x ≥ 10^10`.

On the margins of this double-precision evaluation: they are `≥ 0.029` for
`10^10 ≤ x ≤ 10^30` and `≥ 0.011` beyond. The value 0.011 cannot be improved: it
is `c_0 − (2/3)ε_1 − 1/2`, the distance between the constant 1/2 of the corollary
and the limit of the bound. (For χ_3 the two margins are `0.0731` and `0.0534`; in
variant C they are `0.0287`, `0.0105` for χ_4 and `0.0726`, `0.0529` for χ_3.)
For `x > 10^20` this step uses Lemma 4(b), i.e. BKLNW Corollary 11.2; with the
fallback (b') the limit of L for χ_4 is `0.4821 < 1/2` and the corollary is only
obtained for `x ≤ 10^20`.

*Correction.* Earlier versions of this note applied the grid argument from
`a = 1.514·10^8` on and described its margins as "≥ 0.03". That was wrong near the
lower end: the margin is `1.6·10^-4` at the first grid point `1.514·10^8`,
`7.9·10^-4` at `1.585·10^8` (χ_4), and `3.9·10^-4` at `2.512·10^7` (χ_3); it reaches
0.03 only near `10^10`. The inequality was true, but resting it on a double-
precision evaluation with a margin of `10^-4` was not sound practice. The grid
argument is therefore no longer used below `10^10`; that range is now covered by
the certified computation.

*29 ≤ x ≤ 10^10 (certified computation).* On a gap `[p, p']` between consecutive
race primes with `p ≥ 11`, R is monotone and `log x/√x` decreases, so
`N(x) ≥ min(R(p), R(p')) · log p'/√p'`. `cert_gaps.py` (segmented sieve; exact
integers `D_0`, `D_1`; 64-bit Arb balls for the division, `log` and `√`; every
comparison certain) certifies that `R(p) log p'/√p' > 1/2` and
`R(p') log p'/√p' > 1/2` on all 455 052 501 gaps with `29 ≤ p < p' ≤ 10^10`, and
treats the last piece `[9 999 999 967, 10^10]` in the same way. No failure. The
smallest gap bound is `0.603378720997762 ± 5·10^-16` (the gap starting at 47), so
the smallest certified margin over 1/2 is `0.1033787`. The final state
`D_0 = 5960`, `D_1 = 27 201 020 910 942` agrees with the independent C program.
Run time: 1203 s (`cert_gaps_q4.out`).

*23 ≤ x ≤ 29 (by hand; certified in `cert_small.py`, part D).* On `[23, 29]`,
`R(x) = 2 − 28/x`. On `[23,25]`: `R ≥ 0.7826`, `log x/√x ≥ log 25/5 = 0.6438`,
product `≥ 0.5038`. On `[25,29]`: `R ≥ 0.88`, `log x/√x ≥ 0.6253`, product `≥ 0.550`. ∎

(`N(19) = 0.4977`, so the statement fails at 19; on `[19, 23)`, `R = 1 − 5/x` and
`N = 1/2` is crossed at `x ≈ 19.46`.)

For χ_3: analytic for `x ≥ 10^10` as above (margins `0.0731` on the grid,
`0.0534` in the tail); certified computation (`cert_gaps.py`): the two inequalities
hold on all 455 052 506 gaps with `11 ≤ p < p' ≤ 10^10` and on the last piece, the
smallest gap bound being `0.711388955945845 ± 2·10^-16` (the gap starting at 11),
certified margin `0.2113889`; final state `D_0 = 6264`, `D_1 = 27 948 685 258 617`
as in the C program; run time 1132 s (`cert_gaps_q3.out`). On `[7, 11]`, `R ≡ 1`
and `N = log x/√x ≥ 0.7230`; on `[6, 7]`, `R = 2 − 7/x ≥ 5/6` and
`log x/√x ≥ log 6/√6` (it increases up to `e²`), product `≥ 0.6096`. (`N(5) = 0.4319`;
on `[5, 7)` the value 1/2 is crossed at `x ≈ 5.343`.)

### Proof of Corollary 4

`D_0(x) = D_1(x)/x + R(x)` and `R(x) > 0` for `x > 3`. ∎

### Proof of Theorem 5

Fix `x > p_0` and put `f(t) = t^m − x^{-κ} t^{m'}`, `F = f/log t` on `[p_0, x]`. Then
`Σ_{p≤x} s(p) log p · F(p) = D_m(x) − x^{-κ} D_{m'}(x)`, `F(x) = 0`, and
`F'(x) = f'(x)/log x = −κ x^{m−1}/log x`.

* `f ≥ 0` on `[p_0, x]`.
* `f'(t) = t^{m−1} [m − m'(t/x)^κ] ≤ 0`: if `m' ≥ 0` because `m ≤ 0`; if `m' < 0`
  because `|m'|(t/x)^κ ≤ |m'| < |m|`.
* `f''(t) = t^{m−2} [m(m−1) − m'(m'−1)(t/x)^κ] ≥ 0`: `m(m−1) ≥ 0`; if `0 ≤ m' ≤ 1`
  then `m'(m'−1) ≤ 0`; if `m' < 0` then `m'(m'−1) < m(m−1)` since `u(u−1)` is
  decreasing on `u < 1/2`.
* `1/log t` is positive, decreasing and convex on `t > 1`.

A product of two nonnegative nonincreasing convex functions is convex, so
`F'' ≥ 0`. Lemma 2 and Theorem 1 give
`D_m(x) − x^{-κ}D_{m'}(x) = −F'(x)Θ_1(x) + ∫ Θ_1 F'' ≥ κ x^{m−1} Θ_1(x)/log x > 0`. ∎

(Checked at all primes up to 2^21 for ten pairs, including (−1,0), (−2,1),
(−3,−1), (0,0.01): `verify_small.py`, item 5.)

### Proof of Theorem 6

(a)⟹(b) is Theorem 2; (b)⟹(c) is trivial. (c)⟹(a) is Landau's argument; its
classical prototype for this very character is Hardy–Littlewood and Landau
(1918), see §4. Let `A(x) = xR(x) = Σ_{p≤x} s(p)(x−p)`.

*Step 1 (Mellin transform).* `0 ≤ |A(x)| ≤ x π(x) ≤ x²`, A is continuous, and
for `Re s > 1`, integrating term by term,

    F(s) := ∫_1^∞ A(x) x^{-s-2} dx = Σ_p s(p) ∫_p^∞ (x − p) x^{-s-2} dx = P(s)/(s(s+1)),   P(s) = Σ_p s(p) p^{-s}.

From the Euler product, `log L(s,χ) = Σ_p Σ_k χ(p)^k/(k p^{ks})`, so

    P(s) = −log L(s,χ) + (1/2) log L(2s,χ_0) + H(s),   L(2s,χ_0) = ζ(2s)(1 − ℓ^{-2s}),

with `H(s) = Σ_p Σ_{k≥3} c_k(p) p^{-ks}` (bounded coefficients) analytic in `Re s > 1/3`.

*Step 2 (the real axis).* On the real segment `σ > 1/2`:
`ζ(2σ)(1 − ℓ^{-2σ}) > 0`, and `L(σ,χ) > 0`, because
`L(σ,χ_4) = Σ_{n≥0} ((4n+1)^{-σ} − (4n+3)^{-σ})` and
`L(σ,χ_3) = Σ_{n≥0} ((3n+1)^{-σ} − (3n+2)^{-σ})` are series of positive terms for
every `σ > 0`. Hence `P(s)/(s(s+1))` continues analytically to a neighbourhood of
every real point `σ > 1/2`. **This is where the absence of real zeros is used,
and it cannot be dropped:** a real zero `β > 1/2` of `L(s,χ)` would give P a
logarithmic singularity at β with `P → +∞` as `s → β+`, which is perfectly
consistent with `A ≥ 0`; Landau's theorem would then only give `σ_c = β`. For
`χ_4` and `χ_3` the elementary positivity above settles it; for a general real
character the same proof needs "no real zero in `(1/2, 1)`" as a hypothesis
(compare Zhao's Remark 1).

*Step 3 (Landau).* Suppose `A(x) ≥ 0` for `x ≥ x_0`. By Landau's theorem (for a
function of constant sign near infinity, the abscissa of convergence `σ_c` of
its Mellin transform is a singularity of the transform; Montgomery–Vaughan,
*Multiplicative Number Theory I*, Lemma 15.1) and Step 2, `σ_c ≤ 1/2`. So F, and
therefore `P(s) = s(s+1)F(s)`, is analytic in `Re s > 1/2`.

*Step 4 (conclusion).* In `Re s > 1/2`, `L(2s,χ_0)` is analytic and zero-free and
H is analytic, so `log L(s,χ) = −P(s) + (1/2) log L(2s,χ_0) + H(s)` has a
single-valued analytic continuation there; its exponential is `L(s,χ)`, which
therefore has no zero with `Re s > 1/2`. By the functional equation of the real
primitive character χ (zeros symmetric about `Re s = 1/2`) there is none with
`0 < Re s < 1/2` either. ∎

---------------------------------------------------------------------------

## 3. Cross-check of the analytic bound against the sieve

`bounds.out`, variant B (= Lemma 4(b), BKLNW). "B(x)" is the bound of Theorem 2 with
`t_1` optimised; the last column is the resulting lower bound for `R(x)`.

χ_4:

    x        L(x)      B(x)      N(x) exact   R(x) exact     lower bound for R
    10^4     0.06815   0.06815   0.830478     9.0168         0.74
    10^5     0.23659   0.25661   0.825252     22.66735       7.05
    10^6     0.34322   0.37807   0.878495     63.587553      27.37
    10^7     0.40862   0.45215   0.867385     170.175890     88.71
    10^8     0.44822   0.49582   0.777441     422.047828     269.17
    10^9     0.47217   0.52089   0.710167     1083.682422    794.86
    10^10    0.48673   0.53484   0.746014     3239.897909    2322.78

χ_3:

    x        L(x)      B(x)      N(x) exact   R(x) exact     lower bound for R
    10^6     0.35882   0.39320   0.853793     61.799590      28.46
    10^8     0.48366   0.53415   0.796002     432.124105     289.97
    10^10    0.52739   0.57926   0.798797     3469.131474    2515.72

The bound lies below the truth everywhere, at 43% (10^6), 64% (10^8), 72%
(10^10) of the exact value of R for χ_4. Against the *minimum* of N in each
dyadic range the ratio rises to 0.764 (χ_4, at x = 1 070 765 639) and 0.800
(χ_3). Also `L(2^k) ≤ G_low(k)` in every dyadic range `k = 14..33`, for both
characters and both variants (rigorous consistency test, all passed).

Where the loss comes from at `x = 10^10` (χ_4): `σ_0²` is a worst-case bound for
an oscillating sum (0.156 is subtracted from 2/3, whereas the true value is
`G(10^10) = 0.6485`); odd prime powers bounded trivially (0.019); prime squares
(0.005).

**Slack.** The analytic argument proves `Θ_1(t) > 0` from `t = 4892` (χ_4) and
`t = 8526` (χ_3); the theorem needs the computation only up to `10^4`, where it
is done in exact integer arithmetic; the floating-point computation was carried
to `10^10`. Lowering the analytic threshold further is
possible (bound `θ(y,χ)` for `y ≤ Y_1` by the sieve constant
`|Θ(y)| ≤ 2.154√y` (χ_4), `2.316√y` (χ_3), which replaces `t^{-1/6}` by
`t^{-1/3}`) but pointless.

---------------------------------------------------------------------------

## 4. Relation to prior work

Only sources whose text I read are described, and for each, how much was read.
Where a classical paper is known to me only through a later paper, that is said.

* **A. E. Ingham, *The distribution of prime numbers* (1932), pp. 103–104 and
  Theorem 27.** Not consulted directly; known through Johnston's citations. The
  second-order partial summation (our Lemma 2) and the explicit formula for
  `ψ_1` in the ζ case are from there.
* **G. H. Hardy and J. E. Littlewood, "Contributions to the theory of the Riemann
  zeta-function and the theory of the distribution of primes", Acta Math. 41,
  119–196; E. Landau, "Über einige ältere Vermutungen und Behauptungen in der
  Primzahltheorie" I, II, Math. Z. 1 (1918), 1–24, 213–219.** Not read in the
  original; bibliographic data and statement as given in the introduction and
  reference list of Platt–Trudgian (below), which I did read: Chebyshev's
  conjecture `Σ_{p>2} (−1)^{(p−1)/2} e^{−xp} → −∞` as `x → 0` is equivalent to GRH
  for `χ_4`. This is the classical prototype, for the same character, of the
  equivalence in Theorem 6 (Abel mean and limit `−∞` there; Cesàro mean and sign
  for every x here), and of the Landau mechanism behind (c)⟹(a). The exact form
  of the original statements needs verification against the originals.
* **D. Platt and T. Trudgian, "Fujii's development on Chebyshev's conjecture",
  arXiv:1802.00089.** Read: abstract, introduction, reference list. Fujii showed
  the same equivalence with `e^{−(xp)^α}` for `0 < α < 4.19`; Platt–Trudgian
  enlarge the range of α using computed zeros of `L(s,χ_4)` (the precise upper
  limit is cut off in my text extraction; `../PRIOR-ART.md` gives 20.40442).
* **J. Pintz, "On an assertion of Riemann concerning the distribution of prime
  numbers", Acta Math. Hungar. 58 (1991), 383–387.** Not read; as described by
  Johnston: Pintz asserts that `∫_2^x (π(t) − li(t)) dt < 0` for large x is
  equivalent to RH and likely holds for all `x > 2`, and considers the weights
  `exp(−log² t / y)`. Johnston's theorem "makes explicit a claim of Pintz".
* **D. R. Johnston, "On the average value of π(t) − li(t)", arXiv:2201.06184v2
  (Canad. Math. Bull. 66 (2023) 185–195).** Read in full. Theorem 1.1: RH is
  equivalent to `∫_2^x (π(t) − li(t)) dt < 0` for all `x > 2`; Theorem 1.2: the
  same for `θ(t) − t`. The proof is the template for the present one: Ingham's
  explicit formula for `ψ_1` (his Lemma 2.6, with `Σ 1/γ² = 0.0462`), the bias
  from `ψ − θ ≥ 0.98√x` (Rosser–Schoenfeld), the second-order partial summation
  (his Lemma 2.7, "cf. Ingham pp. 103–104") to pass from ψ to Π and π, Büthe's
  verification `π(x) < li(x)`, `θ(x) < x` for `x ≤ 10^19`, and an Ω-result (his
  Lemma 2.8) for the converse. His §6 lists "the Chebyshev bias for primes in
  arithmetic progressions" as further work. Differences here: the character
  twist; and the kernel term of the partial summation is dropped by sign
  (Θ_1 is one-signed) instead of being estimated.
* **Hamieh–Kadiri–Martin–Ng, "Comparative prime number theory problem list",
  arXiv:2407.03530.** Read Problem #23 (Johnston): with
  `A_1^π(x;q,a) = ∫_0^x (φ(q)π(t;q,a) − li(t)) dt`, "assuming GRH, classify those
  values of q and a such that `A_1^π(x;q,a) < 0` for all sufficiently large x";
  Martin's remark there is that the zero sum outgrows the bias as q grows. That
  problem concerns one class against `li`, not the difference of two classes.
  The link is the identity `A_1^π(x;4,3) − A_1^π(x;4,1) = 2xR(x)`. A consequence
  of Theorem 2 together with Johnston's Theorem 1.1: under RH and GRH(χ_4),
  `∫_2^x (2π(t;4,1) − li(t)) dt = ∫_2^x (π − li) − (x − 2) − xR(x) < 0` for
  *every* `x > 2` (and the same for the class 1 mod 3, with `−(x−3)_+` in place
  of `−(x−2)`). The classes 3 mod 4 and 2 mod 3 are not decided this way (the two
  biases cancel to leading order).
* **Liu–Wang–Shi, Zenodo 10.5281/zenodo.22208026 (31 Aug 2026).** Abstract only
  (the record page). It says they "classify the residue-class pairs for which
  the integrated prime counting error is eventually negative", answering HKMN.
  I could not read the paper, so I do not know whether the sign of
  `∫_0^x (π(t;4,3) − π(t;4,1)) dt` is stated there. From the abstract it concerns
  "eventually", i.e. large x, and the HKMN quantity (one class against li).
* **T. Zhao, "On the mean values of the error terms in Mertens' theorems",
  arXiv:2411.18903v2.** Read the statements of Theorems 1, 3, 4, Corollary 1 and
  Remarks 1, 3. Theorem 3: for a finite explicit set D of fundamental
  discriminants and `i ∈ {1,2}`, GRH for `L(s,χ_d)` is equivalent to
  `∫_2^X E_i(x;d) dx > 0` for all sufficiently large X, where `E_2(x;d)` is the
  error term of `Σ_{p≤x} χ_d(p)/p`; under LI the set D is complete. These are the
  Mertens weights (`p^{-1}`, `p^{-1} log p`), and "sufficiently large X". Remark 1
  is the same Landau argument as in Theorem 6, with the same proviso about real
  zeros.
* **L. C. Washington, "Sums of powers of primes II", arXiv:2209.12845.** Read the
  statement of Theorem 3 and the start of its proof: for `0 < k ≤ 10.32`, RH is
  equivalent to `∫_1^x (π_k(t) − π(t^{k+1})) dt < 0` for all sufficiently large x.
* Suzuki (arXiv:2411.07436) and Fujii's papers: not read by me; see
  `../PRIOR-ART.md`.

**What is and is not new.** Nothing in the method is new:

* the mechanism — a positive bias from prime squares beating an absolutely
  convergent sum over zeros after one integration, and Landau's theorem for the
  converse — is that of Hardy–Littlewood and Landau (1918) for `χ_4`, of Pintz
  and Johnston for `π − li`, and of Zhao and Washington for other weights;
* Lemma 2 is Ingham's second-order partial summation, used by Johnston as his
  Lemma 2.7 with the same positive kernel;
* Theorem 6 is the Landau mechanism applied to one more kernel.

What I did not find in what I read (and my reading is limited: Liu–Wang–Shi,
Suzuki, Fujii, and the originals of Hardy–Littlewood, Landau, Pintz and Ingham
were not read) is:

(a) the statement itself for the Cesàro mean of Chebyshev's race
    `π(t;4,3) − π(t;4,1)`, for **every** `x > 3` (resp. `x > 2` for q = 3), with
    explicit constants; it is the direct character analogue of Johnston's
    theorem and should be presented as such;
(b) its reading as inheritance of leads across weights (Corollary 4), at every x;
(c) the observation that the integrated θ-race `Θ_1` is one-signed, so that the
    positive kernel of Ingham's partial summation can be dropped with no loss —
    and its two uses: Lemma 6 (the same reduction for an arbitrary pair of
    weights) and the cone `m ≤ 0, m < m' ≤ 1` of Theorem 5.

---------------------------------------------------------------------------

## 5. Numerics

Program: `cesaro_scan.c` (segmented sieve adapted from `../sieve/nest_scan.c`;
`gcc -O2 cesaro_scan.c -lm`; `./cesaro_scan q 1e10`; about 65 s per character on
one core). Raw output: `out_q4_1e10.txt`, `out_q3_1e10.txt`. Certified companion:
`cert_small.py` (exact integers and Arb balls; 140 s), output `cert_small.out`.

What is exact and what is floating point:

* `D_0` (int64) and `D_1` (int128) are exact. The test `R(p) > 0` is the exact
  integer test `p·D_0(p) − D_1(p) > 0`. Between consecutive primes
  `R(x) = A − B/x` is monotone and R is continuous, so `min_{[p,p']} R = min(R(p), R(p'))`:
  positivity at all primes is positivity everywhere.
* In the C program, `N` and the gap bounds for `N` are computed from the exact
  `D_0`, `D_1` by a fixed handful of `long double` operations (one division and
  one subtraction for R; `logl`, `sqrtl`, one division, one multiplication for
  the weight): no accumulation, relative error below `10^-17`, but not interval
  arithmetic. Corollary 3 does not rely on them: the whole range
  (`29 ≤ x ≤ 10^10` for χ_4, `11 ≤ x ≤ 10^10` for χ_3) is certified in ball
  arithmetic by `cert_gaps.py`. The smallest gap bound is `0.603379` resp.
  `0.711389` in both computations (line `NGAP` of the C output), i.e. the margin
  over 1/2 is `0.103` resp. `0.211`.
* `Θ`, `W = Σ s(p) p log p`, `θ` are **accumulated** in x87 `long double` (64-bit
  mantissa). Worst-case accumulated rounding at `10^10`: `≤ 0.21` in θ and
  `≤ 6·10^4` in `Θ_1` (which is `≈ 6·10^14` there), i.e. relative `10^-10`. True
  margins of the quantities tested (lines `T1MIN`, `C1B`, `CUP` and the `G_low`
  column): `Θ_1(p) ≥ 1.1756` (χ_4, at p = 7, where `G = 0.0635`) resp. `2.0794`
  (χ_3, at p = 5, `G = 0.1860`), and `G ≥ 0.13` on `[11, 10^10]`, `G ≥ 0.52` on
  `[2^13, 10^10]`; `(y − ϑ_o(y))/√y ≤ 1.98197` resp. `1.99271` for `y ≥ 1427`
  against `c_1 = 2.0712` resp. `2.0819`; `ϑ_o(y) − y ≤ −0.26√y`. (The first
  version of this note said "margin `0.14 t^{3/2}`" for `Θ_1`; that is the gap
  bound for `t ≥ 11` only and is wrong at p = 7.) The sign of `Θ_1(p)` in the C
  program is thus a floating-point decision; it is not what the proof uses.
* **What the proof uses is certified separately** (`cert_small.py`,
  `cert_small.out`): (A) `Θ_1(p) > 0` for all race primes `p ≤ 10^4` as an exact
  inequality between integers; (B) `Θ_1(p) > 0` for `p ≤ 2^22` in 256-bit Arb
  ball arithmetic; (C) Lemma 4(a) for `y ≤ 2^22` in ball arithmetic, including
  the extremum at 1423; (D) `N > 1/2` on `[23, 29]` and on all gaps with
  `29 ≤ p < p' ≤ 2^22` (χ_4), on `[6, 11]` and all gaps with `11 ≤ p < p' ≤ 2^22`
  (χ_3), with exact rational R; (E) the same gap bounds up to `1.6·10^8` (χ_4)
  and `2.6·10^7` (χ_3) in 64-bit balls. All five parts pass for both characters.
  Part E is superseded by `cert_gaps.py`, which carries the same certified gap
  check to `10^10` with a segmented sieve (455 052 501 resp. 455 052 506 gaps, no
  failure; 1203 s and 1132 s, run as two `nice` processes in parallel;
  checkpoints every `5·10^8` in `cert_gaps_q4.out`, `cert_gaps_q3.out`).
* "low" columns are lower bounds over **all real** x in the range, obtained per
  gap as `min(R(p),R(p')) · log p'/√p'` and `min(Θ_1(p),Θ_1(p'))/p'^{3/2}`
  (valid for `p ≥ 11`), each gap being credited to the ranges of both endpoints.
* `verify_small.py` recomputes everything up to `2^21` with numpy (different
  code path), reproduces the dyadic ranges `k = 3..20` of the C output, checks
  the checkpoints, and treats `(p_0, 64]` on a grid of step `10^-4`.

Summary, q = 4, `x ≤ 10^10`:

    primes p ≥ 5 with R(p) ≤ 0:            none (exact)
    primes p ≥ 5 with Θ_1(p) ≤ 0:          none (exact for p ≤ 10^4, Arb for p ≤ 2^22, long double beyond)
    min of Θ_1 over primes p ≥ 5:          1.175573 at p = 7  (G(7) = 0.063475)
    min of R over x ≥ 5:                   2/7 at x = 7      (R → 0 as x → 3+)
    inf of N over [10^4, 10^10]:           0.682694 at x = 1 070 765 639 (a prime)
    inf of N over [11, 10^10]:             N(11) = (6/11) log 11/√11 = 0.394360
                                           (on [11,13), R = 2 − 16/x and N increases; every later gap bound is ≥ 0.4817)
    inf of N over [5, 11]:                 0.2101 at x = 7
    sup of N over [11, 10^10]:             0.982469 at x = 7477
    inf of G over [2^13, 10^10]:           0.5196 ;  sup 0.7527
    c_1 = sup (y − ϑ_o(y))/√y, y ≤ 10^10:  2.071192618 (y → 1423^-);  over y ≥ 1427: 1.981975 (y → 1427^-)
    smallest gap bound for N, gaps p ≥ 29:  0.603379 (gap starting at 47)
    sup (ϑ_o(y) − y)/√y, y ≤ 10^10:        −0.2605  (so ϑ_o(y) < y)
    sup |Θ(y)|/√y, y ≤ 10^10:              2.153717 at y = 126227
    checkpoints: D_0(10^10) = 5960, D_1(10^10) = 27 201 020 910 942, R(10^10) = 3239.89790891

(`D_0`, `D_1` and the N-minima agree with `../sieve/results.txt`, pair 0<1.)

Dyadic ranges, q = 4 (`R_low`, `N_low`, `G_low` are rigorous lower bounds over all
real x in the range; the range k = 33 stops at 10^10):

| k | range | R_low | first prime | N_low | N_min at a prime (prime) | N_max | G_low | L(2^k) |
|---|---|---|---|---|---|---|---|---|
| 3 | [2^3, 2^4) | 0.5455 | 11 | 0.3880 | 0.394360 (11) | 0.5472 | 0.1475 | — |
| 4 | [2^4, 2^5) | 0.7368 | 19 | 0.4817 | 0.497737 (19) | 0.6469 | 0.1411 | — |
| 5 | [2^5, 2^6) | 1.0323 | 47 | 0.6034 | 0.621347 (47) | 0.7059 | 0.2109 | — |
| 6 | [2^6, 2^7) | 1.1940 | 67 | 0.6040 | 0.613355 (67) | 0.7949 | 0.2192 | — |
| 7 | [2^7, 2^8) | 1.6693 | 131 | 0.7059 | 0.715332 (131) | 0.8782 | 0.3070 | — |
| 8 | [2^8, 2^9) | 2.5179 | 257 | 0.7635 | 0.766972 (487) | 0.9242 | 0.3700 | — |
| 9 | [2^9, 2^10) | 2.8566 | 521 | 0.7829 | 0.802744 (521) | 0.9099 | 0.3882 | — |
| 10 | [2^10, 2^11) | 4.0496 | 1451 | 0.7699 | 0.771797 (1471) | 0.9035 | 0.4058 | — |
| 11 | [2^11, 2^12) | 5.3173 | 2053 | 0.8328 | 0.834154 (3923) | 0.9243 | 0.4999 | — |
| 12 | [2^12, 2^13) | 6.5101 | 4099 | 0.8458 | 0.846710 (4099) | 0.9825 | 0.5151 | — |
| 13 | [2^13, 2^14) | 8.9973 | 9871 | 0.8266 | 0.827079 (10639) | 0.9463 | 0.5196 | — |
| 14 | [2^14, 2^15) | 11.9807 | 16411 | 0.8214 | 0.821868 (31883) | 0.9600 | 0.5463 | 0.1106 |
| 15 | [2^15, 2^16) | 14.3282 | 32771 | 0.8217 | 0.821914 (34603) | 0.9515 | 0.5497 | 0.1641 |
| 16 | [2^16, 2^17) | 20.1089 | 65537 | 0.8091 | 0.809134 (104831) | 0.9163 | 0.5671 | 0.2110 |
| 17 | [2^17, 2^18) | 28.1513 | 131101 | 0.8122 | 0.812270 (234511) | 0.9480 | 0.5875 | 0.2519 |
| 18 | [2^18, 2^19) | 34.0247 | 262147 | 0.8236 | 0.823690 (317371) | 0.9047 | 0.6070 | 0.2875 |
| 19 | [2^19, 2^20) | 47.1430 | 674231 | 0.7534 | 0.753464 (721291) | 0.9047 | 0.5525 | 0.3183 |
| 20 | [2^20, 2^21) | 65.8738 | 1048583 | 0.7776 | 0.777576 (2097143) | 0.8948 | 0.5957 | 0.3449 |
| 21 | [2^21, 2^22) | 76.8569 | 2132563 | 0.7581 | 0.758138 (2238323) | 0.8920 | 0.5780 | 0.3679 |
| 22 | [2^22, 2^23) | 111.3006 | 4960811 | 0.7628 | 0.762835 (5175347) | 0.8580 | 0.5947 | 0.3877 |
| 23 | [2^23, 2^24) | 149.9659 | 8388617 | 0.7330 | 0.733019 (15344207) | 0.9028 | 0.5799 | 0.4047 |
| 24 | [2^24, 2^25) | 187.9261 | 16777259 | 0.7632 | 0.763246 (16777259) | 0.8578 | 0.6119 | 0.4193 |
| 25 | [2^25, 2^26) | 276.0561 | 33554467 | 0.7258 | 0.725821 (52302223) | 0.8409 | 0.5878 | 0.4319 |
| 26 | [2^26, 2^27) | 382.2286 | 67108879 | 0.7701 | 0.770054 (134217689) | 0.8669 | 0.6397 | 0.4427 |
| 27 | [2^27, 2^28) | 471.4346 | 148988299 | 0.7219 | 0.721905 (152368063) | 0.8561 | 0.5944 | 0.4520 |
| 28 | [2^28, 2^29) | 682.2928 | 344339183 | 0.7195 | 0.719460 (349362011) | 0.8402 | 0.5985 | 0.4600 |
| 29 | [2^29, 2^30) | 946.5043 | 536870923 | 0.6827 | 0.682694 (1070765639) | 0.8368 | 0.5706 | 0.4668 |
| 30 | [2^30, 2^31) | 1076.2839 | 1073741827 | 0.6830 | 0.683005 (1073741827) | 0.8222 | 0.5710 | 0.4727 |
| 31 | [2^31, 2^32) | 1723.3893 | 2156756983 | 0.7061 | 0.706111 (3417316939) | 0.7995 | 0.6023 | 0.4778 |
| 32 | [2^32, 2^33) | 2298.9618 | 4294967311 | 0.7006 | 0.700573 (7461775091) | 0.8544 | 0.6006 | 0.4822 |
| 33 | [2^33, 10^10] | 2968.2095 | 8589934609 | 0.7323 | 0.732340 (8619119543) | 0.7460 | 0.6337 | 0.4860 |

"first prime" is the prime at which the minimum of R over the primes of the
range is attained (R is roughly increasing, so it is at or near the left end).

Summary, q = 3, `x ≤ 10^10`: no prime `p ≥ 5` with `R(p) ≤ 0` or `Θ_1(p) ≤ 0`;
`min_{x≥5} R = 0.6` at x = 5; `inf N` over `[10^4, 10^10]` is `0.720493` at
`x = 323 405 729`; over `[7, 10^10]` it is `≥ 0.7114`; `sup N = 0.973683` at 3511;
`inf G` over `[2^13, 10^10]` is `0.5411`; `min Θ_1(p) = 2.079442` at p = 5;
`c_1 = 2.081941197` (`y → 1423^-`; `1.992708` over `y ≥ 1427`); smallest gap bound
for N over gaps `p ≥ 11`: `0.711389`; `D_0(10^10) = 6264`,
`D_1(10^10) = 27 948 685 258 617`. The dyadic table is in `out_q3_1e10.txt`.

---------------------------------------------------------------------------

## 6. Table of constants and where each comes from

| constant | value (χ_4 / χ_3) | source | used in |
|---|---|---|---|
| `σ_0²` | 0.155567979923586 / 0.113229969857472 | closed form `log(q/π) − γ + 2L'/L(1,χ)`; `constants.py` (mpmath, 40 digits; `L'/L(1,χ_4) = γ + 2log2 + 3logπ − 4logΓ(1/4)`); agrees with the certified value in `RESEARCH-CONTEXT.md` §5 and with the zero lists | Lemma 3(b) |
| `b = L'/L(0,χ)` | 0.783188785 / 0.948198827 | functional equation; direct evaluation | Step 2 (only its positivity margin) |
| `b_0` | −0.060266914 / 0.143831645 | functional equation; direct evaluation | Step 2 |
| `c_1` | 2.071193 / 2.081942 | extremum at `1423^-`, certified (Arb) for `y ≤ 2^22`, rounded up; for `2^22 < y ≤ 10^10` long-double sieve with margin `0.089√y` (not needed in variant C) | Step 3 |
| `ϑ_o(y) < y`, `y ≤ 10^10` | margin `0.26√y` | Arb for `y ≤ 2^22`, long-double sieve beyond (not needed in variant C); independently Büthe Thm 2 (1.7) to `10^19` | Step 4 |
| `ε_1` | `1/log³(10^10) = 8.191·10^-5` | BKLNW (arXiv:2002.11068v2) Cor. 11.2, row `k = 3, M_3 = 1, X_0 = 89 967 803`; its RH-verification input is covered by Platt–Trudgian | Steps 3, 4 for `y > 10^10` only |
| sign of `Θ_1(p)`, `p ≤ 10^4` | min `1.1756` / `2.0794` | exact integers (`cert_small.py` A) | Step 7 |
| variant C: `Y_1^C`, `ε_1^C`, `c_θ^C` | `2^22`, `8.601·10^-4`, `1.000860` | Arb to `2^22`; BKLNW Cor. 11.2, row `k = 2, M_2 = 0.2, X_0 = 3 594 641` | Remark after Lemma 4 |
| gap bounds for N | min `0.603379` / `0.711389` | Arb to `10^10` for both characters (`cert_gaps.py`) | Corollary 3 |
| margins of the analytic part of Cor. 3 (`x ≥ 10^10`) | `0.0293` grid, `0.0110` tail / `0.0731`, `0.0534` | `bounds.py`, double precision | Corollary 3 |
| `c_θ` | `1 + ε_1` | from the two lines above | Step 4 |
| 3/4, 5/6, 4/5, 8/3 | — | elementary integrals | Steps 3, 4, Thm 2 |
| `t_0` | `10^4` (thresholds 4892 / 8526) | `bounds.py` | Step 6 |

With the classical fallback (b') everything stays true with
`ε_low = 0.0434`, `c_θ = 1.0217`: thresholds 5229 / 9148, `L(10^4) = 0.0622 / 0.0103`,
limit of L equal to `0.4821 / 0.5245`. Theorems 1, 2, 5, 6, 7 and Corollary 4 are
unaffected. Corollary 3 for χ_4 beyond `x = 10^20` is the only statement that
needs (b) (for χ_3 the fallback limit 0.5245 still exceeds 1/2). The fallback
keeps `Y_1 = 10^10`, i.e. it still uses the long-double part of Lemma 4(a); the
variant with certified inputs only is variant C (Remark after Lemma 4).

---------------------------------------------------------------------------

## 7. Hypotheses actually used

* GRH for the single L-function `L(s,χ_q)`; used only in Step 2 (`|t^ρ| = t^{1/2}`)
  and Lemma 3(b).
* No RH for ζ. The θ-bounds are unconditional: computation below `10^10`; above,
  BKLNW Corollary 11.2 (which uses the partial verification of RH up to height
  `2.446·10^12`, rigorous by Platt–Trudgian's verification up to `3·10^12`, and
  explicit zero-free regions), or Rosser–Schoenfeld as fallback.
* No LI, no simplicity of zeros (multiplicities are allowed in Lemma 3).

---------------------------------------------------------------------------

## 8. General pairs of weights (m, m')

For real `m < m'`, `κ = m' − m`, `M = m + 1/2`, `M' = m' + 1/2`, put

    I_{m,m'}(x) = D_m(x) − x^{-κ} D_{m'}(x) = Σ_{p≤x} s(p) p^m (1 − (p/x)^κ),
    Φ_{m,m'}(t) = Σ_{p≤t} s(p) log p · p^m (1 − (p/t)^κ)        (the log-weighted Riesz mean).

**Lemma 6 (no hypothesis).** For `x > p_0`,

    I_{m,m'}(x) = Φ_{m,m'}(x)/log x + ∫_{p_0}^x Φ_{m,m'}(t) w_x(t) dt,
    w_x(t) = (1 + (t/x)^κ)/(t log² t) + 2(1 − (t/x)^κ)/(κ t log³ t) > 0.

In particular, if `Φ_{m,m'} ≥ 0` on `[p_0, x]` then `I_{m,m'}(x) ≥ Φ_{m,m'}(x)/log x`.

*Proof.* With `A(t) = Σ_{p≤t} s(p) p^{m'} log p`: `I(x) = ∫ (t^{-κ} − x^{-κ})/log t dA(t)`
and `dΦ(t) = κ t^{-κ-1} A(t) dt`, `Φ(p_0) = 0`. Two integrations by parts give
`I(x) = ∫ g dΦ = g(x)Φ(x) − ∫ Φ g'` with `g(t) = 1/log t + (1 − (t/x)^κ)/(κ log² t)`,
`g(x) = 1/log x`, `−g' = w_x`. ∎ (Checked to 23 digits for (0,1), (1,2), (0.5,3),
(−1,0.5) at x = 700.3.) For (m,m') = (0,1) this is (2.1), since `Φ_{0,1}(t) = Θ_1(t)/t`.

So for **any** pair the problem reduces, without loss, to the sign of
`Φ_{m,m'}`, and `Φ_{m,m'}` has an absolutely convergent explicit formula. For
`M > 0`, under GRH,

    Φ_{m,m'}(t) = κ t^M [ 1/(2MM') + Σ_ρ t^{iγ}/((ρ+m)(ρ+m')) ] + (lower order),
    Φ_{m,m'}(t) ≥ (κ/2) t^M [ Λ(M,M') − E_{m,m'}(t) ],     Λ(M,M') = 1/(MM') − Σ_{γ>0} 4/(|M+iγ||M'+iγ|),

where the main term `κ/(2MM')` comes from the prime squares
(`∫_0^{√t} y^{2m}(1 − (y²/t)^κ) dy = t^M (1/(2M) − 1/(2M'))`) and the zero sum is
the residue sum of `(−L'/L)(s) t^{s+m} κ/((s+m)(s+m'))`. This is the constant of
Theorem A in `../CORE-NOTES.md`. To make a given pair explicit one needs:

1. **`Λ(M,M') > 0` with a certified value.** Either the closed-form sufficient
   condition `σ_M σ_{M'} < 2` (Cauchy–Schwarz and `prop:varid`), or a certified
   zero list plus a rigorous tail bound for `Σ_{γ>T} 4/γ²`. The uncertified
   values in `general_pairs.out` agree with `../sieve/SUMMARY.md`.
2. **The lower-order terms `E_{m,m'}(t)` explicitly:**
   (i) prime squares: `∫ w dϑ_o` with `w(y) = y^{2m}(1 − (y²/t)^κ)`, which for
   `m > 0` is not monotone, so the two-sided bound `|ϑ_o(y) − y| ≤ c_1√y`
   (`y ≤ 10^10`, sieve; note the upper side is `ϑ_o(y) < y`) and Lemma 4(b) beyond
   are needed; relative size `≍ t^{-1/4}`;
   (ii) odd prime powers `k ≥ 3`: trivially `≪ t^{m+1/3}`, relative `t^{-1/6}`;
   (iii) the residues at `s = −m`, `s = −m'` and at the trivial zeros:
   `(L'/L)(−m,χ) − t^{-κ}(L'/L)(−m',χ) + O(t^{m−1} log t)`, with double poles when
   m or m' is an odd positive integer. Relative size `t^{-M}`: for small `M > 0`
   this is the dominant error and it decays slowly, so `t_0(m,m')` grows like
   `(C/Λ)^{1/M}`.
   None of these loses a logarithm; all are powers of t.
3. **A finite computation of the sign of `Φ_{m,m'}` on `[p_0, t_0(m,m')]`.** If it
   is positive throughout, the conclusion holds for every `x > p_0`. If not (for
   instance `Φ_{0,3}(7) < 0` for q = 4), Lemma 6 still gives an explicit `x_1`:
   the negative part contributes at most `∫ Φ_- w_x ≤ ∫ Φ_-(t) (2/(t log² t) + 2/(κ t log³ t)) dt`,
   a constant independent of x, to be beaten by `Φ(x)/log x`.

For pairs in the cone `m ≤ 0, m < m' ≤ 1` of Theorem 5 none of this is needed: by Lemma 2 (applied
with `F = f_t` for each t) `Φ_{m,m'}(t) ≥ κ t^{m−1} Θ_1(t) ≥ 0`.

Numerics for some pairs up to `2^23` (`general_pairs.out`; Λ from zero lists plus
a smooth tail, not certified):

    q = 4   (m,m')        Λ        primes with Φ ≤ 0     min Φ/(κ t^M) on [10^4, 2^23]   primes with I ≤ 0
            (0,1)         1.0257   none                  0.5213                          none
            (0,1.5)       0.6953   none                  0.3758                          none
            (0,2)         0.4988   none                  0.2867                          none
            (0,3)         0.2786   1 (p = 7)             0.1847                          none
            (0.5,1.5)     0.1966   none                  0.1526                          none
            (1,2)        −0.0312   1 (p = 7)             0.0415                          none
            (2,3)        −0.1697   37447 (last 5234683)  −0.0269                         32950
    q = 3   (0,2)         0.5774   none                  0.3082                          none
            (1,2)         0.0454   none                  0.0739                          none
            (2,3)        −0.1011   9 (last 257)          0.0003                          9 (last 257)

So (0, 1.5), (0, 2), (0.5, 1.5) for q = 4 and (0, 2), (1, 2) for q = 3 are
candidates for an "every x > p_0" theorem by the same method; (0,3) mod 4 needs
the `x_1` variant; and the pairs with `Λ < 0` (e.g. (1,2) mod 4, where no
violation exists below `10^10`) are outside the method — under LI they must
fail for some large x.

---------------------------------------------------------------------------

## 9. What I am not fully sure of

1. **Lemma 3(a), the explicit formula for `ψ_1(x,χ)`.** The contour estimate is
   now written out (§2). The two bounds (i), (ii) for `L'/L` that it uses are
   the standard ones of Davenport Ch. 19 (Ch. 17 for ζ); I quote the chapter
   from memory and have not checked equation numbers. I still have no numbered
   citation for the formula itself in the L-function case (Ingham Thm 27 is the
   ζ case). The formula is confirmed three ways numerically (§2). Only the size
   `O(x)` of the non-zero terms and the sign of `b t − log t − 1 − b_0` matter
   for the theorems.
2. **Citations checked only partially.**
   - BKLNW Corollary 11.2 and Appendix A.1.1, Büthe Theorem 2, Johnston, Zhao,
     Washington, Platt–Trudgian (introduction): read in the arXiv versions; the
     numbering in the published versions was not checked. The arXiv page of
     BKLNW shows no journal reference; "Math. Comp. 90 (2021)" is from the
     referee report and memory.
   - Platt–Trudgian, "The Riemann hypothesis is true up to 3·10^12": not read;
     the height `3 000 175 332 800` is as quoted by BKLNW.
   - Büthe's (1.6) fails as printed (arXiv text) on `(1426.48, 1427)` by my
     computation; it is not used. Worth a look at the published version.
   - Rosser–Schoenfeld (3.15)–(3.16): quoted from Johnston's Lemma 2.2, not from
     the original; the constant 2.05282 is from memory.
   - Hardy–Littlewood, Landau, Pintz, Ingham: originals not read (§4 says through
     whom each is known).
   - Davenport Ch. 12 (17)–(18), Ch. 19 (3), and Montgomery–Vaughan Lemma 15.1
     (Landau): quoted from memory; the numbers need a check against the books.
3. **Floating point.** In variant C (Remark after Lemma 4) every finite
   computation that the proofs of Theorems 1, 2, 5, 6, 7 and Corollaries 3, 4 use
   is exact or ball-certified: `Θ_1(p) > 0` for `p ≤ 10^4` (exact integers),
   Lemma 4(a) for `y ≤ 2^22` (Arb), the gap bounds for N up to `10^10` (Arb). The main text (`Y_1 = 10^10`) uses in addition Lemma 4(a) on
   `(2^22, 10^10]` from the `long double` sieve, with margins `0.089√y` and
   `0.26√y` against a worst-case accumulated error of 0.21; that range is not
   interval-certified. The evaluations of the closed-form function L in
   `bounds.py` are ordinary double precision; their margins are: `L(10^4) = 0.068`
   (χ_4) and `0.018` (χ_3) for Theorem 1, and for Corollary 3 see the next item.
4. **Corollary 3** additionally uses the grid argument of `bounds.py` in double
   precision, now only for `x ≥ 10^10`: margins `≥ 0.029` on `[10^10, 10^30]` and
   `≥ 0.011` beyond for χ_4 (`0.073` and `0.053` for χ_3); and Lemma 4(b) for
   `x > 10^20` (χ_4 only). The earlier statement "margins ≥ 0.03" was wrong near
   the old junction at `1.5·10^8` (true margin `1.6·10^-4`); see the Correction in
   the proof of Corollary 3 and §11.
5. **Theorem 6, direction (c)⟹(a)** is now written in four steps; it relies on
   Landau's theorem as cited (number from memory).
6. **Novelty.** I could not read Liu–Wang–Shi beyond its abstract, nor Suzuki,
   Fujii, or the 1918 originals. The claim "not found" in §4 is limited
   accordingly.
7. `Λ` values in §8 are not certified and the programme of §8 is a plan, not a
   proof, for pairs outside the cone.

---------------------------------------------------------------------------

## 10. Files

    cesaro_scan.c          sieve scan (C); cesaro_scan is the gcc -O2 binary
    out_q4_1e10.txt        raw output, q = 4, X = 10^10        (time_q4.txt: 68.8 s)
    out_q3_1e10.txt        raw output, q = 3, X = 10^10        (time_q3.txt: 69.1 s)
    cert_small.py / .out   exact-integer and Arb certification of the finite checks the proofs use (time_cert.txt: 140 s)
    cert_gaps.py           Arb-certified gap check N > 1/2 up to 10^10, segmented sieve
    cert_gaps_q4.out, cert_gaps_q3.out   its output with checkpoints (time_cert_gaps_q4.txt: 1203 s, time_cert_gaps_q3.txt: 1132 s)
    constants.py / .out    closed-form constants; checks of Lemma 3 against the zero lists
    bounds.py / .out       L(t), thresholds, Theorem 2 bound, comparisons with the sieve, Corollary 3 grid
    verify_small.py        independent recomputation to 2^21; identities; cone check
    verify_small_q4.out, verify_small_q3.out
    general_pairs.py / .out   §8 numerics
    REFEREE-REPORT.md      the referee's report (not written by me)

Reproduce: `gcc -O2 -o cesaro_scan cesaro_scan.c -lm && ./cesaro_scan 4 1e10 > out_q4_1e10.txt && ./cesaro_scan 3 1e10 > out_q3_1e10.txt && python3 cert_small.py && python3 cert_gaps.py 4 1e10 > cert_gaps_q4.out && python3 cert_gaps.py 3 1e10 > cert_gaps_q3.out && python3 constants.py && python3 bounds.py && python3 verify_small.py 4 && python3 verify_small.py 3 && python3 general_pairs.py`

---------------------------------------------------------------------------

## 11. Review history

Referee report: `REFEREE-REPORT.md` (2026-11-24). Verdict: no fatal error; four
repairable items and seven minor ones. Resolution of each:

| item | referee's point | resolution |
|---|---|---|
| REP-1 | Lemma 4(b) cited Dusart arXiv:1002.0442 Thm 5.2, which rests on Gourdon's unpublished zero verification | Replaced by Broadbent–Kadiri–Lumley–Ng–Wilk, arXiv:2002.11068v2, Corollary 11.2 (row `k = 3, M_3 = 1, X_0 = 89 967 803`), read online; it is the same inequality, so no constant changes. BKLNW's own input is "RH up to `2 445 999 556 030`" (Gourdon–Wedeniwski), which their Appendix A.1.1 notes is confirmed by Platt–Trudgian's verification up to `3 000 175 332 800`; this is stated in Lemma 4(b). Not done: I did not read Platt–Trudgian's paper itself, and the journal pages of BKLNW are unchecked (§9.2). |
| REP-2 | contour estimate in Lemma 3(a) only cited | Written out: choice of `T_k`, the bounds (i), (ii), the horizontal and left integrals, the limit, the case `x = 1`. Davenport equation numbers remain unchecked (§9.1). |
| REP-3 | sign of `Θ_1(p)` decided in long double; wrong margin "0.14 t^{3/2}"; N values for Cor. 3 uncertified | New script `cert_small.py`: exact integer test of `Θ_1(p) > 0` for all race primes `p ≤ 10^4` (both characters), Arb certification to `2^22`, Arb certification of Lemma 4(a) to `2^22`, and Arb certification of the `N > 1/2` gap bounds on the hand intervals and on every gap up to `1.6·10^8` (χ_4) / `2.6·10^7` (χ_3); this last range was later extended to `10^10` (see the follow-up below). Margin statement corrected (min `Θ_1 = 1.1756` at 7, `G(7) = 0.0635`; `2.0794` at 5 for χ_3). §5 states plainly that the N values of the C program are `long double` (no accumulation, margin `≥ 0.1`: 0.603379 vs 1/2); Corollary 3 now uses the certified ones. In addition a variant C with certified inputs only was added (Remark after Lemma 4: `Y_1 = 2^22` and BKLNW's row `k = 2`), in which no long-double result is used at all. `cesaro_scan.c` prints the true margins (`T1MIN`, `NGAP`, `C1B`); both scans were rerun. |
| REP-4 | novelty of Lemma 2; priority for the Landau mechanism | Lemma 2 is now attributed to Ingham (pp. 103–104) and Johnston's Lemma 2.7 (status paragraph, "Attribution" after Lemma 2, §4). Hardy–Littlewood and Landau (1918) and Pintz (1991) are cited in §4 and in the proof of Theorem 6; bibliographic data verified online through Platt–Trudgian's and Johnston's reference lists, originals not read (said so). The novelty paragraph is rewritten: what is claimed is only the all-x statement, its reading as lead inheritance, and the one-signedness observation with its uses (Lemma 6, Theorem 5). |
| MINOR-1 | cone written as "m ≤ 0 < m' ≤ 1" in places | "`m ≤ 0, m < m' ≤ 1`" everywhere (status paragraph, Theorem 5 heading, §4, §8). |
| MINOR-2 | wrong range in the quotation of Johnston's Lemma 2.2; RS 2.05282 without equation number; cited bounds do not give `c_1` for small y | (b') now quotes `x ≥ 41` for both inequalities and says the original was not consulted; 2.05282 is marked "from memory"; Lemma 4(a) says explicitly that below 1427 it rests on the computation alone (now Arb-certified to `2^22`). While checking this I found that Büthe's (1.6) fails as printed on `(1426.48, 1427)`; it was only a consistency remark and is now explicitly not used. |
| MINOR-3 | inf of N over `[11, 10^10]` | Stated exactly: `N(11) = (6/11) log 11/√11 = 0.394360` (§1, §5). |
| MINOR-4 | k = 10 line printed twice in `out_q4_1e10.txt` | Could not reproduce: the file (old and regenerated) contains each of the 31 range lines once (`awk '$1 ~ /^[0-9]+$/' out_q4_1e10.txt \| sort \| uniq -d` is empty). No change. |
| MINOR-5 | checks of the explicit formula at `10^3`, `10^4` are vacuous; χ_3 value omitted | `constants.py` now tests `x = 1, 3, 5.5, 10, 20, 50, 200`; table for both characters in §2; the vacuous comparisons were removed. |
| MINOR-6 | "linearly independent" should include simplicity; say why liminf N = liminf G | Corollary 3 parenthetical rewritten. |
| MINOR-7 | Corollary 3 for χ_3 holds below 7 | Now stated for `x ≥ 6` with proof on `[6, 7]` (certified in `cert_small.py`); crossing at `x ≈ 5.343`. |
| §5 of the report | write out the Landau steps and the role of real zeros | Proof of Theorem 6 rewritten in four steps, with the dependence on the absence of real zeros stated. |

**Follow-up (coordinator, after the transcription into the manuscript).**

| point raised | resolution |
|---|---|
| The statement "margins ≥ 0.03" for the double-precision grid evaluation behind Corollary 3 holds only from `10^10` on; near the junction with the certified range the margin is about `10^-4` (χ_4: `1.5·10^-4` at `1.514·10^8`, `7.9·10^-4` at `1.585·10^8`; χ_3: `3.9·10^-4` at `2.512·10^7`). | Confirmed (`bounds.out`, "(record)" lines: `1.6·10^-4`, `7.9·10^-4`; `3.9·10^-4`). The statement is corrected in the proof of Corollary 3 ("Correction"), in §6 and in §9.4. While redoing this I found a second weak point of the same kind: beyond the old grid (`a > 5.5·10^11`) the proof switched to "L alone ≥ 1/2", whose margin starts at 0 there. The analytic part is now: grid of β from `10^10` to `10^30` (margin `≥ 0.0293`, χ_4) and `inf L` beyond `10^30` (margin `0.0110`, which is `c_0 − (2/3)ε_1 − 1/2` and cannot be improved for the constant 1/2). So the honest statement is "≥ 0.029 up to `10^30`, ≥ 0.011 beyond", not "≥ 0.03". |
| Remove the issue at its root: extend the Arb-certified gap check to `10^10` for both characters. | Done: `cert_gaps.py` (segmented sieve, bounded memory, checkpoints every `5·10^8`), run as two processes under `nohup nice -n 10`. χ_4: 455 052 501 gaps with `29 ≤ p < p' ≤ 10^10` plus the last piece up to `10^10`, no failure, smallest gap bound `0.603378720997762`, i.e. smallest certified margin `0.1033787` (gap starting at 47), 1203 s. χ_3: 455 052 506 gaps with `11 ≤ p`, no failure, smallest gap bound `0.711388955945845`, margin `0.2113889` (gap starting at 11), 1132 s. The limit `10^10` was reached (no need to stop earlier). Corollary 3 now rests on the certified computation on `[23, 10^10]` (χ_4) / `[6, 10^10]` (χ_3) and on the analytic bound beyond. |
| "Arb-certified to `2^21`" in the review table versus `2^22` in the lemma and in `cert_small.out`. | Fixed (MINOR-2 row): `2^22`. The remaining occurrences of `2^21` in this note refer to `verify_small.py`, which does stop at `2^21`. |

Not fixed, and why: (1) interval certification of the `long double` sieve itself
on `(2^22, 10^10]` (would need directed rounding in the C code or a much slower
ball-arithmetic sieve); instead the margins are stated and variant C avoids
that range; (2) verification of equation numbers in Davenport,
Montgomery–Vaughan, Rosser–Schoenfeld, of the published versions of the arXiv
papers, and of the 1918 originals and Pintz (no access here); (3) MINOR-4 could
not be reproduced; (4) the closed-form evaluations in `bounds.py` remain double
precision (margins as stated in §9.3–9.4).
