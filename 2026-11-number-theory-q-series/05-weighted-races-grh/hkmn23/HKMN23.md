# HKMN Problem #23: the sign of the integrated prime-counting error, for every x

Working notes, 2026-11-24. Everything in this directory (`monotonicity/hkmn23/`) was produced
for this note. Reviewed twice: `REFEREE-HKMN23.md` (a pass commissioned by the authoring session,
own scripts) and `AUDIT-INDEPENDENT.md` (an outside session, independent code in `audit2/`);
the corrections of both are incorporated (§8, review history). Not refereed externally. Status words used below: **proved** (a written proof, conditional
only on the stated hypotheses and on the certified computations named in it), **certified**
(decided in ball arithmetic, Arb via python-flint 0.9.0, or in exact integer arithmetic),
**numerical** (floating point, no rigor claimed), **cited** (taken from a source I read).

## Summary (15 lines)

1. **Problem.** HKMN #23 (Johnston): for `A_1^π(x;q,a) = ∫ (φ(q)π(t;q,a) − li(t)) dt`, "assuming GRH, classify the (q,a) with `A_1^π < 0` for all sufficiently large x". Liu–Wang–Shi (2026, Zenodo, unrefereed): 18 pairs under GRH ("eventually"); no others under GRH + LI + central non-vanishing.
2. **Theorem 1 (work order T4, complete).** For each of the 18 pairs, under GRH(q) only (RH and GRH for the characters mod q): `∫_2^x (φπ(t;q,a) − li(t)) dt ≤ −li(2)(x−2) + W_1(x)/log x < 0` for **every** `x > 2`.
3. **Lower limit.** The problem list prints the integral from 0 (Johnston: from 2). From 0 it is negative exactly for `x > 2.6111394457` (the 17 pairs with q ≥ 2) and `x > 3.5651763736` (q = 1), positive before.
4. **Method.** One exact identity (Lemma A) reduces everything to the sign of `W_1(t) = φ Σ_{p≤t, p≡a} log p·(t−p) − (t−2)²/2`, with a positive kernel: no loss of `1/log x`, no estimate of the kernel term.
5. **Inputs.** Explicit formula for `ψ_1(x,χ)` (odd, even, imprimitive); `L_1(q)` in closed form *without zeros* (Lemma S; every value printed by LWS is reproduced); prime squares from one certified constant and, beyond `10^12`, from GRH(q) itself plus Chebyshev.
6. **Margins.** `(2/3)N − L_1(q)` ranges from 0.385 (q = 20) to 2.94 (q = 24); the analytic range starts at `X_0` between 58.5 and `8.17·10^5`; `W_1 < 0` is certified in ball arithmetic for all `t ≤ 10^9` (overlap ≥ 1200-fold). Table 1.
7. **By-product.** A short proof of "RH ⟹ `∫_2^x (π − li) < 0` for all x > 2" (half of Johnston's Theorem 1.1) that does not use `π(x) < li(x)` for `x ≤ 10^19`.
8. **Witness principle (Theorem 3, T2).** Under GRH(q) alone: if `F_{q,a}(y_0) > (2/3)N(q,a)` at one real `y_0`, then `A_1^π(x;q,a) > 0` for arbitrarily large x. Certified values of `F` come from primes (exact formula (F)) or from zeros with a closed-form remainder (no completeness needed).
9. **Non-square classes (Theorem 4).** Under GRH(q) and `L(1/2,χ) ≠ 0` (certified for every q ≤ 400; an extra hypothesis for q > 400), every non-square class is excluded: `F_{q,a}` is never constant. This removes LI from all 49 classes LWS leave to it, and their side condition.
10. **The class a = 1 (Theorem 5).** `F_{q,1}(0) = Φ(q)` is a closed form in `(L'/L)(1,χ)`, `(L'/L)(2,χ)`; `Φ(q) > (2/3)N(q,1)` for every q outside the 13 member moduli and {7, 9, 14, 16, 18} (Arb for q ≤ 400, explicit bound for q > 400). So `(q,1)` is classified under GRH alone except at five moduli.
11. **Candidate moduli (Theorem 6a).** Of the 160 classes at the 23 moduli q ≤ 60 that matter under LI: 18 members, 128 excluded under GRH alone (LWS: 67), 14 undecided — exactly the square classes at the knife-edge moduli q = 7, 9, 14, 16, 18. The square classes ≠ 1 at q = 21, 40, 42, 48, 60 are excluded by certified witnesses at `|y_0| ≈ 10^{115}–10^{122}` found by lattice reduction (LLL on the heaviest ordinates).
12. **All q ≤ 60 (Theorem 6b).** <!--S6-START-->1070 of the 1084 pairs with `3 ≤ q ≤ 60` outside the 18 members are excluded under GRH(q) alone; 14 square classes are undecided (the 14 of line 11, and no others). For `61 ≤ q ≤ 120`: all pairs excluded except 18 square classes `a ≠ 1` whose zero-side search was not run or not finished.<!--S6-END-->
13. **Rigor.** Proofs in §§2–4; finite parts decided in Arb ball arithmetic or on integers; floating point only proposes witnesses. Two reviews (`REFEREE-HKMN23.md`, `AUDIT-INDEPENDENT.md`): no fatal error; the audit reproduced Theorem 1 end to end for five pairs with independent code; repairs applied (§8).
14. **Not done.** The 14 knife-edge classes (they need 97–98% of the bound `L_1` for `sup F`; the lattice method reaches about 80%, and thousands of ordinates would have to be aligned, §5); the converse (equivalence with GRH(q)); explicit O-constants in Theorem 3.
15. **Caveats (§6).** Nothing is externally refereed; the certificates trust python-flint/Arb; Theorem 5(b) for q > 400 is a new argument checked only internally; LWS is an unrefereed preprint (not used in any proof).

---------------------------------------------------------------------------

## 0. The problem, the sources, the conventions

**Sources read (pdftotext of the originals).**

* Hamieh–Kadiri–Martin–Ng, *Comparative prime number theory problem list*, arXiv:2407.03530v1.
  Problem #23 (Johnston), verbatim: "It has been shown [21, Theorem 1.1] that RH is equivalent
  to `A_1^π(x) = ∫_0^x (π(t) − li(t)) dt < 0` for all `x > 2`. […] Martin pointed out […] that
  when looking at the analogous `A_1^π(x; q, a) = ∫_0^x (φ(q)π(t; q, a) − li(t)) dt`, the
  corresponding sum over zeros of Dirichlet L-functions grows faster than the bias as q
  increases […]. **Problem: Assuming GRH, classify those values of q and a such that
  `A_1^π(x; q, a) < 0` for all sufficiently large x.**" The integrals are printed with lower
  limit **0**.
* D. R. Johnston, *On the average value of π(t) − li(t)*, arXiv:2201.06184v2 (Canad. Math.
  Bull. 66 (2023)). Theorem 1.1: RH ⟺ `∫_2^x (π(t) − li(t)) dt < 0` for all `x > 2`. Lower
  limit **2**. With lower limit 0 the statement quoted in the problem list is false for small
  x: `∫_0^2 (π − li) = −∫_0^2 li = +0.877257…` (§3 below). This does not affect the problem,
  which is about large x; it matters for "every x" statements, so both limits are treated.
* Liu–Wang–Shi, *Two problems in comparative prime number theory…*, Zenodo
  10.5281/zenodo.22537927 (v1.1, 6 Sep 2026; md5 of the PDF `943c3028…ba3d`), §§1, 7–11 read.
  Not refereed. Their Theorem 1.2: under GRH, `A_1^π(x;q,a) < 0` for all sufficiently large x
  at the 18 pairs
  `(1,1), (2,1), (3,1), (4,1), (5,1), (5,4), (6,1), (8,1), (10,1), (10,9), (12,1), (15,1),
  (15,4), (20,1), (20,9), (24,1), (30,1), (30,19)`; under GRH(q), LI(q) and central
  non-vanishing for every q there are no other pairs. Their Theorem 8.1 (GRH only) is the
  criterion `(2/3)N(q,a) > L_1(q)`; their Theorem 8.2 (the converse) uses LI; their
  Proposition 8.3 removes LI for 67 of the 116 non-square classes at their 23 candidate moduli;
  §8.5: at a square class "GRH alone supplies no number here, and the undecided band at a
  square class is the whole range below L_1(q)". (The numbers 67 and 116 are from their §8.4;
  "49" = 116 − 67. Values of theirs that I compare with: (49), (57)–(59), and
  `L_1(11) = 3.7735558` in §11.1.)

**Notation.** `q ≥ 1`, `a` a reduced residue mod q, `φ = φ(q)`,

    N(q,a)  = #{b mod q : b² ≡ a}                       (0 if a is not a square mod q)
    L_1(q)  = Σ_{χ mod q} Σ_ρ 1/|ρ(ρ+1)|                (ρ over the nontrivial zeros of the primitive
                                                         character χ* inducing χ; χ_0 contributes the zeros of ζ)
    A(x)    = A(x;q,a)   = ∫_2^x (φ π(t;q,a) − li(t)) dt        (Johnston's lower limit)
    A_0(x)  = A_0(x;q,a) = ∫_0^x (φ π(t;q,a) − li(t)) dt        (HKMN's lower limit;  A_0 = A + c_0 for x ≥ 2)
    θ_1(t;q,a) = Σ_{p≤t, p≡a (q)} log p · (t − p),      ψ_1(t;q,a) = Σ_{n≤t, n≡a (q)} Λ(n)(t − n)
    W_1(t)  = W_1(t;q,a) = φ θ_1(t;q,a) − (t − 2)²/2             (t ≥ 2)
    F_{q,a}(y) = − Σ_{χ mod q} χ̄(a) Σ_{ρ of χ*} e^{iγy} / (ρ(ρ+1))      (ρ = 1/2 + iγ under GRH; LWS (48))

`li(x) = PV∫_0^x dt/log t`, `li(2) = 1.04516378011749278…`.
**GRH(q)**: RH, and every zero in `0 < Re s < 1` of `L(s,χ*)` has real part 1/2, for every
primitive character χ* of conductor `f | q`, `f > 1` (equivalently: for every χ mod q, the
zeros of `L(s,χ)` in the open strip are on the line; the zeros of the removed Euler factors
on `Re s = 0` do not count). It *includes RH* (χ_0), as in LWS §8. **NV(q)**: `L(1/2,χ) ≠ 0` for every χ mod q. No
linear-independence hypothesis is used anywhere in this note.

Two trivial reductions. (i) If m is odd and a is odd, then `π(t;2m,a) = π(t;m,a)` unless
`2 ≡ a (mod m)` (which among the 18 pairs happens only for m = 1), `φ(2m) = φ(m)`,
`GRH(2m) = GRH(m)`, `F_{2m,a} = F_{m,a}`: (this holds for every modulus `2m`, m odd, and every odd a with `2 ≢ a (mod m)`, and is
used again in Theorem 6); so `(6,1) = (3,1)`, `(10,1) = (5,1)`,
`(10,9) = (5,4)`, `(30,1) = (15,1)`, `(30,19) = (15,4)` are literally the same functions.
(ii) `π(t;2,1) = π(t) − 1` for `t ≥ 2`. The 18 pairs are therefore 12 distinct cases with
`q ≥ 2`, plus `(1,1)`, which follows from `(2,1)` (Corollary 1.3).

---------------------------------------------------------------------------

## 1. Results

### 1.1 "For every x" at the 18 pairs (work order T4)

**Theorem 1 (proved; finite part certified).** Let (q,a) be one of the 18 pairs of
Liu–Wang–Shi with `q ≥ 2`, and assume GRH(q). Then

    W_1(t;q,a) < 0   for every t > 2,

and consequently, for every `x > 2`,

    ∫_2^x (φ(q) π(t;q,a) − li(t)) dt  ≤  − li(2)·(x − 2) + W_1(x)/log x  <  0.

For `(1,1)`, under RH: `∫_2^x (π(t) − li(t)) dt ≤ −(li(2) − 1)(x − 2) + W_1(x;2,1)/log x < 0`
for every `x > 2` (`li(2) − 1 = 0.04516…`).

So for all 18 pairs the integrated error is negative for **every** `x > 2`, not only for all
sufficiently large x; `x_1(q,a) = 2` with Johnston's lower limit.

**Theorem 2 (the explicit bound; proved).** Same hypotheses, `q ≥ 2`. For every `t ≥ 8`,

    W_1(t) ≤ t^{3/2} · U(t),      U(t) = −m + K(t) t^{-1/4} + r_+(t),      m = (2/3)N(q,a) − L_1(q) > 0,

where `K(t) = (4/5)c_1` for `t ≤ 10^12` and `K(t) = K_B(t)` for `t > 10^12`, and `c_1`,
`K_B`, `r_+` are explicit (§2.4–2.5; `K_B` is decreasing, `r_+(t) ≪ t^{-1/2} log t` is
decreasing). `U` is decreasing on `[8, 10^12]` and on `(10^12, ∞)`, `U(X_0) < 0`,
`U(10^12 +) < 0`, `U(t) → −m`. The thresholds `X_0` and the constants are in Table 1.
In particular `A(x) ≤ −li(2)(x−2) − |U(x)| x^{3/2}/log x` for `x ≥ X_0`, and
`limsup A(x) log x / x^{3/2} ≤ −m`.

**Table 1** (`l1const.out`, `w1_analytic.out`, `raw/w1_cert_*.json`). `m` = margin;
`X_0` = start of the analytic range; "certified" = `W_1 < 0` certified in ball arithmetic for
all `2 < t ≤` that bound; `sup G` = certified upper bound for `sup W_1(t)/t^{3/2}` over the pieces of `W_1`
that start at primes in `[10^8, 10^9)` (the maximum of `M_i/p_{i+1}^{3/2}`; in this decade it
is close to the true supremum; compare `−(2/3)N ± L_1`).

    (q,a)            φ  N  L_1(q)        m = (2/3)N − L_1   c_1      X_0       certified to   overlap    sup G [1e8,1e9]   x_1 (limit 2)   x_1 (limit 0)
    (1,1)            1  1  0.0461176444  0.6205490          —        (from (2,1))             —          —                 2               3.5651763736
    (2,1)            1  1  0.0461176444  0.6205490          2.0712   58.5      1 000 000 007  1.7·10^7   −0.673            2               2.6111394457
    (3,1) = (6,1)    2  2  0.1586735289  1.1746598          4.1639   66.9      1 000 000 009  1.5·10^7   −1.287            2               2.6111394457
    (4,1)            2  2  0.1999123498  1.1334210          4.1424   75.2      1 000 000 009  1.3·10^7   −1.293            2               2.6111394457
    (5,1) = (10,1)   4  2  0.5992348066  0.7340985          7.1660   3730      1 000 000 021  2.7·10^5   −1.215            2               2.6111394457
    (5,4) = (10,9)   4  2  0.5992348066  0.7340985          4.0960   689       1 000 000 009  1.5·10^6   −1.152            2               2.6111394457
    (8,1)            4  4  0.7356310182  1.9310356          8.2848   149       1 000 000 009  6.7·10^6   −2.463            2               2.6111394457
    (12,1)           4  4  0.6326579695  2.0340087          8.9786   158       1 000 000 009  6.3·10^6   −2.408            2               2.6111394457
    (15,1) = (30,1)  8  4  1.9257330053  0.7409337          14.3321  6.05·10^4 1 000 000 021  1.7·10^4   −2.203            2               2.6111394457
    (15,4) = (30,19) 8  4  1.9257330053  0.7409337          9.6368   1.78·10^4 1 000 000 009  5.6·10^4   −2.372            2               2.6111394457
    (20,1)           8  4  2.2817057508  0.3849609          14.3321  8.17·10^5 1 000 000 021  1.2·10^3   −2.353            2               2.6111394457
    (20,9)           8  4  2.2817057508  0.3849609          9.1035   1.75·10^5 1 000 000 009  5.7·10^3   −1.912            2               2.6111394457
    (24,1)           8  8  2.3925822571  2.9407511          17.9572  683       1 000 000 009  1.5·10^6   −4.770            2               2.6111394457

Rigor class, identical for all rows: analytic proof for `t ≥ X_0` (GRH(q), Chebyshev's
`θ(x) < x log 4`, and constants evaluated in Arb: `L_1(q)` and the constants `B`, `C`, `E_R`
of the explicit formulas — values of `L'/L` on `[1,2]` and Taylor coefficients of `L(s,χ)` at
0 and −1 — and the certified prime-square constant `c_1`), ball-arithmetic certificate for
`2 < t ≤ 10^9`; the two ranges overlap by the factor in the column "overlap" (at least 1200).
The thinnest case is `(20,1)`: `m = 0.385`, `X_0 = 8.17·10^5`.

**Corollary 1.3 (HKMN's normalisation; proved).** Let `c_0 = li(4) − 2 li(2) = 0.8772575348…
= −∫_0^2 li(t) dt`, and let `x_0 = 2.6111394456601689…` be the unique zero in `(1.4513…, ∞)`
of `g(x) = li(x²) − x li(x) = −∫_0^x li(t) dt`. Under GRH(q):

* for the 17 pairs with `q ≥ 2`:  `A_0(x;q,a) > 0` for `0 < x < x_0`, `= 0` at `x_0`, and
  `< 0` for **every** `x > x_0`;
* for `(1,1)`: `A_0(x) > 0` for `0 < x < x_1`, `< 0` for every `x > x_1`,
  `x_1 = 3.5651763735993442…`.

So with the lower limit 0 printed in the problem list, `x_1(q,a) = 2.611…` for all pairs
except `(1,1)`, where it is `3.565…`; the statement "for all x > 2" quoted there for q = 1
holds only with Johnston's lower limit 2.

**Corollary 1.4 (a short proof of the "RH ⟹" half of Johnston's Theorem 1.1).** The case
`(2,1)` of Theorem 1 gives `∫_2^x (π − li) < 0` for all `x > 2` under RH. The proof uses no
verification of `π(x) < li(x)` (Johnston's proof uses Büthe's range `x ≤ 10^19`). Besides RH
it uses the classical constants `Σ_ρ 1/|ρ(ρ+1)| = 0.0461176…`, `log 2π`, `(ζ'/ζ)(−1)`, the
computed constant `c_1(2,1) = 2.0712` (primes below `10^6`; attained just below 1423), and,
below the start `X_0 = 58.5` of the analytic range, the inequality `M_i < 0` (§2.5) on the 15
concave pieces of `W_1(·;2,1)` that begin at the primes 3, …, 53.

### 1.2 Exclusion under GRH alone (work order T2)

**Theorem 3 (asymptotics and the witness principle; proved).** Assume GRH(q). Then
`F_{q,a}` is a uniformly almost periodic function with `|F_{q,a}| ≤ L_1(q)`, and

    W_1(t;q,a) = t^{3/2} ( F_{q,a}(log t) − (2/3)N(q,a) ) + O_q(t^{4/3}),
    A(x;q,a)   = (x^{3/2}/log x) ( F_{q,a}(log x) − (2/3)N(q,a) ) + O_q(x^{3/2}/log² x).

(The second line is Theorem 7.1 of LWS; it is re-derived here from Lemma A.) Consequently:
if `F_{q,a}(y_0) > (2/3)N(q,a)` for **one** real `y_0` (of either sign), then
`A(x;q,a) > 0` — hence also `A_0(x;q,a) > 0` — for arbitrarily large x.

**Theorem 4 (non-square classes; proved; NV certified for q ≤ 400).** Let `q ≥ 3`, assume
GRH(q) and NV(q), and let a be a **non-square** mod q. Then `A(x;q,a)` is positive for
arbitrarily large x and negative for arbitrarily large x. NV(q) is certified in Arb for every
`q ≤ 400` (all 29 468 primitive characters of conductor ≤ 400; `t2_summary.out`). So: under
GRH(q) alone, no non-square class of any modulus `q ≤ 400` belongs to the answer of Problem 23.
This contains the 49 non-square classes at the 23 candidate moduli that LWS leave to LI, and
does not need their condition `A(q,a) ≠ 0`.

**Theorem 5 (the class a = 1; proved; q ≤ 400 certified).**
(a) For a primitive character χ of conductor `f > 1` and parity `a_χ ∈ {0,1}`, with no
hypothesis,

    c(χ) := −Re Σ_ρ 1/(ρ(ρ+1)) = Re[(L'/L)(2,χ) − (L'/L)(1,χ)] + δ(a_χ),    δ(0) = log 2,  δ(1) = 1 − log 2,

and `c(ζ) := −Σ_ρ 1/(ρ(ρ+1)) = log 2π + log 2 − 1/2 − (ζ'/ζ)(−1) = 0.045970…`. Under GRH(q),

    F_{q,1}(0) = Φ(q) := Σ_{χ mod q} c(χ*) = c(ζ) + φ(q)/2 − log 2 + Σ_{χ ≠ χ_0} Re[(L'/L)(2,χ*) − (L'/L)(1,χ*)]     (q ≥ 3).

(b) `Φ(q) > (2/3)N(q,1)` for **every** `q ≥ 1` outside the set

    {1, 2, 3, 4, 5, 6, 8, 10, 12, 15, 20, 24, 30}  ∪  {7, 9, 14, 16, 18}

(for `q ≤ 400`: the blocks `c_f` are Arb balls and both inequalities `Φ(q) ≷ (2/3)N(q,1)` are
decided on balls; for `q > 400`, under GRH(q), by the explicit bound
`Φ(q) ≥ φ(q)/2 + 1 − 2 log(4φ(q)(3 log q + 0.11)) − (log q)/3 − 0.04`, §4.4).
(c) Hence, under GRH(q) alone: `A(x;q,1) > 0` for arbitrarily large x for every q outside
that set. Together with Theorem 1 this classifies the pairs `(q,1)` under GRH alone, up to
five moduli:

    (q,1) belongs to the answer of Problem 23   if   q ∈ {1,2,3,4,5,6,8,10,12,15,20,24,30},
    undecided under GRH alone                   if   q ∈ {7, 9, 14, 16, 18},
    does not belong                             if   q is in neither of these two sets.

The witness is `y_0 = 0`, i.e. `x_0 = 1`: no prime and no zero is computed. The five undecided
moduli are exactly the "knife edges" of LWS §11.1 (`L_1(q)` exceeds `(2/3)N` by 1.9–3.3%);
there `Φ(q)/L_1(q) = 0.92, 0.93, 0.84` (q = 7, 9, 16), short of the required 0.98, 0.98, 0.97.

**Theorem 6 (all classes, q ≤ 120; certified witnesses).** Assume GRH(q). Call a pair
*excluded* if `A(x;q,a) > 0` for arbitrarily large x (so it is not in the answer of
Problem 23), and *member* if Theorem 1 applies.

(a) *The 23 candidate moduli of LWS* (`q ∈ {1,…,10, 12, 14, 15, 16, 18, 20, 21, 24, 30, 40, 42,
48, 60}`; 160 reduced classes). 18 classes are members. **128 classes are excluded**: all 116
non-square classes — each with an explicit certified witness `F_{q,a}(log x_0) > 0`,
`1 ≤ x_0 ≤ 10^7`, so that NV(q) is not even needed here; the class 1 at
`q = 21, 40, 42, 48, 60` (Theorem 5; confirmed independently from zeros at q = 21, 40, 48); and
the square classes `a ≠ 1` at these moduli,

    (21,4), (21,16), (42,25), (42,37), (40,9), (48,25), (60,49),

by zero-side witnesses found by lattice reduction at `|y_0| ≈ 10^{115}–10^{122}` (§4.5; certified
`F(y_0) ≥ 3.326, 3.325, 6.764, 7.038, 5.750` against `(2/3)N = 2.667, 2.667, 5.333, 5.333, 5.333`).
The remaining **14 classes are undecided** under GRH alone — exactly the square classes at the
five knife-edge moduli:

    (7,1), (7,2), (7,4), (9,1), (9,4), (9,7), (14,1), (14,9), (14,11), (16,1), (16,9), (18,1), (18,7), (18,13).

For comparison, LWS without LI: 18 members ("eventually"), 67 excluded (their §8.4; their
Proposition 8.3 assumes central non-vanishing too, which they certify at the candidate moduli,
as is done here with NV(q ≤ 400) — "GRH alone" is used in this sense on both sides), and
therefore 160 − 18 − 67 = 75 classes left to LI (my count from their numbers).

(b) <!--6B-START-->*All moduli q ≤ 60* (1084 pairs besides the 18 members; `t2_summary.out`). 1070 pairs are
excluded: every non-square class (explicit prime-side witness, or Theorem 4 — NV(q) is
certified in this range), the class 1 at every modulus outside the 18 of Theorem 5(b), and
every other square class for which a prime-side or zero-side witness is certified, where for
`q = 2m`, m odd, the witness may be the one of `(m, a mod m)` (the reduction (i) of §0:
`F_{2m,a} = F_{m,a}`; used for (22,3), (22,5), (22,15), (38,25), (38,35), (42,25), (42,37)). **14 square classes are
undecided**: the 14 of part (a) and no others.

(c) *Moduli 61 ≤ q ≤ 120.* Witnesses: prime side; Theorem 4 for non-square classes without
prime-side witness; zero side where it was run (`t2_summary.out` lists the pairs); and, for
{fmt(red120)}, the reduction to `q/2` and the zero-side witnesses at `q/2`. Every pair is
excluded except 18 square classes `a ≠ 1` for which the zero-side search was not run or not
finished: (83,69), (91,30), (93,76), (97,24), (99,49), (99,70), (103,14), (103,81), (104,9), (104,25), (105,4), (105,46), (107,40), (108,97), (109,12), (109,81), (112,81), (120,49). None of these is claimed to be hard.<!--6B-END-->

### 1.3 What is new relative to the sources

* **"Every x" versus "eventually".** LWS prove `A_1^π < 0` *for all sufficiently large x* at
  the 18 pairs: their Theorem 8.1 ends with "which is negative for large x", from the
  asymptotic formula (43) of their Theorem 7.1, whose error term `O_q(1/log x)` is not made
  explicit; I found no threshold stated in the sections I read (§§1, 7–11). Theorem 1 gives it for every `x > 2`, with the exact
  crossover `x_0`, `x_1` for the lower limit 0. Johnston's theorem is the case q = 1; the
  cases `(3,1), (4,1), (6,1)` were derivable from `explicit/EXPLICIT.md` plus Johnston
  (triage §6); the other 9 distinct cases, in particular q = 15 and q = 20 where the
  "difference with π" trick loses too much (§5), are not in LWS or in Johnston (I know of no other source; this is a statement about
  the texts I read, not a verified priority claim).
* **Method.** One exact identity (Lemma A) reduces the sign of the integrated π-error to the
  sign of the single function `W_1`, with a positive kernel and with the constant `−li(2)(x−2)`
  of the right sign; nothing is lost in passing from θ to π. In LWS (§7.1.3) and in Johnston
  ((3.1)–(3.6)) the kernel term is estimated, of size `x^{3/2}/log² x`; Johnston's proof under RH
  treats `x > 10^19` analytically and uses Büthe's `π(x) < li(x)` for `x ≤ 10^19` below
  (his Lemma 2.1 and the sentence "it suffices to consider x > 10^19" in his §3).
* **`L_1(q)` without zeros.** `Σ_ρ 1/|ρ(ρ+1)|` is evaluated from logarithmic derivatives of
  completed L-functions on `[1,2]` (Lemma S), 26-digit balls in seconds; LWS use
  Turing-certified zero lists to height 1000 plus a tail sandwich. Every value they print is
  reproduced: `L_1(7) = 1.35817204963` (their enclosure `[1.3581720402, 1.3581720541]`),
  `L_1(9) = 1.36364757066`, `L_1(16) = 2.75357978807`, `L_1(11) = 3.7735558`, zeta block
  `0.0461176444215`.
* **GRH-only exclusions.** LWS need LI for every square class outside the list and for 49
  non-square classes at the candidate moduli. Here, without LI: all non-square classes at
  every modulus `q ≤ 400` (Theorem 4 with the certified NV(q)); at `q > 400`, all non-square
  classes under GRH(q) + NV(q) — *not* under GRH(q) alone; the class 1 for every modulus
  outside five (Theorem 5, GRH(q) alone); and the explicit finite lists of Theorem 6.
  On LWS §8.5 ("GRH alone supplies no number here"): the mechanism — almost periodicity turns
  one value `F(y_0) > δ` into `F > δ/2` at arbitrarily large y — is theirs (it is the last
  step of their proof of Proposition 8.3, with threshold 0). What is added here is its use at
  the threshold `(2/3)N(q,a)` for square classes, together with certified values of `F`: the
  closed form at `y = 0` for `a = 1`, the prime-side formula (F), and the zero side with a
  closed-form remainder.
* An observation on LWS Proposition 8.3: its proof that `F_{q,a} ≢ 0`, as written, uses
  `A(q,a) ≠ 0` (a zero-counting argument); Lemma N shows that `F_{q,a}` is never constant,
  with no condition, which is what removes the 49 exceptions.
* All statements about LWS above refer to an unrefereed preprint, read in the v1.1 PDF;
  none of the theorems of this note depends on them.

---------------------------------------------------------------------------

## 2. Proof of Theorems 1 and 2

### 2.1 Lemma A (the identity; no hypothesis)

*For `q ≥ 1`, a reduced mod q, and `x > 2`:*

    A(x) = ∫_2^x (φ π(t;q,a) − li(t)) dt = − li(2)(x − 2) + W_1(x)/log x + ∫_2^x W_1(t) h_x''(t) dt,

    h_x(t) = (x − t)/log t,     h_x''(t) = (x + t)/(t² log² t) + 2(x − t)/(t² log³ t) > 0   (1 < t ≤ x).

*Proof.* Let ν be the signed measure `φ Σ_{p≡a (q)} log p · δ_p − dt` on `[2, x]`, and
`W(t) = ν([2,t]) = φ θ(t;q,a) − (t − 2)`, so `W_1(t) = ∫_2^t W`, `W_1(2) = 0`. For
`F ∈ C²[2,x]`, two integrations by parts give
`∫_{[2,x]} F dν = F(x)W(x) − F'(x)W_1(x) + ∫_2^x W_1 F''`. Take `F = h_x`: `h_x(x) = 0`,
`h_x'(x) = −1/log x`. The left side is
`φ Σ_{p≤x, p≡a} (x − p) − ∫_2^x (x − t) dt/log t = φ ∫_2^x π(t;q,a) dt − ∫_2^x (li(t) − li(2)) dt`
(Fubini). ∎

This is Ingham's second-order partial summation (Johnston's Lemma 2.7, `EXPLICIT.md` Lemma 2)
with the base point 2 chosen so that the linear term has the helpful sign. Checked numerically
to 28 digits at `x = 50.5, 200.25` for six pairs (`verify_identity.out`).

**Corollary.** If `W_1 ≤ 0` on `[2,x]` then `A(x) ≤ −li(2)(x−2) + W_1(x)/log x < 0`. For
q = 1 apply this to `(2,1)`: `∫_2^x (π − li) = A(x;2,1) + (x − 2)`.

So Theorem 1 is the statement `W_1(t;q,a) < 0` for all `t > 2`, i.e.

    φ(q) Σ_{p≤t, p≡a (q)} log p · (t − p)  <  (t − 2)²/2          (t > 2).

### 2.2 Lemma E (explicit formulas for ψ_1; no hypothesis)

Let χ be primitive mod `f > 1`, `Z_χ(x) = Σ_ρ x^{ρ+1}/(ρ(ρ+1))` (absolutely convergent). For
`x ≥ 1`:

    χ odd :   ψ_1(x,χ) = −Z_χ(x) − b x + log x + 1 + b_0 − T_o(x),            b = (L'/L)(0,χ),  (L'/L)(s,χ) = 1/(s+1) + b_0 + O(s+1)
    χ even:   ψ_1(x,χ) = −Z_χ(x) − x log x + (1 − b) x + c_{−1} − T_e(x),      (L'/L)(s,χ) = 1/s + b + O(s),  c_{−1} = (L'/L)(−1,χ)
    ζ     :   ψ_1(x)   = x²/2 − Z_ζ(x) − x log 2π + (ζ'/ζ)(−1) − T_e(x)        (Ingham, Thm 27; Johnston Lemma 2.6)

    T_e(x) = Σ_{r≥1} x^{1−2r}/(2r(2r−1)) = (x/2)[(1+u)log(1+u) + (1−u)log(1−u)],   T_o(x) = Σ_{n≥1} x^{−2n}/(2n(2n+1)) = 1 − ½ log(1−u²) − atanh(u)/u,   u = 1/x;
    0 < T_e(x) ≤ 0.6/x,   0 < T_o(x) ≤ 0.2/x²   (x ≥ 3);     T_e(1) = log 2,   T_o(1) = 1 − log 2.

*Proof.* `ψ_1(x,χ) = (1/2πi)∫_{(2)} (−L'/L)(s,χ) x^{s+1}/(s(s+1)) ds`; move the line to
`−∞` exactly as in `EXPLICIT.md`, Lemma 3(a) (the contour estimates there, Davenport Ch. 19,
do not use that χ is real or odd). Residues: at ρ, `−x^{ρ+1}/(ρ(ρ+1))`. Odd χ: as in
`EXPLICIT.md`. Even χ: `L(0,χ) = 0` is a simple zero, so at `s = 0` the integrand has a double
pole; with `−L'/L = −1/s − b + …` and `x^{s+1}/(s+1) = x(1 + s(log x − 1) + …)` the residue is
`−x log x + x − bx`; at `s = −1` (simple pole of the kernel only) it is `(L'/L)(−1,χ)`; at
`s = −2r`, `−x^{1−2r}/(2r(2r−1))`. ∎

Functional equation (used in §4): for both parities,
`b = −(L'/L)(1,χ̄) − log(f/π) + γ_E + log 2`, and the second constant
(`b_0` resp. `c_{−1}`) `= −(L'/L)(2,χ̄) − log(f/π) + γ_E + log 2 − 1`.

*Checks* (`verify_formula.out`, numerical): for nine primitive characters (odd/even,
real/complex, conductors 3–24) and for ζ, "direct prime-power sum − formula with the zeros
below height 200" is within the size of the omitted tail at `x = 1, 2.5, 5.5, 10, 30, 100`,
and at `x = 1` it agrees to three digits with the omitted tail of `Σ_ρ 1/|ρ|²` (a wrong
constant term would show as O(1)). The closed form of Theorem 5(a) agrees with
`−Re(b_0 − b + log 2)` resp. `−Re(1 − b + c_{−1} − log 2)` to 12 digits, which tests the
constants `b`, `b_0`, `c_{−1}` of `l1const.json` against values of `L'/L` at 1 and 2.

**Lemma I (imprimitive characters).** For χ mod q induced by χ*, put
`imp_{q,a}(x) = Σ_{p | q} φ(m_p) Σ_{k≥1, p^k ≡ a (m_p)} log p · (x − p^k)_+`,
`m_p = q/p^{v_p(q)}`. Then

    φ ψ_1(x;q,a) = Σ_{χ mod q} χ̄(a) ψ_1(x,χ*) − imp_{q,a}(x),        0 ≤ imp_{q,a}(x) ≤ I_q x log x,   I_q = Σ_{p|q} φ(m_p).

*Proof.* `ψ_1(x,χ*) − ψ_1(x,χ) = Σ_{p|q} Σ_k χ*(p^k) log p (x − p^k)_+`. For `p | q`,
`χ*(p^k) ≠ 0` only if `p ∤ f_χ`, i.e. χ is induced from a character χ' mod `m_p`, and then
`χ*(p^k) = χ'(p^k)`, `χ̄(a) = χ̄'(a)`; orthogonality mod `m_p` gives `φ(m_p) 1_{p^k ≡ a (m_p)}`. ∎

Summing Lemma E over χ with Lemma I: for `x ≥ 1`, with no hypothesis,

    φ ψ_1(x;q,a) = x²/2 − Σ_χ χ̄(a) Z_{χ*}(x) + ℓ_{q,a}(x) − imp_{q,a}(x),                                        (E)
    ℓ_{q,a}(x) = −n_e x log x + B x + n_o log x + C − (1 + n_e) T_e(x) − n_o T_o(x),

    n_e = Σ_{χ even, χ≠χ_0} χ̄(a),   n_o = Σ_{χ odd} χ̄(a)     (integers;  1 + n_e ≥ 0),
    B = −log 2π + Σ_{even, ≠χ_0} Re[χ̄(a)(1 − b_χ)] − Σ_{odd} Re[χ̄(a) b_χ],
    C = (ζ'/ζ)(−1) + Σ_{odd} Re[χ̄(a)(1 + b_{0,χ})] + Σ_{even, ≠χ_0} Re[χ̄(a) c_{−1,χ}].

Under GRH(q), `Σ_χ χ̄(a) Z_{χ*}(x) = −x^{3/2} F_{q,a}(log x)` and `|F_{q,a}| ≤ L_1(q)`:

    F_{q,a}(log x) = x^{−3/2} [ φ ψ_1(x;q,a) − x²/2 − ℓ_{q,a}(x) + imp_{q,a}(x) ]        (x ≥ 1).                 (F)

Every quantity on the right of (F) is computable exactly (a finite prime-power sum and four
constants): this is the "prime-side" evaluation of `F` used in §4.

### 2.3 Lemma S (the zero sums in closed form)

*Let χ be primitive mod `f ≥ 1` and assume GRH(χ). For real `α > 0` put
`V_χ(α) = Σ_ρ α/(α² + γ²) = Σ_ρ Re 1/(1/2 + α − ρ)`. Then (Davenport Ch. 12 (17)–(18); `s = 1/2 + α`)*

    V_χ(α) = Re (L'/L)(s,χ) + ½ log(f/π) + ½ ψ((s + a_χ)/2)                     (f > 1),
    V_ζ(α) = (ζ'/ζ)(s) + 1/(s−1) + 1/s − ½ log π + ½ ψ(s/2)                     (f = 1),

*and*

    S(χ) := Σ_ρ 1/|ρ(ρ+1)| = (1/π) ∫_0^π V_χ(α(θ))/α(θ) dθ,       α(θ) = (5/4 − cos θ)^{1/2} ∈ [1/2, 3/2].

*Moreover the n-point trapezoidal sum `T_n` of this integral satisfies
`T_n (4^n − 1)/(4^n + 1) ≤ S(χ) ≤ T_n`.*

*Proof.* With `A = 1/4 + γ²`: `1/|ρ(ρ+1)| = (A(A+2))^{-1/2} = (1/π)∫_0^π dθ/(A + 1 − cos θ)`
(since `∫_0^π dθ/(c − cos θ) = π(c²−1)^{-1/2}` for `c > 1`), and `A + 1 − cos θ = α(θ)² + γ²`;
sum over ρ (positive terms). For the error: `1/(c − cos θ) = (c²−1)^{-1/2}(1 + 2Σ_{k≥1} r^k cos kθ)`
with `r = c − (c²−1)^{1/2} ≤ 1/2` for `c ≥ 5/4`. So the integrand is
`a_0 + Σ_{k≥1} a_k cos kθ` with `a_0 = S(χ)` and `0 ≤ a_k ≤ 2^{1−k} S(χ)`, and the trapezoidal
sum with n intervals on `[0,π]` equals `a_0 + a_{2n} + a_{4n} + …`. ∎

`l1const.py` uses `n = 64` at 128 bits; only values of `L'/L` at real points of `[1,2]`
enter. Closed-form sandwich used as a second check: `1/(A+1) ≤ (A(A+2))^{-1/2} ≤ ½(1/A + 1/(A+2))`,
i.e. `V(√5/2)/(√5/2) ≤ S ≤ ½[2V(1/2) + (2/3)V(3/2)]`. Table (`l1const.out`):

    q        1,2       3,6       4         5,10      8         12        15,30     20        24
    L_1(q)   0.046118  0.158674  0.199912  0.599235  0.735631  0.632658  1.925733  2.281706  2.392582
    (2/3)N   0.666667  1.333333  1.333333  1.333333  2.666667  2.666667  2.666667  2.666667  5.333333
    margin m 0.620549  1.174660  1.133421  0.734099  1.931036  2.034009  0.740934  0.384961  2.940751

and the ten failing candidate moduli of LWS: `L_1 − (2/3)N` = 0.0248 (q = 7, 14), 0.0303
(9, 18), 0.0869 (16), 1.5655 (21, 42), 3.0768 (40), 3.4048 (48), 1.9604 (60). The 18 pairs
are confirmed as the square classes at the 13 moduli with `m > 0`.

### 2.4 The bound for `W_1` (proof of Theorem 2)

Let (q,a) be as in Theorem 1, `R = {b mod q : b² ≡ a}`, `N = #R`, `Θ(y) = φ Σ_{b∈R} θ(y;q,b)`,
`Θ_1 = ∫_0^· Θ`. Assume GRH(q). Let `t ≥ 8`.

*Step 1 (prime powers).* `ψ(u;q,a) = Σ_{k≥1} Σ_{b^k ≡ a} θ(u^{1/k};q,b)`. Dropping the terms
`k ≥ 3` (they are ≥ 0):

    φ θ_1(t;q,a) ≤ φ ψ_1(t;q,a) − Sq(t),       Sq(t) := ∫_0^t Θ(u^{1/2}) du = ∫_0^{√t} Θ(y) 2y dy.

*Step 2 (zeros).* By (E), `imp ≥ 0`, and `|Σ_χ χ̄(a) Z_{χ*}(t)| ≤ L_1(q) t^{3/2}`:
`φ ψ_1(t;q,a) ≤ t²/2 + L_1 t^{3/2} + ℓ_{q,a}(t)`, and
`ℓ_{q,a}(t) ≤ −n_e t log t + B t + n_o log t + C + max(0,−n_o)·0.2/t²` (the two series enter
with coefficients `−(1+n_e) ≤ 0` and `−n_o`).

*Step 3 (prime squares, `√t ≤ Y_1 = 10^6`; certified constant).* Let

    c_1 = c_1(q,a) ≥ sup_{0<y≤Y_1} (N y − Θ(y))/√y

(certified in ball arithmetic, `w1_cert.py` part B; between consecutive primes of R the
function `(Ny − const)/√y` is increasing, so the supremum is a maximum over left limits at
primes). Then `Sq(t) ≥ ∫_0^{√t} (N y − c_1 √y) 2y dy = (2/3)N t^{3/2} − (4/5) c_1 t^{5/4}`.

*Step 4 (prime squares, `√t > Y_1`; from GRH(q) itself).* Write `Y = √t` and
`Sq(t) = 2YΘ_1(Y) − 2∫_0^Y Θ_1(y) dy`. Upper bound: by (E) for the classes `b ∈ R`
(`Σ_{b∈R} |Σ_χ χ̄(b) Z_{χ*}| ≤ N L_1 y^{3/2}`, `imp ≥ 0`, `θ_1 ≤ ψ_1`),

    Θ_1(y) ≤ N y²/2 + N L_1 y^{3/2} + E_R y log y     (y ≥ 3),        Θ_1(y) ≤ φ log 2 · (y − 2)_+   (y ≤ 3),

where `|Σ_{b∈R} ℓ_{q,b}(y)| ≤ E_R y log y` for `y ≥ 3`,
`E_R = |n_e^R| + |B^R|/log 3 + |n_o^R|/3 + |C^R|/(3 log 3) + 0.6(N + |n_e^R|)/(9 log 3) + 0.2|n_o^R|/(27 log 3)`
and `X^R = Σ_{b∈R} X(q,b)`. Lower bound: by (E), Lemma I, and `θ(v;q,·) ≤ θ(v) ≤ v log 4`
(Chebyshev: `Π_{p≤n} p < 4^n`; Hardy–Wright §22.2, theorem number not re-checked) for the prime
powers `k ≥ 2`,

    Θ_1(Y) ≥ N Y²/2 − λ(Y) Y^{3/2},
    λ(Y) = N[L_1 + (2/3)φ log 4] + Nφ log 4 [(3/4)Y^{−1/6} + (log Y/log 2) Y^{−1/4}] + (N I_q + E_R) Y^{−1/2} log Y.

Hence `Sq(t) ≥ (2/3)N t^{3/2} − K_B(t) t^{5/4}` with

    K_B(t) = 2λ(√t) + (4/5) N L_1 + E_R t^{−1/4} (log t)/2 + 2φ log 2 · t^{−5/4}      (decreasing for t ≥ e^8)

(the last term bounds `2∫_0^3 Θ_1 ≤ φ log 2` with a factor 2 to spare).

*Step 5.* Adding up and using `−(t−2)²/2 = −t²/2 + 2t − 2`:

    W_1(t) ≤ −m t^{3/2} + K(t) t^{5/4} + [ℓ^+(t) + 2t − 2]  ≤  t^{3/2} U(t),
    r_+(t) = max(0,−n_e) t^{−1/2} log t + max(0, B+2) t^{−1/2} + max(0, n_o) t^{−3/2} log t + max(0, C−2) t^{−3/2} + max(0,−n_o)·0.2·t^{−7/2}.

Each term of `r_+` is decreasing for `t ≥ e²`. ∎

Constants (`w1_analytic.out`, `w1_analytic.json`; all evaluated in Arb, upper ends used):

    (q,a)    n_e n_o   B        C        E_R      I_q  K_A=(4/5)c_1  K_B(10^12)  K_B(∞)   U(10^6)   U(10^12 −)  U(10^12 +)  −m
    (2,1)     0   0   −1.8379  +1.9851   2.336    1    1.657         4.06        1.98     −0.5680   −0.6189     −0.6165     −0.6205
    (3,1)     0   1   −2.7861  +3.1289   4.672    1    3.331         16.35       8.28     −1.0693   −1.1713     −1.1583     −1.1747
    (4,1)     0   1   −2.6211  +2.9248   4.672    1    3.314         16.58       8.51     −1.0286   −1.1301     −1.1168     −1.1334
    (5,1)     1   2   −2.1114  +3.0315   4.640    1    5.733         34.03       18.14    −0.5528   −0.7284     −0.7001     −0.7341
    (5,4)     1  −2   +0.4797  −0.0240   8.946    1    3.277         34.21       18.14    −0.6280   −0.7308     −0.6999     −0.7341
    (8,1)     1   2   −1.6810  +2.4792   9.344    1    6.628         69.60       37.81    −1.7211   −1.9244     −1.8614     −1.9310
    (12,1)    1   2   −2.0227  +2.9142   9.344    4    7.183         68.78       36.66    −1.8069   −2.0268     −1.9652     −2.0340
    (15,1)    3   4   +0.3635  +1.3443   9.281    6    11.466        144.34      80.72    −0.3760   −0.7295     −0.5966     −0.7409
    (15,4)   −1   0   −3.8325  +3.7149   17.892   6    7.709         144.70      80.72    −0.4833   −0.7332     −0.5962     −0.7409
    (20,1)    3   4   +1.1433  +0.3025   9.281    6    11.466        148.33      84.70    −0.0192   −0.3735     −0.2366     −0.3850
    (20,9)   −1   0   −3.8152  +3.6891   17.892   6    7.283         148.68      84.70    −0.1408   −0.3777     −0.2362     −0.3850
    (24,1)    3   4   +1.3550  +0.0285   18.687   6    14.366        299.14      171.89   −2.4831   −2.9264     −2.6416     −2.9408

`X_0` is the point (rounded up to 3 digits) from which `U < 0` on `[X_0, 10^12]`; `U(X_0) < 0`
is re-verified in Arb at the printed `X_0`. Consistency test (passed for every decade above
`X_0`, all cases): the certified bound `max M_i/p_{i+1}^{3/2}` over the pieces starting at
primes in `[10^k, 10^{k+1})` is below `U(10^k)`.
The prime-square step uses only `Y_1 = 10^6`; `(4/5)c_1 t^{-1/4}` is where `X_0` comes from.

### 2.5 The finite range (certified)

On `[2, p_1]` (`p_1` the first prime ≡ a mod q), `W_1(t) = −(t−2)²/2 < 0` for `t > 2`. On
`[p_i, p_{i+1}]` (consecutive primes ≡ a mod q), with `ϑ_i = Σ_{j≤i} log p_j`,
`Λ_i = Σ_{j≤i} p_j log p_j`: `W_1(t) = φ(ϑ_i t − Λ_i) − (t−2)²/2`, a concave quadratic whose
maximum over all real t is

    M_i = (φ ϑ_i)²/2 + 2φ ϑ_i − φ Λ_i.

`w1_cert.py` part A certifies `M_i < 0` (128-bit balls; the upper end of the ball is negative)
for **every** prime `p_i ≡ a (mod q)`, `p_i ≤ 10^9`, in all 12 cases: 50 847 533 primes for
`(2,1)` down to 6 354 932 for `(24,1)`, no failure, no undecided comparison
(`raw/w1_cert_*.json`, `raw/w1_cert_*.log`). Hence `W_1(t) < 0` for `2 < t ≤ 10^9`. The primes
come from an integer sieve (numpy); no floating-point comparison decides anything.
The least margins occur at the first primes: e.g. `sup W_1/t^{3/2} ≤ −0.0443` for `(2,1)`
(at `p = 3`), `−0.0955` for `(4,1)` (`p = 5`), `−0.1752` for `(15,4)` (`p = 19`).

**Proof of Theorem 1.** `W_1 < 0` on `(2, 10^9]` by §2.5 and on `[X_0, ∞)` by Theorem 2, and
`X_0 < 10^9`. Then the Corollary of Lemma A. ∎

Independent numerical check (`verify_identity.out`, mpmath, 30 digits): on the grid
`p_1 ≤ x ≤ 2000`, step 0.37, for the 12 cases, `A(x) ≤ −li(2)(x−2) + W_1(x)/log x` and
`W_1 < 0` hold at every grid point.

---------------------------------------------------------------------------

## 3. Lower limit 0 versus 2 (proof of Corollary 1.3)

`∫_0^x li(t) dt = x li(x) − li(x²)` for `x > 0` (differentiate; both sides vanish at 0+ and are
continuous at 1). So `A_0(x) = g(x) + φ Σ_{p≤x, p≡a} (x − p)` with `g(x) = li(x²) − x li(x)`,
`g' = −li`: g increases on `(0, μ)`, decreases on `(μ, ∞)` (`μ = 1.45136923…`, the zero of
li), `g(0+) = 0`, `g(2) = c_0 = 0.87725753480406…`, `g(3) = −0.7695…`; so g has exactly one
zero `x_0 = 2.61113944566016890…` in `(μ, ∞)` (`x1.out`, certified by sign changes in Arb).

*q ≥ 2.* The first prime `p_1 ≡ a (mod q)` is ≥ 3 for all 17 pairs, so `A_0 = g` on `(0, p_1]`:
positive on `(0, x_0)`, negative on `(x_0, p_1]`. For `x ≥ p_1 ≥ 3`, by Theorem 1,
`A_0(x) = c_0 + A(x) ≤ c_0 − li(2)(x − 2) ≤ 0.8773 − 1.0451 < 0`.

*q = 1.* `A_0(x) = g(x) + Σ_{p≤x}(x − p)`, `A_0' = π − li` off the primes, and `π(p) < li(p)`
at each prime `p ≤ 31` (checked), so `A_0` is decreasing on `[2, 31]`; `A_0(3) = +0.2305`,
`A_0(5^-) = −1.6600`: one zero `x_1 = 3.56517637359934422…` in `(3,5)`. For `x ≥ 22`,
`A_0(x) ≤ c_0 − (li(2) − 1)(x − 2) < 0` by Theorem 1 (`c_0/(li(2)−1) + 2 = 21.42`). ∎

---------------------------------------------------------------------------

## 4. Proofs for §1.2 (GRH-only exclusion)

### 4.1 Proof of Theorem 3

*Almost periodicity.* `F_{q,a}(y) = Σ_j c_j e^{iγ_j y}` with `Σ|c_j| = L_1(q) < ∞` (Lemma S).
Given `ε > 0`, cut the series at J terms with tail `< ε/3`; by Dirichlet's simultaneous
approximation theorem there are arbitrarily large τ with `|e^{iγ_j τ} − 1| < ε/(3L_1)` for
`j ≤ J`, and then `|F(y+τ) − F(y)| < ε` for all real y. So if `F(y_0) > c`, then
`F > c + η` on an unbounded set of y, `η = (F(y_0) − c)/2`.

*Asymptotics.* `W_1(t) = [φψ_1(t;q,a) − t²/2] − Sq(t) − P_{≥3}(t) + 2t − 2` with
`0 ≤ P_{≥3}(t) ≤ φ log 4 [(3/4)t^{4/3} + (log t/log 2) t^{5/4}]` (prime powers `k ≥ 3`,
Chebyshev). By (F), `φψ_1 − t²/2 = t^{3/2} F_{q,a}(log t) + O_q(t log t)`. For the squares:
`Sq ≡ 0` if a is not a square; otherwise `Sq(t) = (2/3)N t^{3/2} + O_q(t^{5/4})` — the lower
bound is Step 4 of §2.4 (valid for all `Y ≥ 3` with `λ* = sup_{Y≥3} λ(Y)`), the upper bound is
the same computation with the roles of the two bounds for `Θ_1` exchanged. Finally
`0 < h_x''(t) ≤ 4.9 x/(t² log² t)` on `[2,x]`, so `|∫_2^x W_1 h_x''| ≪_q x^{3/2}/log² x`, and
Lemma A gives the formula for `A(x)`. The last assertion follows: at `x = e^y` with
`F(y) > (2/3)N + η`, `A(x) ≥ (x^{3/2}/log x)(η − O(1/log x)) > 0` for large y. ∎

### 4.2 Lemma N and the proof of Theorem 4

**Lemma N.** *Assume GRH(q). For every reduced a, `F_{q,a}` is not constant.*

*Proof.* If `F_{q,a} ≡ c`, then by (F), for all `x ≥ 1`,
`φψ_1(x;q,a) + imp_{q,a}(x) = x²/2 + c x^{3/2} + ℓ_{q,a}(x)`. The left side is piecewise
linear in x (its second derivative vanishes on every interval free of prime powers); the
second derivative of the right side is `1 + (3/4)c x^{−1/2} + ℓ''(x) → 1` as `x → ∞`. ∎

*Proof of Theorem 4.* The Bohr mean of `F_{q,a}` is its coefficient at frequency 0, namely
`−(4/3) Σ_χ χ̄(a) m_χ(1/2) = 0` under NV(q). A uniformly almost periodic real function of mean
0 that is not identically zero takes positive and negative values (if `F ≤ 0` and
`F(y_1) < 0`, then `F < F(y_1)/2` on a relatively dense set of intervals of fixed length, and
the mean is negative). By Lemma N, `sup F > 0 > inf F`. For a non-square a, `N(q,a) = 0`; apply
Theorem 3 to both signs. ∎

(Without NV(q) the same argument shows: if `−(4/3)Σ_χ χ̄(a) m_χ(1/2) ≥ (2/3)N(q,a)`, then
`A > 0` for arbitrarily large x; central zeros are part of the bias, as LWS Theorem 11.1 says.)

*Certificate of NV.* `t2_closed.py` evaluates `L(1/2,χ)` in 96-bit balls for all 29 468
primitive characters of conductor `3 ≤ f ≤ 400`; every ball excludes 0; the smallest certified
lower bound of `|L(1/2,χ)|` is `1.40·10^{-4}`, at conductor 101. (`ζ(1/2) ≠ 0`.)

Remark (why scans are not enough). The character `19.3` (Conrey label; order 18, odd) has a
zero at height `γ = −0.018956…` (and `19.13`, its conjugate, at `+0.018956…`);
`|L(1/2,χ)| = 0.01914`. So `F_{19,a}` contains the term `−(8/3)Re(χ̄(a)e^{iγy})` up to a factor
`1 + O(γ)`, whose phase moves by at most 0.31 rad over the whole scan range
`0 ≤ y ≤ log 10^7 = 16.1`: for `a = 14, 15, 12` the prime-side scan finds `F_{19,a}(log x) < 0`
on all of `1 ≤ x ≤ 10^7`. Theorem 4 covers these classes; an explicit witness needs
`|y|` of the order of `π/γ ≈ 166` (the zero side reaches it at once).

### 4.3 Proof of Theorem 5(a) and of (c)

Zeros of `L(s,χ)` are symmetric under `ρ ↦ 1 − ρ̄`, so `Σ_ρ Re(1/ρ) = Σ_ρ Re 1/(1−ρ)` and
`Σ_ρ Re 1/(ρ+1) = Σ_ρ Re 1/(2−ρ)`; by Davenport Ch. 12 (17)–(18) at `s = 1` and `s = 2`
(the terms `½ log(f/π)` cancel):

    −Re Σ_ρ 1/(ρ(ρ+1)) = Σ_ρ Re 1/(2−ρ) − Σ_ρ Re 1/(1−ρ) = Re[(L'/L)(2,χ) − (L'/L)(1,χ)] + ½[ψ(1 + a_χ/2) − ψ((1+a_χ)/2)],

and `½[ψ(1) − ψ(1/2)] = log 2`, `½[ψ(3/2) − ψ(1)] = 1 − log 2`. For ζ, put `x = 1` in Lemma E:
`0 = 1/2 − Z_ζ(1) − log 2π + (ζ'/ζ)(−1) − log 2`. Under GRH(q), `F_{q,1}(0) = −Σ_χ Z_{χ*}(1)`,
which is real, hence equal to `Σ_χ c(χ*)`. For `q ≥ 3` there are `φ/2` odd and `φ/2 − 1` even
non-principal characters, and `(φ/2)(1 − log 2) + (φ/2 − 1) log 2 = φ/2 − log 2`. Part (c) is
Theorem 3 with `y_0 = 0` (`F_{q,1}(0) = Φ(q) > (2/3)N(q,1)`) and Theorem 1. ∎

Why `y = 0` is special: for `ρ = 1/2 + iγ`,
`−1/(ρ(ρ+1)) = (γ² − 3/4 + 2iγ)/((1/4+γ²)(9/4+γ²))`, whose real part divided by the weight
`1/|ρ(ρ+1)|` is `(γ² − 3/4)/((1/4+γ²)(9/4+γ²))^{1/2}`: negative for `|γ| < 0.866`, equal to
0.63 at `γ = 2`, 0.95 at `γ = 6`, and tending to 1. So at `y = 0` all high zeros of all
characters are in phase, with no appeal to linear independence; the loss comes only from low
zeros. This is why `Φ(q)/L_1(q)` is 0.98–0.997 for q ≤ 4 and decreases to 0.64–0.75 at
q = 40, 48, 60 (table in `t2_summary.out`). For `a ≠ 1` the phases `χ̄(a)` destroy the
alignment, and no closed-form witness of this kind is available in general.

### 4.4 Proof of Theorem 5(b)

*q ≤ 400 (certified).* `t2_closed.py` computes `c_f = Σ_{χ primitive mod f} c(χ)` in 96-bit
balls for every conductor `f ≤ 400`; `Φ(q) = c(ζ) + Σ_{f | q, f > 1} c_f`. Result
(`t2_summary.out`): `Φ(q) > (2/3)N(q,1)` for all `q ≤ 400` except the 18 moduli listed. The
smallest relative excess is 2.9% at `q = 60` (`Φ = 5.48583` against `5.33333`); then
`q = 48` (`5.56748`), `q = 40` (`5.75620`), `q = 21, 42` (`3.62950` against `2.66667`).

*q > 400 (proved, under GRH(q)).* Let `T_2 = Σ_{χ≠χ_0} Re (L'/L)(2,χ*)`,
`T_1 = −Σ_{χ≠χ_0} Re (L'/L)(1,χ*)`, so `Φ(q) = c(ζ) + φ/2 − log 2 + T_2 + T_1`.

(i) `Σ_{χ mod q} χ*(n) = φ 1_{n≡1 (q)}` for `(n,q) = 1`, and `= φ(m_p) 1_{p^k ≡ 1 (m_p)}` for
`n = p^k`, `p | q` (as in Lemma I). Hence
`T_2 = −(ζ'/ζ)(2) − φ Σ_{n≡1 (q), n>1} Λ(n)/n² − Σ_{p|q} φ(m_p) Σ_{k: p^k≡1 (m_p)} log p/p^{2k}`.
The first sum is `≤ q^{−2} Σ_{j≥1} (log q + log(j+1))/j² ≤ q^{−2}(ζ(2) log q + 1.81)`; in the
second, `p^k ≥ m_p + 1` and the admissible k are the multiples of the order of p, so the inner
sum is `≤ log p/(m_p(m_p+2))` (`≤ log p/(p²−1)` if `m_p = 1`), and the term is `≤ (log p)/3`. So

    T_2 ≥ −(ζ'/ζ)(2) − (φ/q²)(ζ(2) log q + 1.81) − (log q)/3.

(ii) For real `s > 1`, `−Σ_{χ≠χ_0}(L'/L)(s,χ) = s∫_1^∞ φψ(x;q,1) x^{−s−1} dx + (ζ'/ζ)(s) + Σ_{p|q} log p/(p^s − 1)`.
With `E_1(x) = φψ_1(x;q,1) − x²/2` (`E_1(1) = −1/2`, `E_1(x) = O_q(x^{3/2})` under GRH(q)),
two integrations by parts and `s → 1+` give

    −Σ_{χ≠χ_0} (L'/L)(1,χ) = 1 + γ_E + Σ_{p|q} log p/(p−1) + I,       I := 1/2 + 2∫_1^∞ E_1(x) x^{−3} dx.

Passing to primitive characters, `Σ_{χ≠χ_0}[(L'/L)(1,χ) − (L'/L)(1,χ*)] = Σ_{p|q}[φ(m_p) log p/(p^{o_p} − 1) − log p/(p−1)]`
(`o_p` = order of p mod `m_p`), so

    T_1 = 1 + γ_E + I + Σ_{p|q} φ(m_p) log p/(p^{o_p} − 1)  ≥  1 + γ_E + I.

(Checked numerically for q = 3, 4, 5, 7, 12, 15, 16, 21, 40, 60, 97 with the integral
truncated at `10^7`: agreement to 3–6 digits, within the truncation error; `t2_check_identity.out`.)

(iii) Lower bound for I. On `[1, X_1]` use `E_1(x) ≥ −x²/2`; on `[X_1, ∞)` use (E):
`E_1(x) ≥ −L_1(q)x^{3/2} − |ℓ_{q,1}(x)| − imp_{q,1}(x)`. So
`I ≥ 1/2 − log X_1 − 4L_1 X_1^{−1/2} − R`, where
`R = 2[(n_e + ωφ)(log X_1 + 1) + |B|]/X_1 + [n_o(log X_1 + 1/2) + |C|]/X_1² + 0.4(1+n_e)/X_1³ + 0.1 n_o/X_1⁴`
(`n_e = φ/2 − 1`, `n_o = φ/2`, `imp ≤ ω(q)φ x log x`). A priori bounds under GRH. Since
`(1/2)/(1/4+γ²) ≤ 3·(3/2)/(9/4+γ²)`, for every `χ ≠ χ_0`:
`0 ≤ V_{χ*}(1/2) ≤ 3V_{χ*}(3/2) ≤ 3[0.59 + ½ log(f/π)]` (Lemma S at `s = 2`,
`|(L'/L)(2,χ*)| ≤ −(ζ'/ζ)(2) = 0.570`, `½ψ((2+a)/2) ≤ 0.019`). Hence, by Lemma S at `s = 1`,
`−½ log(f/π) ≤ Re (L'/L)(1,χ*) ≤ log f + 0.91` (the two digamma terms taken with the parity of
χ: `log f + 0.681` for even, `log f + 0.909` for odd χ; combining their worst cases blindly
gives `log f + 1.61`), so in any case `|Re (L'/L)(1,χ*)| ≤ log q + 1.62`; and
`L_1(q) ≤ 2V_ζ(1/2) + Σ_{χ≠χ_0} 2V_{χ*}(1/2) ≤ 0.047 + (φ−1)(3 log q + 0.106) ≤ L^+ := φ(3 log q + 0.11)`
(for ζ the term is `2V_ζ(1/2) = 0.0462`, not bounded by the displayed inequality). Inserting
`|Re (L'/L)(1,χ̄*)| ≤ log q + 1.62`, `|(L'/L)(2,χ̄*)| ≤ 0.570`, `|log(f/π)| ≤ log q` in the
formulas of §2.2 for `b` and for the second constant gives `|Re b_χ| ≤ 2 log q + 2.9`,
`|Re b_{0,χ}|, |Re c_{−1,χ}| ≤ log q + 0.9`, and then, from the definitions of B and C with
`a = 1`: `|B| ≤ φ(2 log q + 4.4)`, `|C| ≤ φ(log q + 2.4)`.
Take `X_1 = 16(L^+)²`: `I ≥ −1/2 − 2 log(4L^+) − R`.

(iv) Collecting, and using `c(ζ) − log 2 − (ζ'/ζ)(2) + 1/2 + γ_E = log 2π − (ζ'/ζ)(−1) − (ζ'/ζ)(2) + γ_E = 1`
(functional equation of ζ at s = 2):

    Φ(q) ≥ LB(q) := φ/2 + 1 − 2 log(4φ(3 log q + 0.11)) − (log q)/3 − (φ/q²)(ζ(2) log q + 1.81) − R(q)      (q ≥ 3, GRH(q)).

`t2_assemble.py` (floating point, exact φ, ω, N from a sieve): `Φ(q) − LB(q) ≥ 6.01` for all
`3 ≤ q ≤ 400` (consistency of the derivation with the certified values);
`LB(q) − (2/3)N(q,1) ≥ 18.6` for all `400 < q ≤ 10^6` — this inequality is then **decided in Arb**
for each of these q (64-bit balls, exact integer φ, ω, N; no failure; smallest certified
difference 18.6056) — (minimum at `q = 420`, where
`LB/((2/3)N) = 2.74`), with `R(q) ≤ 3.9·10^{−4}` there. For `q > 10^6`: `N(q,1) ≤ φ(q)/4` as
soon as `q > 168` (the index of the 2-torsion in `(Z/q)^*` is `Π r(p^e)`, `r = 1` for
`2, 4, 8, 3`, `r = 2` for 5 and 16, `r = 3` for 7 and 9, otherwise ≥ 4), `φ(q) ≥ (q/2)^{1/2} ≥ 707`,
`log q ≤ 2 log φ + log 2`; the last two terms of `LB` are `< 10^{−3}` (with `ω(q) ≤ log q/log 2`,
`log X_1 ≤ 4 log q` and `X_1 ≥ 144 φ² log² q`: `R(q) ≤ 13 log² q/(72 φ log² q) + … < 3·10^{−4}`, and
`(φ/q²)(ζ(2) log q + 1.81) < 3·10^{−5}`); and
`g(φ) = φ/3 + 1 − 2 log(4φ(6 log φ + 2.2)) − (2 log φ + 0.7)/3 − 0.04` is increasing for `φ ≥ 20`
with `g(707) = 208.7 > 0`; so `LB(q) > φ/6 ≥ (2/3)N(q,1)`. ∎

### 4.5 Witnesses for the remaining classes (Theorem 6)

Two certified evaluations of `F_{q,a}` are used (the second with two ways of *proposing* the
point `y_0`); both only need GRH(q) to be *meaningful*
(they are identities about zeros on the line) and neither needs a complete list of zeros.

* **Prime side** (`t2_witness.py`): formula (F) at a rational `x_0 ≥ 1`: a prime-power sum up
  to `x_0` and the four constants `n_e, n_o, B, C`, all in 128-bit balls. A floating-point scan
  of `1 ≤ x ≤ 10^7` proposes `x_0`; the certificate is the Arb evaluation at `x_0`.
* **Zero side** (`t2_zeros.py`): brackets of zeros of Hardy's `Z(t,χ)` with certified opposite
  signs at the ends (each bracket contains a zero; one zero per bracket is counted), giving a
  finite set `Z` of zeros; then for all real y

      F_{q,a}(y) ≥ F_Z(y) − r,      r = Σ_χ [ S(χ*) − Σ_{ρ∈Z_χ} 1/|ρ(ρ+1)| ],

  because the total weight `S(χ*)` is known exactly (Lemma S). No completeness of `Z` is
  needed (a missed zero is part of `r`); this replaces Turing's method. A floating-point scan
  over `|y| ≤ Y_max` proposes `y_0`; the certificate is the ball evaluation of `F_Z(y_0) − r`.

**Results** (`t2_summary.out`, `raw/t2_witness_120_*.jsonl`, `raw/t2_zeros_*.jsonl`).

*Prime side, q ≤ 120, X = 10^7.* Of the 4368 pairs outside the 18 members, 4266 receive a
certified witness: 2767 of the 2788 non-square classes (the other 21 — e.g. `(19,12)`,
`(19,14)`, `(19,15)`, see the Remark in §4.2 — are covered by Theorem 4), and 1499 of the 1580
square classes, among them `(q,1)` for all 102 moduli `q ≤ 120` outside the 18 of Theorem 5(b)
(at `x_0 = 1`, in agreement with the closed form). 81 square classes get no prime-side
witness.

*Zero side at the candidate moduli* (zeros to height 400; `|y| ≤ 5000`, step 0.05):

    (q,a)     (2/3)N    remainder r   best y_0        F_Z(y_0)    certified F(y_0) ≥     verdict
    (21,1)    2.6667    0.0699        −0.00094        3.5631      3.4932                 excluded   (closed form: Φ(21) = 3.6295)
    (40,1)    5.3333    0.0987        −0.00117        5.6654      5.5668                 excluded   (Φ(40) = 5.7562)
    (48,1)    5.3333    0.0991        −0.00118        5.4764      5.3773                 excluded   (Φ(48) = 5.5675)
    (60,1)    5.3333    0.0968        −0.00109        5.3957      5.2988                 (zero side inconclusive; closed form Φ(60) = 5.4858 > 5.3333 decides)
    (21,4)    2.6667    0.0699        −3221.68        1.9710      1.9011                 undecided
    (21,16)   2.6667    0.0699        +4876.70        1.8935      1.8236                 undecided
    (40,9)    5.3333    0.0987        +4237.55        3.9646      3.8659                 undecided
    (48,25)   5.3333    0.0991        −2895.47        4.2410      4.1419                 undecided
    (60,49)   5.3333    0.0968        −4219.91        3.2497      3.1529                 undecided

For `a = 1` the scan finds its maximum at `y ≈ 0`, and `F_Z(≈0)` equals `Φ(q)` up to the
remainder — an independent confirmation of Theorem 5(a) from actual zeros. For the square
classes `a ≠ 1` at these moduli the largest value found in `|y| ≤ 5000` is 61–80% of the
threshold (45–49% of `L_1`): exhaustive search is hopeless there. They are settled by the next
method.

*Zero side with lattice reduction* (`t2_refine.py`, `t2_lll.py`, `t2_lll_check.py`; added after
the independent audit). The sum over the located zeros is `F_Z(y) = Σ_j A_j cos(γ_j y + θ_j)`.
For the K heaviest terms (K = 40–100) one wants `γ_j y + θ_j ≡ 0 (mod 2π)` simultaneously — an
inhomogeneous simultaneous approximation problem, solved approximately by LLL on the lattice
spanned by the rows `(…, s_j·2π·2^ν, …)`, `(s_j γ_j 2^ν)_j ⊕ ε`, `(s_j θ_j 2^ν)_j ⊕ M` with weights
`s_j = (A_j/A_1)^{1/2}` (the method of Odlyzko–te Riele for the Mertens conjecture, with an
embedded inhomogeneous row). The reduced vector containing the last row gives `y_0` with all K
phases small (weighted alignment 0.98–0.9998), at `|y_0|` up to `2^{420}`. This is only a
*proposal*. The certificate is as before: `F(y_0) ≥ F_Z(y_0) − r`, with `F_Z(y_0)` evaluated in
ball arithmetic from **all** located zeros. Because `|y_0| ≈ 10^{120}`, every zero bracket is
first refined to width `2^{-480}` (`t2_refine.py`: secant proposal; the narrow bracket is accepted
only with certified opposite signs of Hardy's Z at its ends, and it lies inside the first-stage
bracket, so one zero per bracket and disjointness are preserved); `y_0` is an exact dyadic
rational, recorded in full in `raw/t2_zeros_lll_*.jsonl`.

<!--LL-START-->    (q,a)        (2/3)N    remainder r   y_0 (approx.)    F_Z(y_0)    certified F(y_0) ≥   / L_1(q)    parameters
    (21,4)    2.6667    0.0699        +4.8870e+114        3.3960      3.3260        0.786       K = 80, nu = 400, T = 400, 4070 cosines
    (21,16)   2.6667    0.0699        +7.0013e+119        3.3948      3.3248        0.786       K = 80, nu = 420, T = 400, 4070 cosines
    (28,9)    2.6667    0.0914        +1.6946e+53        3.7720      3.6806        0.735       K = 40, nu = 200, T = 300, 2970 cosines
    (40,9)    5.3333    0.0987        -2.2098e+115        6.8628      6.7641        0.804       K = 100, nu = 400, T = 400, 5852 cosines
    (48,25)   5.3333    0.0991        -1.4827e+119        7.1375      7.0384        0.805       K = 100, nu = 420, T = 400, 5891 cosines
    (60,49)   5.3333    0.0968        -8.8072e+121        5.8470      5.7502        0.788       K = 100, nu = 420, T = 400, 5707 cosines<!--LL-END-->

Independent re-evaluation (`t2_lll_check.py`, mpmath at 1100 bits, phases reduced mod 2π
explicitly, character values from exact exponents, nesting and disjointness of the refined
brackets re-checked): agreement with the Arb balls to 15 digits for every record
(`t2_lll_check.txt`). With these, **every square class at the candidate moduli outside the five
knife-edge moduli is decided**. What the method cannot do is discussed in §5.

*Zero side at other moduli q ≤ 60* (`run_zeros_queue.sh`, `run_zeros_queue2.sh`; zeros to
height 100–300, `|y| ≤ 2000–3000`; a second pass `|y| ≤ 40000` for q = 28, 36 and a third pass
`|y| ≤ 300000` for `(28,9)`, in which the floating-point proposal step uses only the zeros below
height 60 resp. 40 and the certificate all located zeros).
<!--6Q-START-->Certified zero-side witnesses for square classes `a ≠ 1` (all have `F(y_0) ≥` the last column `> (2/3)N`):

    (q,a)        (2/3)N    remainder r   y_0                F_Z(y_0)    certified F(y_0) ≥
    (11,3)    1.3333    0.0745         -1727.662581      1.8617      1.7872     T = 300, 4800 zeros
    (11,4)    1.3333    0.0745          1123.433616      1.9180      1.8435     T = 300, 4800 zeros
    (11,5)    1.3333    0.0745           666.831337      1.9057      1.8311     T = 300, 4800 zeros
    (11,9)    1.3333    0.0745           331.848918      2.0560      1.9815     T = 300, 4800 zeros
    (19,4)    1.3333    0.1460          2344.968149      5.8330      5.6870     T = 300, 9708 zeros
    (19,5)    1.3333    0.1460         -1998.537938      5.6053      5.4593     T = 300, 9708 zeros
    (19,6)    1.3333    0.1460          -523.271176      5.9242      5.7782     T = 300, 9708 zeros
    (19,7)    1.3333    0.1460            74.455951      5.3872      5.2412     T = 300, 9708 zeros
    (19,9)    1.3333    0.1460          -427.489307      6.4145      6.2684     T = 300, 9708 zeros
    (19,11)   1.3333    0.1460           650.059869      5.6595      5.5135     T = 300, 9708 zeros
    (19,16)   1.3333    0.1460         -2154.939497      5.4833      5.3373     T = 300, 9708 zeros
    (19,17)   1.3333    0.1460         -1579.892490      5.8279      5.6819     T = 300, 9708 zeros
    (28,25)   2.6667    0.0914        -16237.118232      2.7831      2.6916     T = 300, 5940 zeros
    (32,9)    2.6667    0.1297         -1750.181289      4.6931      4.5634     T = 300, 8614 zeros
    (32,17)   2.6667    0.1297           -71.838513      4.9794      4.8497     T = 300, 8614 zeros
    (32,25)   2.6667    0.1297         -1124.025263      5.2554      5.1257     T = 300, 8614 zeros
    (36,13)   2.6667    0.0918        -34300.870586      3.0251      2.9333     T = 300, 5968 zeros
    (36,25)   2.6667    0.0918          1083.021998      2.8377      2.7459     T = 300, 5968 zeros
    (39,22)   2.6667    0.3605         -1740.935884      6.4328      6.0723     T = 150, 5810 zeros
    (44,25)   2.6667    0.2983          1206.627875      6.1831      5.8848     T = 150, 4796 zeros
    (45,34)   2.6667    0.3580         -1836.168698      6.4876      6.1297     T = 150, 5762 zeros
    (49,44)   1.3333    0.9806         -1813.588373     12.3625     11.3819     T = 100, 7134 zeros
    (53,28)   1.3333    1.2683         -1231.992855     15.0224     13.7541     T = 100, 9370 zeros
    (56,25)   5.3333    0.5180          -316.983702      6.4128      5.8949     T = 100, 3646 zeros
    (57,16)   2.6667    0.8132          1842.986187     10.7430      9.9297     T = 100, 5844 zeros
    (57,25)   2.6667    0.8132           815.823671     10.7491      9.9359     T = 100, 5844 zeros
    (59,26)   1.3333    1.4350          1880.465557     20.0721     18.6371     T = 100, 10666 zeros
    (64,41)   2.6667    0.7372           400.554890     10.3392      9.6021     T = 100, 5328 zeros
    (68,33)   2.6667    0.7250         -1231.605711     11.7410     11.0159     T = 100, 5212 zeros
    (72,25)   5.3333    0.5198         -1726.774644      6.5661      6.0463     T = 100, 3668 zeros
    (72,49)   5.3333    0.5198          -992.355517      6.9922      6.4724     T = 100, 3668 zeros
    (75,49)   2.6667    0.9078          -633.554354     12.3966     11.4888     T = 100, 6542 zeros
    (80,49)   5.3333    0.7184          1754.869415     10.4442      9.7257     T = 100, 5150 zeros
    (84,25)   5.3333    0.6122          -703.387041      6.8046      6.1924     T = 80, 2694 zeros
    (84,37)   5.3333    0.6122          -108.480467      6.6877      6.0755     T = 80, 2694 zeros
    (87,28)   2.6667    1.6351           682.089103     17.4815     15.8463     T = 80, 7620 zeros
    (88,9)    5.3333    1.1280          -459.868067     13.5246     12.3966     T = 80, 5186 zeros
    (88,25)   5.3333    1.1280          -976.878023     13.0761     11.9482     T = 80, 5186 zeros

Searched without success (best certified lower bound for `sup F`): (7,2): 1.1002 against 1.3333 (|y| ≤ 83336511594824747871295987749137745124807957229143826707037072720628762964425915919205952792204709792456443932180480); (7,4): 1.1063 against 1.3333 (|y| ≤ 60165159375606463455064398619173677520945666495509634445131601329945691598674205644490317487290907174857020752789504); (16,9): 2.2742 against 2.6667 (|y| ≤ 1067663835875696536634630744567372292504908343981177635246257176915911967806802070766814401758794571396673488027648).
For `q = 11, 28, 32` the class 1 was searched as well: the maximum found is at `y ≈ 0` and equals `Φ(q)` up to the remainder (2.9147, 3.8555, 4.8940 against Φ = 2.9848, 3.9407, 5.0116; remainders 0.0745, 0.0914, 0.1297). For `q = 19` a larger value than `Φ(19) = 4.694` occurs at `y_0 = 2826.3` (`F ≥ 5.553`).<!--6Q-END-->

---------------------------------------------------------------------------

## 5. Remarks

**Why not "difference with π" (the route of triage §6).** Writing
`φπ(t;q,a) − li = (π − li) − #{p | q, p ≤ t} + Σ_{χ≠χ_0} χ̄(a)π(t,χ)` and using Johnston's
theorem for the first term needs the stronger inequality `(2/3)(N − 1) > L_1(q) − L_1(1)`
(the bias of the difference is `N − 1`, its zero sum omits ζ). It holds comfortably for
q = 3, 4, 8, 12, 24, barely for q = 5 (`0.553 < 0.667`) and q = 15 (`1.880 < 2`), and it
**fails** for q = 20 (`2.236 > 2`). So that route cannot give `(20,1)`, `(20,9)`; Lemma A
treats `li` directly and keeps the full margin `m = (2/3)N − L_1(q)`.

**Equivalences (not proved here).** For `a = 1`, Landau's theorem applied to the Mellin
transform of `A(x;q,1)` (as in `EXPLICIT.md`, Theorem 6) would give the converse
"`A(x;q,1) < 0` for all large x ⟹ GRH(q)", *provided* `Π_χ L(σ,χ) ≠ 0` for real
`σ ∈ (1/2, 1)` (a real zero there is compatible with one-signedness and cannot be excluded by
the method). For the moduli in question this is a finite verification that I have not
carried out or located in a source I read; the equivalence is therefore not claimed. For
`a ≠ 1` the zeros enter through `Π_χ L(s,χ)^{χ̄(a)}` and can in principle cancel.

**The five knife-edge moduli (why they are out of reach of witnesses).** For q = 7, 9, 14,
16, 18 a GRH-only exclusion needs a point where `F_{q,a}` exceeds 98.2%, 97.8%, 96.8% of its
sup-norm bound `L_1(q)`. The lattice method reaches about 80% (it aligns the 40–100 heaviest
ordinates; the others contribute nothing on average). To reach 98% one must have *all* zeros in
phase up to a height beyond which the remaining weight is below the slack `L_1 − (2/3)N`
(0.025 for q = 7, 0.030 for q = 9, 0.087 for q = 16): the zeros above height 400 alone carry
0.035 (q = 7) resp. 0.049 (q = 16), so several thousand ordinates (height ≳ 1000 for q = 7,
≳ 400–600 for q = 16) would have to be aligned simultaneously — a lattice of dimension in the
thousands with entries of ~10^4 bits, far beyond LLL. The other definite outcome, a certified
upper bound `sup F < (2/3)N` (membership), would contradict LI for these moduli and cannot be
expected; no such bound follows from finitely many zeros without excluding all linear
relations. For `a = 1` the best certified value remains the closed form `Φ(q)` (92–93% of `L_1`
for q = 7, 9; 84% for q = 16). For `a ≠ 1` the lattice method gives the certified lower bounds
`sup F_{7,2} ≥ 1.1002`, `sup F_{7,4} ≥ 1.1063` (threshold 1.3333, `L_1 = 1.3582`) and
`sup F_{16,9} ≥ 2.2742` (threshold 2.6667, `L_1 = 2.7536`), i.e. 81–83% of `L_1`
(`raw/t2_zeros_lll_7_*.jsonl`, `raw/t2_zeros_lll_16_9.jsonl`; zeros to height 400, where already
`(2/3)N + r > L_1 − r` for q = 7: no point whatever can be a witness with these zero lists). These 14 classes are the genuine residue of the problem
"assuming GRH" at the candidate moduli; LI decides them at once.

---------------------------------------------------------------------------

## 6. Gaps and caveats

1. **Nothing here is externally refereed**, and LWS is a Zenodo preprint. The internal referee
   report (`REFEREE-HKMN23.md`) re-derived Lemmas A, E, I, S, N and Theorems 1–5 and reproduced
   with its own code (mpmath/numpy, not ball arithmetic): Lemma A at four points; `L_1(20)`;
   `S(χ)` and `c(χ)` from actual zeros for two characters; B, C, E_R, `U(t)` and `X_0` for four
   pairs; the finite range up to `10^6 ≥ X_0` and `c_1` for five pairs; `Φ(q)` for seven moduli;
   the scan `400 < q ≤ 10^6` of Theorem 5(b). It did **not** re-run: the certificate between
   `10^6` and `10^9`, seven of the twelve finite cases, NV for the 29 468 characters, the
   witnesses of Theorem 6 (two re-evaluated), and it did not check the quotations from the sources. No statement of
   LWS is used in a proof: their list of 18 pairs is re-derived (Lemma S table), and their
   Theorem 9.6 (cut-off q ≤ 60 under LI) is neither used nor re-verified.
2. **Hypotheses.** Theorem 1 for (q,a) uses GRH(q), which includes RH. I do not see how to
   weaken it: the principal character is part of the zero sum (LWS §7.3).
3. **Certified computations rest on software**: python-flint 0.9.0 / Arb for `L(s,χ)` and its
   Taylor coefficients (Hurwitz zeta), `digamma`, `li`, `log`, Hardy `Z`; numpy for the
   integer sieve. Ball arithmetic is rigorous by design, not formally verified; there is no
   second *ball-arithmetic* implementation of the certificates (see item 1 for what was
   reproduced in floating point). Cross-checks actually made: published values of
   `L_1` (LWS, five numbers); the sandwich of Lemma S; the explicit formulas against zeros
   (`verify_formula.out`); Lemma A to 28 digits; `|F| ≤ L_1` in every evaluation; the imaginary
   parts of zero-side sums vanish; `c_1(2,1) = 2.0711926…` equals the constant of
   `EXPLICIT.md` Lemma 4 for q = 4 (same set of primes).
4. **A bug found during this work** (recorded so that it is not repeated): the primitive
   character inducing the Conrey character `q.n` is *not* always `f.(n mod f)` (it fails for
   q = 16, 27, 32, 48, …). `t2_witness.py` now identifies it by its values (validated for all
   q ≤ 130). `w1_analytic.py` (Theorem 2) always asserted value equality and was never
   affected; `l1const.py` and `t2_closed.py` enumerate primitive characters directly.
5. **Theorem 5(b) for q > 400** is a new argument written at the level of these notes: the
   constants are crude and each step is elementary, but it has had no independent check
   beyond the numerical verification of the key identity and the consistency test
   `LB ≤ Φ` for q ≤ 400. If it were wrong, Theorem 5 would still stand for `q ≤ 400`.
6. **The O-terms of Theorem 3** are not made explicit (they are not needed: only "arbitrarily
   large x" is claimed in §1.2). Theorem 4 gives no explicit x at which `A > 0`.
7. **Residual set.** Under GRH alone the following remain undecided: the 14 square classes
   at q = 7, 9, 14, 16, 18 (Theorem 6(a); §5 explains why no witness is to be expected), and
   whatever Theorem 6(c) lists for `61 ≤ q ≤ 120`. Square classes `a ≠ 1` for moduli beyond 120
   are not covered by any general statement here (only by witnesses, modulus by modulus).
8. **Floating point** is used only (i) to *propose* witnesses that are then certified, (ii) for
   printed margins and tables (the inequality `LB(q) > (2/3)N` on `400 < q ≤ 10^6` is decided
   in Arb since the independent audit), (iii) in the purely numerical checks. The step from 3-digit `X_0` to
   the claim "U < 0 on [X_0, 10^12]" uses the monotonicity of U and an Arb evaluation at `X_0`.
9. **Not done:** the converse/equivalence (§5); explicit constants in Theorem 3; a Lean
   formalisation of Lemma A (it is a two-line integration by parts and would be a natural
   small target next to `lean-race`).

---------------------------------------------------------------------------

## 7. Files and reproduction

    l1const.py            S(χ), L_1(q), constants b, b_0, c_{−1}  (Arb, no zeros)      -> l1const.out, l1const.json     (10 s)
    verify_formula.py     numerical check of Lemma E with zeros to height 200          -> verify_formula.out            (5 min)
    w1_cert.py, run_cert.sh   certificate W_1 < 0 for t ≤ 10^9, constant c_1           -> raw/w1_cert_*.json, *.log     (13 core-min)
    w1_analytic.py        Theorem 2: constants, X_0, overlap, consistency              -> w1_analytic.out, .json
    x1.py                 lower-limit constants c_0, x_0, x_1                           -> x1.out
    verify_identity.py    numerical check of Lemma A and of Theorem 1 on a grid        -> verify_identity.out           (3 min)
    t2_closed.py          Φ(q) blocks c_f and NV for conductors ≤ 400                  -> raw/t2_closed_400_*.jsonl
    t2_witness.py         prime-side witnesses, q ≤ 120, X = 10^7                      -> raw/t2_witness_120_*.jsonl
    t2_zeros.py, run_zeros_queue.sh, run_zeros_queue2.sh   zero-side witnesses          -> raw/t2_zeros_*.jsonl, raw/zeros_cache/
                          (certified brackets; the first queue was stopped after q = 56, its jobs for q ≥ 64 were not run)
    t2_refine.py, run_refine.sh   zero brackets refined to width 2^-480 (certified)       -> raw/zeros_hp/
    t2_lll.py             lattice-reduction search + ball certificate at |y_0| ~ 10^50..10^122  -> raw/t2_zeros_lll_*.jsonl
    t2_lll_check.py       independent re-evaluation of the lattice witnesses (mpmath)  -> t2_lll_check.txt
    fill_doc.py           regenerates the result-dependent passages of this file
    AUDIT-INDEPENDENT.md, audit2/   independent audit by an outside session (its own code)
    t2_check_identity.py  numerical check of the identity for T_1                      -> t2_check_identity.out
    t2_assemble.py        tables of §1.2, bound of Theorem 5(b)                        -> t2_summary.out, t2_summary.json
    REFEREE-HKMN23.md     internal referee report (independent agent, own scripts)
    SHA256SUMS, make_sha.sh   sha256 of every output and data file
    PROGRESS.md           log, open ends

Reproduce: `pip install python-flint mpmath numpy`; `python3 l1const.py 128 > l1const.out`;
`./run_cert.sh 1e9 1e6`; `python3 w1_analytic.py`; `python3 x1.py`; then the `t2_*` scripts in
the order above (`t2_closed.py 400 i 4` and `t2_witness.py 3 120 1e7 i 4` for i = 0..3).

---------------------------------------------------------------------------

## 8. Review history

* **Internal referee** (`REFEREE-HKMN23.md`, same day, commissioned by the authoring session;
  own scripts, mpmath/numpy). No fatal error. Repaired: placeholders and Theorem 6(b) as a fixed
  list (R1); the ζ term in the a priori bound, the small terms for `q > 10^6`, `Φ(q)` compared on
  Arb balls (R2); page references and softer wording about LWS (R3); the low zero mod 19 is at
  0.018956, not 0.015 (M1); definition of GRH(q), wording of Corollary 1.4, rigor label, the
  `sup G` column, outward rounding of reconstructed balls, bracket disjointness asserted on
  loading, stale docstrings (M2–M11).
* **Independent audit** (`AUDIT-INDEPENDENT.md`, scripts in `audit2/`, an outside session with
  its own character construction and L-function evaluation). No fatal error; the sources were
  re-fetched and every quotation found accurate; Theorem 1 was reproduced end to end for
  (2,1), (5,4), (12,1), (20,1), (20,9) (independent constants, `X_0`, and a ball-arithmetic
  certificate to ≥ 10^7); NV(q ≤ 400) and all blocks `c_f` were recomputed; the 5 940 brackets of
  the tightest zero-side certificate (28,25) were re-certified. Repaired after it:
  - R1: §1.3 said that *all* non-square classes are excluded under GRH alone; Theorem 4 needs
    NV(q), certified only for `q ≤ 400`. Reworded (§1.3, summary line 9). No theorem changed.
  - M1: the inequality `LB(q) > (2/3)N(q,1)` on `400 < q ≤ 10^6` (Theorem 5(b)) was a
    floating-point scan; it is now decided in Arb for every q in the range (`t2_assemble.py`).
  - M2: Theorem 6(b),(c) now say that five pairs at `q ≤ 60` and seven at `61 ≤ q ≤ 120` are
    excluded through the reduction `(2m,a) = (m,a)` and zero-side witnesses at `q/2`, and list them.
  - M3: the almost-periodicity mechanism of the witness principle is attributed to LWS
    (proof of their Proposition 8.3); what is new is the threshold `(2/3)N` with certified values.
  - M4: "GRH alone" in the comparison with LWS means "GRH plus a certified central
    non-vanishing" on both sides. M5: parity tracked in the bound for `Re (L'/L)(1,χ*)`.
    M6: "new" replaced by "not in LWS or in Johnston".
* **After the audit: lattice reduction** (requested by the coordinating session: "attack the
  knife-edge cases … report only definite outcomes"). New scripts `t2_refine.py`, `t2_lll.py`,
  `t2_lll_check.py`. Definite outcomes: `(28,9)`, `(21,4)`, `(21,16)`, `(40,9)`, `(48,25)`,
  `(60,49)` (hence also `(42,25)`, `(42,37)`) are **excluded** under GRH(q) alone, with certified
  witnesses at `|y_0|` between `10^{53}` and `10^{122}`; Theorem 6(a) went from 21 to 14 undecided
  classes. No definite outcome for the five knife-edge moduli (§5 gives the quantitative reason).
  Also found: the `*.out` result files had never been committed (a repository-wide ignore rule);
  they are now force-added.

