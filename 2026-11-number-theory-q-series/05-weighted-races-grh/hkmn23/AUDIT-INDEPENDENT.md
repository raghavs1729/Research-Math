# Independent audit of `monotonicity/hkmn23/HKMN23.md`

Date: 2026-11-24. First outside pass (the earlier `REFEREE-HKMN23.md` was written by the same
session as the note). Nothing in `monotonicity/hkmn23/` was edited except that this file and the
directory `audit2/` were added. All my code is in `audit2/` and shares no code with the
author's scripts: I wrote my own Dirichlet-character construction (CRT + primitive roots,
`achars.py`), I evaluate `L(s,χ)` through the Hurwitz zeta function (mpmath, or flint's
`acb.zeta(s,a)`), never through flint's `dirichlet_char`, and `(L'/L)(1,χ)` through digamma and
Stieltjes constants `γ_1(k/f)`. CPU: at most two cores, `nice -n 10`. Nothing committed.

**Verdict in one line.** No FATAL error. The mathematics of Theorems 1–5 and of Lemmas A, E, I,
S, N re-derives, and every number I recomputed agrees with the note. I found one REPAIRABLE
overstatement: §1.3 says all non-square classes are excluded "under GRH alone", but Theorem 4
also needs NV(q), which is certified only for q ≤ 400. There are also several MINOR wording and
status points. The note does not claim to have solved Problem 23, and it correctly lists the
residue: five knife-edge moduli for a = 1, square classes a ≠ 1 at the computed moduli, and
everything about square classes a ≠ 1 for q > 120.

---------------------------------------------------------------------------

## 0. Sources (fetched and read; not taken from the note)

* HKMN arXiv:2407.03530v1, Problem #23: the wording quoted in §0 of the note is exact,
  including lower limit 0 and "Assuming GRH, classify those values of q and a such that
  `A_1^π(x;q,a) < 0` for all sufficiently large x".
* Johnston arXiv:2201.06184v2: Theorem 1.1 has lower limit 2; Lemma 2.1 is Büthe's
  `π(x) < li(x)` for `2 ≤ x ≤ 10^19`; §3: "it suffices to consider x > 10^19 in light of Lemma 2.1".
  The note quotes this correctly.
* Liu–Wang–Shi, Zenodo 10.5281/zenodo.22537927 v1.1, file `P1-cpnts-problems.pdf`, md5
  `943c30287804f3218cbaf7879d21ba3d` (the same file the note read). I checked against the text:
  Theorem 1.2 lists exactly the 18 pairs printed in the note. Theorem 8.1 is GRH only and ends
  "which is negative for large x", from `O_q(1/log x)`; I found no explicit threshold anywhere
  (grep). Theorem 8.2 uses LI. Proposition 8.3 assumes GRH(q) + central non-vanishing and
  `A(q,a) ≠ 0`, and its proof uses a zero-counting argument. §8.4 has "67 of the 116 non-square
  classes". §8.5 has "So GRH alone supplies no number here, and the undecided band at a square
  class is the whole range below L_1(q)". §11.1 gives the enclosures of `L_1(7), L_1(9), L_1(16)`
  and margins −1.86%, −2.27%, −3.26%; `L_1(11) = 3.7735558`; (49) is the zeta block 0.0461176444.
  Every quotation and number the note attributes to LWS is accurate.

## 1. Lemma A, the reduction, the lower limit

*Re-derived by hand.* `∫_2^x li = li(2)(x−2) + ∫_2^x (x−s)/log s ds`. With
`ν = φΣ log p δ_p − dt` and `h_x = (x−t)/log t` (so `h_x(p) log p = x − p`):
`∫_{[2,x]} h_x dν = h_x(x)W(x) − h_x'(x)W_1(x) + ∫_2^x W_1 h_x''`, `h_x(x) = 0`,
`h_x'(x) = −1/log x`. My computation of `h_x''` matches the note:
`(x+t)/(t² log² t) + 2(x−t)/(t² log³ t) > 0` on `(1,x]`. The boundary terms are correct. The base
point 2 gives the constant `−li(2)(x−2)`. An atom at t = 2 (only when 2 ≡ a, which no member
has) lies inside `[2,x]` with `W(2−) = 0`. So A(x) ≤ −li(2)(x−2) + W_1(x)/log x once
`W_1 ≤ 0` on `[2,x]`. For q = 1, `π(t;2,1) = π(t) − 1` and `li(2) − 1 > 0`. All correct.

The reductions `(2m,a) = (m,a)` (m, a odd, `2 ≢ a mod m`) are correct: same primes, same φ,
same primitive characters, the same N. They are also used, silently, in Theorem 6 (see M2).

*Lower limit (Corollary 1.3).* I recomputed it by direct quadrature of `li`, not through the
antiderivative (`audit2/ax1.py`, `ax1.out`, 30 digits): `∫_0^2 li = −0.877257534804065308…`,
`x_0 = 2.61113944566016890518…`, `x_1 = 3.56517637359934422071…`, `A_0(3) = 0.2304720`,
`A_0(5−) = −1.6599835`, and `π(p) < li(p)` at every prime ≤ 31. The case analysis in §3 is right:
for q ≥ 2, `A_0 = g` on `(0,p_1]` with `p_1 ≥ 3 > x_0`, and `c_0 − li(2) < 0`; for q = 1, `A_0` is
decreasing on `[2,31]` and `c_0/(li(2)−1) + 2 = 21.4 < 31`. **Correct.**

## 2. Explicit formula, imprimitive characters, Lemma S

*Lemma E.* Re-derived the residues. Odd χ: the residue at 0 is `−bx`, and at the double pole −1 it
is `log x + 1 + b_0`. Even χ: the residue at the double pole 0 is `−x log x + x − bx`, and at −1 it is
`(L'/L)(−1,χ)`. Trivial zeros give `T_e`, `T_o`. ζ gives `x²/2 − x log 2π + (ζ'/ζ)(−1) − T_e`. The
functional-equation forms of b, b_0, c_{−1} also re-derive, with
`−½ψ(1/2) − ½ψ(1) = γ + log 2` and `−½ψ(3/2) = γ/2 + log 2 − 1`.
*Numerics.* I computed b, b_0/c_{−1} by the functional-equation route (L'/L at 1 via Stieltjes
constants, at 2 via Hurwitz derivatives) for 3.2, 4.3, 5.2, 5.4, 12.11, 15.2, 15.8, 15.14 (odd,
even, complex, real). They agree with the author's Taylor-coefficient values in `l1const.json`
to every printed digit (13–15) (`al1.out`).

*Lemma I (imprimitive, primes dividing q).* Correct: `ψ_1(x,χ*) − ψ_1(x,χ) = Σ_{p|q}Σ_k χ*(p^k)…`;
`χ*(p^k) ≠ 0` exactly when the conductor divides `m_p`. The characters mod q of conductor dividing
`m_p` are in bijection with the characters mod `m_p`, and orthogonality gives
`φ(m_p)1_{p^k≡a (m_p)}`. Sign `φψ_1 = Σχ̄(a)ψ_1(χ*) − imp`, `imp ≥ 0`, `imp ≤ I_q x log x`.
`1 + n_e ≥ 0` because `Σ_{χ even}χ̄(a) ≥ 0`. A numerical end-to-end test of formula (E)/(F)
with my constants: the prime-side `F(log x)` must satisfy `|F| ≤ L_1` for all x ≥ 1, and a wrong
constant would violate this at small x. It holds on a 4000-point grid of `[1, 2·10^6]` for (2,1),
(5,4), (12,1), (20,1), (20,9). The maximum `|F|` occurs at x = 1 and equals Φ(q) for a = 1
(0.045971, 0.608550, 2.054172 against `L_1` = 0.046118, 0.632658, 2.281706) (`aw1_*.out`).

*Lemma S.* Re-derived. Davenport (17)–(18), with `Re B(χ) = −Σ Re 1/ρ`, gives V. For ζ,
`ψ(s/2+1) = ψ(s/2) + 2/s`. The identity `(1/π)∫_0^π dθ/(c − cos θ) = (c²−1)^{−1/2}` holds, and
`A + 1 − cos θ = α² + γ²`. For the trapezoid enclosure, the Poisson-kernel coefficients satisfy
`0 ≤ a_k ≤ 2^{1−k}S` (r ≤ 1/2 for c ≥ 5/4), and the n-interval rule aliases exactly `k ∈ 2nℤ`. So
`S ≤ T_n ≤ S(4^n+1)/(4^n−1)`. **Correct.**

*Independent L_1(q)* (`al1.py`, `al1.out`: mpmath, own characters, Gauss–Legendre instead of the
trapezoid, Laurent expansion at s = 1 for ζ, Stieltjes route at s = 1 for χ). Agreement with
`l1const.json`:

    q        L_1 (mine)                         |mine − author|
    1        0.04611764442150902382709084       6e−39
    5        0.5992348065748709826568775        2e−35     (even real character 5.4: S = 0.155171980676144658)
    8        0.7356310181755837921263302        6e−36
    12       0.6326579695105418934017192        2e−35     (imprimitive: ζ, 3.2, 4.3 induced; 12.11 primitive even)
    15       1.925733005290319645975207         4e−35     (imprimitive: ζ, 3.2, 5.x induced; 15.x primitive)
    20       2.281705750772886159371835         2e−35

In every case the closed-form sandwich `V(√5/2)/(√5/2) ≤ S ≤ ½[2V(1/2) + (2/3)V(3/2)]` contains the
value. Also: `L_1(7) = 1.358172`, `L_1(16) = 2.753580`, `L_1(21) = 4.232181`, `L_1(28) = 5.007635`,
`L_1(60) = 7.293766` (`ascan.out`). These are inside the LWS enclosures where LWS print them.

## 3. Theorem 1 / Theorem 2 for selected pairs

I checked five pairs: the tightest pair (20,1); (20,9), with the same modulus and `n_e = −1`; (12,1),
imprimitive with two primes dividing q and `I_q = 4`; (5,4), with a ≠ 1 and `n_o = −2`; and (2,1),
which gives Johnston's case. Script `aw1.py`; outputs `aw1_{q}_{a}.out/.json`.
Constants: my characters and my functional-equation constants. `c_1`: my numpy sieve, with the
maximising point re-evaluated in Arb. U and X_0: my re-implementation of §2.4 as I re-derived it.
Certificate: Arb, 128 bits, `M_i < 0` on every piece.

    (q,a)   m            n_e n_o  B          C          E_R        I_q  c_1 (at)              X_0 (exact root)  U(1e12−)  U(1e12+)  certified W_1<0 up to    pieces   failures
    (20,1)  0.384960916   3   4   1.1433438  0.3024841  9.2809590  6    14.332099224 (379−)   8.1603·10^5       −0.37349  −0.23663  20 000 081               158 711  0
    (20,9)  0.384960916  −1   0  −3.8151759  3.6891204  17.891555  6     9.103522435 (37−)    1.7495·10^5       −0.37765  −0.23625  20 000 069               158 711  0
    (12,1)  2.034008697   1   2  −2.0227008  2.9141512  9.3435267  4     8.978588064 (11−)    157.24            −2.02683  −1.96523  20 000 077               317 441  0
    (5,4)   0.734098527   1  −2   0.4797282 −0.0240460  8.9457776  1     4.095950439 (37−)    688.82            −0.73082  −0.69988  20 000 059               317 584  0
    (2,1)   0.620549022   0   0  −1.8378771  1.9850537  2.3358817  1     2.071192618 (1423−)  58.419            −0.61889  −0.61649  10 000 019               664 578  0

All of these agree with Table 1, the constants table of §2.4, and `w1_analytic.json`. The note's
X_0 values (8.17·10^5, 1.75·10^5, 158, 689, 58.5) are my exact roots rounded up. Since my own
certificate reaches at least 2·10^7, it overlaps the analytic range by a factor of 24 for (20,1)
and more for the others. **So for these five pairs Theorem 1 is reproduced end to end without
the author's 10^9 certificate.** (My U uses `K_B(t)` with the term `2φ log 2 t^{−5/4}`. My
`K_B(10^30) = 85.70` against the note's `K_B(∞) = 84.70`. This is no discrepancy: the
`Y^{−1/6}` term has not yet decayed at 10^30, and the limit is
`2N(L_1 + (2/3)φ log 4) + (4/5)NL_1 = 84.70`.)

*Where the maxima of W_1 are (the continuum between primes).* `W_1' = W = φθ(t;q,a) − (t − 2)`.
It jumps upward by `φ log p` at each prime of the class. So W_1 is continuous, concave on each
piece, and has convex kinks at the primes; a prime is never a local maximum. On a piece
`[p_i, p_{i+1}]`:
* If `t* = 2 + φϑ_i` lies in the piece, the maximum is interior and equals `M_i`.
* If `t* < p_i`, W_1 is decreasing on the piece and the maximum is at the left end.
* If `t* > p_{i+1}`, W_1 is increasing and the maximum is at the right end.

I computed the true maximum of every piece. Interior maxima are rare: 52 of 158 711 pieces for
(20,1), 0 for (20,9) and (12,1), 54 for (5,4), 1 for (2,1). That one case is on `[3,5]`, at
t = 3.10, where the true supremum of `W_1/t^{3/2}` over `(p_1, 10^7]` is −0.0908. A check of W_1
only at the primes would therefore not be a proof in general. The note's test `M_i < 0` (the
maximum over all real t of the concave quadratic) is the right one, and conservative. True
suprema of `W_1/t^{3/2}` over `(p_1, 2·10^7]`: −1.678 (20,1, at t = 30661), −2.153 (20,9),
−1.291 (12,1), −1.188 (5,4). Per-decade bounds `max M_i/p_{i+1}^{3/2}` agree with the
`per_decade_sup` of the author's JSON to four digits (e.g. (20,1): −0.7402 at 101, −1.8059 at
1901, −1.648 at 30661, −2.1895 at 304301, −2.0893 at 4232821).

*Theorem 2 derivation.* Steps 1–5 re-derived. The upper bound for `Θ_1` (y ≥ 3) gives exactly
E_R. The lower bound gives λ(Y). The factor N on the prime-power term is crude but valid, since
the sets `{p^k ≡ b}` for distinct b are disjoint. In
`Sq = 2YΘ_1(Y) − 2∫_0^YΘ_1`, the integral of the upper bound over `[3,Y]` is at most that over
`[0,Y]` because `∫_0^3 y log y > 0`; the part `y ≤ 3` costs at most `φ log 2`. So `K_B` is valid,
with the factor 2 to spare in its last term noted by the internal referee. `r_+` is right, and
every term of `U` is decreasing on `[8,10^12]`. `K_B` is decreasing for `t ≥ e^8`. **Correct.**

## 4. Theorems 3 and 4

*Theorem 3.* The witness principle is correct. `F` is a real absolutely convergent sum
(`Σ|c| = L_1`), hence uniformly almost periodic, and Dirichlet's theorem gives arbitrarily large
positive almost periods, so `F > c + η` on an unbounded set. The error terms re-derive: `W_1`
carries `O(t^{4/3})` from the powers k ≥ 3, `Sq = (2/3)Nt^{3/2} + O(t^{5/4})`, and
`h_x'' ≤ 4.89 x/(t² log² t)`. Together these give
`A(x) log x / x^{3/2} = F(log x) − (2/3)N + O_q(1/log x)`, with an O-constant uniform in x. On
the unbounded set where `F > (2/3)N + η`, this is positive for large x. y_0 may be negative.
*Attribution:* LWS use the same almost-periodicity step at the end of their proof of
Proposition 8.3 (with threshold 0). The note's new contribution is to apply it with threshold
(2/3)N and certified values of F. See M3.

*Theorem 4, the threshold for a non-square class.* `N(q,a) = 0`, so the threshold is 0, and
"A > 0 for arbitrarily large x" follows from one `y_0` with `F(y_0) > 0`. Lemma N is correct: if
`F ≡ c`, then `φψ_1 + imp` would equal `x²/2 + c x^{3/2} + ℓ(x)` for all x ≥ 1. The left side is
piecewise linear; the right side has second derivative → 1. No `A(q,a) ≠ 0` is needed. Under
NV(q), the Bohr mean `−(4/3)Σχ̄(a)m_χ(1/2)` is 0. A non-zero real u.a.p. function of mean 0
takes both signs, so A is positive and negative for arbitrarily large x. **Correct.** This
strictly improves LWS Prop. 8.3 (which needs `A(q,a) ≠ 0`, failing for 49 candidate classes,
e.g. 15: 2, 8, 14).

*NV certificate, re-done completely* (`anv.py`, `anv.out`, `anv_all.out`). I computed
`L(1/2,χ) = f^{−1/2}Σχ(k)ζ(1/2,k/f)` in 96-bit balls with my own characters for **every**
primitive character of conductor `3 ≤ f ≤ 400`: 29 468 characters, every ball excludes 0, and
the minimum `|L(1/2,χ)| ≥ 1.4008·10^{−4}` is at f = 101. Both numbers are identical to the
note's. Sample moduli: 19 (min 0.019140, the Remark's character), 101, 163, 315, 400. The low
zero of the Remark: my Hardy-type function `ε^{−1/2}Λ(1/2+it,χ)` changes sign at
`γ = −0.018956399` for the character with `χ(2) = e(13/18)`, i.e. Conrey 19.3 (2 is a primitive
root mod 19 and `3 ≡ 2^{13}`), and at `+0.018956399` for its conjugate (`alowzero.out`).
Correct as stated.

## 5. Theorem 5

*Closed form (a).* Re-derived. Zeros are symmetric under `ρ ↦ 1 − ρ̄`, and Davenport (17)–(18) at
s = 1, 2 gives `½[ψ(1) − ψ(1/2)] = log 2` and `½[ψ(3/2) − ψ(1)] = 1 − log 2`. Lemma E at x = 1
gives `c(ζ) = 0.045970522563879642`. Counting `φ/2` odd and `φ/2 − 1` even non-principal
characters gives `φ/2 − log 2`.

*Φ(q), recomputed* (`aphi.py`, `aphi.out`: mpmath, own characters, L'/L(1) via Stieltjes
constants):

    q     φ   N   Φ(q)            (2/3)N     Φ − (2/3)N    relative
    7     6   2   1.25422876396   1.33333    −0.0791       −5.9%     (undecided: closed form fails)
    9     6   2   1.27376460969   1.33333    −0.0596       −4.5%     (undecided)
    14    6   2   1.25422876396   1.33333    −0.0791                 (= q 7)
    16    8   4   2.3257354978    2.66667    −0.3409       −12.8%    (undecided)
    18    6   2   1.27376460969   1.33333    −0.0596                 (= q 9)
    21   12   4   3.62949751457   2.66667    +0.9628       +36%
    40   16   8   5.75619743027   5.33333    +0.4229       +7.9%
    48   16   8   5.56747850632   5.33333    +0.2341       +4.4%     (near boundary)
    60   16   8   5.48583306819   5.33333    +0.1525       +2.9%     (nearest to the boundary)
    84   24   8   7.52847887442   5.33333    +2.1951
    120  32  16   12.8489305723   10.6667    +2.1823
    168  48  16   20.172231135    10.6667    +9.5056

All agree with the note and `t2_summary.out`. The ratios Φ/L_1 = 0.923, 0.934, 0.845 at q = 7, 9,
16 and the required 0.982, 0.978, 0.968 are right.

*Full scan q ≤ 400* (`aphi_all.py`, `aphi_all.out`, mpmath, about 20 core-minutes). I computed
every block `c_f` (3 ≤ f ≤ 400) independently. They agree with the author's Arb blocks in
`raw/t2_closed_400_*.jsonl` to `7·10^{−23}`. The set of q ≤ 400 with `Φ(q) ≤ (2/3)N(q,1)` is
exactly {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 15, 16, 18, 20, 24, 30}. That is the 13 member
moduli plus 7, 9, 14, 16, 18, as claimed. The smallest relative excess among the others is
2.86% at q = 60, as claimed. (This scan is floating point at 25 digits; the author's comparison
is in Arb. The two agree.)

*The explicit bound for q > 400, refereed line by line.*
- (i) `T_2`: `Σ_χ χ*(n)` is `φ1_{n≡1}` for `(n,q) = 1` and `φ(m_p)1_{p^k≡1(m_p)}` for `p | q`.
  The bound `n = jq + 1`, `Λ(n)/n² ≤ (log q + log(j+1))/(jq)²` is right, and
  `Σ log(j+1)/j² = 1.8007 ≤ 1.81`. In the imprimitive term, `p^{o} ≥ m_p + 1` and
  `φ(m_p) ≤ m_p`, so each term is `≤ log p/(m_p+2) ≤ (log p)/3`, and `log p/(p²−1) ≤ (log p)/3`
  when `m_p = 1`. **Correct.**
- (ii) The Dirichlet series identity is right: `Σ_{χ≠χ_0}χ(n) = φ1_{n≡1} − 1_{(n,q)=1}`. At
  `s = 1+ε`, `s(s+1)/(2(s−1)) = 1/ε + 3/2 + O(ε)` cancels `ζ'/ζ = −1/ε + γ`. Under GRH
  `∫E_1 x^{−3}` converges, giving `1 + γ + Σ log p/(p−1) + I`. The passage from χ to χ* uses
  `(L'/L)(s,χ) − (L'/L)(s,χ*) = Σ_p Σ_k χ*(p^k) log p p^{−ks}`, with the sign checked.
  **Correct.**
- (iii) The six tail integrals give R exactly as printed. `(1/2)/(1/4+γ²) ≤ 3(3/2)/(9/4+γ²)`
  holds. With the parity tracked: for a = 0, `½ψ(1) = −0.2886` at s = 2 and `−½ψ(1/2) = 0.9818`
  at s = 1; for a = 1, `½ψ(3/2) = 0.0182` and `−½ψ(1) = 0.2886`. This gives
  `Re(L'/L)(1,χ*) ≤ log f + 0.681` (a = 0) and `≤ log f + 0.909` (a = 1). So the note's
  "log f + 0.91" is right *only with parity tracked*: the parity-blind combination
  `0.59 + 0.98` would give `log f + 1.61`. The final bound `|Re(L'/L)(1,χ*)| ≤ log q + 1.62` is
  safe either way. `L_1 ≤ φ(3 log q + 0.11)` holds, with the ζ term treated separately as the note
  now says. `|Re b| ≤ 2 log q + 2.9`, the second constant `≤ log q + 0.84`,
  `|B| ≤ φ(2 log q + 4.4)` and `|C| ≤ φ(log q + 2.4)` (φ ≥ 2) all check. The choice
  `X_1 = 16(L^+)²` gives `−log X_1 − 4L_1X_1^{−1/2} ≥ −2 log(4L^+) − 1`. **Correct.**
- (iv) `log 2π − (ζ'/ζ)(−1) − (ζ'/ζ)(2) + γ = 1.0000000` (functional equation), so the
  constant is exactly 1. **Correct.**
- 400 < q ≤ 10^6 (`alb.py`, `alb.out`, my own sieve, float64): min of `LB − (2/3)N` is
  `18.6056` at q = 420; max R = `3.94·10^{−4}`; max of the `(φ/q²)` term is 0.029. Identical
  to the note. The exceptions to `N(q,1) ≤ φ/4` are exactly the 30 moduli
  1–10, 12, 14, 15, 16, 18, 20, 21, 24, 28, 30, 36, 40, 42, 48, 56, 60, 72, 84, 120, 168, all
  ≤ 168. `g(707) = 208.67`, and g is increasing on `[20,10^5]`. For q = 9 699 690, 38 798 760,
  primorial-type moduli up to 6·10^17 and 1 000 003, `LB/((2/3)N)` is huge and `R < 10^{−7}`.
  **The argument is correct.** Its only non-rigorous step is the floating-point evaluation on
  `(400, 10^6]`. With an absolute margin ≥ 18.6, that step cannot fail, but by the note's own
  status words it is "numerical", not "proved" (M1).

## 6. Theorem 6 spot checks

*Excluded pairs (prime side; `awit.py`, `awit.out`).* My own constants and prime-power sums, 40
digits, at the author's x_0:

    (q,a)    type                                  x_0            F (mine)              author's ball
    (15,2)   non-square, A(q,a)=0 (LWS need LI)   1              0.0273386978614986    0.0273386978614986 ± 1.7e−17
    (15,14)  non-square, A(q,a)=0                 1              0.0162518591075299    0.0162518591075299 ± 9.8e−18
    (16,7)   non-square, smallest margin seen     1              0.00141555265776567   0.00141555265776567 ± 3.4e−19
    (60,7)   non-square                           16.477568976   2.00023472793523      2.00023472793523 ± 1.8e−15
    (13,3)   square, a ≠ 1                        15994.496927   1.58079537669282      1.58079537669282 ± 3.5e−15   (> 4/3)
    (42,5)   non-square                           1              0.160109528423505     agrees
    (7,3)    non-square                           2 600 104      0.436891089793774     agrees

The constants n_e, n_o, B, C agree to all 12 printed digits.

*Zero side, the tightest certificate (28,25)* (`azero.py`, `azero_28_25.out`). I took the
author's 5 940 cached brackets. I matched each cache file to one of my characters, then
re-certified **every** bracket. For each, my real function `ε^{−1/2}Λ(1/2+it,χ)` (Hurwitz route,
root number from the Gauss sum) has opposite certified signs at the two ends, and the brackets of
each character are pairwise disjoint. I recomputed the remainder with my own `S(χ)` and `F_Z(y_0)`
in Arb with my χ̄(a): `L_1(28) = 5.00763547317`, `r = 0.091410642`, `F(y_0) ≥ 2.69164 > 8/3`.
This matches the note (2.6916). The zero-side witness method itself is sound: missed zeros go
into r, and the disjointness is now asserted on loading (`_load`), which resolves the internal
referee's M8.

*Undecided pairs* (`ascan.py`, `ascan.out`). All are square classes with `L_1(q) > (2/3)N`. For
the knife-edge moduli, Φ < (2/3)N. My float scan of the prime-side F on `[1, 2·10^6]` finds no
value near the threshold. These maxima agree with the author's scan maxima where the ranges
overlap:

    (7,2)   L_1 − (2/3)N = +0.0248   max F = 0.58667 (x = 946 285)   [author 0.586674 at 946 283]
    (16,9)  +0.0869                  0.90661 (186.25)               [0.906601 at 186.42]
    (21,4)  +1.5655                  1.38285 (500 671)              [1.38285 at 500 696]
    (28,9)  +2.3410                  1.36737 (1007.85)              [1.36737 at 1008.25]
    (60,49) +1.9604                  1.82768 (4025, range 2e6)      [1.83505 at 9.89e6, range 1e7]

*Counts* (my own enumeration plus the author's raw witness files):
* There are 160 classes at the 23 candidate moduli, 44 of them squares and 116 non-squares.
  All 116 non-squares have a certified prime-side witness in `raw/t2_witness_120_*.jsonl`.
* For 3 ≤ q ≤ 60 there are 1084 pairs outside the members. For 61 ≤ q ≤ 120 there are 3284.
* For q ≤ 120, 4266 of the 4368 pairs have prime-side witnesses (2767 of the 2788 non-squares
  and 1499 of the 1580 squares). The 21 non-squares without a witness are exactly those named in
  §4.5.
* Applying Theorem 4, Theorem 5, the prime-side witnesses and the zero-side witnesses, and then
  the reduction (2m,a) = (m,a), leaves 22 undecided pairs for q ≤ 60 and 29 for 61 ≤ q ≤ 120.
  These are exactly the note's lists.
* Without the reduction there would be 5 + 7 more: (22,3), (22,5), (22,15), (38,25), (38,35),
  and (78,61), (90,79), (98,93), (106,81), (114,25), (114,73), (118,85). See M2.

## 7. Wording, novelty, status

* **No claim to have solved Problem 23.** The note claims "for every x" at the 18 pairs, and
  under GRH alone it claims a classification of a = 1 except at five moduli. It names an
  explicit undecided list (21 + 1 classes at q ≤ 60, 29 at 61–120). §6.7 says square classes
  a ≠ 1 beyond the computed range are not covered. That is accurate.
* **Every x versus eventually.** This is correctly attributed: LWS give no threshold, and Johnston
  covers q = 1. "New" for the other cases is a priority claim I cannot fully verify (I know of no
  such result). It rests on the LWS text, which I checked.
* **Status words.** "proved", "certified", "numerical" and "cited" are used accurately throughout,
  with two exceptions: R1 below, and M1 (the float step in Theorem 5(b)). The Φ assembly for
  q ≤ 400 is now done in Arb (`t2_assemble.py`, `Phi_arb`), which resolves the internal
  referee's R2(c).

---------------------------------------------------------------------------

## Findings

### FATAL

None.

### REPAIRABLE

**R1. §1.3, fourth bullet ("GRH-only exclusions"): "Here: all non-square classes
(Theorem 4), the class 1 for every modulus outside five (Theorem 5), and the explicit finite list
of Theorem 6 are excluded under GRH alone."** This overstates Theorem 4. Theorem 4 assumes GRH(q)
**and** NV(q), and NV(q) is certified only for q ≤ 400. For q > 400 nothing in the note excludes
a non-square class under GRH alone. The note's own Theorem 4, summary item 9 and §6 are
correctly worded; only this sentence is not. The same unqualified "under GRH alone" also appears
implicitly in summary item 12 / Theorem 6(b) ("every non-square class (… or Theorem 4)"). That
use is fine, since q ≤ 60 ≤ 400. *Fix:* "all non-square classes at every modulus q ≤ 400
(Theorem 4 with the certified NV(q)); at q > 400, all non-square classes under GRH(q) + NV(q)".

### MINOR

**M1. Theorem 5(b), status of 400 < q ≤ 10^6.** The bound LB(q) is evaluated in float64
(`t2_assemble.py`). The margin is ≥ 18.6, so the result is not in doubt (reproduced), but by §0's
definitions this range is "numerical", not "proved". Either run the scan in Arb or integers
(seconds), or label it. §6.8 already discloses this; the theorem label should match.

**M2. Theorem 6(c) and summary item 12, "(prime side only)".** Seven of the pairs counted as
excluded at 61 ≤ q ≤ 120 have no prime-side witness: (78,61), (90,79), (98,93), (106,81),
(114,25), (114,73), (118,85). They are excluded only through the reduction (2m,a) = (m,a) of §0(i)
and the zero-side witnesses at q = 39, 45, 49, 53, 57, 59. The same holds at q ≤ 60 for (22,3),
(22,5), (22,15), (38,25), (38,35). The conclusion is correct, but Theorem 6(b),(c) should say
that the 2m → m reduction and zero-side witnesses at q/2 are used. §0(i) is stated only for the
18 pairs; it holds whenever `2 ≢ a (mod m)`.

**M3. Attribution of the witness principle.** The almost-periodicity step "F(y_0) > δ ⇒ F > δ/2
at arbitrarily large y" is already the last step of LWS's proof of Prop. 8.3 (threshold 0). The
note's §1.3 presents "one certified value of F above (2/3)N is such a number" as an observation
against LWS §8.5. It should say that the mechanism is theirs and that what is new is the
threshold (2/3)N together with certified evaluations (closed form at y = 0, prime side, zero side).

**M4. "LWS under GRH alone: … 67 excluded".** LWS's Prop. 8.3 assumes central non-vanishing as
well. They discharge it computationally at the candidate moduli, as the note does with NV(q ≤ 400),
so the comparison is fair but should name the convention.

**M5. §4.4 (iii), "Re(L'/L)(1,χ*) ≤ log f + 0.91".** This is true only when the digamma terms are
taken with the character's parity (the a = 1 case gives 0.909, a = 0 gives 0.681). A reader
combining the worst cases of the two digamma terms gets 1.61. One clause would fix it. It has no
effect, because 1.62 is used afterwards.

**M6. §1.3 "the other 9 distinct cases … are new".** This is a priority claim for which the
only checked evidence is the LWS text (no threshold) and Johnston (q = 1). Phrase it as "not in
LWS or Johnston".

---------------------------------------------------------------------------

## What is established, and under which hypotheses

1. **Theorem 1 (proved; finite part certified)**, under GRH(q) (which includes RH): for each of
   the 18 LWS pairs, `∫_2^x (φ(q)π(t;q,a) − li(t)) dt < 0` for **every** x > 2, with the explicit
   bound of Theorem 2. For (2,1), (5,4), (12,1), (20,1) and (20,9) I reproduced this end to end:
   independent constants, independent X_0, and an independent ball-arithmetic certificate of
   `W_1 < 0` to ≥ 10^7, which overlaps X_0 by a factor ≥ 24. With the HKMN lower limit 0 the
   integral is negative exactly for `x > 2.6111394457` (q ≥ 2) and `x > 3.5651763736` (q = 1).
2. **Theorem 3 (proved), GRH(q):** one certified value `F_{q,a}(y_0) > (2/3)N(q,a)` implies
   `A > 0` for arbitrarily large x.
3. **Theorem 4 (proved), GRH(q) + NV(q):** every non-square class changes sign infinitely often.
   NV(q) is certified for all q ≤ 400 (independently re-verified, 29 468 characters). So the
   exclusion of non-square classes is under GRH alone **for q ≤ 400 only** (corrected statement,
   R1).
4. **Theorem 5 (proved; q ≤ 400 certified, 400 < q ≤ 10^6 by a float scan with margin ≥ 18.6,
   q > 10^6 analytic), GRH(q):** (q,1) is in the answer iff q ∈ {1,2,3,4,5,6,8,10,12,15,20,24,30};
   it is undecided under GRH alone at q ∈ {7,9,14,16,18}; and it is excluded for every other q.
5. **Theorem 6 (certified witnesses), GRH(q):** at the 23 candidate moduli there are 18 members,
   121 excluded and 21 undecided square classes. For 3 ≤ q ≤ 60 there are 1062 excluded and 22
   undecided pairs. For 61 ≤ q ≤ 120, all pairs are excluded except 29 square classes a ≠ 1. Both
   ranges use, besides prime-side witnesses, Theorem 4, Theorem 5, zero-side witnesses and the
   2m → m reduction (M2).
6. **Not established:** Problem 23 under GRH alone is **not** solved. Open are the classes
   at the five knife-edge moduli, the 7 square classes a ≠ 1 at q = 21, 40, 42, 48, 60, the class
   (28,9), the 29 classes at 61–120, all square classes a ≠ 1 at q > 120, and (without NV) the
   non-square classes at q > 400. The complete classification still holds only under GRH + LI +
   central non-vanishing (LWS).

## Files

`audit2/achars.py` (characters), `al1.py`/`al1.out` (L_1, constants), `aw1.py` + `aw1_*.out/.json`
(Theorems 1–2), `anv.py`/`anv.out`/`anv_all.out` (NV), `alowzero.py`/`.out`, `aphi.py`/`aphi.out`,
`aphi_all.py`/`aphi_all.out` (Φ), `alb.py`/`alb.out` (q > 400 bound), `awit.py`/`awit.out`,
`azero.py`/`azero_28_25.out` (Theorem 6), `ascan.py`/`ascan.out` (undecided), `ax1.py`/`ax1.out`
(lower limit).
