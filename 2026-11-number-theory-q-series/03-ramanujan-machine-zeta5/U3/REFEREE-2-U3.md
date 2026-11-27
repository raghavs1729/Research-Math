# Referee report 2 on U3/THEOREM.md (independent, adversarial)

Date: 2026-11-26. Object: `U3/THEOREM.md` as of commit bba4603. The numbered results of `U2/THEOREM.md` were taken as
given after I read their statements (U2 was refereed separately, `U2/REFEREE-2-U2.md`). `U3/AUDIT.md` was not read
before the verdict below was fixed. My scripts are in
`/tmp/claude-0/-home-user-weighted-prime-races/7b9224b4-6348-502d-a364-b106f526c180/scratchpad/u3ref/`
(`ex1.py`, `ser.py`, `b5.py`, `cf2.py`, `rhs.py`, `keyc.py`, `gg.py`, `bil.py`, `tay.py`, `rem2b.py`). None of the
document's scripts were used to produce my numbers; I wrote the local series, the Frobenius basis, the connection
matrix and the quadratures myself.

## Verdict

**PROOF COMPLETE.** I found no fatal error and no gap. Every step is either derived in the text or reduces to a
finite exact identity, and I re-derived each of those identities independently (by hand or with my own exact code).
No numerically observed identity is used in the proof. The findings below are minor: wording, a notation clash,
one script assertion that checks less than the text claims, and some interpretive remarks that are stated more
firmly than they are proved. None of them affects the theorem.

## 1. Reading of the statement (numerical)

* **Continued fraction.** I computed exact integer convergents `p_n, q_n` (`p_n = a_n p_{n-1} − 64 n^10 p_{n-2}`,
  `(p_{−1},p_0) = (1,23)`, `(q_{−1},q_0) = (0,1)`) up to n = 6000. I stored `q_n/p_n` every 200 steps and fitted
  `ℓ + Σ_{j≥1} c_j n^{−j}` with 8, 12, 16 and 20 nodes. Results 12, 16 and 20 agree to 45 digits:
  `1/PCF = 0.043780903856351970838937321518603279588971872…`,
  `PCF = 22.841008565767985615909667462437328325…`. The raw error behaves like `n^{−6}`, not `c/n`: the tail is
  `Σ 1/(8(m+1)^5 w_m w_{m+1})` and `w^{(3)}` grows.
* **Right-hand side** (mpmath `nsum`, 50 digits; `S_0` and `Z_0` also via `hyper`, difference 0):
  `S_0 = 1.11863638716418706834961925…`, `S_1 = −0.70900136746192943753111978…`,
  `Z_0 = 1.04335518933431988860760740…`, `K = −7.78349425139529619…`, `1 + S_1K = 6.51850806787133178…`.
  The formula minus the fraction is `9.7·10^{−48}`, which is at the precision of the extrapolation.
* **Sign and factor tests.** Replacing `−4` by `+4` changes the value by exactly 8. Replacing `1 + S_1K` by `1 − S_1K` is
  off by 49. Both alternatives fail, as they should.

The statement as written in §0 (all three displayed forms) is what is proved. I checked the algebra of the second
and third forms by hand, given `S_0 = 16L(f,2)/π^2`, `S_1 = −2L(f,1)`; that L-value input comes from `../Z`, which
is outside this referee's scope.

## 2. Line-by-line

### §1 (notation, (1.2))
Correct. `A_{k+1}(u) = κ(u)A_k(u+1)` gives both relations in (1.2): the `k=0` term gives `1`, resp. `1/(u+1)`. `Z` is
analytic for `Re u > −1/2` (the terms are `O(k^{−3})` locally uniformly).

### §2, Lemma 2.1
Correct. I confirmed `ã_3(n) = a_n/8` and `ã_3(n) − (n+1)^5 − n^5 = 5(2n+1)(4n^2+4n+3)/8` (sympy). The induction
and the Casoratian argument are complete. `v^{(3)}_1/w^{(3)}_1 = 1/23 > 0`. The identification
`p_{n−1} = 8^n n!^5 w^{(3)}_n` is correct: `p_{−1} = 1 = w_0`, `p_0 = a_0 = 23 = 8·23/8`, `q_0 = 1 = 8 v_1`, and the
substitution `U_{n−1} = 8^n n!^5 c_n` turns `U_n = a_nU_{n−1} − 64n^{10}U_{n−2}` into `R_3`.

### §2, Lemma 2.2 (gauge map)
Correct, and proved for all n.
* **Operator identity.** I derived it independently as an explicit certificate. For arbitrary `X_{n−1..n+2}`:

      R_3(G X)_n = P_1(n−1) · R_2(X)_n + 32 (n+1)^5 · R_2(X)_{n+1},

  where `R_j(X)_n := (n+1)^5X_{n+1} + n^5X_{n−1} − ã_j(n)X_n`. Here `P_1(n−1) = −(2n−1)(16n^4−8n^3+4n^2−2n+1)`. The
  remainder is exactly 0 (`ex1.py`). So `G` maps solutions of `R_2` to solutions of `R_3` at every n where `R_2` holds
  at n and n+1.
* **Initial values.** I checked them by hand. `G(w)_0 = −11 + 32·(3/8) = 1`. `G(w)_1 = −615·(3/8) + 1024·(467/2048) = 23/8`.
  `G(v)_0 = 32/8 = 4 = v^{(3)}_0 + 4w^{(3)}_0`. `G(v)_1 = −615/8 + 1024·(177/2048) = 93/8 = 1/8 + 4·23/8`.
  So `G(w) = w^{(3)}` and `G(v) = v^{(3)} + 4w^{(3)}` hold for all n ≥ 0.
* **Injectivity** is not needed. `G` is only applied forwards, to the identity (3.1), and linearity of `G` is all
  that Prop. 4.1 uses.
* So the `−4` (it comes from `G(v) = v^{(3)} + 4w^{(3)}`) and the factor 192 (it comes from (3.1) with `v_1 = 1/8`,
  `g_*(1) = ã(0)g_*(0) − 24μ_4^*`) are derived, not fitted.

### §3, the solution r_*, Prop. 3.1
Correct.
* I recomputed the local series with my own code in the variable s. The pivot structure is triangular, with lowest
  degree m−3 in `L4(s^m)`. Results: `r_* = 1 − s^2/16 − (7/96)s^3 − (169/2304)s^4 …` and
  `μ_* = 1 + s/2 + (3/8)s^2 + (5/16)s^3 + (839/3072)s^4 …`.
* `(θ^iμ_*(1))_{i=0..4} = (1, −1/2, 1/4, −1/8, 7/128)`, and `(θ^i r_*(1))_{i=0..2} = (1, 0, −1/8)`.
* Lemma 3 of U2 gives the right side `F_2 + F_1/2 − u(u+1/2)F_0 = −1/8 − u(u+1/2)`. The identity for `R_*` follows from
  `(u+1/2)^4 − 1/16 = u(u+1)(u^2+u+1/2)`, which I checked by hand.
* The periodicity and decay argument is complete: `|r_* − 1| ≤ Cs^2`, `0 < Z ≤ Y ≤ u+2` for `u > 0`, and `1/Γ-ratio ~ u^2`.
* The Frobenius coordinates follow from U2 §5, since `Y_*` is analytic at 0.

### §3, Prop. 3.2
Correct.
* **(i)** I recomputed the boundary polynomial myself and got `B_5(u) = (u + 1/2)(M_0/2 + M_1)` with
  `M_i = θ^iμ(1)`. The coefficients of `M_2, M_3, M_4` cancel identically. For `μ_ρ = νθρ − ρθν`: `M_0 = ρ(1)` and
  `M_1 = −ρ(1)/2`, using `ν(1) = 0`, `θν(1) = −1`, `θ^2ν(1) = 1/2`. Hence `B_5 ≡ 0` for every ρ holomorphic at 1, as
  claimed. `L5μ_* = 0` follows from U2 Lemma 1, and the boundary terms at 0 vanish for `u > 0`.
* **(ii)** The limit `u → 0+` is exactly as in U2: `μ_* = μ_4^* log^4 x + O(1+|log x|^3)` near 0, and `μ_*` is bounded on
  `[1/2,1]`. Then (3.1) follows by uniqueness for a second-order recurrence with `n^5 ≠ 0` for `n ≥ 1`.
* **(iii)** The log^4 coefficient `(c_2(ν)c_3(ρ) − c_3(ν)c_2(ρ))/12` with `c_j(r_*) = 2(−1)^jη*_{3−j}` gives
  `(N_0η*_1 − N_1η*_0)/6 = −(1/6)[ε^1](N(ε)Y_*(−ε))`, which I checked by hand. `g_*(0) = Σ c_j(r_*)ξ_j` uses only U2's
  integral identity for `ν` and `y_j`. Correct.

### §3, Lemma 3.3 (asymptotics)
Correct. The error terms are proved by the elementary bound `∫|m|x^n ≤ Cβ_j(n) + O(2^{−n})`; no Watson lemma is
needed or used.
* **(i)** Correct, from `μ = −s^2 − (3/2)s^3 + …` and `μ_* = 1 + s/2 + …`.
* **(ii)** I computed the exact rational function. It is

      Σ_{j≤3} m_j G(β_j)_n = (6712n^4 + 23412n^3 + 31402n^2 + 18297n + 4015) / (8(n+1)(n+2)(n+3)(n+4)(n+5)) ~ 839/(8n),

  so the polynomial part is indeed `O(1/n)`. The `ρ_4` part is `O(c(n)β_4 + n^5β_5) = O(1/n)`. The claim needs the
  exact value `m_3 = 5/16`: with general `m_3` the limit is `180 − 576 m_3`. I reproduced `m_3 = 5/16` exactly (see
  Finding 2).

### §4, Prop. 4.1
Correct. `G(g_*)_n = (g_*(0) − 768μ_4^*)w^{(3)}_n − 192μ_4^* v^{(3)}_n` follows from linearity, (3.1) and Lemma 2.2.
Dividing by `w^{(3)}_n ≥ 1` and using `G(g_*)_n → 0` gives the limit. This is a derivation of the limit of
`p_n/q_n`, not a fit.

### §5, Prop. 5.1 (Legendre-type relation)
Correct. I re-derived the bilinear computation: `C[g,g_*] = −192(g(0)μ_4^* − μ_4g_*(0))C[w,v]`, with
`C[w,v]_n = 1/(8(n+1)^5)`. The asymptotic side, `(−2/n^3)(−1/n^2) − (6/n^4)(1/n) + O(n^{−6}) = −4/n^5 + O(n^{−6})`,
uses only Lemma 3.3(i), whose error terms are proved. Then `−24X = −4`, so `X = 1/6`. The value 1/6 is derived, not
observed. The asymptotics come from the integral representation `g_*(n) = ∫_0^1 μ_* x^n dx` (the definition of
`g_*`), not from any unproved representation.

### §6, Prop. 6.1 (bilateral sum)
Correct.
* **Reflection.** `γ(z) = tan^4(πz)γ(−z−1/2)` is correct.
* **Poles.** Each term `γ(u+k)/(u+k+1)` has poles only at `u ∈ 1/2 + Z`: the zero of order 4 of `1/Γ(z+1)^4` at
  `z = −1` absorbs `1/(z+1)`. In a period strip the only pole is `u = 1/2`, of order ≤ 4.
* **Decomposition (6.2).** I checked it by hand: `m = j+1`, `u − m + 1 = −(v + j − 1/2)`.
* **Expansion at u = 0.** Only the `j = 0` term contributes at order `u^3`, giving `16π^2u^3` (`N(1/2) = π^2/16`).
* **Growth.** The bound is uniform in `Re u ∈ [0,1]`, `|Im u| ≥ 1`: Stirling for `k ≥ 0`, the reflected form with
  `|tan| ≤ coth π` for `k < 0`, and `Σ_k |u+k|^{−3} = O(|Im u|^{−2})`. This is the same argument accepted in U2 Prop. 8.
* **Liouville step.** `sec^4`, `sec^2 tan`, `sec^2` and `tan` have triangular principal parts `v^{−4..−1}`. After
  subtracting them, the remainder is entire, periodic and bounded, hence constant. Its limits `∓iπc` at `±i∞` force
  `c = 0` and constant 0. The span of `sec^4, sec^2 tan, sec^2` is exactly `{(1+t^2)(e_0+e_1t+e_2t^2)}`.
* **(6.1).** I expanded it by hand and it is correct.
* **e_2.** `e_2 = −π^2 Σ A_j/(j−1/2)` holds by uniform convergence (terms `O(j^{−3})` near `v = 0`).
* **Telescoping certificate** (verified exactly, `ex1.py`):

      H(k) = A_k k^3/(k − 1/2):   H(k+1) − H(k) = A_k[(k+1/2)^3/(k+1) − k^3/(k−1/2)] = −(1/8)A_k(1/(k+1) + 1/(k−1/2)).

  With `H(0) = 0` and `H(N) = A_N N^3/(N−1/2) → 1/π^2` (`A_N ~ 1/(π^2N^2)`, standard), summing gives
  `Z_0 + Σ A_k/(k−1/2) = −8/π^2`. That is equivalent to `e_2 = e_0 + 8` and to the stated by-product
  `Σ A_k/(k+1) = 2Σ A_k/(1−2k) − 8/π^2`.

### §7 (end of proof)
Correct. I redid the coefficient extraction by hand, independently of `u3_algebra.py`.
* The `ε^3` coefficient of `ζ(3,1/2+ε)(1 − 4π^2ε^2/3)` is 0.
* The `ζ(3)` terms cancel exactly because `ω_3 − (4π^2/3)ω_1 = −π^2/2` (`ω_i = −τ_i/16 + ζ_i/32`), and this cancels
  the `−7ζ(3)/π^2` from the `u^3/2` part of `Y_*`.
* `g_*(0) = 2ω_2 + 744ζ(5)ω_1/π^4 = (93/4)ζ(5)E_1 + (π^2/16)(e_0+e_2−4τ_0)`, and `192μ_4^* = 32ω_1/π^4 = E_1`.
* The Legendre substitution gives (7.2) and the final formula.
* The signs are correct: `2S_0 − Z_0 ≥ 1 > 4/π^2`, so `K < 0`, so `S_1K > 0`, so `E_1 > 0`.
* The by-product `e_1 = 8(1+4S_1+π^2S_1Z_0)/(πS_0)` follows from (7.2) and `τ_1 = 8π^2S_1`; I checked it by hand.

### (d) Interchanges and convergence
Every series and interchange is justified: `Z`, `Y`, `Z~` converge locally uniformly with `O(k^{−2})` or `O(k^{−3})`
terms, the bilateral sum converges, and Taylor coefficients follow by Weierstrass/Cauchy. The Fubini step for
`g_*(0) = Σ c_jξ_j` is just linearity of a finite sum plus U2's justified identity. No gap.

### (e) Hidden numerical inputs
None. Every "numerically" parenthesis in §3–§6 is a cross-check only. The exact machine checks the proof relies on
(v6 symbolic, u3_exact (a)–(f), u3_algebra) are all finite identities, and I reproduced each one independently.

## 3. Findings

No fatal findings, no gaps. Minor:

1. **Notation clash.** In Prop. 3.1 the proof uses `G(u) = Γ(u)^4/Γ(u+1/2)^4`, the same letter as the gauge map `G` of
   Lemma 2.2, which is used again in §3.3 and §4. Rename one of them.
2. **Lemma 3.3(ii): text and script.**
   * The parenthesis says "the terms of order n^2 and n cancel exactly". For `O(1/n)` the constant term must cancel
     as well; it does, but only because `m_3 = 5/16` (the limit is `180 − 576 m_3`). The text should say
     "orders `n^2`, `n` and `1`".
   * `u3_exact.py` (d) asserts only `deg num ≤ deg den`, which is `O(1)`. That is weaker than the `O(1/n)` the text
     claims, and `O(1)` would not suffice for Prop. 4.1, since only `w^{(3)}_n ≥ 1` is proved. The script prints the
     limit (0), so the information is there, but the assertion should be `dn < dd`. The true exact rational function
     is given in §2 above; it is `~839/(8n)`.
3. **Lemma 2.2 certificate.** The proof says only "polynomial identity … symbolic". The explicit certificate
   `R_3(GX)_n = P_1(n−1)R_2(X)_n + 32(n+1)^5R_2(X)_{n+1}` is short and makes the lemma checkable by hand. Suggest adding it.
4. **Remark 3 overstates.** It says `G(g_{r_*})` "is the minimal solution of the U3 recurrence". Minimality is not
   proved: no growth lower bound for a second solution beyond `w^{(3)} ≥ 1` is established, and none is needed.
   Numerically, `G(g_*)_n` decays much faster than the proved `O(1/n)`, roughly like `n^{−5}` for n ≤ 60, and
   `G(g)_n/n → 128` is confirmed. Rephrase as "it is `O(1/n)`, while `G(g)_n ~ 128n`".
5. **Interpretive claims about Z_0.** §0 says "`Z_0` plays the role of a quasi-period of `f`", and Remark 2 says
   "(7.2) is the corresponding Legendre relation between periods and quasi-periods". Neither is proved; they are
   heuristic interpretations. Mark them as such in §0 as well, not only in Remark 2. The disclaimer "Nothing is
   claimed or proved about its transcendence or independence" is present and adequate. The integral formula of
   Remark 2 is labelled "numerically checked". I confirmed it: the quadrature, truncated at `10^{−14}` from both
   endpoints, agrees to `1.5·10^{−13}`, consistent with the truncation.
6. **Status paragraph cites AUDIT.md.** It quotes another auditor's agreement (`7.7·10^{−35}`). That is fine as a
   record but is not part of the proof.
7. **Convergence rate.** Not an error in the document, only a note: `q_n/p_n − ℓ_3 = O(n^{−6})` empirically, much
   faster than `c/n`.

## 4. What I recomputed (all independent of the document's scripts)

* **Fraction.** 45 digits from exact convergents to n = 6000 with Richardson extrapolation. Agreement with the
  closed form: `9.7·10^{−48}`.
* **Key constants, by a method different from the document's.**
  * Local series at x = 1 (exact Fractions and 40-digit floats), and a Frobenius basis at x = 0 built from
    ε-power series of `A_k(ε)`.
  * Connection coefficients obtained numerically by matching values and three θ-derivatives at x = 1/2, rather than
    by using Prop. 3.1. Consistency at x = 0.45 and 0.55: about 1e-40.
  * `μ_4, μ_4^*` from the log^4 coefficient, and `g(0), g_*(0), g(n), g_*(n)` by direct tanh-sinh quadrature of
    `μ x^n`. The values:

        μ_4    = −0.09578247700703834186578197919229825377
        μ_4^*  =  0.00249257671064393388038748052867741286
        g(0)   = −7.50087288837386745900494753468994137063
        g_*(0) =  1.93525134794861177029108768866756514715

* **Relation tests at 39–41 digits.**
  * `g(0)μ_4^* − μ_4 g_*(0) = 0.1666…6667`, i.e. 1/6 to 40 digits.
  * `μ_4 − 4S_1/(3π^2) ≈ 2e-41`, `g(0)/(192μ_4) − 1/U2 ≈ 1e-39`, `192μ_4^* − E_1 ≈ 5e-41`.
  * `g_*(0) − [(93/4)ζ(5)E_1 + π^2K/8] ≈ 8e-40`.
  * `g_*(0)/(192μ_4^*) − 4 − (fraction) ≈ 1.6e-39`.
  * (3.1) for g and g_* at n = 1, 2, 3, 5, with g(n) and g_*(n) from quadrature: residuals of about 1e-39 to 1e-40.
  * Asymptotics: `n·g_*(n) → 1`, `n^3 g(n) → −2`, `n^5 C[g,g_*]_n → −4`, `G(g)_n/n → 128`, and `G(g_*)_n` decays.
* **Bilateral formula** at the complex points 0.23+0.31i and −0.17+0.05i, and at the real point 0.3: residual
  about 1e-29 to 1e-30. Tests designed to fail do fail: `e_2 = e_0 − 8` gives a residual of order 5–90, and flipping the
  sign of `e_1` gives a residual of order 20–170.
* **Telescoping identity:** `Z_0 − 2Σ A_k/(1−2k) + 8/π^2 = −1.3·10^{−51}`.
* **Taylor coefficients of Z/N at 0.** I computed them with a 48-point trapezoidal Cauchy integral on |u| = 0.15
  (`tay.py`), with Z summed directly. The four relations (6.1) hold to 6e-21, 1e-20, 3e-20 and 6e-20, at the
  precision of the run. As a test designed to fail, replacing `−16π^2` by `+16π^2` in `ζ_3` is off by 316.

## 5. What I could not check

* `../Z/THEOREM.md` (`S_0 = 16L(f,2)/π^2`, `S_1 = −2L(f,1)`), which is needed only for the L-value forms in §0, and
  the U2 results themselves, which were taken as given per instructions.
* The literature claim of novelty ("not found after that search").
* The PSLQ claims about `Z_0`. They are non-mathematical and disclaimed; I did not rerun them.
