# Second referee report on TAU-RACE.md (independent, adversarial)

Referee: an AI agent (Claude), 2026-11-26. I read TAU-RACE.md in full (lines 1-658), RESULTS.md, the
programs `cert_constants.py`, `cert_small.py`, `bounds.py`, `sym2_phi.py`, `sym2_verify.py`,
`tau_table.py`, `arbutil.py`, and the data files `constants.json` and `bounds.txt`. I did **not** read
REFEREE-REPORT.md or AUDIT-NUMERICS.md, before or after reaching this verdict. My own scripts are in
`/tmp/claude-0/-home-user-weighted-prime-races/7b9224b4-6348-502d-a364-b106f526c180/scratchpad/tauref/`.
No human has checked anything here.

## Verdict

| | verdict |
|---|---|
| **Theorem 1** (a), (b) | **COMPLETE MODULO LISTED GAPS.** The gaps are G1-G3. None is a known error. All are citation or presentation gaps that a standard reference closes. |
| **Theorem 2** | **COMPLETE MODULO LISTED GAPS** (the same ones, inherited through Theorems 1 and 4). The deduction from Theorem 1 itself is correct. |
| **Theorem 3** | **COMPLETE MODULO LISTED GAPS.** (c)⟹(a) and (e)⟹(a) (Landau) are correct as written and need neither Theorem 4 nor G1. (a)⟹(b) and (a)⟹(d) inherit Theorems 1 and 2. |
| **Theorem 4** | **COMPLETE MODULO LISTED GAPS** as a computer-assisted proof. The method is rigorous: certified sign changes plus a completeness test from an explicit-formula positivity argument (Lemma V3), with no Turing-method or argument-principle count needed. It rests on (i) the functional equation of `Λ_S` with exactly the stated `γ_S` and the Hadamard product (F3, F4: cited, standard, not consulted by the authors or by me); (ii) FLINT/Arb; (iii) the authors' two programs. I reviewed both programs line by line and reran them: the output is byte-identical. |

I found **no fatal error**. I found **no wrong step** in the mathematics. The three repairs that were never
re-refereed (Lemma 0, the "first term + integral" n-tail in `sym2_phi.py`, and the gap rule in `cert_small.py`) are correct (findings 4-6).

## Findings

### Fatal
None.

### Gaps

**G1. (F3)/(F4) are load-bearing and were taken on trust. They are cited correctly, but no one has read them** (§0 lines 39-52; used in Lemmas 0, 2-5, V0-V3).
Theorem 4, Lemma 3 (`a_{S,0} ≥ 0`, `K_S`), Lemma 4 and Lemma 5 all use:
- the functional equation `Λ_S(s) = Λ_S(1−s)`, with `γ_S = Γ_R(s+1)Γ_C(s+11)`, conductor 1, sign +1;
- the fact that `Λ_S` is entire of order 1;
- the Hadamard identity `Re Λ_S'/Λ_S(s) = Σ_ρ Re 1/(s−ρ)`.

Theorem 4 also uses the functional equation directly: the table only stores `k ≥ 0`, and `Φ(e^{−kh}) = e^{kh}Φ(e^{kh})` in Lemma V1 is the functional equation.

The citations are the right ones: Shimura, Proc. LMS 31 (1975); Gelbart–Jacquet, Ann. ENS 11 (1978). The gamma factor is the standard one; in the duplicated form `π^{-3s/2}Γ((s+1)/2)Γ((s+11)/2)Γ((s+12)/2)` it matches Iwaniec–Kowalski §5.12 as I recall it. I confirmed it numerically with PARI 2.15.4, using a constructor different from the authors': analytic normalisation, `Vga = [1,11,12]`, exact rational `a_S(n)`.
- `lfuncheckfeq` gives −111 bits.
- The same data with the wrong `Vga = [0,11,12]` gives −2 bits.

This is evidence, not proof.

*To close:* quote the theorem, with its gamma factor, from one source actually consulted (Shimura's main theorem for level 1, or IK §5.12 together with Gelbart–Jacquet). Order 1, boundedness in strips and the Hadamard form then follow, because the Φ-representation of Lemma V0(c) already gives them once the functional equation is known.

**G2. Lemma V0(c): the contour shift that derives `Φ(1/y) = yΦ(y)` from (F3) is asserted, not justified** (lines 414-416).
Moving Mellin inversion from `Re s = 2` to `Re s = −1` needs `Λ_S(σ+iT) → 0` uniformly in `−1 ≤ σ ≤ 2` as `|T| → ∞`, and integrability on interior lines. "Bounded in vertical strips" alone gives neither.

The standard repair is short:
1. `L_S` is polynomially bounded on `Re s = 2` (trivially) and on `Re s = −1` (by the functional equation and Stirling).
2. Phragmén–Lindelöf then gives a polynomial bound for `L_S` in the strip, using finite order.
3. Stirling gives `γ_S ≪ e^{−3π|T|/4}·poly`, which supplies the decay.

Lemma V2 cannot be used here, because it presupposes V0(c). *To close:* add these three lines.

**G3. The citation (F6) could not be verified, and the wording of the fallback in Remark 4.3 is inaccurate** (line 56; Remark 4.3, lines 391-394).
I could not check that `|θ(y) − y| < 0.2y/log² y` for `y ≥ 3 594 641` is BKLNW (Math. Comp. 2021) Corollary 11.2, row k = 2. The same inequality is, as far as I recall, Dusart's (arXiv:1002.0442, Ramanujan J. 2018); citing that would be safer.

The input is not load-bearing for the signs, because of the Rosser–Schoenfeld fallback. But Remark 4.3 says "Theorems 1–3 hold as stated except for the numerical constants in Theorem 2". That is backwards:
- Theorem 1(a) states `L_B` explicitly, with `ε_1` and `c_θ`. So its bound changes under the fallback (e.g. `L_B(10^6) ≥ 0.0230`).
- Theorem 2's constants survive. With `c_θ = 1.01624`, my estimate is `L_B(9·10^6) ≈ 0.1562 > 0.15` and `L_A(9·10^6) ≈ 0.3408 > 1/3`; the change is `2Δc_θ E_odd(9·10^6) ≈ 0.0027`.

*To close:* state precisely which constants change.

### Minor

**4. Lemma 0 (re-checked; correct).**
- The choice of left lines (`mU' ∈ ℕ` for Δ, `mU' ∈ ℕ + 1/2` for `L_S`) keeps distance ≥ 1/2 from trivial zeros in w, and distance ≥ 1/4 from the digamma poles in the half-arguments.
- The horizontal integrals are `≪ log² T_k/T_k²`, uniformly in U', since `∫ t^{σ+1}dσ` converges for t > 1.
- The left integral is `≪ t^{1−U'} log U'/U'`.
- Residues at `ρ/m` are `−(μ/m)t^{ρ/m+1}/((ρ/m)(ρ/m+1))`, and `|ρ/m + 1| ≥ |ρ/m|` because `Re ρ ≥ 0`.
- The order of limits (U' → ∞, then k → ∞) is legitimate.

Cosmetic: the bound `|s(s+1)| ≥ c(U'² + y²)` on the left line needs U' ≥ 2, which is harmless. The reference to IK Prop. 5.7(2) concerns zeros near s; on the horizontal sides the trivial zeros are at distance ≥ T_k, so there is no issue.

**5. The n-tail in `sym2_phi.py` (re-checked; correct).**
`n ↦ n³e^{−βn^{2/3}}` is decreasing for n ≥ 1 when `βn^{2/3} > 9/2`, and the code asserts `β > 4.5`. So `Σ_{n>N} ≤ f(N+1) + ∫_{N+1}^∞ f`, and `∫_{N+1}^∞ n³e^{−βn^{2/3}} dn = (3/2)Γ(6, β(N+1)^{2/3})/β^6` is exact.

The other ingredients also check out:
- `|a_S(n)| ≤ d_3(n) ≤ n²`.
- The u-integral tails `u_lo^{10}e^{−πx²/u_lo²}` and `Γ(10, 2πu_hi)/(2π)^{10}` are valid: `u^9e^{−πx²/u²}` is increasing.
- The recurrence `a_S(p^k) = c a_{k−1} − c a_{k−2} + a_{k−3}`, with `c = τ(p)²/p^{11} − 1`, is right: `e_1 = e_2 = c`, `e_3 = 1`.

I reran `sym2_phi.py 640 160 4`: sha256 `ada15ff1be6ce72b…`, identical to the stored file; maximum radius `2.906·10^{-84}`.

**6. The gap rule in `cert_small.py` (re-checked; correct).**
- For each gap `[q, p]`, Θ_1 is linear and R is monotone and continuous. So `min(T1(q), T1(p))/p^{3/2}` and `min(R(q), R(p))·log p/√p` are valid lower bounds on all of `[q, p]`. The second needs `log x/√x` decreasing, i.e. x > e². That fails for the single gap `[7, 11]`, but that gap is not used for a = 11 (condition p > a).
- Taking every gap with p > a covers `[a, 9 999 991]`.

My own recomputation, with an independent τ table, gives the same minima: G gap-min 0.581460 on [100, ·] and 0.625770 on [10⁴, ·]; N gap-min 0.171897 on [11, ·] and 0.795106 on [100, ·].

**7. The "checks at t = 1" of Lemma 2 (line 199) are an algebraic tautology.**
With `a_0 = log 2π − ψ(11/2) − Λ'/Λ(1)`, `a_{−1} = 2 log 2π − ψ(9/2) − ψ(15/2) − L'/L(2)` and `Σ1/(ρ(ρ+1)) = Λ'/Λ(1) − Λ'/Λ(2)`, the identity reduces to `ψ(11/2) − ψ(9/2) = 2/9`. That holds whatever the values of `L'/L`. So it tests only the trivial-zero bookkeeping. The direct-sum check (ii) is the meaningful one. I did not work out whether the analogous check for `J_S` is also a tautology; it probably is. §7 presents (i) as a second implementation, which overstates it.

**8. Lemma 4/V3 Euler-product tail.**
`Σ_{p>P} log p/(p²−1) ≤ 1.01(log P + 1)/P` is stated without proof. It is true, and very loose: it bounds the sum over all integers n > P by an integral. Add one line.

**9. Naming collision in RESULTS.md.**
Item 1 says "T2 is known (Sheth 2025, Cor. B), T3/T4 dropped", meaning work-order items. Theorems 2-4 of the note carry the same labels. A reader will take this to mean "Theorem 2 is known". Rename.

**10. Literature.**
Hiary–Ireland–Kyi, "A method for verifying the generalized Riemann hypothesis" (arXiv:2408.00187, 2024), is not cited. It is an explicit-formula verification of GRH up to a height, in interval arithmetic, with degree-2 examples including `L(s, Δ)` (window η = 84, 46 zeros). It contains no degree-3 or `Sym²Δ` example, so it does not anticipate Theorem 4. It is the natural reference for the method class of Lemma V3, alongside the Büthe/Weil–Barner method that LITERATURE.md already quotes.

### Checked and found correct (no finding)
- **Lemma 1.** The even prime powers satisfy `Σ_{p^j≤√u} c_{2j}(p) log p = ψ_S(√u) − ψ(√u)`, since `Λ_S(p^j) = (c_{2j}+1) log p`. The signs in (2.1) are right. The bias enters through `+J_ζ ≈ +(2/3)t^{3/2}`, so it is positive, and `λ(p)² − 2 = λ(p²) − 1` has mean −1.
- **Lemma 2.** Residues:
  - at ρ: `−t^{ρ+1}/(ρ(ρ+1))`;
  - at 0: `−a_0 t`;
  - at −1: `+a_{−1}`;
  - at the trivial zeros `−11/2 − n`, which are simple: `−t^{−9/2−n}/((n+11/2)(n+9/2))`, with sum `2/9` at t = 1.

  The other steps check out:
  - Under GRH, `Re 1/(1−ρ) = (1/2)/|ρ|²`, so `σ_0² = 2Λ'/Λ(1)` with `Λ(s) = (2π)^{−s}Γ(s+11/2)L(s)` and ρ the zeros of L, on the line `Re = 1/2`.
  - The inequality is applied in the right direction: `Σ1/|ρ(ρ+1)| < Σ1/|ρ|²`.
  - `a_{−1} − 2/9 > 0.0596`.
- **Lemma 3.** The trivial zeros of `L_S` are at −1, −3, …, −9 (simple), −11, −13, … (double) and −12, −14, … (simple). That gives the residue `2t^{1/2}` at s = −1/2, and `(L_S'/L_S)(−2)` at s = −1 (regular, since `γ_S(−2)` is finite). The subtracted trivial sum is positive.
  - The functional-equation expressions for `a_{S,0}` and `(L_S'/L_S)(−2)` are right.
  - The comparison `Re 1/(1−ρ) ≤ 3 Re 1/(2−ρ)` holds both on the line and for `|γ| > 2`.
  - Hence `a_{S,0} ≥ 0.1577`.
- **Lemma 4.** `25(g+9/4)² ≤ 81(g+1/4)(g+25/4) ⇔ 0 ≤ 56g² + 414g`, so `C_S = (12/5)·Λ_S'/Λ_S(2)`.
- **Lemma 5.** This is truly unconditional and explicit, given F3/F4.
  - `a/(a²+u²) ≥ 1/(2+u²)` for a ∈ [1, 2].
  - The convexity step gives `Σ_{|γ|>H} 2/γ² ≤ 2c∫_{H'}^∞ F T^{−2}`, with the factor 4 · c/2 = 2c correct.
  - Integration by parts gives `≤ 2/(aH'²)` and the term `6 log ζ(2)/H'²`.
  - The duplication-formula form of `γ_S'/γ_S` is right.
  - Binet's formula with `|u²+z²| ≥ 2xy` gives `Re ψ(z) ≤ log|z| + 1/(24xy)`; the factors 1/2 and `y = T/2` combine to `1/(24 x_j T)`.

  I recomputed `T(60, 51.45) = 0.1846043`; the optimum over H' is 0.1846031 at H' = 51.50. The true value of `Σ_{|γ|>60} 2/γ²` is about 0.104 (PARI: `Σ1/|ρ|² = 0.1984216`, of which the zeros below 60 give 0.1463768).
- **Lemma 6 and the assembly of §3.** The `E_odd` count of odd k ≥ 7 is right. The monotonicity threshold is 526.2. The assembly gives exactly `L_A` and `L_B`.
- **Lemma 7 / Theorem 2.**
  - `h''` is right and positive.
  - The closed forms for `A_1` and `A_2` are right. I checked by direct mpmath quadrature of the kernel over prime gaps: `A_1(97) = 0.67197227368526`, `A_2(97) = 10.32644664916`. I also checked the identity `xR(x) = Θ_1(x)/log x + ∫_2^xΘ_1h''` at x = 1500, where it agrees to 40 digits.
  - The coverage of all real x is complete: R is monotone and continuous on each `[p, p')`, x* lies in (5, 7), and R < 0 on (2, 5].
- **Theorem 1, finite range.** Θ_1 is piecewise linear and continuous with nodes at the primes. t* lies in (3, 5). The prime checks plus linearity, together with the analytic range from `10^6`, cover every real t.
- **Theorem 3.**
  - The Mellin identity `∫A x^{−s−2} = P(s)/(s(s+1))` holds.
  - The decomposition `P = log L − (1/2)log L_S(2s) + (1/2)log ζ(2s) + H` holds, with `H = −Σ_{k≥3 odd}`, analytic for `Re s > 1/3`.
  - `L(σ, Δ) > 0` for every real `σ > −11/2`, via the theta integral, so real zeros are excluded.
  - `L_S(2s)` and `ζ(2s)` are zero-free and finite for `Re s > 1/2`.

  So Landau's theorem applies. An off-line zero does not need a separate size comparison against the bias: Landau handles it. The `Θ_1` version is also correct.
- **Theorem 4.**
  - Lemma V0: (a) is right, and the AM–GM bound in (b) is right.
  - Lemma V1 is Poisson summation, with `f̂(ξ) = Λ_S(s−2πiξ)`.
  - Lemma V2: the ray rotation and the bound for `M_θ` check out, and the code reproduces `M_θ = 109.73`.
  - The bounds for `E_tail` and `E_alias` check out.
  - Lemma V3: `|f'| ≤ f/a ≤ f`, the bound `1/(2−β) + 1/(1+β) ≥ 4/3`, and the e^{1/16} Lipschitz factor on a grid of step 1/8 all check out. Multiple zeros, or extra zeros inside a sign-change interval, are caught.
  - In the code, the found-zero terms are lower-bounded over the whole isolating interval (`arb.union`), which is the right direction, and the worst-case tracking is sound.

## What I recomputed (independent code)

1. **τ(n), n ≤ 10^7.** I computed it from `Δ = (E_4³ − E_6²)/1728`, with σ_3 and σ_5 by sieve and exact FLINT polynomial products. This is a different formula from the authors' `η^24`.
   - It agrees with `tau_table.py` at all 10^7 values (0 mismatches).
   - From it, in long double: `Θ_1(p) > 0` at every prime 5 ≤ p ≤ 9 999 991, and `R(p) > 0` at every prime 7 ≤ p ≤ 9 999 991. The minimum of R is 0.1608 at p = 7.
   - The float minimum of G over primes ≥ 10⁴ is 0.625775 at 9 132 733. The infima at primes are 0.6117 (≥ 100) and 0.1747 (N, ≥ 11).
   - The gap minima are as in finding 6. `D_0 < 0` at p = 2 and 17.
2. **Thresholds.** In mpmath at 40 digits: `t* = 4.266791956347541089…` and `x* = 5.518054631888274981…`. Θ_1 and R vanish there to 10⁻⁴¹. `R(3) = −0.17678`, `R(5) = −0.07870`.
3. **L(s, Δ), PARI `lfunetaquo([1,24])`.**
   - `σ_0² = 2Λ'/Λ(1) = 0.11225214047190852674558342625…`, all digits as certified.
   - `a_0 = 0.1706578475916…`, `a_{−1} = 0.2817995196…`, `a_{−1} − 2/9 = 0.05958`, `L(1, Δ) = 0.8393455120319…`.
   - As an independent route, I summed over the 404 zeros with 0 < γ ≤ 400: `2Σ1/|ρ|² = 0.104047`. A Riemann–von Mangoldt tail estimate of 0.008202 brings this to 0.112249, against 0.112252, consistent. Also `2Σ_{γ≤400} 1/|ρ(ρ+1)| = 0.103654 < 2Σ1/|ρ|²`.
4. **Sym²Δ, PARI, my own `lfuncreate`.**
   - Values:
     - `L_S(2) = 0.80587520944866875961…`, `L_S(1) = 0.63179294572788320301…` and `L_S(1/2) = 0.50554937522271351500…`, identical to the Arb values;
     - `Λ_S'/Λ_S(2) = 0.2941512513`, inside the certified `[0.2941460, 0.2941565]`;
     - `a_{S,0} = 0.9410337`;
     - `(L_S'/L_S)(−2) = −0.2282527`, by the functional equation plus an Euler product at 3;
     - `K_S` bound `0.3417553`.
   - `lfunzeros` finds 41 zeros in (0, 60] and 55 in (0, 72]. The first is 7.07022414028743, inside the authors' interval. The last below 60 is 58.96287; the next is 60.62122.
   - Riemann–von Mangoldt: `N(T) − θ_S(T)/π ∈ [−0.50, 0.60]` for T = 10, 20, …, 72, with `θ_S(T) = arg γ_S(1/2+iT)`. There is no unit jump, so no missing or extra zero is indicated. This is a float sanity check, not a rigorous count.
   - `Rem(T)` with the 55 zeros up to 72 removed: 0.06835 (T = 0), 0.11766 (T = 50), 0.18486 (T = 60). This is consistent with the certified bound 0.185062.
5. **Theorem 4 programs.** I reran `sym2_phi.py` and `sym2_verify.py` in a scratch copy (repository untouched). Both outputs are byte-identical (sha256 `ada15ff1…`, `1e6a096b…`).
   - Error terms: `E_tail = 3.17·10⁻¹⁰⁰`, `E_alias = 4.38·10⁻⁹⁵`.
   - The radius of `Z(t)` is `≤ 1.21·10⁻⁸³`.
   - `max Rem` is `≤ 0.185062` on the grid, and 0.196998 after the Lipschitz factor; verdict PASS.
6. **Analytic bounds.** In mpmath, from the certified constants:
   - `L_B(10^6) = 0.027442`, `L_B(9·10^6) = 0.158966`, `L_B(10^7) = 0.163794`, `L_B(10^{14}+) = 0.361228`;
   - `lim L_A = 0.553901`, `lim L_B = 0.369297`;
   - `L_B` vanishes at t ≈ 7.0509·10^5;
   - a monotonicity scan at 200 001 log-spaced points on `[527, 10^{14}]` found no violation.

## What I could not check
- I did not consult the originals of Shimura 1975, Gelbart–Jacquet 1978, Deligne, Hecke, Iwaniec–Kowalski (Theorem 5.6 and Proposition 5.7 are quoted from memory), BKLNW, Dusart or Rosser–Schoenfeld (G1, G3).
- I wrote no second *rigorous* implementation of the Φ table or of the completeness test. My Sym² checks are PARI floating point. Theorem 4's rigour therefore rests on my line-by-line reading of the two programs, on Arb, and on F3/F4.
- I did not rerun `cert_constants.py`, `cert_small.py` or `selftest.py` themselves. Their outputs were reproduced independently (items 1-3 and 6 above) in floating point or via PARI, not in balls.
- I did not reverify the direct-sum checks (ii) of Lemma 2, or the `J_S` check against truncated zero sums.
- Novelty: four web searches found no explicit all-x sign theorem for sums of τ(p), and no GRH equivalence of the kind in Theorem 3. They turned up Koyama–Kurokawa (arXiv:2203.12791), arXiv:2509.04187, and Hiary–Ireland–Kyi (arXiv:2408.00187, finding 10). I did not read Sheth or Hayani beyond what LITERATURE.md records.
