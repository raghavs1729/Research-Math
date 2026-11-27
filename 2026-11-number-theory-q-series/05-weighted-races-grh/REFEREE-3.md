# Referee report 3: `BRIDGE-ATTEMPT.md` (Lemma D, Theorem C, Corollary C′, Proposition N, §4–§6)

Date 2026-11-24. I re-derived every proof from scratch and re-implemented the criterion
independently: my own code, my own closed forms through `mpmath.dirichlet`, my own zeros for
conductors 163 and 232 through python-flint `hardy_z`, and my own Monte Carlo. I did not edit the
document and committed nothing. Referee scripts were kept in the session scratchpad. The criterion
script (`crit.py`, ~60 lines) is reproduced in the appendix so that the numbers below can be
regenerated.

**Verdict: no FATAL error.** There are two REPAIRABLE items: a sign condition missing from the
"K < 1" corollary of Theorem C, and a point-value table that does not show Theorem C's quantities.
There are also several MINOR items. Every proved statement checks out after the stated repairs. The
two numerical claims about real characters are reproduced independently. The conductor-163 claim is
robust. The conductor-232 claim is real, but it is fragile with respect to how the tail is modelled.

---------------------------------------------------------------------------

## 1. Lemma D (derivative identity for all M > 0): correct

Re-derivation.
* `dS/dM = Σ c_j' cos θ_j` with `c_j' = c_j(1−ε_j)/M`, so `dS/dM = (S − Z)/M`. Formally this gives
  `δ' = f(−1)·E[(S−Z)/M | S=−1]`, hence `−Mδ' = f(−1)(1 + E[Z|S=−1]) = f(1) + E[Z;S=−1]`, because f is
  even. Signs and conventions agree with the document.
* Fourier form: `E[X e^{itX}] = i c J_1(ct)` for `X = c cos θ`. Hence
  `E[X_j;S=s] = (1/π)∫_0^∞ c_j J_1(c_j t) Φ^{(j)}(t) sin(ts) dt`. At `s = −1` this is exactly the
  document's `R(M)` with its minus sign.
* Gil–Pelaez gives `P(S > −1) = 1/2 + (1/π)∫_0^∞ Φ sin t / t dt`. The sign is correct.
* Algebra of step (5): `∂_MΦ = (t/M)Φ' + (t/M)Σ ε_j c_j J_1 Φ^{(j)}` is correct. Integration by
  parts gives `∫Φ' sin t = −∫Φ cos t = −π f(1)`. This yields `−Mδ' = f(1) + R`.
* Independent sign check: for χ_4 at M = 2, a central difference of my own Gil–Pelaez δ gives
  `−Mδ' = 0.198276`. The document's `f(1) − |E[Z;S=−1]| = 0.20823 − 0.0100 = 0.1982` agrees.

Justification.
* **Why 9 zeros.** In a leave-one-out product at least 8 of the 9 lowest factors survive. By (J2)
  each is `≤ (8/π)(κt)^{−1/2}` with `κ = c_9(M_a) ≤ c_i(M)` (`c_i` increases in M and decreases in
  γ). Their product is `≤ C t^{−4}`. The weight `t²` (from `c_j' t |J_1(c_j t)| ≤ c_j² t²/(2M)`)
  leaves `t^{−2}`, which is integrable. With fewer zeros this particular argument fails; the number
  9 is a convenient sufficient choice, not a necessity.
* Domination is uniform on `[M_a, M_b]`, so differentiation under the integral is justified, δ is
  C¹, and the series in R converges absolutely:
  `Σ ε_j c_j ∫|J_1||Φ^{(j)}| ≤ (Σc_j²/2)∫ t·min(1, Ct^{−4}) < ∞`. Local uniform convergence of
  `P_N` and of the termwise derivatives is correct.
* **No absolute-continuity or conditional-law argument is used.** f is defined by Fourier inversion
  of the L¹ function Φ, and R is defined as a Fourier series. The probabilistic reading
  `R = E[Z;S=−1]` is only interpretive, as the document says. Finitely many large amplitudes cause
  no problem.
* MINOR: the hypotheses GRH and `L(1/2) ≠ 0` are stronger than needed. The proof uses only "≥ 9
  positive ordinates, `Σγ^{−2} < ∞`".

## 2. Theorem C and its ingredients: correct, with one repair to the "K < 1" corollary

| ingredient | check |
|---|---|
| `\|J_1(z)\| ≤ z/2` | paper Lemma J0env(i); grid check on [0,200]: max(\|J_1\| − z/2) = 0 |
| `\|sin t\| ≤ t` | trivial |
| `−log J_0 = Σ (λ_n/n) z^{2n}` with `λ_n = Σ_k j_{0,k}^{−2n} > 0` | Hadamard product, valid for `\|z\| < j_{0,1}`; `λ_1 = 1/4` (numerically 0.249995 from 2·10⁴ zeros) |
| Lemma J, `\|J_0\| ≤ b` | **true**. 4·10⁶-point grid on [0,200]: `max(\|J_0\| − b) = 1.1e−16` (rounding, at z → 0). Analytically: z ≤ z_* < j_{0,1} by the λ_n; on `[z_*, 2√2]` J_0 decreases from `J_0(z_*) = 0.3416 ≤ ρ` to `J_0(2√2) = −0.19655 > −0.197`; beyond 2√2 the bound `√P(2√2) < 0.44587` holds; for z ≥ z_far, (J2) applies. Large z: `sup √z\|J_0\| = 0.79788 ≈ √(2/π) ≪ 8/π`. |
| `z_* = 1.79715`, `z_far = 32.5995` | recomputed. The summary's "1.797" is slightly below z_*, so `exp(−min(z,1.797)²/4) ≥ b`; it is still a valid envelope (also grid-checked). The true `sup_{z≥j01}\|J_0\| = 0.402759` is confirmed. |
| Lemma L | correct, by positivity of the λ_n. Grid check for z_1 ∈ {0.5,…,2.3}: min(J_0 − lower bound) ≥ −3e−16. |
| `Σ ε_j c_j² = S_4/16` | `ε c² = (M²/(M²+γ²))·16M²/(M²+γ²) = 16M⁴/(M²+γ²)² = c⁴/16`. Correct. |
| closed form `S_4 = 64Mh − 64M²h'` | differentiating `Σ(M²+γ²)^{−1} = h/(2M)` gives `Σ(M²+γ²)^{−2} = h/(4M³) − h'/(4M²)`. Correct. My independent mpmath values at M = 2 are `σ² = 2.387846`, `S_4 = 3.179067`, identical to the document's. |
| n ≥ 9 | needed so that `B¹` has ≥ 8 factors (decay `t^{−4}`) and `∫t²B¹ < ∞`. In case (a) nothing else supplies decay. |

The proof steps are all correct:
* `|Φ| ≤ B` in both cases. In case (b), `Π_{i>n} b(c_i t) = exp(−v t²/2)` for `t ≤ z_*/c_n`, and it
  is `≤ exp(−v z_*²/(2c_n²))` beyond.
* `|Φ^{(j)}| ≤ B¹` for every j: for j ≤ n because `c_1 ≥ c_j`; for j > n via
  `c_j·min(t, z_*/c_n) ≤ c_1 t`. This uses that γ_1 is the lowest ordinate overall.
* The bound on R and the lower bound `MAIN − TAILF ≤ f(1)` both hold, the latter because
  `c_j t ≤ z_1` and `cos ≥ 0` on `[0, t_1]`.

**REPAIRABLE (R1).** "In particular `δ'(M) < 0` whenever `K := NUM/(MAIN − TAILF) < 1`" is false as
written. If `MAIN < TAILF`, then K is negative and so `< 1`, but nothing follows.
*Fix:* state the conclusion as "`δ'(M) < 0` whenever `MAIN − TAILF − NUM > 0`" (equivalently
`MAIN > TAILF` and `K < 1`). The scripts already set `K = ∞` when `MAIN ≤ TAILF`, so no numerical
conclusion is affected.

### Corollary C′ (cell form): correct; proof terse (MINOR)

Monotonicity in M:
* increasing: `c_i(M)`, `σ²(M)`, `S_4(M)`, `ε_j(M)`, `v_n(M)`;
* decreasing: `t_1`'s cap `z_1/c_1(M)` and the breakpoint `z_*/c_n(M)`.

The decreasing breakpoint is the only place where plugging `M_a` into Theorem C's B is not obviously
monotone. It is fine for two independent reasons:
1. The document's chain `Σ_{i>n} min(c_i(M)t, z_*)² ≥ 2v_n(M_a)·min(t, z_*/c_n(M_a))²` bounds the
   true `|Φ|` directly.
2. `v_n/c_n² = ½Σ_{i>n}(M²+γ_n²)/(M²+γ_i²)` is itself increasing in M, so in fact
   `B(M) ≤ B_cell` pointwise.

`B¹` is a product of factors each nonincreasing in M, so `B¹(M) ≤ B¹_cell`. MAIN is monotone
because its integrand decreases in σ² and S_4, and `t_1 ≤ z_1/c_1(M_b)` works for the whole cell.

MINOR: the proof should say explicitly that the leave-one-out bound for j > n (case (b)) also passes
to the cell, for example via (2) above. The statement is true.

Remark 1's heuristic `K ≈ S_4/(32σ²) = Σε_jc_j²/Σc_j²` and its relation to the Gaussian heuristic are
correct algebra.

## 3. Independent re-implementation (χ_4, `prototype/explore/chi4_zeros.txt`)

Setup of my code (appendix):
* Theorem C exactly as stated: `B¹ = B/b(c_1 t)` for all j, closed-form σ² and S_4 via
  `mp.dirichlet`, adaptive quadrature to t = 10⁵, and z_1 maximised over 40 values.
* True `−Mδ'`: central differences (step 10⁻³M) of my own Gil–Pelaez δ, using all 511 zeros, a
  Gaussian factor for the variance above the list (closed form minus list), and Simpson's rule on
  4·10⁵ points.

| M | K (case a, 511) | MAIN | TAILF | NUM | LB = MAIN−TAILF−NUM | K (case b, n=12) | LB (b) | true −Mδ' | certified δ' cell (`cert/results_q4.txt`) |
|---|---|---|---|---|---|---|---|---|---|
| 2 | 0.1056 | 0.20674 | 7.51e−3 | 0.02104 | 0.17819 | 0.0995 | 0.18006 | 0.198276 | (outside cert range) |
| 4 | 0.2341 | 0.12430 | 1.57e−3 | 0.02873 | 0.09400 | 0.2208 | 0.09579 | 0.108809 | δ'(4) = −0.027202 ∈ [−0.027661, −0.026577] ✓ |
| 8 | 0.3399 | 0.07381 | 2.04e−4 | 0.02502 | 0.04858 | 0.3190 | 0.05015 | 0.056613 | −0.0070767 ∈ [−0.0072123, −0.0069296] ✓ |
| 17 | 0.3770 | 0.04348 | 1.30e−5 | 0.01639 | 0.02708 | 0.3470 | 0.02839 | 0.030244 | −0.0017790 ∈ [−0.0018175, −0.0017462] ✓ |
| 50 | 0.4365 | 0.02135 | 5.3e−7 | 0.00932 | 0.01203 | 0.3711 | 0.01343 | 0.013681 | (outside cert range) |

Results:
* **The lower bound is below the truth at every point** (ratios LB/truth: 0.90, 0.86, 0.86, 0.90,
  0.88).
* My true values lie inside the certified derivative enclosures at M = 4, 8, 17.
* K stays below 0.44 in case (a) and below 0.38 in case (b), consistent with the document's
  cell-wise maxima (0.4666 at the top cell 59–60 in case (a); 0.3889 in case (b) on [1.5, 60]).
* MAIN and TAILF agree with the document's table to all printed digits (M = 2: 0.20674, 7.5e−3;
  M = 8: 0.07380, 2.0e−4).

**REPAIRABLE (R2), presentation only.** The "Point values" table of §3 and the sentences drawn from
it ("NUM overestimates |E[Z;S=−1]| by a factor between 1.2 and 2.0", the K column) are **not
Theorem C's quantities**. `bridge/envelope.py` computes something different:
* a per-j leave-one-out NUM, `Σ_j ε_jc_j² ∫t² B/b(c_j t)`, which is sharper than Theorem C's
  `(S_4/16)∫t² B/b(c_1 t)`. It is valid, but it is not stated or proved in the document;
* tail sums of σ² and S_4 from the zero-density integral instead of the closed forms;
* the one-line envelope without the far piece, truncated at t = 60 with no remainder term (harmless
  with 511 zeros, since ρ^511 ≈ 0, but not a bound).

Theorem C as stated gives NUM 6–8 % larger (M = 2: 0.02104 vs 0.01985; M = 8: 0.02502 vs 0.02311)
and K = 0.106 / 0.340 instead of 0.100 / 0.314. The cell runs (`verify_envelope.py`) do use Theorem C
proper (`B¹` without factor 1, closed forms, one-sided sums, a valid far remainder; I checked the
remainder formula and the direction of every Riemann sum). So the headline claim "K_cell ≤ 0.47" is
not affected.
*Fix:* recompute the point table with Theorem C's NUM, or state the per-j variant as a lemma (its
proof is one line: `|Φ^{(j)}| ≤ B/b(c_j t)` for j ≤ n, and `≤ B/b(c_J t)` for j ≥ J).

MINOR: even in the document's own table, the factor NUM/|E[Z;S=−1]| at M = 1.5 is 0.01597/0.0067 = 2.4.
That is outside the stated range "1.2–2.0".

## 4. Proposition N: correct; example reproduced

* `P(c_1cos θ_1 > −a) = 1/2 + arcsin(a/c_1)/π` holds for `|a| < c_1`. The symmetrisation of
  `∂G/∂c_2` is correct; so are `g' ≥ u c_1^{−3}` on `[1−c_2, 1+c_2] ⊂ (0, c_1)`, which gives
  `∂G/∂c_2 ≥ c_2/(2c_1³)`, the bound on `|∂G/∂c_1|`, and the chain rule.
* The limiting claim (M < γ_2/√15 ⇔ c_2 < 1; c_1 → 4 > 1 + c_2; c_1' → 0) is correct.
* Example `{0.03, 10}`, M = 1: `c_1 = 3.99820`, `c_2 = 0.39801`. The two sides of (N) are
  `1.2270e−3 > 3.356e−4`. By my own exact quadrature over θ_2, `δ'(1) = +3.615e−4`, the same as
  `two_zero.out`. Monte Carlo (4·10⁶ samples) gives `δ(1) = 0.58072`, `δ(2) = 0.58137`.
  `δ' > 0` on all of `[0.8, 4]`, including M > γ_2/√15 = 2.58: (N) is sufficient, not necessary.
* **Well-definedness.** `δ = P(S > −1)` is always defined. Here −1 lies strictly inside
  `(−(c_1−c_2), c_1−c_2)`, where the conditional probability is a smooth function of (c_1, c_2)
  uniformly in θ_2 (`|1 + c_2x| ≤ 1.398 < c_1`). So δ is C^∞ near the example even though the law of
  S has unbounded (logarithmic) density at `±c_1 ± c_2`. The proof does not use Lemma D or
  Theorem C, which would need ≥ 9 zeros.
* MINOR: Proposition N proves `δ' > 0` pointwise. "Increasing on an interval" needs the one-line
  remark that (N) is an open condition, uniform for M in compacts of `(0, γ_2/√15)` once γ_1 is
  small (`c_1' ≤ 4γ_1²/M³`).

## 5. Real characters of conductor 163 and 232: plausible, independently reproduced

My own zeros (python-flint `hardy_z`, 64-bit precision, step 0.01, bisection), up to height 40:
* q = 163 (odd, Conrey index 162): 38 zeros, `γ_1 = 0.20290`, `γ_2 = 2.36853`, `γ_3 = 4.05507`.
  The smooth count is ≈ 37.8.
* q = 232 (odd, index 115): 41 zeros, `γ_1 = 0.31068`, `γ_2 = 2.49222`. The smooth count is ≈ 40.1.

Both match the document. Monte Carlo with common random phases: 2·10⁶ samples for 163, 8·10⁶ for
232. Zeros above 40 enter as an independent Gaussian whose variance is my own closed-form σ²
(`mp.dirichlet` with the Legendre/Kronecker values) minus the list part.

| | δ(0.6) | δ(0.9) | δ(1.2) | δ(2.0) | rise |
|---|---|---|---|---|---|
| 163, Gaussian tail | 0.59029 | 0.59404 | 0.59601 | 0.58676 | δ(1.2)−δ(0.6) = **+0.00572 ± 0.00021** |
| 163, list only (zeros ≤ 40) | 0.58949 | 0.59333 | 0.59622 | 0.58865 | +0.0067 |
| document | 0.59033 | 0.59414 | 0.59603 | 0.58664 | +0.0057 |

| | δ(0.72) | δ(0.86) | δ(1.0) | rise |
|---|---|---|---|---|
| 232, Gaussian tail | 0.59746 | 0.59789 | 0.59827 | δ(1.0)−δ(0.72) = **+0.00080 ± 0.00006** |
| 232, list only | 0.59674 | 0.59736 | 0.59795 | +0.0012 |

Assessment:
* **163 is robust.** The rise is ~27 standard errors, and changing the tail model from "dropped" to
  "Gaussian" moves it by only 0.001.
* **232 is real but fragile.** The rise is ~13 standard errors, but the tail treatment moves it by
  0.0004, half the signal. A certificate for 232 therefore needs rigorous tail control at the 10⁻⁴
  level, not a surrogate.

The method's weak points are the following.
1. Sign changes on a 0.01 grid with float midpoints do not prove isolation or completeness. Agreement
   with the *smooth* Riemann–von Mangoldt count "to < 1" is not a proof, because S(T) is unbounded.
2. The Gaussian surrogate is not a bound.
3. Everything is conditional on GRH(χ), and on LI for the race interpretation.

What a certification needs:
* Arb-certified sign changes of Z, plus Turing's method (or an argument-principle count) for
  `N(T) = #list` at some T ≈ 40–80, and enclosures of the ordinates.
* Rigorous δ(M) at two values of M only, for example `δ(1.2) − δ(0.6) > 0` for 163 and
  `δ(1.0) − δ(0.72) > 0` for 232. No derivative is needed. Gil–Pelaez in ball arithmetic: the
  t-integral truncates easily because σ² ≈ 8–12.
* A two-sided bound on the factor from the zeros above T. Lemma L gives
  `exp(−v t²/2 − β S_4^{tail} t⁴) ≤ Π_{γ>T} J_0(c t) ≤ exp(−v t²/2)` for `t ≤ z_1/c_{max,tail}`,
  where `S_4^{tail} ~ Σ c⁴ ≲ 10⁻³`. Beyond that range the bound comes from b and the list factors.
  This should give errors ≪ 10⁻⁴, enough for both characters, with 232 the tighter one.

I could not certify these claims here. My checks raise their plausibility to high for 163 and
moderately high for 232.

## 6. "Theorem 4 loses ~10⁷ in its tail at M = 17.3": plausible

The paper's Proposition criterion-eval reports a worst relative margin `(MAIN−TAIL)/MAIN = 0.105` at
the left edge. Hence the paper's TAIL bound is ≈ 0.895 × its MAIN bound. With MAIN bound / total =
0.277 (the document's `thm4_loss_q4.out`, which I spot-checked through its own model; its values of
−δ' at M = 4, 8 match my independent ones to 4 digits), that is ≈ 0.25 of −δ'.

The true value of the paper's TAIL integral, `(1/π)∫_{3/5}^∞ |∂_MΦ| |sin t|/t`, is ≈ 2.9e−8 of −δ'.
The ratio is ≈ 8.5·10⁶, so "~10⁷" is correct in order of magnitude. Mixing a point MAIN ratio with
a cell margin affects this only at the 10–20 % level. The stated cause (only the ~12 ordinates
≤ 2M_a ≈ 34.6 are kept, so the Gaussian factor of the higher zeros, e^{−v t²/2} at t ≥ 0.6, is
lost) is consistent.

**MINOR (wrong side remark in §6(c)).** "With exact pieces the same split would work down to M ≈ 5
but not to M_0" is contradicted by the same model. Running `thm4_loss` at smaller M gives:

| M | A+B | \|C\| |
|---|---|---|
| 3.689 | 0.617 | 0.39 |
| 3 | 0.449 | 0.56 |
| 2 | 0.205 | 0.80 |

So the exact split still works at M_0 and fails only between M = 3 and 3.69. The sentence should be
corrected; it does not affect any result.

## 7. Other MINOR items

* Summary item 1 says "the one-line envelope `exp(−min(z,1.797)²/4)`", but Theorem C uses the
  three-piece b. The far piece `(8/π)z^{−1/2}` is what gives decay, and the one-line envelope alone
  does not decay.
* §3 states "K_cell ≤ 0.39 on [M_0, 20]". I did not rerun the full cell scan (load limits). My point
  values (0.377 at M = 17, case (a)) are consistent with it.
* The numerical sections §4.2 and §5 were not independently re-checked beyond the spot checks above.

---------------------------------------------------------------------------

## Certification list

**Proved as written** (I re-derived each):
* Lemma D, for every M > 0 (δ is C¹; the identity holds; the integrals and series converge
  absolutely).
* Lemma J; Lemma L.
* Theorem C's inequality `−Mδ'(M) ≥ MAIN − TAILF − NUM`, in cases (a) and (b), including
  `Σε_jc_j² = S_4/16`.
* The closed forms for σ² and S_4 (given the variance identity).
* Corollary C′. The statement is true; the proof is terse (MINOR, §2).
* Proposition N, including the limiting claim and the numerical example `{0.03, 10}`, M = 1.

**Proved after repair:**
* "δ' < 0 whenever K < 1". Repair (R1): require `MAIN − TAILF − NUM > 0`.

**Not proved (numerical evidence only; independently corroborated where stated):**
* §3: the cell criterion on [1.5, 400] for q = 3, 4. Floating point, not Arb. My independent point
  values agree and are lower bounds below the true and certified values. The point-value table
  needs relabelling (R2).
* §4.2.
* §4.3: δ increasing for conductors 163 and 232. Reproduced by independent zeros and Monte Carlo;
  163 robust, 232 tail-sensitive. Not certified.
* §5.
* §6(c): the ~10⁷ loss is plausible; the "M ≈ 5" remark is wrong.

---------------------------------------------------------------------------

## Appendix: referee script `crit.py` (independent implementation)

```python
import numpy as np, mpmath as mp, sys
from scipy.special import j0
from scipy.integrate import quad, simpson
mp.mp.dps=30
g_all=np.loadtxt('prototype/explore/chi4_zeros.txt')
chi=[0,1,0,-1]
def h(M):   # (Lambda'/Lambda)(M+1/2), Lambda(s)=(4/pi)^{s/2} Gamma((s+1)/2) L(s,chi4)
    s=mp.mpf(M)+mp.mpf(1)/2
    return mp.log(4/mp.pi)/2+mp.digamma((s+1)/2)/2+mp.dirichlet(s,chi,1)/mp.dirichlet(s,chi)
def sig2S4(M):
    hh=h(M); hp=mp.diff(h,M); return float(4*M*hh), float(64*M*hh-64*M*M*hp)
rho=0.446; zs=2*np.sqrt(np.log(1/rho)); zfar=(8/(np.pi*rho))**2
def logb(z):
    z=np.asarray(z,float)
    return np.where(z<=zs,-z*z/4,np.where(z<=zfar,np.log(rho),np.log(8/np.pi)-0.5*np.log(np.maximum(z,1e-300))))
def beta(z1): return (-np.log(j0(z1))-z1*z1/4)/z1**4
def criterion(M,n=None):
    g=g_all if n is None else g_all[:n]
    c=4*M/np.sqrt(M*M+g*g); s2,S4=sig2S4(M)
    v=0.0 if n is None else s2-np.sum(c**2)/2; cn=c[-1]
    logB=lambda t: logb(c*t).sum()-v/2*min(t,zs/cn)**2
    B=lambda t: np.exp(logB(t)); B1=lambda t: np.exp(logB(t)-logb(c[0]*t))
    br=np.concatenate([[0],np.geomspace(1e-2,1e5,200)])
    num=S4/(32*np.pi)*sum(quad(lambda t:t*t*B1(t),a,b,limit=200,epsabs=0,epsrel=1e-10)[0] for a,b in zip(br[:-1],br[1:]))
    best=None
    for z1 in np.linspace(0.4,2.35,40):
        t1=min(np.pi/2,z1/c[0])
        main=quad(lambda t:np.exp(-s2*t*t/2-beta(z1)*S4*t**4)*np.cos(t),0,t1,epsrel=1e-12)[0]/np.pi
        brt=np.concatenate([[t1],br[br>t1]])
        tf=sum(quad(B,a,b,limit=200,epsabs=0,epsrel=1e-10)[0] for a,b in zip(brt[:-1],brt[1:]))/np.pi
        if best is None or main-tf>best[0]-best[1]: best=(main,tf,z1)
    main,tf,z1=best
    return dict(K=num/(main-tf),MAIN=main,TAILF=tf,NUM=num,LB=main-tf-num)
def delta(M,N=400001):     # Gil-Pelaez with 511 zeros + Gaussian for the closed-form remainder variance
    c=4*M/np.sqrt(M*M+g_all**2); s2,_=sig2S4(M); vt=s2-np.sum(c**2)/2
    T=max(40.0,60/np.sqrt(s2)); t=np.linspace(0,T,N); t[0]=1e-300
    Phi=np.prod([j0(ci*t) for ci in c],axis=0)*np.exp(-vt*t*t/2)   # (log-sum form used in the run)
    y=Phi*np.sin(t)/t; y[0]=1.0
    return 0.5+simpson(y,x=np.linspace(0,T,N))/np.pi
# -M delta'(M) = -M (delta(M+h)-delta(M-h))/(2h), h = 1e-3 M
```
