# The analytic leg `M ≥ 300`: referee report and ball-arithmetic certificate

Date: 2026-11-24. Object: `paper/weighted-races.tex`, Section "Monotonicity in the
weight for large m" (l. 3621–4700): `lem:phi`, `lem:J0env`, `lem:binet`, `lem:Nplus`,
`lem:Nall`, `lem:Nlow`, `lem:Dprop`, `lem:diffint`, `lem:negstruct`, `lem:tailexp`,
the proof of `thm:monotone-large` (`eq:split`, `eq:mainlb`, `eq:tailexp`),
`lem:analytic-leg`, `rem:monotone-rigor`, and what they import from earlier sections
(`prop:varid`/`eq:varidq`, `lem:S4mom`/`eq:S4id`, `lem:J0facts` (J1)–(J3), `lem:invert`).
The paper was not edited. Scripts and outputs are in `monotonicity/analytic-leg/`.

## 0. Verdict

**(a).** The companion's chain, read from first principles, proves

> **Proposition (analytic leg).** Assume GRH(χ_q), q ∈ {3, 4}. Then δ_q is
> differentiable on (299, ∞) and δ_q′(M) < 0 for all M ≥ 300.

No FATAL finding. The one place where the paper itself says the argument was "checked
in float64 only" (the extension from the point M = 300 to the half-line) is a correct
hand proof with three unstated monotonicity facts; I wrote it out in full (§4, Step 7)
and certified every numerical fact it uses in Arb ball arithmetic
(`analytic-leg/certify_analytic_leg.py`, 90 checks, all PASS, 1.3 s). §4 is a complete,
self-contained proof of the Proposition, so the result does not depend on trusting the
paper's wording anywhere. Three REPAIRABLE items (one of them outside the Proposition)
and a list of MINOR items are in §3.

Quantitative form that comes out of the proof: for M ≥ 300,

    −δ_q′(M) ≥ (1 − R_q) · MAIN_an(M),   R_4 = 4.6×10⁻⁴,  R_3 = 3×10⁻³,

with `MAIN_an` the explicit function of §4. The bound sits between 3.6 % and 4.5 % below
the true |δ′| at M = 300, 310, 320, 350 (§6), which is the loss the constants predict.

## 1. Method

1. Every lemma of the chain re-derived by hand: hypotheses, each inequality, each
   constant (§2).
2. Every numerical constant re-evaluated in Arb (python-flint 0.9.0, 50 digits), where
   possible by a route different from the paper's (sections [A], [A′] of the script).
3. The final step proved analytically for all M ≥ 300 (§4 Step 7) with its finitely
   many numerical inputs certified ([B], [C]); then confirmed independently by a cell
   sweep of [300, 10⁷] that does not use the hand monotonicity argument ([D]).
4. The explicit closed-form bounds on D, σ², S₄ and the zero-count bounds compared with
   Arb evaluations of the exact closed forms and with the certified zero lists
   ([E], [F]; sanity only).
5. The resulting bound for δ′ compared with the direct certificate
   (`cert/results_q*_stretch.txt`) and with seven new sharp cells computed with the
   same prover ([G], §6).
6. A float64 test of `eq:Abound`, the zone products and the pointwise main-term
   inequalities against the actual zeros at M = 40, 100, 300
   (`analytic-leg/sanity_float.py`; not part of the proof).

## 2. Lemma-by-lemma referee

Notation: χ real primitive odd mod q ∈ {3,4}; γ the positive ordinates with
multiplicity; a_γ = 2M/√(M²+γ²); ℓ = L_M = log(qM/2π); N(T) = #{γ ≤ T}.

| Item | What I checked | Result |
|---|---|---|
| `prop:varid`, `eq:varidq` | root number +1 (τ(χ₄)=2i, τ(χ₃)=i√3), Hadamard product for Ξ(z)=ξ(½+z), termwise log-derivative; GRH used only to make γ real; L(½,χ)≠0 needed for Ξ(0)≠0 | correct; L(½,χ₄)=0.66769…, L(½,χ₃)=0.48087… certified > 0 (also elementary: group the series in pairs) |
| `lem:S4mom` (`eq:S4id`) | termwise algebra 16·8M²/(M²+γ²) − 64M²·2(γ²−M²)/(M²+γ²)² = 256M⁴/(M²+γ²)² | correct |
| `lem:J0facts` (J1)–(J3) | (J2) proof via the integral representation, split at u = 1 − 1/z; (J3) alternating series | correct. Confirmed: \|J₀(z)\|√z < 8/π on [1,400] (bisection); sup of (−log J₀ − z²/4)/z⁴ on (0,1] is its value at 1 = 0.01762 ≤ 1/17 (positivity of the Rayleigh coefficients) |
| `lem:phi` | λ₁=¼, λ₂=1/32, λ₃=1/192 (re-expanded); φ ≥ z/2; φ ≤ 0.517z on (0,½] (sup = φ(½)/½ = 0.51631); J₀(2.4) = 0.0025077 > 0 | correct; wording slip in (iii), see m1 |
| `lem:J0env` | P′ = −2J₁²/z; √P(2√2) = 0.445855 < 0.446; J₀(1.0733) = 0.73209 < 0.733; 12/(5√5) = 1.07331 > 1.0733 | correct. Independent bisection: sup_{[1.0733,64]}\|J₀\| < 0.733, sup_{[1.697,64]}\|J₀\| < 0.447, and (8/π)/8 = 0.318 beyond |
| `lem:binet` | Binet's formula for ψ; ∫₀^∞ t/(e^{2πt}−1) dt = 1/24; the ψ′ formula | correct |
| `lem:Nplus` | each γ ≤ T contributes ≥ 4 to 4T ξ′/ξ(T+½); ψ(x) ≤ log x; \|L′/L(s)\| ≤ −ζ′/ζ(2) = 0.569961 < 0.57 for s ≥ 2 ⇔ T ≥ 3/2 | correct; GRH-conditional and labelled so |
| `lem:Nall` | rests on `lem:rvm-joint`(b), i.e. on [MV07, Thm 14.5] / [BMOR21, Cor. 1.2] | not needed for q ∈ {3,4}: see R2. The BMOR quotation is accurate (I read the arXiv text: conductor q > 1, T ≥ 5/7, bound min{0.247 log qT + 6.894, 0.298 log qT + 4.358}) |
| `lem:Nlow` | S₄/8 − σ² = Σ 8u(4u−1) ≤ 24·N(√3Y); explicit bounds re-derived (§4 Step 2) | correct; GRH-conditional and labelled so. No zero-count theorem from the literature is used |
| `lem:Dprop` (i)–(iii) | D = Σ 4a a′; closed form; (iii) re-derived: \|D−2ℓ−2\| ≤ 5.92/M² + (2.28+3.56M)2^{−(M−3/2)} ≤ 8/M | correct (very generous: 0.02/M would do at M ≥ 300) |
| `lem:diffint` | (a) head/tail product rule; (b) \|J₁(z)\| ≤ min(1, z/2), a′ ≤ 2/M, count of a_γt > 1; (c) dominated convergence | (a), (b) correct; (c) correct, with an unnecessary dependence, R2 |
| `lem:negstruct` | 4·(3/5) = 2.4 < j₀,₁; −∂_MΦ = ΦB ≥ 0; Φ ≥ e^{−βt²}; B ≥ Dt²/2; sin t/t ≥ 1 − 1/384 > 0.9973; the closed form of ∫₀^{1/8} t²e^{−βt²}dt | correct |
| `lem:tailexp` | cell bounds, drop-the-smallest-factor, κ-product, I₁, I₂ | correct (not used by the Proposition) |
| proof of `thm:monotone-large` | 0.189, 0.18, 0.226, 1.133, 0.90, 0.996, 0.9998, 0.003 | correct except one arithmetic slip in the sandwich step, R3 |
| `lem:analytic-leg`, bounds | zone products, 64D₊+512ℓ₀+1165, c₃⁹, I₁′, I₂′, β₊ | correct |
| `lem:analytic-leg`, M = 300 | ratio 4.5013×10⁻⁴ (q=4), 2.9098×10⁻³ (q=3); weakened variant 1.822×10⁻³, 1.180×10⁻² | all four printed claims certified |
| `lem:analytic-leg`, [300,∞) | log-derivative argument | correct; exposition gap and float64 status, R1 |

## 3. Findings

### FATAL
None.

### REPAIRABLE

**R1. The extension from M = 300 to [300, ∞) (last paragraph of `lem:analytic-leg`;
`rem:monotone-rigor`, "third").** The statement and the constants (−0.15, −0.09,
0.001, −0.007, hence −0.14 and −0.08) are right. What was missing:

* the supporting script (`paper/v3/verify_monotonicity.py`, item [5]) only takes a
  float64 central difference of log(T_i/MAIN_an) at the single point M = 300, which
  says nothing about M > 300; `verify_monotonicity_arb.py` certifies only the value
  at M = 300;
* the written proof uses, without saying so, that (i) D₊ is increasing (D₊′ = 2/M −
  8/M² > 0), so the prefactors 64D₊+512ℓ₀+1165 and D₊/2+4ℓ₀ are increasing and the two
  positive correction terms are largest at M = 300; (ii) β₊ ≥ Mℓ needs ℓ > 1;
  (iii) the bound −(3/2)(1.03ℓ+1)/(Mℓ) is increasing in M.

*Repair:* §4 Step 7 (full statement and proof) plus sections [B], [C] of the
certificate; section [D] re-proves the same thing on [300, 10⁷] by an independent
route (exact log-derivatives enclosed on 1047 cells; they are ≤ −0.1470 and ≤ −0.0861
for q = 3, ≤ −0.1587 and ≤ −0.0926 for q = 4).

**R2. `lem:diffint`(c) depends on `lem:Nall`.** The dominating function on
t ∈ [1/8, 16] is written with Ñ_χ of `lem:Nall`, i.e. with the unspecified absolute
constant C_R′ of `lem:rvm-joint` (literature input, hypothesis "(H)"). For the
Proposition this is the only place where anything outside the self-contained chain
enters, and it is unnecessary. *Repair:* on the compact range [1/8, 16] the constant
(64/M_a)·N_χ(32M_b) < ∞ dominates (zeros of an entire function are discrete); for
t ≥ 16 use N⁺ of `lem:Nplus` (2M_at ≥ 3/2 trivially). Written out in §4 Step 3.

**R3. Arithmetic slip in the sandwich step of `thm:monotone-large` (l. 4391–4396).
Outside the Proposition.** The text says "On [1/8, 3/5] every window argument is
≥ (4/√5)·(1/8) = 0.447, so each window factor is ≤ J₀(0.447) ≤ 0.9510". But
(4/√5)/8 = 1/(2√5) = 0.2236 (0.447 is 1/√5; `lem:bigt` has the right value). The
correct factor bound is J₀(1/(2√5)) = 0.98754 ≤ 0.9876. *Repaired statement:* "…is
≥ 1/(2√5), so each window factor is ≤ J₀(1/(2√5)) ≤ 0.9876, and the integrand is
≪_q log M · 0.9876^{N₀}". The conclusion (superexponentially small, hence
θ(m) ≤ 1.14 for M ≥ M₂(χ)) is unchanged, since N₀ ≥ (M/4π) log M; only the
(unspecified) threshold M₂(χ) grows. Nothing in `lem:analytic-leg` or in
`cor:monotone-explicit` uses this step.

### MINOR

* **m1.** `lem:phi`(iii): the text bounds "the alternating remainder after c₆" but the
  displayed sum stops at c₅ and subtracts 10⁻⁶. Correct justification: the remainder
  after c₅ lies in [0, c₆]. The claim J₀(2.4) > 0.0024 is true (0.0025077).
* **m2.** `lem:analytic-leg` states no hypotheses; it inherits GRH(χ) from the
  section's standing assumption. Since it is cited on its own elsewhere
  (`monotone-coupling.tex`, REFEREE-1), it should begin "Assume GRH(χ_q)".
* **m3.** `lem:analytic-leg`, proof: "for M ≥ 300, q ≥ 3". The inputs
  (`lem:Dprop`(iii), `lem:Nlow`) are proved for q ∈ {3,4} only; write q ∈ {3,4}.
* **m4.** "the remaining ≥ n₂ − 10 keep 0.733" presupposes N(2M) ≥ 10; true
  (n₂(300) ≥ 66.7).
* **m5.** `lem:Nlow`: 8u(4u−1) ≤ 24, and the proof divides by 32; dividing by 24
  would give 4/3 more zeros. Not needed.
* **m6.** `lem:tailexp`: "p = n′/2 − 1 ≥ 2"; with n ≥ 12 in fact p ≥ 4.5. Harmless.
* **m7.** `verify_monotonicity.py`: stale docstring for item [5] ("M ≥ 250",
  "0.61^{n1}"); item [2] prints "claim ≤ 1.1" for Σ Λ(n) log n / n², the paper uses
  0.89 (true value 0.884482).
* **m8.** `rem:monotone-rigor`: "can be done by hand" and "still checked in float64
  only" can both be replaced by a reference to the Arb run of this audit.
* **m9.** `lem:Nplus` loses a factor π against the true count and `lem:Nlow` a factor
  ≈ 7 (n₁(300) = 28.8 against N(300) = 203 for χ₄). This is why the ratio at M = 300
  is only 10⁻³–10⁻⁴ rather than astronomically small. Not an error.
* **m10.** Arb's `bessel_j` on a wide input ball is valid but very pessimistic
  ([6 ± 0.03] ↦ [± 1.06]); the script therefore encloses sup\|J₀\| on an interval by
  \|J₀(mid)\| + radius (mean value theorem, \|J₁\| ≤ 1). Relevant to anyone reusing
  wide-ball Bessel evaluations.

Is anything GRH-conditional that is not labelled so? No. `lem:Nplus`, `lem:Nlow`,
`lem:Dprop` carry the hypothesis; `lem:diffint`–`lem:analytic-leg` inherit it from the
section preamble (see m2). `lem:window` is unconditional and labelled so, and is not
used for M ≥ 300.

## 4. The Proposition, with a complete proof

**Setting.** q ∈ {3,4}, χ = χ_q the real primitive (odd) character mod q. Assume
GRH(χ). Let γ run over the positive ordinates of L(s,χ) with multiplicity. For M > 0,

    a_γ(M) = 2M/√(M²+γ²),   Φ_M(t) = Π_γ J₀(2a_γ(M)t),
    δ_q(M) = P(1 + Σ_γ 2a_γ(M) cos θ_γ > 0),   θ_γ i.i.d. uniform on [0,2π).

(M = m + ½ in the paper; ′ = d/dM.) Write ℓ = log(qM/2π), N(T) = #{γ ≤ T}.

**Proposition.** δ_q is differentiable on (299, ∞) and δ_q′(M) < 0 for every M ≥ 300.
More precisely −δ_q′(M) ≥ MAIN_an(M) − TAIL_an(M) ≥ (1 − R_q)MAIN_an(M) > 0 with
R_4 = 4.6×10⁻⁴, R_3 = 3×10⁻³ and

    MAIN_an(M) = (0.9973/π)(D₋/2) G(β₊),      G(β) = ∫₀^{1/8} t² e^{−βt²} dt
                                                   = (√π/(4β^{3/2})) erf(√β/8) − e^{−β/64}/(16β),
    D_± = 2ℓ + 2 ± 8/M,     β₊ = Mℓ + 1.55 + (32M(ℓ−1) + 100)/1088,
    TAIL_an(M) = T₁(M) + T₂(M),
    T₁ = (1/(0.447π)) P₁ r₁^{n₁} r₂^{n₂},      P₁ = 64D₊ + 512ℓ₀ + 1165,
    T₂ = (c₃⁹/π) P₂ r₂^{n₂−10},                P₂ = (D₊/2)I₁′ + 4(ℓ₀I₁′ + I₂′),
    r₁ = 0.447/0.733,  r₂ = 0.733,  c₃ = 4·5^{1/4}/π,
    ℓ₀ = log(q(2M+2)/2π) + 1.14,   I₁′ = 16^{−5/2}/2.5,   I₂′ = 16^{−5/2}(2.5 log 16 + 1)/6.25,
    n_i = Y_i(log(qY_i/2π) − 2)/16 − 1/2,   Y_i = iM/√3   (i = 1, 2).

**Step 0. Inputs.**

* (I1) L(½,χ) > 0 (pairing consecutive terms of the series; values 0.6677, 0.4809).
* (I2) ξ(s) = (q/π)^{(s+1)/2}Γ((s+1)/2)L(s,χ) is entire of order 1 and ξ(s) = ξ(1−s)
  (root number +1). With GRH and (I1), Hadamard gives ξ(½+z)/ξ(½) = Π_γ(1 + z²/γ²),
  and termwise differentiation (dominated by const/γ², Σγ⁻² < ∞) gives, for M > 0,

      σ²(M) := Σ 2a_γ² = Σ 8M²/(M²+γ²)     = 4M (log ξ)′(M+½),
      D(M)  := Σ 4a_γa_γ′ = Σ 16Mγ²/(M²+γ²)² = 4(log ξ)′(M+½) + 4M(log ξ)″(M+½) = (σ²)′,
      S₄(M) := Σ (2a_γ)⁴ = Σ 256M⁴/(M²+γ²)² = 16σ² − 64M²(log ξ)″(M+½),

  with (log ξ)′(s) = ½log(q/π) + ½ψ((s+1)/2) + L′/L(s), (log ξ)″(s) = ¼ψ′((s+1)/2) + (L′/L)′(s),
  and a_γ′ = 2γ²/(M²+γ²)^{3/2} ∈ (0, 2/M].
* (I3) Bessel. (B1) \|J₀\| ≤ 1, \|J₁(z)\| ≤ min(1, z/2). (B2) \|J₀(z)\| ≤ (8/π)z^{−1/2}.
  (B3) for 0 ≤ z ≤ 1: J₀ ≥ ¾ and log J₀(z) ≥ −z²/4 − z⁴/17. (B4) j₀,₁ > 2.4, and on
  (0, j₀,₁): φ := J₁/J₀ = Σ 2λ_n z^{2n−1} with all λ_n > 0, so φ(z) ≥ z/2; φ(z) ≤ 0.517z
  on (0, ½]. (B5) P = J₀²+J₁² is decreasing, J₀ is decreasing on [0, 2√2], hence
  sup_{z≥z₀}\|J₀\| ≤ max(J₀(z₀), √P(2√2)) for z₀ ≤ 2√2; with √P(2√2) < 0.446,
  J₀(1.0733) < 0.733, J₀(1.697) < 0.410:

      sup_{z ≥ 1.0733} |J₀| ≤ 0.733,      sup_{z ≥ 1.697} |J₀| ≤ 0.447.

* (I4) Binet: log x − 1/(2x) − 1/(12x²) ≤ ψ(x) ≤ log x − 1/(2x),
  1/x + 1/(2x²) ≤ ψ′(x) ≤ 1/x + 1/(2x²) + 1/(6x³) (x > 0).
* (I5) For s ≥ 2: \|L′/L(s)\| ≤ 0.57·2^{−(s−2)}, \|(L′/L)′(s)\| ≤ 0.89·2^{−(s−2)}
  (n^{−s} ≤ n^{−2}2^{−(s−2)}; Σ Λ(n)n⁻² = 0.569961, Σ Λ(n) log n·n⁻² = 0.884482).
* (I6) If Φ_M ∈ L¹ then δ_q(M) = ½ + (1/π)∫₀^∞ Φ_M(t) (sin t/t) dt (`lem:invert`).

**Step 1. Explicit bounds for M ≥ 20.** Put x = M/2 + 3/4 = (M/2)(1+u), u = 3/(2M),
E_L = (2.28 + 3.56M)·2^{−(M−3/2)} (≤ 2×10⁻⁴ for M ≥ 20).

* (E1) 2Mℓ ≤ σ² ≤ 2Mℓ + 1.4. Upper: 2Mψ(x) ≤ 2M log(M/2) + 2M log(1+u) − 2/(1+u)
  ≤ 2M log(M/2) + 3 − 2 + 3/M, plus 2.28M·2^{−(M−3/2)}. Lower: 2M log(1+u) ≥ 3 − 9/(4M),
  −M/x ≥ −2, −M/(6x²) ≥ −2/(3M); total ≥ 2Mℓ + 1 − 3/M − E_L > 2Mℓ.
* (E2) \|D − 2ℓ − 2\| ≤ 5.92/M² + E_L ≤ 8/M, i.e. D₋ ≤ D ≤ D₊. (D = 2log(q/π) + 2ψ(x)
  + Mψ′(x) ± E_L; upper: 2u − (2/M)/(1+u) + [2 − 2u + 2u²] + (2/M)/(1+u)² + 4/(3M²);
  lower: 2u − u² − (2/M)/(1+u) − 2/(3M²) + [2 − 2u] + (2/M)/(1+u)², and
  1/(1+u) − 1/(1+u)² ≤ u.)
* (E3) 32M(ℓ−1) − 34 ≤ S₄ ≤ 32M(ℓ−1) + 100. (16M²ψ′(x) ≤ 32M + 32 + 64/(3M) and
  ≥ 32M − 48 + 32 − 96/M; combine with (E1).)
* (E4) Hence β := σ²/2 + S₄/1088 ≤ β₊ (the paper's weaker σ² ≤ 2Mℓ + 3.1 is used in β₊).

**Step 2. Zero counts (GRH).**

* (C1) N(T) ≤ N⁺(T) := (T/2) log(q(T+2)/2π) + 0.57T for T ≥ 3/2: each γ ≤ T
  contributes 8T²/(T²+γ²) ≥ 4 to σ²(T) = 4T(log ξ)′(T+½), and
  (log ξ)′(T+½) ≤ ½log(q/π) + ½log((T+3/2)/2) + 0.57.
* (C2) N(√3·Y) ≥ Y(log(qY/2π) − 2)/16 − ½ for Y ≥ 20: with u_γ = Y²/(Y²+γ²),
  S₄(Y)/8 − σ²(Y) = Σ 8u_γ(4u_γ − 1) ≤ 24·#{u_γ ≥ ¼} = 24·N(√3Y) ≤ 32·N(√3Y)
  (terms with u_γ < ¼ are negative, the others are ≤ 24), and by (E1), (E3) the left
  side is ≥ 4Yℓ_Y − 4Y − 4.25 − 2Yℓ_Y − 1.4 ≥ 2Y(ℓ_Y − 2) − 6; divide by 32.
  Consequently, for M ≥ 35: **N(M) ≥ n₁(M), N(2M) ≥ n₂(M).**

**Step 3. Differentiability (M₀ ≥ 300 fixed, M ∈ [M_a, M_b] := [M₀−1, M₀+1]).**
By (C2), S := {γ ≤ 2M_a} has ≥ n₂(299) ≥ 66 elements, and a_γ(M) ≥ 2/√5 for γ ∈ S.

* Pointwise: split Φ = (finite head)·(tail with all arguments ≤ ½). The tail is
  exp Σ log J₀(2a_γt), termwise M-derivatives −φ(2a_γt)·2a_γ′t are bounded by
  0.517·4a_γa_γ′t² ≤ 0.517·16M_b t²/γ², summable uniformly on [M_a,M_b]. So
  M ↦ Φ_M(t) is differentiable, and

      ∂_MΦ_M(t) = Σ_γ J₀′(2a_γt)·2a_γ′t·Π_{γ′≠γ} J₀(2a_γ′t),

  whence, by (B1) and for every finite set S′ of ordinates,

      |∂_MΦ_M(t)| ≤ A(t)·sup_γ Π_{γ′∈S′∖{γ}} |J₀(2a_γ′t)|,
      A(t) := Σ_γ min(1, a_γt)·2a_γ′t ≤ D t²/2 + (4t/M)·N(2Mt)          (*)

  (terms with a_γt ≤ 1 give ≤ Σ 2a_γa_γ′t² = Dt²/2; terms with a_γt > 1 have γ < 2Mt
  and are ≤ 2a_γ′t ≤ 4t/M). For 4t < j₀,₁ all factors are positive and
  ∂_MΦ = −Φ·B, B(t) = Σ φ(2a_γt)·2a_γ′t ≥ 0.
* Integrable dominant for ∂_MΦ·(sin t)/t, uniform in M ∈ [M_a, M_b], with
  D_max = max D: on [0, 1/8], 0.517·D_max t² (as B ≤ 0.517Dt², \|Φ\| ≤ 1); on [1/8, 16],
  the constant 128D_max + (64/M_a)N(32M_b); on [16, ∞), by (*) with S′ = S, (B2) and (C1),
  [D_max t/2 + (4/M_a)N⁺(2M_bt)]·(c₃/√t)⁷ ≪ t^{−5/2} log t.
  The same bounds give Φ_M ∈ L¹. By (I6) and dominated convergence,

      δ_q′(M) = (1/π) ∫₀^∞ ∂_MΦ_M(t) (sin t/t) dt        (M₀−1 < M < M₀+1).

**Step 4. Split.** For 0 ≤ t ≤ 3/5 all arguments are ≤ 4t ≤ 2.4 < j₀,₁, so
−∂_MΦ = ΦB ≥ 0 there, and sin t/t > 0. Dropping (1/8, 3/5]:

    −δ_q′(M) ≥ MAIN − TAIL,   MAIN = (1/π)∫₀^{1/8} Φ B (sin t/t) dt,
                              TAIL = (1/π)∫_{3/5}^∞ |∂_MΦ| (|sin t|/t) dt.

**Step 5. MAIN ≥ MAIN_an(M).** On [0, 1/8] all arguments are ≤ ½. By (B3),
Φ ≥ exp(−σ²t²/2 − S₄t⁴/17) ≥ e^{−βt²} (t⁴ ≤ t²/64, 64·17 = 1088); by (B4),
B ≥ Σ a_γt·2a_γ′t = Dt²/2; sin t/t ≥ 1 − 1/384 > 0.9973. So
MAIN ≥ (0.9973/π)(D/2)G(β) ≥ (0.9973/π)(D₋/2)G(β₊) by (E2), (E4) and G decreasing.

**Step 6. TAIL ≤ TAIL_an(M).** Take S′ = {γ ≤ 2M} in (*); it has N(2M) elements, N(M)
of them with γ ≤ M.

* 3/5 ≤ t ≤ 16. For γ ≤ 2M: 2a_γt ≥ (4/√5)(3/5) = 1.07331 > 1.0733, factor ≤ 0.733.
  For γ ≤ M: 2a_γt ≥ 2√2·(3/5) = 1.69706 > 1.697, factor ≤ 0.447. The product over S′
  minus one factor is ≤ 0.447^{N(M)}·0.733^{N(2M)−N(M)}/0.447 = r₁^{N(M)} r₂^{N(2M)}/0.447
  ≤ r₁^{n₁} r₂^{n₂}/0.447 (Step 2; r₁, r₂ < 1). And
  ∫_{3/5}^{16} A(t) min(1, 1/t) dt ≤ ∫_{3/5}^{16} [D₊t/2 + (4/M)N⁺(2Mt)] dt, with
  (4/M)N⁺(2Mt) ≤ 4tℓ₀ for t ≤ 1 and ≤ 4t(ℓ₀ + log t) for t ≥ 1 (2Mt+2 ≤ (2M+2)t), so the
  integral is ≤ (D₊/2 + 4ℓ₀)·127.82 + 4·291.142 ≤ P₁. This part of TAIL is ≤ T₁.
* t ≥ 16. Of the ≥ N(2M) − 1 remaining factors, bound nine by (B2):
  \|J₀(2a_γt)\| ≤ (8/π)(4t/√5)^{−1/2} = c₃t^{−1/2}, and the other ≥ N(2M) − 10 by 0.733.
  With \|sin t\|/t ≤ 1/t and A(t)/t ≤ D₊t/2 + 4t(ℓ₀ + log t):
  ∫₁₆^∞ [(D₊/2 + 4ℓ₀)t^{−7/2} + 4t^{−7/2} log t] dt = (D₊/2 + 4ℓ₀)I₁′ + 4I₂′ = P₂.
  This part is ≤ (c₃⁹/π)P₂ r₂^{N(2M)−10} ≤ T₂.

So far: **−δ_q′(M) ≥ MAIN_an(M) − TAIL_an(M) for every M ≥ 300.**

**Step 7. TAIL_an(M) < MAIN_an(M) for every M ≥ 300.** All functions involved are
smooth and positive on [300, ∞). Elementary derivatives:

    D₊′ = 2/M − 8/M² ∈ (0, 2/M],   D₋′ = 2/M + 8/M² > 0,   ℓ₀′ = 1/(M+1) ∈ (0, 1/M),
    n_i′ = i(ℓ_{Y_i} − 1)/(16√3),   ℓ_{Y_i} = log(qY_i/2π) increasing in M,
    β₊′ = (ℓ + 1) + ℓ/34 = (35/34)ℓ + 1 > 0,
    G′(β)/G(β) = −H(β)/G(β),   H(β) = ∫₀^{1/8} t⁴e^{−βt²}dt,   H/G ≤ 3/(2β).

(The last: H(T)/G(T) for the integrals truncated at T is increasing in T, because
(H/G)′ = T²e^{−βT²}(T²G − H)/G² ≥ 0, and its limit at T = ∞ is 3/(2β).)
Since ℓ_{Y₁}(300) > 1 and ℓ(300) > 1 (certified: ≥ 4.41 and ≥ 4.96), n_i′ > 0 and
increasing, P₁, P₂, D₊/2 + 4ℓ₀ are increasing, and β₊ ≥ Mℓ. Therefore, for M ≥ 300,

    (log T₁)′ = P₁′/P₁ + n₁′ log r₁ + n₂′ log r₂
              ≤ 640/(300·P₁(300)) + n₁′(300) log r₁ + n₂′(300) log r₂            =: −κ₁,
    (log T₂)′ = P₂′/P₂ + n₂′ log r₂
              ≤ 5/(300·(D₊/2 + 4ℓ₀)(300)) + n₂′(300) log r₂                        =: −κ₂,
    (log MAIN_an)′ = D₋′/D₋ − β₊′·H(β₊)/G(β₊)
              ≥ −(3/(2M))·(35/34 + 1/ℓ) ≥ −(3/600)(35/34 + 1/ℓ(300))               =: −κ₃

(P₁′ = 64D₊′ + 512ℓ₀′ ≤ 640/M; P₂′ = (D₊′/2 + 4ℓ₀′)I₁′ ≤ 5I₁′/M and P₂ ≥ (D₊/2+4ℓ₀)I₁′).
Certified values:

| q | κ₁ | κ₂ | κ₃ | κ₁ − κ₃ | κ₂ − κ₃ |
|---|---|---|---|---|---|
| 3 | 0.152651 | 0.091592 | 0.006154 | 0.146497 | 0.085438 |
| 4 | 0.164247 | 0.098062 | 0.006099 | 0.158148 | 0.091963 |

Hence each of T₁/MAIN_an and T₂/MAIN_an is strictly decreasing on [300, ∞) (indeed
≤ its value at 300 times e^{−(κ_i−κ₃)(M−300)}), so

    TAIL_an(M)/MAIN_an(M) ≤ TAIL_an(300)/MAIN_an(300) =: R_q(300),

and the certified values are

| q | n₁(300) | n₂(300) | MAIN_an(300) | T₁(300) | T₂(300) | R_q(300) |
|---|---|---|---|---|---|---|
| 4 | 28.7597 | 73.0264 | 1.352415×10⁻⁵ | 3.7325×10⁻¹³ | 6.0873×10⁻⁹ | 4.5013×10⁻⁴ |
| 3 | 25.6454 | 66.7979 | 1.404338×10⁻⁵ | 1.1657×10⁻¹¹ | 4.0852×10⁻⁸ | 2.9098×10⁻³ |

Both are < 1. With Step 6: −δ_q′(M) ≥ (1 − R_q(300))·MAIN_an(M) > 0 for all M ≥ 300. ∎

The finitely many numerical facts Step 7 needs are exactly: R_q(300) < 1;
κ₁ > κ₃ and κ₂ > κ₃; ℓ_{Y₁}(300) > 1. They involve only π, log, exp, √, erf.

## 5. The certificate

`monotonicity/analytic-leg/certify_analytic_leg.py` (python-flint 0.9.0, 50 digits;
every decimal constant enters as an exact decimal ball; a check passes only if the
ball comparison is certain). Output: `certify_analytic_leg.out`, 90 PASS, 0 FAIL.

| Section | Content | Role |
|---|---|---|
| [A] | J₀(2.4) > 0.0024; J₁ > 0 on (0, 2√2]; √P(2√2) < 0.446; J₀(1.0733) < 0.733; J₀(1.697) < 0.410; thresholds 12/(5√5), 6√2/5; brute-force suprema; (J2), (J3), φ ≤ 0.517z; c₃ < 1.9041; 0.9973; 0.57; 0.89; 127.82; 291.142; L(½,χ) > 0 | constants of Steps 0, 5, 6 |
| [A′] | 0.189, 0.18, 0.226, 1.133, 0.90, 0.9998, 0.003 and the corrected sandwich constant (R3) | `thm:monotone-large`; not needed for the Proposition |
| [B] | R_q(300) < 5×10⁻⁴ resp. 3×10⁻³; the script's weakened variant; n₂(299) − 10 ≥ 8 | Step 7, Step 3 |
| [C] | κ₁, κ₂, κ₃; the paper's displayed constants −0.15, −0.09, 0.001, 0.007, 4.4, 5.1, 27.8, 13.9; exact log-derivatives at 300 | Step 7 |
| [D] | 1047 cells of [300, 10⁷]: ratio < 1 by monotone enclosures; exact (log T_i/MAIN_an)′ enclosed on each cell ≤ −0.14, −0.08; drop of log(T_i/MAIN_an) across each cell ≥ (κ_i − κ₃)·width | independent confirmation of Step 7 |
| [E] | D₋ ≤ D ≤ D₊, 2Mℓ ≤ σ² ≤ 2Mℓ+1.4, the S₄ bounds and (C2) against Arb values of the exact closed forms at 17 weights in [10, 10⁵] (max \|M(D − 2ℓ − 2)\| = 0.117) | sanity of Step 1 |
| [F] | n₁ ≤ N(M), n₂ ≤ N(2M) on [35, γ_max/2] (slack ≥ 10 resp. 26), N ≤ N⁺, against the certified zero lists | sanity of Step 2; not used |
| [G] | comparison with the direct certificate | §6 |

Steps 0–6 are hand proofs (re-derived above); the script certifies their numerical
constants, it does not re-prove the lemmas. Step 7 is a hand proof whose only
numerical inputs are certified in [B], [C]; [D] verifies its conclusion on [300, 10⁷]
without using it.

Reproduce: `python3 monotonicity/analytic-leg/certify_analytic_leg.py` (1.3 s);
`sh monotonicity/analytic-leg/sharp_cells.sh` (≈ 2.5 min, regenerates
`sharp_cells.txt`); `python3 monotonicity/analytic-leg/sanity_float.py`.

## 6. Cross-check against the direct certificate of δ′

The analytic claim is δ_q′(M) ≤ U(M) := −(MAIN_an − TAIL_an)(M). The direct
certificate gives lo ≤ δ_q′ ≤ hi on cells. A contradiction would be U < lo.

| q | M | stretch cell enclosure [lo, hi] | sharp cell enclosure (new, tol 0.003) | U(M) | sup_cell MAIN_an / \|hi_sharp\| |
|---|---|---|---|---|---|
| 4 | 300 | [−1.457079, −1.345098]×10⁻⁵ | [−1.412512, −1.404063]×10⁻⁵ | −1.351806×10⁻⁵ | 0.9632 |
| 4 | 310 | [−1.388129, −1.281489]×10⁻⁵ | [−1.339273, −1.331260]×10⁻⁵ | −1.282117×10⁻⁵ | 0.9633 |
| 4 | 320 | [−1.326022, −1.224122]×10⁻⁵ | [−1.272020, −1.264410]×10⁻⁵ | −1.217905×10⁻⁵ | 0.9633 |
| 3 | 300 | [−1.514954, −1.398567]×10⁻⁵ | [−1.466235, −1.457465]×10⁻⁵ | −1.400252×10⁻⁵ | 0.9636 |
| 3 | 310 | [−1.442886, −1.332072]×10⁻⁵ | [−1.389814, −1.381499]×10⁻⁵ | −1.329571×10⁻⁵ | 0.9636 |
| 3 | 320 | [−1.374321, −1.268728]×10⁻⁵ | [−1.319656, −1.311763]×10⁻⁵ | −1.263404×10⁻⁵ | 0.9636 |
| 3 | 350 | [−1.189839, −1.098446]×10⁻⁵ | [−1.140353, −1.133531]×10⁻⁵ | −1.092353×10⁻⁵ | 0.9637 |

* **Consistent everywhere.** On all 7 (q=4) and 16 (q=3) stretch cells meeting
  [300, 320] resp. [300, 350], sup_cell(MAIN_an − TAIL_an) ≤ 0.933·\|lo\|: the analytic
  bound never claims more negativity than the certificate allows.
* The stretch cells have enclosure radius ≈ 4 %, about the same as the gap between U
  and the truth, so they cannot by themselves show hi ≤ U (U/hi ranges over
  0.994–1.005 at the tabulated points). I therefore ran the same prover
  (`cert/certify_derivative.py cell --tol 0.003`) on seven tiny cells
  ([M, M+0.001], or [M−0.001, M] at the range ends). There
  **hi ≤ −sup_cell MAIN_an ≤ U** in every case (last column < 1): the direct
  certificate independently implies the analytic inequality at these weights, with
  δ′ more negative than U by at least 3.6 % and at most 4.5 % (U/lo_sharp ≥ 0.955).
* That gap is the expected one: MAIN_an/[D/(2√(2π)σ³)] = 0.9597 at (q, M) = (4, 300),
  from the factor 0.9973, D₋/D and (β₊/(σ²/2))^{−3/2}.

## 7. Hypotheses of the conclusion

* **GRH(χ_q)**, and nothing else that is unproved. It is used to make the ordinates
  real (the definition of a_γ, δ_q, Φ) and, through positivity of the terms
  8M²/(M²+γ²), in the variance identity and both zero-count bounds.
* **L(½, χ_q) ≠ 0**: needed (Ξ(0) ≠ 0 in the Hadamard product), but unconditional for
  q = 3, 4.
* **Zero lists: not used at all.** No ordinate, no completeness statement, no
  simplicity of zeros (multiplicities are allowed), no explicit zero-count theorem
  from the literature. The counts N(M) ≥ n₁, N(2M) ≥ n₂, N(T) ≤ N⁺(T) are derived
  from the closed forms of σ² and S₄.
* **LI: not used.** δ_q is the probability P(1 + Σ 2a_γ cos θ_γ > 0) with independent
  phases by definition; LI is only what identifies it with the logarithmic density of
  the race.
* Classical inputs taken from the literature: functional equation and root number of
  real primitive characters, Hadamard factorisation (order 1), reality of the zeros
  of J₀ and its product formula, Binet's formula for ψ, Fourier inversion for an
  integrable characteristic function.
* Software: correctness of Arb for π, log, exp, √, erf (Step 7) and J₀, J₁, ζ, ψ at
  points (constants that also have the paper's alternating-series proofs).

## 8. What I did not do

* I did not audit `cert/certify_derivative.py` or the zero-list certification; they
  enter only the cross-check of §6, not the proof.
* The cell sweep [D] stops at 10⁷; beyond that the Proposition rests on Step 7 alone
  (which covers all M ≥ 300 anyway).
* `sanity_float.py` is vacuous for `eq:Abound` at M = 300 (both sides underflow); it
  is informative at M = 40 and 100.
* The paper is unchanged. Suggested edits: R1 (state the three monotonicity facts,
  cite the Arb run), R2 (drop `lem:Nall` from `lem:diffint`(c) for odd χ), R3
  (0.447/0.9510 → 1/(2√5)/0.9876 in the sandwich step), m1–m3, m8.
* I made no commit. The repository has an automatic "Work in progress" committer that
  picked up some files of `monotonicity/analytic-leg/` while I was working (commit
  `082a69a`); that was not my action.
