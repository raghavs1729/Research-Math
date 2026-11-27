# Monotonicity of the weighted Chebyshev bias: the coupling mechanism

Working notes (2026-11-24). Notation is that of `paper/weighted-races.tex`:
`M = m + 1/2 > 0`, positive ordinates `γ` of `L(s,χ)` (χ real primitive, GRH,
`L(1/2,χ) ≠ 0`), and under LI

    X_M = 1 + Σ_{γ>0} c_γ(M) cos θ_γ ,   c_γ(M) = 2a_γ = 4M / sqrt(M² + γ²),

θ_γ i.i.d. uniform, `δ(m) = P(X_M > 0)`, `σ_M² = Σ c_γ²/2 = 4M (log Λ)'(m+1)`.
The paper's Conjecture `conj:monotone` says `δ_q` is strictly decreasing on all
of `(−1/2, ∞)` for q = 3, 4; `thm:monotone-large` proves it for
`M ≥ 17.317` (q=4) and `M ≥ 19.119` (q=3).  The paper says the family "is not
a dilation" and that near `m = −1/2` only the Chernoff envelope is available.

## The observation

Divide by M.  For a FIXED phase vector θ,

    φ(M) := X_M / M = 1/M + Σ_γ 4 cos θ_γ / sqrt(M² + γ²).

The bias term `1/M` decreases in M at rate `1/M²`.  Each oscillating
coefficient `4/sqrt(M²+γ²)` also DEcreases in M, at rate
`4M/(M²+γ²)^{3/2}`, and these rates are ABSOLUTELY SUMMABLE over the zeros
(they are O(γ^{-3})), although the coefficients themselves are not summable.
So for 0 < M < M',

    φ(M) − φ(M') = (1/M − 1/M') + Σ_γ cos θ_γ · Δ_γ ,
    Δ_γ := 4/sqrt(M²+γ²) − 4/sqrt(M'²+γ²) > 0,
    Σ_γ Δ_γ = ∫_M^{M'} Q(u) u^{-2} du ,
    Q(u) := Σ_{γ>0} 4u³ / (u² + γ²)^{3/2}      (finite, continuous, strictly increasing, Q(0+) = 0, Q(∞) = ∞).

Hence, pathwise (for every θ for which the series converges),

    X_M / M − X_{M'} / M'  ≥  η(M,M') := ∫_M^{M'} (1 − Q(u)) u^{-2} du .

Let `M_0 = M_0(χ)` be the unique root of `Q(M_0) = 1`.  For `M < M' ≤ M_0`,
`η > 0`, so `{X_{M'} > 0} ⊂ {X_M > 0}` SURELY, and
`δ(M) − δ(M') ≥ P(−η M' < X_{M'} ≤ 0) > 0`.

**Theorem B (coupling monotonicity; GRH + NCZ + LI).** For every real primitive
χ, `m ↦ δ_χ(m)` is strictly decreasing on `(−1/2, M_0(χ) − 1/2]`.
Numerically `M_0(χ_4) ≈ 3.689` (m ≤ 3.19) and `M_0(χ_3) ≈ 4.484` (m ≤ 3.98).
Q<1 is sharp for the pathwise statement (Σ c_γ = ∞).
Heuristically `Q_χ(M) ≈ (2/π) M log q` for large q, so `M_0 ≍ 1/log q`: the
theorem lives exactly at the "critical corner" scale of `rem:joint-corner`.

Monte Carlo (511 zeros of χ_4, 4e5 samples): zero nesting violations for all
tested pairs with M' ≤ 3.6; violations appear at once beyond (1 at 3.6→5,
1382 at 5→8).

## LI-free version (same x, GRH only)

In the true race the phases are `γ log x − arg(M + iγ)`, so the M-derivative
of a term has modulus `4/(M²+γ²)` (still summable):

    E_M(x)/M − E_{M'}(x)/M' = (1/M − 1/M') − Σ_{γ>0} 4 Re[ x^{iγ} (1/(M+iγ) − 1/(M'+iγ)) ] + o(1),

an ABSOLUTELY CONVERGENT explicit formula (this is the exact identity
`D_m(x) − x^{-κ} D_{m'}(x) = κ ∫_2^x t^{-κ-1} D_{m'}(t) dt`, κ = m' − m, in
disguise: a Riesz-type mean of the heavier race).  So if

    Λ(M,M') := 1/(M M') − Σ_{γ>0} 4 / ( |M+iγ| |M'+iγ| )  >  0

(by Cauchy–Schwarz it suffices that `σ_M σ_{M'} < 2`), then for all
sufficiently large x:  `D_{m'}(x) ≥ 0  ⟹  D_m(x) > 0`.

**Theorem A (eventual nesting of lead sets; GRH only, no LI).**  Leads are
inherited downward in the weight.  Consequences without LI: upper and lower
logarithmic densities are monotone on `σ_M² ≤ 2` (q=4: M ≤ 1.824, m ≤ 1.32;
q=3: M ≤ 2.128, m ≤ 1.63); e.g. for q = 4, `σ_0 σ_1 = 0.461 < 2`, so every
lead of class 1 mod 4 in Chebyshev's classical race (x large) is also a lead
in the race weighted by p.  Under LI the condition is sharp (Kronecker).

## What remains for the full conjecture (q = 3, 4)

The compact middle range `M ∈ [M_0, M_1]`:  [3.689, 17.317] for q=4,
[4.484, 19.119] for q=3.  There `δ'` is of size 1e-1 … 2e-3, so a certified
(Arb ball arithmetic) evaluation of the differentiated Gil–Pelaez integral
over M-cells closes it.  Exact identity for the derivative, with
`Z = Σ ε_γ c_γ cos θ_γ`, `ε_γ = M²/(M²+γ²)`:

    δ'(M) = −(f_{S}(−1)/M) · E[ 1 + Z | S = −1 ],   S = X_M − 1,  |Z| ≤ Q(M).

Reduction worth knowing: by Abel summation, monotonicity for ALL M and ALL χ
would follow from `E[T_k | S = −1] ≤ 0` for every tail sum
`T_k = Σ_{γ>γ_k} c_γ cos θ_γ` (true whenever the head sum is symmetric
unimodal, by Anderson's inequality).
