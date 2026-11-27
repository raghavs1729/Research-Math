# Theorem U1 — closed form of an "unidentified" ζ(5)-level continued fraction of the Ramanujan Machine

## 0. Statement, source, status

**Source (read 2026-11-24).** R. Elimelech, O. David, C. De la Cruz Mengual, R. Kalisch, W. Berndt, M. Shalyt,
M. Silberstein, Y. Hadad, I. Kaminer, *Algorithm-assisted discovery of an intrinsic order among mathematical
constants*, arXiv:2308.11829v2 (PNAS 121 (2024) e2321440121), §3.6 "Results for ζ(5) and the quest to prove its
irrationality", Table 1 ("Sporadic results suspected as related to ζ(5), found via the Distributed Factorial
Reduction algorithm for polynomial continued fractions of degree 5"), fourth row:

    a_n = 8(n^5 + (n+1)^5) − 15(n^3 + (n+1)^3) + 9(2n+1),     b_n = −64 n^10,     Formula: 1.20426...

with the comment: *"Three of the discovered formulas in the table did not fit any combination of zetas (in
numerical tests using PSLQ) and are thus presented with their numerical values."*
The convention of the paper (and of this file) is

    PCF(a,b) = a_0 + b_1/(a_1 + b_2/(a_2 + b_3/(a_3 + ...))) = lim p_n/q_n,
    u_n = a_n u_{n−1} + b_n u_{n−2}  (n ≥ 1),   (p_{−1}, p_0) = (1, a_0),   (q_{−1}, q_0) = (0, 1).

(The printed "1.20426..." is not a truncation of the true value 1.2042534...: the fraction converges like `c/n`,
and 1.20426 is what convergents of depth about 4·10^4–1.1·10^5 give; see `audit/a5_source_value.py`.)

**Openness / novelty (searched 2026-11-24..26).** No closed form is given in the paper (both arXiv versions; the
row is numerical in Table 1 and in App. C Table 4, whose caption says "we did not find a closed-form expression
for them"). Two independent searches (`../STATUS-MAP-RM.md` §2.1 and `AUDIT-2.md`) found no identification of
the constant, of the continued fraction or of the limit of the recurrence in: Cohen's databases (arXiv:2409.06086,
2607.06581) and Cohen arXiv:2212.01095 (which has the degree-3 analogues `28ζ(3) ± π^3` from `ψ''(1/4), ψ''(3/4)`
but only Euler-type fractions in degree 5), Yamamoto arXiv:2403.09729, about twelve later papers of the Ramanujan
Machine group (up to arXiv:2610.27799 and the Ramanujan Challenge papers), Shvets's 2026 papers, Chamberland–Straub,
Golyshev, Almkvist's tables of 5th-order equations, the project's site and GitHub, `openai/math`, and the list of
papers citing arXiv:2308.11829. Not reachable: the PNAS supplement. **Prior art on the sequence:** OEIS A397472
(Shvets, June 2026) lists `64^n w_n = 1, 16, 496, 19456, ...` with the binomial formula and the same recurrence; it
does not mention the continued fraction, the second solution, the limit, ζ(5) or π^5. So what is new here is the
link between the Ramanujan Machine fraction and that sequence, and the evaluation of the limit. "New" means "not
found after these searches".

> **Theorem U1.** The continued fraction above converges, and
>
>     PCF(a,b) = 96 / (372 ζ(5) − π^5) = 1.2042534211559321151983839420481208612545...
>
> Equivalently `1/PCF = (372 ζ(5) − π^5)/96 = (31/8) ζ(5) − π^5/96`.
> Moreover, with `w_n = Σ_{k=0}^{n} C(2k,k)^3 C(2n−2k,n−k)^3 / 64^n = [x^n] 3F2(1/2,1/2,1/2;1,1;x)^2`
> (`= [x^n] 2F1(1/4,1/4;1;x)^4`), the numerators are `p_{n−1} = 8^n n!^5 w_n`, and
>
>     Σ_{n≥0} 1/((n+1)^5 w_n w_{n+1}) = (372 ζ(5) − π^5)/12.
>
> (The equality `[x^n] 2F1(1/4,1/4;1;x)^4 = Σ_k C(2k,k)^3 C(2n−2k,n−k)^3/64^n` is Clausen's formula
> `2F1(a,b;a+b+1/2;x)^2 = 3F2(2a,2b,a+b;2a+2b,a+b+1/2;x)` with `a = b = 1/4`.) Convergence is only polynomial:
> `p_n/q_n − PCF ~ c/n`, `c = 4608 π^9/(Γ(1/4)^8 (372ζ(5) − π^5)^2) = 0.72394...` (Remark 2; not needed for the
> theorem).

So the constant *is* a ζ(5)-constant, but with a π^5 term that a PSLQ search over "combinations of zetas" does not
contain. The proof identifies the generating functions of p_n and q_n with a modular form and an Eichler integral
of a weight-6 Eisenstein series, and evaluates the limit by a Mellin transform (the mechanism of Ramanujan's
formula for ζ(2k+1)).

**Verification status.** Proof below: complete, classical ingredients cited in §1. Refereed twice by independent
agents: `audit/` (independent scripts: own Richardson evaluation of the fraction, agreement with the closed form to
~498 digits; independent numerical re-verification of Lemma A, (C1)–(C8), Lemma B, (4.3), (5.1)) and `AUDIT-2.md`
(line-by-line proof-read: "complete after listed fixes" — the fixes are incorporated in this version). Machine checks:
`verify/v1_operator_identity.py` (Lemma A, exact symbolic), `verify/v2_series_checks.py` (exact rational
arithmetic: recurrence, `p_{n−1} = 8^n n!^5 w_n` and `q_{n−1} = 8^n n!^5 v_n` for n ≤ 40, the theta/Eisenstein
identity to O(q^300) and as a polynomial identity), `verify/v3_numerics.py` (the continued fraction evaluated by an
independent method agrees with `96/(372ζ(5) − π^5)` to 300 digits; moments `M(4) = 1/4`, `M(5) = 186ζ(5)/π^5` to
50 digits). Not formalized in Lean (the proof uses modular forms and analytic continuation of Dirichlet series).

---

## 1. Notation and the classical facts used

`θ = x d/dx`. For `t > 0` put `q = e^{−πt}` (i.e. `τ = it`, `q = e^{iπτ}`) and

    θ_2 = 2 Σ_{n≥0} q^{(n+1/2)^2},   θ_3 = 1 + 2 Σ_{n≥1} q^{n^2},   θ_4 = 1 + 2 Σ_{n≥1} (−1)^n q^{n^2},
    λ = θ_2^4/θ_3^4,   E_6(q) = 1 − 504 Σ_{n≥1} σ_5(n) q^n.

All of the following are classical (Jacobi); references: Whittaker–Watson, *Modern Analysis*, ch. 21–22;
NIST DLMF ch. 15, 20, 23; Zagier, *Elliptic modular forms and their applications* (in "The 1-2-3 of modular forms"),
§3.1.

* (C1) `θ_3^4 = θ_2^4 + θ_4^4`; hence `1 − λ = θ_4^4/θ_3^4`.
* (C2) `t ↦ λ(it)` is a decreasing bijection `(0,∞) → (0,1)` with `λ(i) = 1/2`, and
  `2F1(1/2,1/2;1;λ(it)) = θ_3(it)^2` (i.e. `K(k) = (π/2) θ_3^2`, `k^2 = λ`).
* (C3) `q dλ/dq = λ(1−λ) θ_3^4`, i.e. `dλ(it)/dt = −π λ(1−λ) θ_3^4`.
* (C4) (imaginary transformation) `θ_3(i/t)^2 = t θ_3(it)^2`, `θ_2(i/t)^2 = t θ_4(it)^2`, `θ_4(i/t)^2 = t θ_2(it)^2`.
* (C5) with `q` replaced by `q^2` on the left (`E_6(q^2)` is "E_6 at τ"):
  `E_6(q^2) = (1/2)(θ_2^4 + θ_3^4)(θ_3^4 + θ_4^4)(θ_4^4 − θ_2^4)` (Zagier, loc. cit., section on theta series, where `E_4` and `E_6` are
  expressed through the Thetanullwerte; classical, cf. also DLMF §23.6(i) with §20.7; checked against the q-series to `O(q^60)` in `v2_series_checks.py`); and the Landen/Gauss transformations
  `θ_3(2τ)^2 = (θ_3^2 + θ_4^2)/2`, `θ_4(2τ)^2 = θ_3 θ_4`, `θ_2(2τ)^2 = (θ_3^2 − θ_4^2)/2`,
  `θ_3(τ/2)^2 = θ_3^2 + θ_2^2`, `θ_4(τ/2)^2 = θ_3^2 − θ_2^2`, `θ_2(τ/2)^2 = 2 θ_2 θ_3`.
* (C6) (quadratic transformation, DLMF 15.8.18) `2F1(1/4,1/4;1;4z(1−z)) = 2F1(1/2,1/2;1;z)` for `0 ≤ z ≤ 1/2`.
* (C7) `Σ_{n≥1} σ_5(n) n^{−s} = ζ(s) ζ(s−5)` (Re s > 6); `ζ(0) = −1/2`, `ζ(−1) = −1/12`, `ζ(4) = π^4/90`.
* (C8) (Abelian lemma) If `α_n > 0`, `Σ α_n x^n` and `Σ β_n x^n` converge for `|x| < 1`, `Σ α_n x^n` tends to `+∞` as
  `x → 1−`, and `β_n/α_n → ℓ ∈ [0,∞]`, then `Σ β_n x^n / Σ α_n x^n → ℓ` as `x → 1−`. (Given ε, split the series at an index beyond
  which `|β_n − ℓ α_n| ≤ ε α_n`; the finite part is bounded and the denominator tends to infinity. Same for ℓ = ∞.)

Constants of the problem:

    ã(n) := a_n/8 = 2T^5 + (5/4) T^3 + (1/16) T,   T = n + 1/2,
    L := θ^5 − 2x (θ+1/2) ((θ+1/2)^4 + (5/8)(θ+1/2)^2 + 1/32) + x^2 (θ+1)^5.

(The first identity is a direct expansion.) For a power series, `L Σ c_n x^n = Σ_n (n^5 c_n − ã(n−1) c_{n−1} +
(n−1)^5 c_{n−2}) x^n`, since `2(θ+1/2)((θ+1/2)^4 + ...)` multiplies `x^{n−1}` by `ã(n−1)`.

    F := θ_2^4 θ_4^4 (θ_4^4 − θ_2^4) = Σ_{n≥1} a°_n q^n = 16 q − 512 q^2 + 3904 q^3 − ...,
    Φ := (1/2) Σ_{n≥1} a°_n q^n / n^5.

## 2. Two lemmas

**Lemma A (operator identity).** Let `u` be a solution on an interval `I ⊂ (0,1)` of
`x(1−x) u'' + (1 − 3x/2) u' − u/16 = 0` without zeros, let `T` be a function on `I` with
`T' = 1/(x √(1−x) u^2)`, and let `h` be five times differentiable on `T(I)`. Then

    L[ u^4 · h(T) ] = h^{(5)}(T) / ( √(1−x) · u^6 ).

*Proof.* This is a computer-algebra identity: a finite computation in the differential ring generated by
`x, √(1−x), u, u', h^{(k)}(T)`, in which `u, u', h^{(k)}` are free symbols, `u''` is eliminated by the hypergeometric
equation and `d/dx h^{(k)}(T) = h^{(k+1)}(T) T'`. Done exactly by `verify/v1_operator_identity.py` (sympy; output `0`);
re-verified independently by numerical Taylor expansion at six points in `audit/a2_lemmaA_independent.py`.
(Conceptually it is Bol's identity: with `W = 1/(x√(1−x))` the Wronskian of the hypergeometric equation and
`T' = W/u^2`, the functions `u^4 T^j`, `j ≤ 4`, span the kernel of the symmetric fourth power `L`, whose leading
coefficient is `x^5(1−x)^2`, and `L[u^4 h(T)] = x^5(1−x)^2 T'^5 u^4 h^{(5)}(T)`.) ∎

(In particular `L[u^4 T^j] = 0` for `j ≤ 4`: `L` is the symmetric fourth power of the hypergeometric operator of
`2F1(1/4,1/4;1;x)`.)

**Lemma B (Eisenstein).** `a°_n = 16 ( σ_5(n) − 65 σ_5(n/2) + 64 σ_5(n/4) )` (with `σ_5(r) = 0` for
`r ∉ Z`). Equivalently `63 F = −2 E_6(q) + 130 E_6(q^2) − 128 E_6(q^4)`. Moreover `F(i/t) = −t^6 F(it)`.

*Proof.* Put `A = θ_3^4, B = θ_4^4, C = θ_2^4` (at τ), so `A = B + C` by (C1). By (C5),

    E_6(q^2) = (1/2)(C + A)(A + B)(B − C),
    E_6(q^4) = (1/64)(A + B)(36AB − (A + B)^2),          [θ_3(2τ)^4 = (θ_3^2+θ_4^2)^2/4 etc.; θ_3^2 θ_4^2 = (AB)^{1/2} cancels]
    E_6(q)   = (A + C)((A + C)^2 − 36AC),                [θ_3(τ/2)^4 = (θ_3^2+θ_2^2)^2 etc.]

and `63 C B (B − C) = −2 E_6(q) + 130 E_6(q^2) − 128 E_6(q^4)` is then a polynomial identity in `B, C`
after substituting `A = B + C` (checked by `verify/v2_series_checks.py`, which also confirms the three theta
expressions and the final expansion of `F` against the q-series to `O(q^60)` and `O(q^300)`). Comparing coefficients
gives `a°_n`. (Alternative proof not using the Landen formulas: read as power series in `q = e^{2πiτ'}` (`τ = 2τ'`),
`Σ a°_n q^n`, `E_6(q)`, `E_6(q^2)`, `E_6(q^4)` are the expansions of modular forms in `M_6(Γ_0(4))` — the theta
functions with nome `e^{iπτ}` are modular on `Γ(2)`, which is conjugate to `Γ_0(4)` by `τ ↦ 2τ` — and the Sturm bound
for `M_6(Γ_0(4))` is `6·[SL_2(Z):Γ_0(4)]/12 = 3`, so agreement of the coefficients of `q^0, ..., q^3` suffices; the
script checks them up to `q^300`.) The transformation law follows from (C4): `θ_2^4, θ_4^4` are interchanged up to the factor `t^2`,
so `F(i/t) = t^6 θ_4^4 θ_2^4 (θ_2^4 − θ_4^4) = −t^6 F(it)`. ∎

## 3. The generating functions of p_n and q_n

Let `u(x) = 2F1(1/4,1/4;1;x) = Σ u_n x^n`, `u_n = ((1/4)_n/n!)^2 > 0`, and

    P(x) := u(x)^4 = Σ_{n≥0} w_n x^n      (|x| < 1),   w_0 = 1,  w_n > 0.

(By Clausen's formula (DLMF 16.12.2) `u^2 = 3F2(1/2,1/2,1/2;1,1;x) = Σ C(2k,k)^3 (x/64)^k`, which gives the binomial-sum form of
`w_n` in the statement; this is not used in the proof. Checked exactly for n ≤ 40.)

For `x ∈ (0,1)` let `t = t(x) ∈ (1,∞)` be defined by `x = 4λ(it)(1 − λ(it))`. By (C2), (C3), `t ↦ x` is a
real-analytic decreasing bijection `(1,∞) → (0,1)` with `dx/dt ≠ 0`, so `t(x)` is real-analytic on `(0,1)` (inverse
function theorem), and `√(1−x) = 1 − 2λ`. Consequently `P̃ := θ_3(it(x))^8` and `Q̃ := θ_3(it(x))^8 Φ(it(x))` below
are real-analytic on `(0,1)`. By (C6) and (C2),

    u(x) = θ_3(it)^2,     P(x) = θ_3(it)^8     (0 < x < 1).                                        (3.1)

Put `T(x) := −π t(x)`. By (C3), `dx/dt = 4(1−2λ)·(−π λ(1−λ) θ_3^4) = −π x √(1−x) u^2`, hence
`T' = 1/(x √(1−x) u^2)`: the hypothesis of Lemma A holds with this `u` and `T`.

* `h ≡ 1` gives `L P = 0` on `(0,1)`, hence (power series) `n^5 w_n = ã(n−1) w_{n−1} − (n−1)^5 w_{n−2}` for all
  `n ≥ 1` (with `w_{−1} = 0`).
* Take `h(T) := Φ` evaluated at `t = −T/π`, i.e. `h(T) = (1/2) Σ a°_n e^{nT}/n^5`. Then `h^{(5)}(T) = (1/2) Σ a°_n
  e^{nT} = F/2`, and Lemma A gives, for `Q̃(x) := θ_3(it)^8 Φ(it)`,

      L Q̃ = F / (2 (1−2λ) θ_3^{12}) = θ_2^4 θ_4^4 / (2 θ_3^8) = λ(1−λ)/2 = x/8        (0 < x < 1).        (3.2)

`Q̃` is the restriction to a small interval `(0,ε)` of a function analytic at `x = 0`: `λ = (1 − √(1−x))/2` is
analytic at 0; `λ = θ_2^4/θ_3^4 = 16q − 128q^2 + ...` is a convergent power series in `q` with non-zero linear
coefficient, hence locally invertible at `q = 0`; for large `t` the pair `(q, λ) = (e^{−πt}, λ(it))` lies on the graph
of this local inverse (both are small and `λ = λ(q)`), so `q = e^{−πt}` is an analytic function of `λ` near 0 with
`q = λ/16 + O(λ^2)`; and `θ_3^8 Φ = 8q + O(q^2)` is a convergent power series in `q` without constant term. So `Q̃ = Σ_{n≥1} v_n x^n` near 0, with `v_1 = 1/8` and, by (3.2),

    n^5 v_n = ã(n−1) v_{n−1} − (n−1)^5 v_{n−2}  (n ≥ 2),      v_0 = 0,  v_1 = 1/8.

**Proposition 1.** `p_{n−1} = 8^n n!^5 w_n` and `q_{n−1} = 8^n n!^5 v_n` for all `n ≥ 0`.

*Proof.* For `n = 0, 1`: `p_{−1} = 1 = w_0`, `p_0 = a_0 = 8 ã(0) = 8 w_1`; `q_{−1} = 0 = v_0`, `q_0 = 1 = 8 v_1`.
If `c_n` satisfies `(n+1)^5 c_{n+1} = ã(n) c_n − n^5 c_{n−1}` then `U_{n−1} := 8^n n!^5 c_n` satisfies
`U_n = 8^{n+1} n!^5 (ã(n) c_n − n^5 c_{n−1}) = a_n U_{n−1} − 64 n^{10} U_{n−2}`. ∎
(Checked exactly for n ≤ 40 in `v2_series_checks.py`, parts (a) and (c).)

## 4. Convergence, and reduction to a boundary value

**Casoratian.** `D_n := w_n v_{n+1} − w_{n+1} v_n` satisfies `(n+1)^5 D_n = n^5 D_{n−1}` for `n ≥ 1` (multiply the two
recurrences crosswise), and `D_0 = 1/8`. Hence `D_n = 1/(8 (n+1)^5)` and

    v_{n+1}/w_{n+1} − v_n/w_n = 1/(8 (n+1)^5 w_n w_{n+1}) > 0      (n ≥ 0).                         (4.1)

So `v_n > 0` for `n ≥ 1`, `q_n ≠ 0`, and `q_n/p_n = v_{n+1}/w_{n+1}` increases to
`ℓ := Σ_{n≥0} 1/(8(n+1)^5 w_n w_{n+1}) ∈ (0,∞]`.
*Lower bound for w_n.* `u^3` has non-negative coefficients and constant term 1, so `w_n ≥ u_n`; and
`(1/4)_n/n! = Π_{j≤n}(1 − 3/(4j)) ≥ (1/4) Π_{2≤j≤n}(1 − 1/j)^{3/4} = n^{−3/4}/4` (Bernoulli), so
`w_n ≥ n^{−3/2}/16` for `n ≥ 1`. Therefore `1/(w_n w_{n+1}) ≤ 256 (n+1)^3`, the series for ℓ converges
(`ℓ ≤ 1/(8 w_1) + 32 ζ(2)`), the continued fraction **converges** to `1/ℓ ∈ (0,∞)`, and `Σ v_n x^n` has radius of
convergence ≥ 1. Define `Q(x) := Σ v_n x^n` on `|x| < 1`; `Q = Q̃` near 0, hence on `(0,1)` (both are real-analytic
there). Also `Σ n w_n = ∞`, i.e. `P'(x) → +∞` as `x → 1−`.

Since `(n+1) v_{n+1} / ((n+1) w_{n+1}) → ℓ`, the Abelian lemma (C8) applied to `P'` and `Q'` gives

    ℓ = lim_{x→1−} Q'(x)/P'(x).                                                                   (4.2)

**The limit.** Write `g(t) = θ_3(it)^8`, `φ(t) = Φ(it) = (1/2) Σ a°_n e^{−πnt}/n^5`; both are real-analytic on
`(0,∞)`, and `P = g`, `Q = gφ` as functions of `t ∈ (1,∞)`. Since `dx/dt ≠ 0` there, and `g_t = P'(x) dx/dt ≠ 0`
(`P' > 0` on `(0,1)` because all `w_n > 0`),

    Q'(x)/P'(x) = (gφ)_t / g_t = φ + g φ_t / g_t.

By (C4), `g(1/t) = t^4 g(t)`; differentiating at `t = 1` gives `−g_t(1) = 4 g(1) + g_t(1)`, i.e.
`g_t(1) = −2 g(1) ≠ 0`. Letting `t → 1+` (`x → 1−`) in (4.2):

    ℓ = φ(1) − φ_t(1)/2 = Σ_{n≥1} a°_n e^{−πn} (2 + πn) / (4 n^5).                                   (4.3)

## 5. Evaluation

For `n ≥ 1`, `∫_1^∞ (t−1)^3 (t+1) e^{−πnt} dt = e^{−πn} (24/(πn)^5 + 12/(πn)^4) = 12 e^{−πn} (2 + πn)/(πn)^5`.
As `|a°_n| ≤ 2080 σ_5(n)`, termwise integration is legitimate and (4.3) becomes

    ℓ = (π^5/48) ∫_1^∞ (t−1)^3 (t+1) F(it) dt = (π^5/48) ( M(5) − 2 M(4) ),                          (5.1)
    M(s) := ∫_1^∞ (t^{s−1} − t^{5−s}) F(it) dt,

because `(t−1)^3(t+1) = (t^4 − 1) − 2(t^3 − t)`. `M` is an entire function of `s` (`F(it) = O(e^{−πt})`).
For every `s`, Lemma B (`F(i/t) = −t^6 F(it)`) gives `−∫_1^∞ t^{5−s} F(it) dt = ∫_0^1 t^{s−1} F(it) dt` (substitute
`t → 1/t`), so `M(s) = ∫_0^∞ t^{s−1} F(it) dt`; and for `Re s > 6` the series may be integrated termwise:

    M(s) = ∫_0^∞ t^{s−1} F(it) dt = Γ(s) π^{−s} Σ a°_n n^{−s}
         = 16 Γ(s) π^{−s} ζ(s) ζ(s−5) (1 − 65·2^{−s} + 64·4^{−s})          (Re s > 6),

by Lemma B and (C7) (absolute convergence, `Σ |a°_n| n^{−Re s} < ∞` for `Re s > 6`, justifies the interchange).
Both sides are meromorphic in `s`, so the identity holds wherever the right-hand side is regular, in particular at
`s = 4, 5` (there `Γ(s)`, `ζ(s)` and `ζ(s−5)` are all finite; in fact the right-hand side is entire: the pole of
`ζ(s−5)` at `s = 6` is cancelled by `1 − 65·2^{−6} + 64·4^{−6} = 0`, that of `ζ(s)` at `s = 1` by `ζ(−4) = 0`, and
those of `Γ(s)` at `s = 0, −1, −2, ...` by the zeros of `1 − 65·2^{−s} + 64·4^{−s}` at `s = 0`, of `ζ(s−5)` at
`s = −1, −3, ...` and of `ζ(s)` at `s = −2, −4, ...`):

    M(5) = 16 · 24 · π^{−5} ζ(5) ζ(0) (1 − 65/32 + 1/16) = 16·24·(−1/2)(−31/32) ζ(5)/π^5 = 186 ζ(5)/π^5,
    M(4) = 16 · 6 · π^{−4} ζ(4) ζ(−1) (1 − 65/16 + 1/4) = 16·6·(1/90)(−1/12)(−45/16) = 1/4.

Hence `ℓ = (π^5/48)(186 ζ(5)/π^5 − 1/2) = (372 ζ(5) − π^5)/96`, and the continued fraction equals `1/ℓ`. The
series identity in the statement is `8ℓ`, by (4.1). ∎

(Numerical confirmation of (4.3), (5.1), `M(4)`, `M(5)`: `verify/v3_numerics.py`.)

## 6. Remarks

1. **Why PSLQ over zeta values failed (the displayed series identity is (4.3), which is proved; the rest of this
   remark is commentary).** The weight-6 form `F` is a pure Eisenstein series for the theta group, and
   its Eichler integral has period polynomial with coefficients in `Q ζ(5) + Q π^5` — as in Ramanujan's formula
   for ζ(5). (4.3) is such a formula:
   `Σ_{n≥1} (σ_5(n) − 65σ_5(n/2) + 64σ_5(n/4)) (2 + πn) e^{−πn}/n^5 = (372 ζ(5) − π^5)/384`.
2. **Convergence rate.** `p_n/q_n − 1/ℓ ~ c/n` with `c = 1/(8 C^2 ℓ^2)`, `C = Γ(1/4)^4/(2 π^{9/2})`: by (4.1),
   `ℓ − v_n/w_n = Σ_{k≥n} 1/(8(k+1)^5 w_k w_{k+1})`, and `w_n ~ C n^{−3/2}` by singularity analysis of `P = u^4` at
   `x = 1` (`u = A − (2√π/Γ(1/4)^2)√(1−x) + O(1−x)`, `A = √π/Γ(3/4)^2`, the only singularity on `|x| = 1`; so
   `w_n ~ 4A^3 · (2√π/Γ(1/4)^2) · n^{−3/2}/(2√π) = C n^{−3/2}`). Numerically `c = 0.72394...`, in agreement with
   `n (p_n/q_n − 1/ℓ) → 0.7239` (`verify/v3_numerics.py`, `audit/a4_proof_chain.py`). The two formal solutions of the
   recurrence behave like `n^{−3/2}` and `n^{−5/2}`.
3. **The other two unidentified constants of Table 1** are *not* of this type; see `../U23-STRUCTURE.md`.
4. **How it was found.** (i) Petkovšek's algorithm shows the recurrence has no hypergeometric solution, so the
   fraction is not of Euler type; (ii) its numerator sequence was recognised as the coefficients of
   `3F2(1/2,1/2,1/2;1,1;x)^2` from the local exponents `{0,1,2,1/2,3/2}` at `x = 1` of the associated 5th-order
   ODE (an arithmetic progression, the signature of a symmetric 4th power); (iii) variation of constants gave
   `PCF = 3/∫_0^1 K(z)^3 (2K(1−z) − K(z))(1−2z) dz` (K the complete elliptic integral with parameter z), verified to
   60 digits (this integral formula is a numerically verified remark, not used in the proof); (iv) the integrand is a weight-6 Eisenstein series, which predicted the basis `{ζ(5), π^5}`; PSLQ
   then returned `(96, −372, 1)` at 60 and at 300 digits. Scripts: `../code/z5_*.py`, `../code/u1_*.py`.
