# Theorem U2 — the second "unidentified" ζ(5)-level continued fraction of the Ramanujan Machine

## 0. Statement, source, status

**Source (read 2026-11-24).** arXiv:2308.11829v2 (Elimelech et al., PNAS 121 (2024) e2321440121), §3.6, Table 1,
fifth row (also App. C, Table 4):

    a_n = 8(n^5 + (n+1)^5) − 12(n^3 + (n+1)^3) + 7(2n+1),     b_n = −64 n^10,     Formula: 2.45174...

("Three of the discovered formulas in the table did not fit any combination of zetas (in numerical tests using
PSLQ) and are thus presented with their numerical values"; Table 4: "we did not find a closed-form expression for
them"). Convention: `PCF = a_0 + b_1/(a_1 + b_2/(a_2 + ...)) = lim p_n/q_n`, `u_n = a_n u_{n−1} + b_n u_{n−2}`,
`(p_{−1},p_0) = (1,a_0)`, `(q_{−1},q_0) = (0,1)`.

Notation: `A_k = C(2k,k)^4/256^k = ((1/2)_k/k!)^4`, `H_m` harmonic numbers,

    S_0 := Σ_{k≥0} A_k = 4F3(1/2,1/2,1/2,1/2; 1,1,1; 1) = 1.11863638716...,
    S_1 := Σ_{k≥0} A_k (H_{2k} − H_k − log 2) = −0.70900136746...      (every term is ≤ 0; terms are O(k^{−3})).

> **Theorem U2.** The continued fraction converges, and
>
>     1/PCF = (93/4) ζ(5) + (π^6/64) · S_0/S_1 ,       i.e.   PCF = 64 S_1 / (1488 ζ(5) S_1 + π^6 S_0) = 2.4517460646287829170...
>
> Moreover `p_{n−1} = 8^n n!^5 w_n` with
> `w_n = Σ_{k=0}^n A_k A_{n−k} (1 + 4(2k−n)(h_k − h_{n−k}))`, `h_k = H_{2k} − H_k`, and
> `Σ_{n≥0} 1/((n+1)^5 w_n w_{n+1}) = 8/PCF`.

> **Corollary (L-value form).** Let `f = η(2τ)^4 η(4τ)^4` be the weight-4 newform of level 8. One has
> `S_0 = 16 L(f,2)/π^2` (proved: Rogers–Wan–Zucker 2015; Zagier; see Zudilin, SIGMA 14 (2018) 086, eq. (4)).
> One also has `8 S_1 = −16 L(f,1)` (the case `(r,t) = (1/2,1/2)` of the relation "L(f,1)/F_1(1) ∈ Q" that Zudilin,
> loc. cit. eq. (6), states as suggested by numerical calculation). This is **proved**: it follows from published
> formulas (Borwein–Straub–Wan 2013, `W_4(−1) = −(4/π) S_1 = (4/π^3)∫_0^1 K'(k)^2 dk`; Rogers–Wan–Zucker 2015, (31)),
> and a self-contained proof is in `../Z/THEOREM.md`. Hence, unconditionally,
>
>     PCF = 8 L(f,1) / (186 ζ(5) L(f,1) − π^4 L(f,2)),      1/PCF = (93/4) ζ(5) − (π^4/8) L(f,2)/L(f,1).
>
> This last identity is confirmed numerically to 299 digits (`verify/v4_numerics.py`, part 5).

So this constant was not found to be a combination of zeta values (and we do not prove that it is not one): besides
ζ(5) and π its closed form involves the two periods of the rigid Calabi–Yau modular form attached to the conifold
point of the family `4F3(1/2,1/2,1/2,1/2;1,1,1;x)`.

**Openness / novelty.** As for U1 (`../U1/THEOREM.md` §0 and `../U1/AUDIT-2.md`): the row is numerical in both arXiv
versions; no identification found in the literature searched there. "New" means "not found after that search".

**Status of the proof.** The theorem (hypergeometric form) is proved below. The Corollary is unconditional since
2026-11-26: the identity for `S_1` (which an earlier version of this file called unproved — wrongly, it was already
derivable from Borwein–Straub–Wan and Rogers–Wan–Zucker) is proved in `../Z/THEOREM.md`, audited in `../Z/AUDIT.md`.
(Zudilin, loc. cit., Remark 2, says none of the relations (6) seem to be proved;
Bönisch–Klemm–Scheidegger–Zagier, arXiv:2203.09426, §4, prove for this family that
the conifold periods are combinations of the periods and quasi-periods of `f` only up to multiplicative constants,
the exact rationals being numerical). Machine checks, all passing:
`verify/v1_wedge_identity.py` (Lemma 1, exact), `v2_local_and_boundary.py` (local series at x = 1, Lemma 3 and the
vanishing of the boundary polynomial in Prop. 5, exact), `v3_series_checks.py` (Prop. 2, exact, n ≤ 60),
`v5_final_algebra.py` (§7, exact), `v6_u3_gauge.py` (§9, exact), `v4_numerics.py` (every analytic statement of §4–§7 re-checked numerically to
15–50 digits, and the final formula against the 300-digit value of the fraction). Not formalized in Lean.
Independent adversarial audit: `AUDIT.md` (line-by-line; verdict "complete after listed fixes", no fatal error; the
fixes — a citation, two justifications of interchanges, a uniformity remark — are incorporated in this version; the
auditor re-derived Lemma 3, the local data, the Casoratian, Prop. 8's principal parts and §7 by hand and checked the
value to 18 digits by brute-force summation of `S_0`, `S_1`).

---

## 1. Notation

`θ = x d/dx`, `s = 1 − x`,

    L4 := θ^4 − x (θ + 1/2)^4,
    L5 := θ^5 − 2x (θ+1/2)((θ+1/2)^4 + (θ+1/2)^2 + 1/16) + x^2 (θ+1)^5,
    ã(u) := 2(u+1/2)((u+1/2)^4 + (u+1/2)^2 + 1/16),          ã(n) = a_n/8   (direct expansion; checked in v1).

`L5 Σ c_n x^n = Σ (n^5 c_n − ã(n−1) c_{n−1} + (n−1)^5 c_{n−2}) x^n`.

    N(ε) := (Γ(1+ε)/Γ(1/2+ε))^4,        Y(ε) := Σ_{k≥0} ((1/2+ε)_k/(1+ε)_k)^4 = Σ_i η_i ε^i,
    Ξ(ε) := 2 ζ(3, 1/2+ε) N(ε)          (ζ(s,a) the Hurwitz zeta function).

`Y` is analytic at least for `Re ε > −1/2` (the terms are `O(k^{−2})` locally uniformly), `Y(0) = S_0`, and
`Y'(0) = Σ A_k λ_k`, `λ_k := 4 Σ_{j<k} (1/(j+1/2) − 1/(j+1)) = 8(H_{2k} − H_k)`; `N'(0)/N(0) = 4(ψ(1) − ψ(1/2)) = 8 log 2`.
Hence

    D := Y'(0) − 8 log 2 · Y(0) = 8 S_1 < 0,                                                        (1.1)

since `H_{2k} − H_k = Σ_{j=k+1}^{2k} 1/j < log 2`.

**Frobenius basis at 0.** `A_k(ε) := ((1/2+ε)_k/(1+ε)_k)^4`, `Φ(x,ε) := Σ_k A_k(ε) x^{k+ε}`; then `L4 Φ = ε^4 x^ε`,
so `y_j := [ε^j] Φ` (`j = 0..3`) are solutions of `L4 y = 0` on `(0,1)`,

    y_j = Σ_{i≤j} (log^i x / i!) φ_{j−i}(x),   φ_m(x) = Σ_k [ε^m]A_k(ε) x^k,   φ_0(0) = 1,  φ_m(0) = 0 (m ≥ 1),      (1.2)

and they form a basis of the solution space on `(0,1)`. In particular `y_0 = Σ A_k x^k`, `y_1 = y_0 log x + Σ A_k λ_k x^k`.

## 2. The operator L5 and the sequences w_n, v_n

**Lemma 1.** If `f, h` are solutions of `L4 y = 0` on an interval, then `m := f θh − h θf` satisfies `L5 m = 0`.
*Proof.* Exact computation in the differential ring `Q(x)[f, θf, θ^2f, θ^3f, h, ..., θ^3h]`, with `θ^4` eliminated
by `L4`: `verify/v1_wedge_identity.py` (output `0`). ∎ (`L5` is the exterior square of `L4`.)

**Proposition 2.** Let `ψ_0 := y_0 θy_1 − y_1 θy_0`. Then `ψ_0 = Σ_{n≥0} w_n x^n` (a power series, radius 1) with

    w_n = Σ_{k=0}^n A_k A_{n−k} (1 + (2k−n) λ_k) = Σ_{k=0}^n A_k A_{n−k} (1 + (2k−n)(λ_k − λ_{n−k})/2) ≥ Σ_k A_k A_{n−k} ≥ A_n ≥ 1/(16 n^2),

`(n+1)^5 w_{n+1} + n^5 w_{n−1} = ã(n) w_n` for `n ≥ 0`, and `p_{n−1} = 8^n n!^5 w_n`. If `v_0 = 0`, `v_1 = 1/8` and
`(n+1)^5 v_{n+1} + n^5 v_{n−1} = ã(n) v_n` for `n ≥ 1`, then `q_{n−1} = 8^n n!^5 v_n`, and

    v_{n+1}/w_{n+1} − v_n/w_n = 1/(8 (n+1)^5 w_n w_{n+1}) > 0     (n ≥ 0).                              (2.1)

*Proof.* With `g_1 := Σ A_k λ_k x^k`: `θy_1 = θy_0 · log x + y_0 + θg_1`, so `ψ_0 = y_0^2 + y_0 θg_1 − g_1 θy_0`, whose
`x^n`-coefficient is the first sum; the second follows by averaging with `k ↔ n−k`. `λ_k` is increasing, so
`(2k−n)(λ_k − λ_{n−k}) ≥ 0`; and `(1/2)_n/n! = C(2n,n)/4^n ≥ 1/(2√n)`. By Lemma 1, `L5 ψ_0 = 0`, which is the
recurrence. The statements on `p, q` follow as in Proposition 1 of `../U1/THEOREM.md` (`p_{−1} = 1 = w_0`,
`p_0 = a_0 = 8ã(0) = 8w_1`, `q_{−1} = 0 = v_0`, `q_0 = 1 = 8v_1`, and `U_{n−1} := 8^n n!^5 c_n` turns the recurrence of
`c` into that of the continued fraction). (2.1): the Casoratian `D_n = w_n v_{n+1} − w_{n+1} v_n` satisfies
`(n+1)^5 D_n = n^5 D_{n−1}`, `D_0 = 1/8`. (All checked exactly for `n ≤ 60`: `v3_series_checks.py`.) ∎

Consequently `v_n > 0` for `n ≥ 1`, `q_n ≠ 0`, and `q_n/p_n = v_{n+1}/w_{n+1}` increases to a limit `ℓ ∈ (0, ∞]`; the
continued fraction converges (to `1/ℓ`, possibly 0) and `8ℓ = Σ 1/((n+1)^5 w_n w_{n+1})`.

## 3. Two solutions of L4 that are holomorphic at x = 1, and their Mellin transforms

In the variable `s = 1 − x`, `θ = −(1−s) d/ds`, and the coefficient of `s^p` in `L4 Σ c_m s^m` equals
`(p+1)(p+2)^2(p+3) c_{p+3} + (a linear form in c_0, ..., c_{p+2})` (`v2`, part (a)). Hence for arbitrary `c_0, c_1, c_2`
there is a unique formal solution `Σ c_m s^m`, and it converges for `|s| < 1`: `x = 1` is a regular singular point
of `L4`, and the nearest other singularity is `x = 0`. Precisely: with `δ = s d/ds` one has `θ = −((1−s)/s) δ`, and
`s^3 · L4 = s^3 (s θ^4 − (1−s)(2θ^3 + (3/2)θ^2 + (1/2)θ + 1/16))` is of the form `(1−s)^4 δ^4 + Σ_{i≤3} P_i(s) δ^i` with
polynomials `P_i`; so `L4 y = 0` is equivalent to a first-order system `s w' = A(s) w` for `w = (y, δy, δ^2y, δ^3y)`
with `A` holomorphic in `|s| < 1` (a singularity of the first kind), and every formal power-series solution of
such a system converges in the disc where `A` is holomorphic (Coddington–Levinson, *Theory of Ordinary
Differential Equations*, ch. 4, §3, Theorem 3.1). Define

    ν  := s + (3/4) s^2 + (29/48) s^3 + ...,          r := s^2 + (5/4) s^3 + (377/288) s^4 + ...,
    μ  := ν θr − r θν = −s^2 − (3/2) s^3 − ...        (a solution of L5 by Lemma 1).

They are real-analytic solutions on `(0,1]`; near `x = 0` each is a combination of `y_0..y_3`, hence `O(|log x|^3)`
(resp. `μ = O(log^4 x)`), and the same holds for their θ-derivatives.

**Lemma 3 (boundary terms).** Let `F` be a solution of `L4` on `(0,1)`, holomorphic at `x = 1`, `F_i := θ^iF(1)`, and
`F̂(u) := ∫_0^1 F(x) x^{u−1} dx` (`u > 0`). Then

    u^4 F̂(u) − (u+1/2)^4 F̂(u+1) = F_2 + F_1/2 − u(u+1/2) F_0 .

*Proof.* `∫_0^1 x^c θG dx = G(1) − (c+1) ∫_0^1 x^c G dx` for `c > −1` and `G = O(|log x|^k)` at 0. Apply it four
times to each term of `0 = ∫ x^{u−1} L4F dx`; the boundary polynomial is computed symbolically in `v2`, part (c). ∎

For `ν`: `(F_0,F_1,F_2) = (0, −1, 1/2)`, right side `0`. For `r`: `(0, 0, 2)`, right side `2`. (`v2`, part (c).)

**Proposition 4.** For `u > 0`:  `ν̂(u) = Γ(u)^4/Γ(u+1/2)^4 = N(u)/u^4`  and  `r̂(u) = 2Y(u)/u^4`.

*Proof.* `G(u) := Γ(u)^4/Γ(u+1/2)^4` satisfies `u^4 G(u) = (u+1/2)^4 G(u+1)` and `G(u) ~ u^{−2}` (`u → +∞`).
(i) By Lemma 3, `π_1(u) := ν̂(u)/G(u)` has period 1. As `u → +∞`, `ν̂(u) = ∫_{1/2}^1 ((1−x) + O((1−x)^2)) x^{u−1} dx +
O(2^{−u} u^{−1}) = u^{−2}(1 + o(1))`, so `π_1(u) → 1`; a 1-periodic function with a limit at `+∞` is constant.
(ii) `R_0(u) := 2Y(u)/u^4` satisfies `u^4 R_0(u) − (u+1/2)^4 R_0(u+1) = 2Y(u) − 2((u+1/2)/(u+1))^4 Y(u+1) = 2`,
because `((u+1/2)/(u+1))^4 ((u+3/2)_k/(u+2)_k)^4 = ((u+1/2)_{k+1}/(u+1)_{k+1})^4`. So `π_2 := (r̂ − R_0)/G` has
period 1. Now `r̂(u) = O(u^{−3})` (as `r = O((1−x)^2)` near 1), and `0 < Y(u) ≤ u + 2`: indeed
`(u+1/2)_k/(u+1)_k = Π_{j<k}(1 − (1/2)/(u+1+j)) ≤ exp(−(1/2) Σ_{j<k} 1/(u+1+j)) ≤ ((u+1)/(u+1+k))^{1/2}`, so
`Y(u) ≤ (u+1)^2 Σ_{k≥0} (u+1+k)^{−2} ≤ 1 + (u+1)`. Hence `π_2(u) = O(u^{−3}) · u^2 → 0`, so `π_2 ≡ 0`. ∎

## 4. The minimal solution of the recurrence

Put `g(u) := ∫_0^1 μ(x) x^u dx` for `u > −1`.

**Proposition 5.** (i) `u^5 g(u−1) + (u+1)^5 g(u+1) = ã(u) g(u)` for all `u > 0`.
(ii) Let `μ_4 := lim_{x→0} μ(x)/log^4 x`. Then `g(1) = ã(0) g(0) − 24 μ_4`.
(iii) `g(n) = g(0) w_n − 192 μ_4 v_n` for all integers `n ≥ 0`, and `g(n) = O(n^{−3})`.

*Proof.* (i) Integrate `0 = x^{u−1} L5μ` over `(0,1)` by parts as in Lemma 3, for the three terms of `L5`:

    0 = −u^5 g(u−1) + ã(u) g(u) − (u+1)^5 g(u+1) + B_5(u),

where `B_5(u)` is an explicit polynomial in `u` and in `θ^iμ(1)`, `i = 0..4`. From the local series,
`(θ^iμ(1))_{i=0..4} = (0, 0, −2, 3, −9/4)`, and with these values `B_5(u) ≡ 0` (`v2`, part (d), exact).
(ii) By (1.2), `μ` is, near 0, a polynomial of degree ≤ 4 in `log x` with constant coefficients plus terms
`O(x |log x|^5)`; so `μ = μ_4 log^4 x + O(1 + |log x|^3)`, and `∫_0^1 log^i x · x^{u−1} dx = (−1)^i i!/u^{i+1}`; hence
`u^5 g(u−1) → 24 μ_4` as `u → 0+`. Let `u → 0+` in (i) (`g` is continuous on `(−1,∞)`).
(iii) By (i) at `u = n ≥ 1`, `(g(n))_{n≥0}` satisfies the recurrence of Prop. 2 for `n ≥ 1`, as do `w` and `v`; the
right side has the same values at `n = 0` (`g(0)`) and `n = 1` (`g(0) ã(0) − 192 μ_4/8`, by (ii)). Finally
`|μ(x)| ≤ C(1−x)^2` on `[1/2,1]` gives `g(n) = O(n^{−3})`. ∎

**Corollary 6.** If `μ_4 ≠ 0`, then `ℓ = lim v_n/w_n = g(0)/(192 μ_4)`, in particular `ℓ < ∞`, `g(0) ≠ 0`, and

    PCF = 192 μ_4 / g(0).

*Proof.* `v_n/w_n = (g(0) − g(n)/w_n)/(192 μ_4)` and `g(n)/w_n = O(n^{−3} · n^2) → 0` by Prop. 2. As `ℓ > 0`,
`g(0) ≠ 0`. ∎

## 5. Frobenius coordinates; evaluation of μ_4 and g(0)

For `u > 0`, `ŷ_j(u) = ∫_0^1 y_j x^{u−1} dx = [ε^j] Σ_k A_k(ε)/(k+ε+u)`, and the term `k = 0` is
`[ε^j] 1/(ε+u) = (−1)^j u^{−j−1}`, the others being analytic at `u = 0`. Therefore, if `ρ = Σ_j c_j y_j`, then
`ρ̂(u) = Σ_j (−1)^j c_j u^{−j−1} + (a function analytic at u = 0)`. Comparing with Proposition 4 (`N` and `Y` are
analytic at 0, `N = Σ N_i u^i`):

    ν = Σ_j (−1)^j N_{3−j} y_j,          r = 2 Σ_j (−1)^j η_{3−j} y_j.                                  (5.1)

(Checked numerically at `x = 1/2` against the local series, values and three derivatives, to 49 digits: `v4`, (1).)

**Proposition 7.** `μ_4 = (N_0 η_1 − N_1 η_0)/6 = N(0) D/6 ≠ 0`, and `g(0) = −2 [ε^3] (Ξ(ε) Y(−ε))`. Hence

    PCF = 16 · [ε^1](N(ε) Y(−ε)) / [ε^3](Ξ(ε) Y(−ε)).                                               (5.2)

*Proof.* By (1.2) the highest power of `log x` in `y_j θy_k` is `log^{j+k−1} x` with coefficient `1/(j!(k−1)!)`. In
`μ = Σ_{j,k} c_j(ν) c_k(r) (y_j θy_k − y_k θy_j)` the power `log^4` comes only from `(j,k) = (2,3), (3,2)`, with
coefficient `±(1/4 − 1/6) = ±1/12`: `μ_4 = (c_2(ν) c_3(r) − c_3(ν) c_2(r))/12 = (N_1(−2η_0) + N_0 (2η_1))/12`.
With `N_1 = 8 log 2 · N_0` and (1.1), `μ_4 = N(0) D/6`, and `N(0) = π^{−2}`, `D < 0`.
Next, for small `ε`, using `∫_0^1 x^c θν dx = −(c+1) ν̂(c+1)` (`ν(1) = 0`) and Proposition 4,

    ∫_0^1 (ν θΦ − Φ θν) dx = Σ_k A_k(ε) (2k+2ε+1) ν̂(k+ε+1) = N(ε) Σ_k (2k+2ε+1)/(k+ε+1/2)^4 = Ξ(ε),

since `A_k(ε) Γ(k+ε+1)^4/Γ(k+ε+3/2)^4 = N(ε)/(k+ε+1/2)^4`. Taking `[ε^k]` (uniform convergence),
`∫_0^1 (ν θy_j − y_j θν) dx = ξ_j := [ε^j]Ξ`, and by (5.1)
`g(0) = Σ_j c_j(r) ξ_j = 2 Σ_j (−1)^j η_{3−j} ξ_j = −2 [ε^3](Ξ(ε) Y(−ε))`. (Justification of the two interchanges:
`|A_k(ε)| ≤ C k^{−2}` for complex `|ε| ≤ 1/4`, `k ≥ 1`. On `[1/2,1]`: `|ν(x)| ≤ C(1−x)`, `|θν| ≤ C`, so
`∫_{1/2}^1 (|ν| (k+1) + |θν|) x^{k−1/4} dx ≤ C(k+1) ∫_0^1 (1−x) x^{k−1/4} dx + C/(k+3/4) ≤ C_1/k`. On `(0,1/2]`:
`|ν|, |θν| ≤ C(1 + |log x|^3)`, so `∫_0^{1/2} (|ν|(k+1) + |θν|) x^{k−1/4} dx ≤ C (k+2) 2^{−k} ∫_0^{1/2} (1+|log x|^3) x^{−1/4} dx
≤ C_2 (k+2) 2^{−k}`. Hence `Σ_k |A_k(ε)| ∫_0^1 (|ν| |k+ε| + |θν|) x^{k+Re ε} dx ≤ C_3 Σ_k k^{−3} < ∞` uniformly in
`|ε| ≤ 1/4`: Fubini applies, the integral is analytic in `ε`, and coefficient extraction commutes with it by
Cauchy's formula.) Finally (5.2) is Corollary 6 with
`192 μ_4 = 32 (N_0η_1 − N_1η_0) = −32 [ε^1](N(ε)Y(−ε))`. ∎

(Numerically: `g(0)`, `g(1) − ã(0)g(0) + 24μ_4 = 0`, the functional equation of `g`, and `192μ_4/g(0) =
2.451746064628782917`: `v4`, (3).)

## 6. A bilateral sum: two relations between the Taylor coefficients of Y

Write `(Y/N)(u) = Σ_i τ_i u^i` near `u = 0`; `τ_0 = π^2 Y(0)`, `τ_1 = π^2 D`.

**Proposition 8.** `τ_2 = 2π^2 τ_0` and `τ_3 = (4π^2/3) τ_1`. More precisely,

    Y(u)/N(u) + tan^4(πu) · Y(1/2−u)/N(1/2−u) = sec^2(πu) ( τ_0 sec^2(πu) + (τ_1/π) tan(πu) )      (−1/2 < Re u < 1, u ≠ 1/2).

*Proof.* Let `T(u) := Σ_{k∈Z} (Γ(u+k+1/2)/Γ(u+k+1))^4`. For `k → +∞` the terms are `O(k^{−2})`; for `k = −m`, the
reflection formulas `Γ(1/2+z)Γ(1/2−z) = π/cos πz`, `Γ(1+z)Γ(−z) = −π/sin πz` give
`Γ(z+1/2)/Γ(z+1) = −tan(πz) Γ(−z)/Γ(1/2−z)`, so the term is `tan^4(πu) (Γ(m−u)/Γ(m+1/2−u))^4 = O(m^{−2})`. Hence `T` is a
meromorphic function of period 1, with poles at most at `1/2 + Z`, and

    T(u) = Y(u)/N(u) + tan^4(πu) · Y(1/2−u)/N(1/2−u)

(the part `k ≥ 0` is `(Γ(u+1/2)/Γ(u+1))^4 Y(u)`; in the part `m ≥ 1` put `m = k+1` and `Γ(k+1−u) = Γ((1/2−u)+k+1/2)`).
*At `u = 0`:* `T(u) = Σ_{i≤3} τ_i u^i + O(u^4)`.
*At `u = 1/2`,* with `v = 1/2 − u`: `tan^4(πu) = cot^4(πv) = (πv)^{−4}(1 − (4π^2/3) v^2 + O(v^4))`, so the principal
part of `T` is `π^{−4} [τ_0 v^{−4} + τ_1 v^{−3} + (τ_2 − (4π^2/3)τ_0) v^{−2} + (τ_3 − (4π^2/3)τ_1) v^{−1}]`.
*At `±i∞`:* uniformly for `Re u ∈ [0,1]`, `T(u) → 0` as `|Im u| → ∞` (for `|Im u| ≥ 1` each term is `O(1/|u+k|^2)`
uniformly, by Stirling for `k ≥ 0` and by the reflected form with `|tan πu| ≤ coth π` for `k < 0`, so `|T| ≤ C Σ_k ((Re u + k)^2 + (Im u)^2)^{−1} =
O(1/|Im u|)`).
The functions `S_2 = π^2 sec^2 πu = v^{−2} + π^2/3 + O(v^2)`, `S_3 = π^3 sec^2 πu tan πu = v^{−3} + O(v)`,
`S_4 = π^4 sec^4 πu = v^{−4} + (2π^2/3) v^{−2} + O(1)` have period 1 and tend to 0 at `±i∞` (`v5`). Put
`T_0 := π^{−4} [τ_0 S_4 + τ_1 S_3 + (τ_2 − 2π^2 τ_0) S_2]` and `c := π^{−4}(τ_3 − (4π^2/3)τ_1)`. Then
`E := T − T_0 − c π tan(πu)` is entire and 1-periodic and bounded on the strip, hence constant; its limits at
`+i∞` and `−i∞` are `−iπc` and `+iπc`, so `c = 0` and `E = 0`. Thus `τ_3 = (4π^2/3) τ_1` and `T = T_0`. Evaluating
at `u = 0`: `τ_0 = T(0) = T_0(0) = τ_0 + π^{−2}(τ_2 − 2π^2 τ_0)`, i.e. `τ_2 = 2π^2 τ_0`; then `T_0` is the right side
of the displayed formula. ∎ (Checked numerically, including at a complex point: `v4`, (4).)

## 7. End of the proof

`N(ε)N(−ε) = (Γ(1+ε)Γ(1−ε))^4/(Γ(1/2+ε)Γ(1/2−ε))^4 = (ε cot πε)^4 = π^{−4}(1 − (4π^2/3)ε^2 + O(ε^4))`, and
`ζ(3,1/2+ε) = Σ_m C(−3,m) ζ(3+m,1/2) ε^m = 7ζ(3) − (π^4/2) ε + 186 ζ(5) ε^2 − (2π^6/3) ε^3 + O(ε^4)`
(`ζ(s,1/2) = (2^s − 1)ζ(s)`). Writing `N(ε)Y(−ε) = N(ε)N(−ε) · (Y/N)(−ε)` and using Proposition 8
(`v5_final_algebra.py`, exact):

    [ε^1](N(ε)Y(−ε)) = −τ_1/π^4,
    [ε^3](Ξ(ε)Y(−ε)) = (2/π^4) { 7ζ(3)(−τ_3 + (4π^2/3)τ_1) − (π^4/2)(τ_2 − (4π^2/3)τ_0) − 186ζ(5) τ_1 − (2π^6/3) τ_0 }
                     = −(2/π^4) ( π^6 τ_0 + 186 ζ(5) τ_1 ).

(The ζ(3)-term vanishes by the second relation of Prop. 8.) By (5.2),

    PCF = 8 τ_1/(π^6 τ_0 + 186 ζ(5) τ_1) = 8D/(π^6 Y(0) + 186 ζ(5) D),      1/PCF = (93/4) ζ(5) + (π^6/8) Y(0)/D,

and `Y(0) = S_0`, `D = 8 S_1`. The series identity of the statement is `8ℓ`, by (2.1). ∎

## 8. Remarks

1. **L-values.** With `S_0 = 16L(f,2)/π^2` and `8S_1 = −16L(f,1)`: `1/PCF = (93/4)ζ(5) − (π^4/8) L(f,2)/L(f,1)`.
   Note `186 = 372/2`, the same number as in U1 (`ζ(5,1/2) = 31ζ(5)`).
2. **Asymptotics (numerical observation, not proved here).** `w_n ~ c_w/n` with `c_w = 16 L(f,1)/π^2` to 300 digits
   (`../code/u2_pslq.py`); also `c_w = −D/π^2` numerically. No asymptotic analysis of `w_n` is carried out in the
   proof, and none is needed.
2a. **Relation to Zudilin's paper (credit).** The method of Proposition 8 — periodicity in `ε` of the bilateral
   hypergeometric sum and the finiteness of its Fourier expansion — is that of Zudilin, SIGMA 14 (2018) 086, §3
   (Lemma 3 and Remark 5, after Golyshev–Mellit), and its second relation `τ_3 = (4π^2/3) τ_1`, i.e. the rational
   proportionality of `F_1(1)` and `F_3(1)/π^2`, is stated there as rigorously established. The first relation
   `τ_2 = 2π^2 τ_0` is the "extra rational relation `F_0(1) = F_2(1)/(2π^2)`" of his Remark 2 for `(r,t) = (1/2,1/2)`.
   Proposition 8 is therefore not new; it is proved here in full to keep the argument self-contained. In
   Zudilin's normalisation `F_j(1) = τ_j/π^2`; combined with the proved `F_0(1) = S_0 = 16 L(f,2)/π^2`, the first
   relation gives **`L(f,2)/F_2(1) = 1/32`**, one of the relations (6) of that paper; the other two,
   `L(f,1)/F_1(1) = −1/16` and `L(f,3)/F_3(1) = −3/256`, are proved in `../Z/THEOREM.md`.
3. **Structure behind the proof.** `ν` is the vanishing period at the conifold `x = 1` (its Mellin transform is the
   Meijer-type `Γ(u)^4/Γ(u+1/2)^4`), `r` is the unique solution vanishing to second order there, and `g(u)`, the
   Mellin transform of the minor `μ`, is the minimal solution of the recurrence continued to real `u`; the pole
   of order 5 of `g` at `u = −1` produces the inhomogeneity `−24μ_4` that separates `p_n` from `q_n`. (Heuristic
   comment, not proved here: the relations of Prop. 8 express the invariance of `r` under the involution
   `y(x) ↦ x^{−1/2} y(1/x)` of `L4`.)
4. **How it was found.** Structure of the recurrence (`../code/z5_struct.py`: `ψ_0` test), asymptotic constant
   identified by PSLQ, then the Mellin analysis above predicted the form `ζ(5) + π^4 L(f,2)/L(f,1)`, which PSLQ
   confirmed at 300 digits (`../code/u23_closed.py`: relation `(8, −186, 1)`).

## 9. The third constant (U3) — SUPERSEDED: the formula below is now proved in `../U3/THEOREM.md` (2026-11-26)

(This section is kept as the record of how the formula was found. Where it says "not proved" or "found by PSLQ",
see `../U3/THEOREM.md`: the Level-A formula is Prop. 4.1 there, and the reduction of the Taylor coefficients of `Z`
to `Z_0` is Prop. 5.1 and 6.1.)

Row 6 of the same table: `a_n = 8(n^5+(n+1)^5) + 20(n^3+(n+1)^3) − 5(2n+1)`, `b_n = −64n^10`, "22.8410...".

* **Exact (proved; `verify/v6_u3_gauge.py`):** if `X_n` satisfies the U2 recurrence
  `(n+1)^5 X_{n+1} + n^5 X_{n−1} = ã(n) X_n`, then `−(2n+1)(16n^4+56n^3+76n^2+46n+11) X_n + 32 (n+1)^5 X_{n+1}`
  satisfies the U3 recurrence; this map sends `w` to the U3 numerator sequence and `v − 4w` to the U3 denominator
  sequence. (Not proved, numerical: the map exchanges dominant and minimal asymptotic behaviour, so that
  `1/U3 + 4 = (coefficient of n^{−3} in v_n)/(coefficient of n^{−3} in w_n)` in the expansions in powers of `n + 1/2`.)
* **Numerical identification (60 digits; `../code/u3_connect.py`, `u3_formula.py`, `u3_z0.py`).** With one more
  hypergeometric constant

      Z_0 := Σ_{k≥0} A_k/(k+1) = ∫_0^1 4F3(1/2,1/2,1/2,1/2;1,1,1;x) dx = 5F4(1/2,1/2,1/2,1/2,1; 1,1,1,2; 1) = 1.04335518933431988860...,

  one has, with `K := 4 + π^2 (Z_0 − 2 S_0)`,

      1/U3 = (93/4) ζ(5) − 4 + (π^6/64) · S_0 K / (1 + S_1 K)
           = 1/U2 − 4 + π^4 L(f,2) / ( 8 L(f,1) (64 L(f,1)L(f,2) − 2π^2 L(f,1) Z_0 − 8 L(f,1) + 1) )

  (second form under `S_0 = 16L(f,2)/π^2`, `8S_1 = −16L(f,1)`); `1/U3 − formula = 1.8·10^{−60}` at 60-digit precision.
  PSLQ finds no expression of `Z_0` through `L(f,1), L(f,2), π` (three bases, 60 digits): it plays the role of a
  quasi-period of `f` (cf. Bönisch–Klemm–Scheidegger–Zagier, arXiv:2203.09426).
  *Derivation (not a proof).* The argument of §3–§5 applies to every solution `ρ` of `L4` holomorphic at `x = 1`:
  `g_ρ(u) = ∫_0^1 (ν θρ − ρ θν) x^u dx` solves the recurrence, and `g_ρ(n) = g_ρ(0) w_n − 192 μ_4(ρ) v_n`. U2 uses the
  `ρ` with `g_ρ(u) = O(u^{−3})` (`ρ = r`); U3, through the gauge map above, uses the `ρ` whose `g_ρ` has no `u^{−3}`
  term in its expansion in powers of `u + 1/2`, namely `r_* = 1 − s^2/16 + O(s^3)`, whose Mellin transform is
  `1/u − (2Y(u) − Z(u))/(16u^4)`, `Z(u) := Σ_k ((u+1/2)_k/(u+1)_k)^4/(u+k+1)`. This gives (checked to 30 digits)
  `1/U3 + 4 = [ε^3](Ξ(ε)Y_*(−ε)) / (16 [ε^1](N(ε)Y_*(−ε)))`, `Y_*(u) = u^3/2 − Y(u)/16 + Z(u)/32`; the reduction of the
  Taylor coefficients of `Z` to `Z(0)` alone (equivalently the relation `ŷ_1 = (1/2 − 4L(f,1)) ν + b r − 16 L(f,1) r_0`
  with `b = (π^2 L(f,1) Z_0 + 4L(f,1) − 1/2)/(32 L(f,2))`) was found by PSLQ and is not proved.
  **Status: U3 is identified numerically in terms of ζ(5), π, S_0, S_1 and the new constant Z_0; not proved.**
