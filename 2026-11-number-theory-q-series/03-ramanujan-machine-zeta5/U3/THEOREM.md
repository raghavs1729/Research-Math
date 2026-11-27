# Theorem U3 — the third "unidentified" ζ(5)-level continued fraction of the Ramanujan Machine

## 0. Statement, source, status

**Source (read 2026-11-24).** arXiv:2308.11829v2 (Elimelech et al., PNAS 121 (2024) e2321440121), §3.6, Table 1,
sixth row (also App. C, Table 4):

    a_n = 8(n^5 + (n+1)^5) + 20(n^3 + (n+1)^3) − 5(2n+1),     b_n = −64 n^10,     Formula: 22.8410...

("Three of the discovered formulas in the table did not fit any combination of zetas (in numerical tests using
PSLQ) and are thus presented with their numerical values"). Convention as in `../U2/THEOREM.md`:
`PCF = a_0 + b_1/(a_1 + b_2/(a_2 + ...)) = lim p_n/q_n`.

Notation: `A_k = C(2k,k)^4/256^k = ((1/2)_k/k!)^4`,

    S_0 := Σ_{k≥0} A_k,      S_1 := Σ_{k≥0} A_k (H_{2k} − H_k − log 2),      Z_0 := Σ_{k≥0} A_k/(k+1),
    S_0 = 1.1186363871641870683...,   S_1 = −0.7090013674619294375...,   Z_0 = 1.0433551893343198886...
    (Z_0 = 5F4(1/2,1/2,1/2,1/2,1; 1,1,1,2; 1) = ∫_0^1 4F3(1/2,1/2,1/2,1/2;1,1,1;x) dx).

> **Theorem U3.** The continued fraction converges, and with `K := 4 + π^2 (Z_0 − 2 S_0)` (`K < 0`),
>
>     1/PCF = (93/4) ζ(5) − 4 + (π^6/64) · S_0 K / (1 + S_1 K),        PCF = 22.841008565767985615909667462437...
>
> Equivalently, with `f = η(2τ)^4 η(4τ)^4` (`S_0 = 16 L(f,2)/π^2`, `S_1 = −2 L(f,1)`, see `../Z/THEOREM.md`):
> `1/PCF = (93/4) ζ(5) − 4 + (π^4/4) L(f,2) K / (1 − 2 L(f,1) K)`, `K = 4 + π^2 Z_0 − 32 L(f,2)`;
> and `1/U3 = 1/U2 − 4 − π^6 S_0/(64 S_1 (1 + S_1 K)) = 1/U2 − 4 + π^4 L(f,2)/(8 L(f,1)(1 − 2 L(f,1) K))`, `U2` the
> constant of `../U2/THEOREM.md`.

Besides `ζ(5)`, `π` and the two periods `S_0, S_1` of `f` that occur in U2, the value involves **one** further
number, `Z_0`. PSLQ at 60 digits found no integer relation (coefficients up to `10^6`, three bases built from
`Z_0, L(f,1), L(f,2)`, powers of `π` and 1); by analogy with the literature (an interpretation, not a statement
proved here) `Z_0` plays the role of a quasi-period of `f` (Remark 2). Nothing is
claimed or proved about its transcendence or independence — it may well be expressible by known constants in a
way those searches did not cover.

**Openness / novelty.** As for U1, U2: the row is numerical in both arXiv versions and no identification was found
in the literature searched for U1/U2 (`../U1/THEOREM.md` §0, `../U1/AUDIT-2.md`). "New" means "not found after that
search". The method is that of U2; what is specific to U3 is §2 (gauge map), §5 (a Legendre-type relation obtained
from the Casoratian) and §6 (the bilateral sum for `Z` and the contiguity relation `e_2 = e_0 + 8`).

**Status.** Proved below, building on `../U2/THEOREM.md` (§1–§6 there; all notation of that file is used).
Machine checks, all passing: `verify/u3_exact.py` (every algebraic identity, exact), `verify/u3_algebra.py` (§7,
exact), `../U2/verify/v6_u3_gauge.py` (§2, exact), `verify/u3_numerics.py` (every analytic statement re-checked
numerically to 13–50 digits; final formula against the 48-digit value of the fraction; run completed, all
assertions pass). Independent adversarial audit: `AUDIT.md` (verdict "complete after listed fixes"; no gap in the
proof; the auditor recomputed the fraction directly from the recurrence, N = 128000 with Richardson extrapolation,
and found agreement with the closed form to 7.7·10^{−35}, re-derived §7 by hand, and checked Prop. 5.1 and 6.1
numerically with his own code; his fixes — a sign in the third displayed form, three digits in Remark 1, the
wording on PSLQ — are incorporated). A second, independent referee engaged by the parent session
(`REFEREE-2-U3.md` on the parent branch) reports "proof complete, no gap, no hidden numerical input" (fraction
recomputed to 45 digits, agreement 9.7·10^{−48}); its minor points (notation clash, the exact rational function
in Lemma 3.3, the explicit certificate in Lemma 2.2, labelling of interpretations) are incorporated.
Not formalized in Lean.

---

## 1. Notation recalled from U2

`θ = x d/dx`, `s = 1 − x`, `L4 = θ^4 − x(θ+1/2)^4`; `ã(u) = 2(u+1/2)((u+1/2)^4 + (u+1/2)^2 + 1/16)`;
`N(ε) = (Γ(1+ε)/Γ(1/2+ε))^4`, `A_k(ε) = ((1/2+ε)_k/(1+ε)_k)^4`, `Y(ε) = Σ_k A_k(ε)`, `Ξ(ε) = 2ζ(3,1/2+ε)N(ε)`;
the Frobenius basis `y_0..y_3` at `x = 0`; the solutions `ν = s + (3/4)s^2 + ...`, `r = s^2 + (5/4)s^3 + ...` of `L4`
holomorphic at `x = 1`; `μ = νθr − rθν`, `g(u) = ∫_0^1 μ x^u dx`; the sequences `w_n, v_n` (`R_2`:
`(n+1)^5 X_{n+1} + n^5 X_{n−1} = ã(n) X_n`, `w_0 = 1, w_1 = ã(0)`, `v_0 = 0, v_1 = 1/8`) with Casoratian

    w_n v_{n+1} − w_{n+1} v_n = 1/(8 (n+1)^5)        (n ≥ 0);                                             (1.1)

`μ_4 = lim_{x→0} μ/log^4 x = N(0)D/6 = 4S_1/(3π^2)`, `g(n) = g(0) w_n − 192 μ_4 v_n` (U2, Prop. 5, 7), and
`g(0)/(192 μ_4) = 1/U2 = (93/4)ζ(5) + π^6 S_0/(64 S_1)` (U2, §7). `(Y/N)(u) = Σ τ_i u^i`, `τ_0 = π^2 S_0`,
`τ_1 = 8π^2 S_1`, `τ_2 = 2π^2 τ_0`, `τ_3 = (4π^2/3) τ_1` (U2, Prop. 8).

New:

    Z(u) := Σ_{k≥0} A_k(u)/(u+k+1)    (analytic for Re u > −1/2;  Z(0) = Z_0),
    Y_*(u) := u^3/2 − Y(u)/16 + Z(u)/32 = Σ_i η*_i u^i,        (Z/N)(u) = Σ_i ζ_i u^i   near u = 0.

Since `A_{k+1}(u) = κ(u) A_k(u+1)` with `κ(u) := ((u+1/2)/(u+1))^4`,

    Y(u) − κ(u) Y(u+1) = 1,          Z(u) − κ(u) Z(u+1) = 1/(u+1).                                         (1.2)

## 2. The U3 sequences and the gauge map

Let `ã_3(n) := a_n/8 = 2T(T^4 + 5T^2 + 25/16)`, `T = n + 1/2`, and `R_3`: `(n+1)^5 X_{n+1} + n^5 X_{n−1} = ã_3(n) X_n`.
Define `w^{(3)}, v^{(3)}` by `R_3` for `n ≥ 1` and `w^{(3)}_0 = 1, w^{(3)}_1 = ã_3(0) = 23/8`, `v^{(3)}_0 = 0, v^{(3)}_1 = 1/8`.
As in U2, Prop. 2: `p_{n−1} = 8^n n!^5 w^{(3)}_n`, `q_{n−1} = 8^n n!^5 v^{(3)}_n`.

**Lemma 2.1.** `w^{(3)}_n` is increasing, `w^{(3)}_n ≥ 1`; `v^{(3)}_n/w^{(3)}_n` is increasing; hence the continued fraction
converges, to `1/ℓ_3` where `ℓ_3 := lim v^{(3)}_n/w^{(3)}_n ∈ (0, ∞]`.

*Proof.* `ã_3(n) − (n+1)^5 − n^5 = 5(2n+1)(4n^2+4n+3)/8 > 0` (`u3_exact.py` (f)). If `w^{(3)}_n ≥ w^{(3)}_{n−1} > 0`, then
`(n+1)^5 w^{(3)}_{n+1} = ã_3(n) w^{(3)}_n − n^5 w^{(3)}_{n−1} ≥ (ã_3(n) − n^5) w^{(3)}_n > (n+1)^5 w^{(3)}_n`; and `w^{(3)}_1 > w^{(3)}_0`.
The Casoratian `D_n = w^{(3)}_n v^{(3)}_{n+1} − w^{(3)}_{n+1} v^{(3)}_n` satisfies `(n+1)^5 D_n = n^5 D_{n−1}`, `D_0 = 1/8`, so
`v^{(3)}_{n+1}/w^{(3)}_{n+1} − v^{(3)}_n/w^{(3)}_n = 1/(8(n+1)^5 w^{(3)}_n w^{(3)}_{n+1}) > 0`, and `v^{(3)}_1/w^{(3)}_1 > 0`. ∎

**Lemma 2.2 (gauge map).** For a sequence `X` put

    G(X)_n := P_1(n) X_n + 32 (n+1)^5 X_{n+1},        P_1(n) := −(2n+1)(16n^4 + 56n^3 + 76n^2 + 46n + 11).

If `X` satisfies `R_2` at `n` and `n+1` (`n ≥ 1`), then `G(X)` satisfies `R_3` at `n`. Moreover `G(w) = w^{(3)}` and
`G(v) = v^{(3)} + 4 w^{(3)}`.

*Proof.* Write `R_j(X)_n := (n+1)^5 X_{n+1} + n^5 X_{n−1} − ã_j(n) X_n` (`ã_2 := ã`). For an arbitrary sequence `X`,

    R_3(G(X))_n = P_1(n−1) · R_2(X)_n + 32 (n+1)^5 · R_2(X)_{n+1}

(an identity of linear forms in `X_{n−1}, X_n, X_{n+1}, X_{n+2}` with polynomial coefficients; `u3_exact.py` (g),
and, in the equivalent eliminated form, `v6_u3_gauge.py`). This gives the first statement. Both sides of each of the two equalities satisfy `R_3` for `n ≥ 1` and agree at
`n = 0, 1` (same script; it also checks them for all `n < 40`). ∎

## 3. A third solution holomorphic at x = 1

By U2 §3, for arbitrary `c_0, c_1, c_2` there is a unique solution `Σ c_m s^m` of `L4`, convergent for `|s| < 1`. Let

    r_* := 1 − (1/16) s^2 − (7/96) s^3 − ...        (c_0, c_1, c_2) = (1, 0, −1/16),
    μ_* := ν θr_* − r_* θν = 1 + (1/2) s + (3/8) s^2 + (5/16) s^3 + ...,        g_*(u) := ∫_0^1 μ_*(x) x^u dx   (u > −1).

(`u3_exact.py` (a), (c).) As in U2, `r_*`, `μ_*` and their θ-derivatives are `O(|log x|^3)`, resp. `O(log^4 x)`, at 0.

**Proposition 3.1.** For `u > 0`: `r̂_*(u) := ∫_0^1 r_* x^{u−1} dx = 2 Y_*(u)/u^4 = 1/u − (2Y(u) − Z(u))/(16 u^4)`.
Consequently `r_* = 2 Σ_j (−1)^j η*_{3−j} y_j`.

*Proof.* `(θ^i r_*(1))_{i=0,1,2} = (1, 0, −1/8)`, so by U2, Lemma 3,
`u^4 r̂_*(u) − (u+1/2)^4 r̂_*(u+1) = −1/8 − u(u+1/2)`. By (1.2), `R_*(u) := 2Y_*(u)/u^4` satisfies the same relation:
`u^4 R_*(u) − (u+1/2)^4 R_*(u+1) = u^3 − (u+1/2)^4/(u+1) − 1/8 + 1/(16(u+1)) = −1/8 − u(u+1/2)`
(`u3_exact.py` (b)). Hence `(r̂_* − R_*)/Γ_4`, `Γ_4(u) := Γ(u)^4/Γ(u+1/2)^4` (called `G` in U2; here `G` denotes the gauge map), has period 1. As `u → +∞`:
`r̂_*(u) − 1/u = ∫_0^1 (r_* − 1) x^{u−1} dx = O(u^{−3})` (`|r_* − 1| ≤ C(1−x)^2` on `[1/2,1]`; the part `(0,1/2)` is
`O(2^{−u})`), and `R_*(u) − 1/u = −(2Y − Z)/(16u^4) = O(u^{−3})` because `0 < Z(u) ≤ Y(u) ≤ u + 2` (U2, proof of
Prop. 4). With `1/Γ_4(u) ~ u^2`, the periodic function tends to 0, so it vanishes. The Frobenius coordinates follow
as in U2 (5.1), `Y_*` being analytic at `u = 0`. ∎ (Numerically: `u3_numerics.py` (1), (2).)

**Proposition 3.2.** (i) `u^5 g_*(u−1) + (u+1)^5 g_*(u+1) = ã(u) g_*(u)` for `u > 0`.
(ii) With `μ_4^* := lim_{x→0} μ_*(x)/log^4 x`: `g_*(1) = ã(0) g_*(0) − 24 μ_4^*`, and for all integers `n ≥ 0`

    g_*(n) = g_*(0) w_n − 192 μ_4^* v_n.                                                                  (3.1)

(iii) `μ_4^* = (N_0 η*_1 − N_1 η*_0)/6 = −(1/6) [ε^1](N(ε) Y_*(−ε))` and `g_*(0) = −2 [ε^3](Ξ(ε) Y_*(−ε))`.

*Proof.* Word for word the proofs of U2, Prop. 5 and Prop. 7, with `r` replaced by `r_*`. (That `μ_*(1) = 1 ≠ 0`,
unlike `μ(1) = 0`, changes nothing: the boundary terms at `x = 1` are collected in `B_5`, which is checked to
vanish for the actual values `θ^iμ_*(1)`; the limit `u → 0+` concerns only `x → 0`; and the integral identity of
U2, Prop. 7 involves only `ν` and the `y_j`.) (i) needs only that
`L5 μ_* = 0` (U2, Lemma 1), the growth at 0, and the vanishing of the boundary polynomial `B_5(u)`, which holds
for `(θ^i μ_*(1))_{i=0..4} = (1, −1/2, 1/4, −1/8, 7/128)` — in fact `B_5 ≡ 0` for `νθρ − ρθν` with *every* solution `ρ`
holomorphic at 1 (`u3_exact.py` (c)). (ii) is U2, Prop. 5 (ii), (iii) (the bound `g(n) = O(n^{−3})` there is not
used for (3.1)). (iii) is the computation of U2, Prop. 7, which uses only the Frobenius coordinates
`c_j(ρ) = 2(−1)^j [u^{3−j}](·)` of `ρ` and the identity `∫_0^1 (ν θy_j − y_j θν) dx = [ε^j] Ξ`. ∎
(Numerically: `u3_numerics.py` (3).)

**Lemma 3.3 (asymptotics).** As `n → ∞`:
(i) `g(n) = −2/n^3 + O(n^{−4})`, `g(n+1) − g(n) = 6/n^4 + O(n^{−5})`;
`g_*(n) = 1/n + O(n^{−2})`, `g_*(n+1) − g_*(n) = −1/n^2 + O(n^{−3})`.
(ii) `G(g_*)_n = O(1/n)`.

*Proof.* Let `β_j(n) := ∫_0^1 s^j x^n dx = j! n!/(n+j+1)!`. If `m` is continuous on `(0,1]`, `O(log^4 x)` at 0, and
`|m(x)| ≤ C s^j` on `[1/2, 1]`, then `∫_0^1 |m| x^n dx ≤ C β_j(n) + O(2^{−n}) = O(n^{−j−1})`.
(i) `μ = −s^2 + O(s^3)` and `g(n+1) − g(n) = −∫ μ s x^n dx`, `μ s = −s^3 + O(s^4)`; `μ_* = 1 + O(s)`,
`g_*(n+1) − g_*(n) = −∫ μ_* s x^n dx`. (ii) For `X(n) = ∫_0^1 m(x) x^n dx`,

    G(X)_n = ∫_0^1 m(x) x^n ( c(n) − 32 (n+1)^5 s ) dx,     c(n) := P_1(n) + 32(n+1)^5 = 32n^4 + 112n^3 + 152n^2 + 92n + 21.

Write `μ_* = m_0 + m_1 s + m_2 s^2 + m_3 s^3 + ρ_4`, `(m_0..m_3) = (1, 1/2, 3/8, 5/16)`, `|ρ_4| ≤ C s^4` on `[1/2,1]`.
The contribution of `ρ_4` is `O(n^4 β_4(n) + n^5 β_5(n)) = O(1/n)`. The contribution of the polynomial part is the
rational function

    Σ_{j≤3} m_j G(β_j)_n = (6712n^4 + 23412n^3 + 31402n^2 + 18297n + 4015) / (8(n+1)(n+2)(n+3)(n+4)(n+5))  ~ 839/n

(`u3_exact.py` (d), which asserts the strict degree drop: individually `G(β_0)_n ~ 16n^2`, `G(β_1)_n ~ −32n^2`,
`G(β_2)_n ~ −128n`, `G(β_3)_n → −576`, and the terms of order `n^2`, `n` **and the constant term** cancel exactly for
these `m_j`; the constant term cancels only because `m_3 = 5/16`). An `O(1)` bound would not suffice for Prop. 4.1. ∎

(Up to scaling and adding multiples of `ν` — which does not change `μ_*` up to scaling, as `νθν − νθν = 0` — the
choice `c_2 = −1/16` in `r_*` is the one for which this cancellation occurs: replacing `r_*` by `r_* + λ r` adds
`λ G(g)_n = 128 λ n + O(1)`.)

## 4. The value of the continued fraction in terms of g_*(0) and μ_4^*

**Proposition 4.1.** If `μ_4^* ≠ 0`, then `ℓ_3 = g_*(0)/(192 μ_4^*) − 4`; in particular `ℓ_3 < ∞`.

*Proof.* Apply the linear map `G` to (3.1) and use Lemma 2.2:

    G(g_*)_n = g_*(0) w^{(3)}_n − 192 μ_4^* (v^{(3)}_n + 4 w^{(3)}_n) = (g_*(0) − 768 μ_4^*) w^{(3)}_n − 192 μ_4^* v^{(3)}_n.

Divide by `w^{(3)}_n ≥ 1` and let `n → ∞`: the left side tends to 0 by Lemma 3.3 (ii), so
`192 μ_4^* · v^{(3)}_n/w^{(3)}_n → g_*(0) − 768 μ_4^*`. ∎

## 5. A Legendre-type relation

**Proposition 5.1.** `g(0) μ_4^* − μ_4 g_*(0) = 1/6`.

*Proof.* For sequences put `C[X,Y]_n := X_n Y_{n+1} − X_{n+1} Y_n` (bilinear, alternating). By `g = g(0) w − 192 μ_4 v`,
(3.1) and (1.1),

    C[g, g_*]_n = −192 (g(0) μ_4^* − μ_4 g_*(0)) · C[w,v]_n = −24 (g(0) μ_4^* − μ_4 g_*(0)) / (n+1)^5.

On the other hand `C[g,g_*]_n = g(n)(g_*(n+1) − g_*(n)) − (g(n+1) − g(n)) g_*(n)`, which by Lemma 3.3 (i) equals
`(−2/n^3)(−1/n^2) − (6/n^4)(1/n) + O(n^{−6}) = −4/n^5 + O(n^{−6})`. Multiply by `(n+1)^5` and let `n → ∞`. ∎
(Numerically `0.16666666666666666667`: `u3_numerics.py` (5).)

## 6. The bilateral sum for Z

Let `γ(z) := (Γ(z+1/2)/Γ(z+1))^4` (zeros of order 4 at `z = −1, −2, ...`; poles of order 4 at `z = −1/2, −3/2, ...`;
`γ(z) = tan^4(πz) γ(−z−1/2)` by the reflection formulas, as in U2, Prop. 8; `γ(u+k) = A_k(u)/N(u)`), and

    T_Z(u) := Σ_{k∈Z} γ(u+k)/(u+k+1).

**Proposition 6.1.** With `t := tan(πu)`,

    T_Z(u) = (1 + t^2)(e_0 + e_1 t + e_2 t^2),        e_0 = π^2 Z_0,     e_2 = e_0 + 8,

for a constant `e_1`, and

    ζ_0 = e_0,     ζ_1 = π e_1,     ζ_2 = π^2 (e_0 + e_2),     ζ_3 = (4π^3/3) e_1 − 16π^2.                     (6.1)

*Proof.* *Convergence, poles.* For `k ≥ 0` the terms are `O(k^{−3})`; for `k = −m ≤ −1` the term is
`tan^4(πu) γ(m−u−1/2)/(u−m+1) = O(m^{−3})`, locally uniformly. The simple pole of `1/(u+k+1)` at `u = −k−1` is
cancelled by the zero of `γ(u+k)`. So `T_Z` is meromorphic, of period 1, with poles at most at `1/2 + Z`, of order
`≤ 4`.
*Decomposition.* The part `k ≥ 0` is `Z(u)/N(u)`. In the part `k = −m`, put `m = j+1`, `v := 1/2 − u`:
`γ(m−u−1/2) = γ(v+j) = A_j(v)/N(v)` and `u − m + 1 = −(v + j − 1/2)`. Hence

    T_Z(u) = Z(u)/N(u) − tan^4(πu) · Z~(v)/N(v),        Z~(v) := Σ_{j≥0} A_j(v)/(v+j−1/2).                (6.2)

*At `u = 0`.* In (6.2) only the term `j = 0` of `Z~` is singular at `v = 1/2`: it contributes
`−tan^4(πu)/((v−1/2) N(v)) = tan^4(πu)/(u N(1/2−u)) = 16π^2 u^3 + O(u^4)` (`N(1/2) = Γ(3/2)^4 = π^2/16`); the rest is
`O(u^4)`. So `T_Z(u) = ζ_0 + ζ_1 u + ζ_2 u^2 + (ζ_3 + 16π^2) u^3 + O(u^4)`.
*At `±i∞`.* Uniformly for `Re u ∈ [0,1]`, each term is `O(|u+k|^{−3})` for `|Im u| ≥ 1` (Stirling for `k ≥ 0`; the
reflected form with `|tan πu| ≤ coth π` for `k < 0`), so `T_Z(u) = O(|Im u|^{−2}) → 0`.
*Liouville.* Exactly as in U2, Prop. 8: subtracting from `T_Z` the combination of `sec^2`, `sec^2 tan`, `sec^4` (in
`πu`) with the same principal part at `u = 1/2` up to order `v^{−2}`, and `c·π tan(πu)` for the `v^{−1}` term, leaves
an entire bounded 1-periodic function, with limits `∓ iπc` at `±i∞`; so `c = 0` and `T_Z` is a combination of
`sec^2 = 1+t^2`, `sec^2 tan = (1+t^2)t`, `sec^4 = (1+t^2)^2`, i.e. of the stated form.
*`e_0` and (6.1).* `t = πu + (πu)^3/3 + O(u^5)` gives
`(1+t^2)(e_0+e_1t+e_2t^2) = e_0 + πe_1 u + π^2(e_0+e_2)u^2 + (4π^3/3) e_1 u^3 + O(u^4)`; compare with the expansion at
`u = 0`. `ζ_0 = Z(0)/N(0) = π^2 Z_0`.
*`e_2`.* Let `u → 1/2`, i.e. `v → 0`, `t → ∞`: `T_Z(u)/t^4 → e_2`. In (6.2), `Z/N` is analytic at `u = 1/2`, and
`Z~(v)/N(v) → Z~(0)/N(0) = π^2 Σ_j A_j/(j − 1/2)` (uniform convergence near `v = 0`). So `e_2 = −π^2 Σ_j A_j/(j−1/2)`.
Finally, with `H(k) := A_k k^3/(k−1/2)` and `A_{k+1} = A_k ((k+1/2)/(k+1))^4`,

    H(k+1) − H(k) = A_k ( (k+1/2)^3/(k+1) − k^3/(k−1/2) ) = −(1/8) A_k ( 1/(k+1) + 1/(k−1/2) )

(`u3_exact.py` (e)), `H(0) = 0`, and `H(k) → 1/π^2` (`A_k ~ 1/(π^2 k^2)`). Summing over `k ≥ 0`:

    Z_0 + Σ_k A_k/(k−1/2) = −8/π^2,        i.e.        e_2 = π^2 Z_0 + 8.  ∎

(Numerically, including at a complex point: `u3_numerics.py` (4).)

## 7. End of the proof

As in U2 §7, `N(ε)N(−ε) = (ε cot πε)^4 = π^{−4}(1 − (4π^2/3)ε^2 + O(ε^4))`, `N(0) = π^{−2}`, and
`ζ(3,1/2+ε) = 7ζ(3) − (π^4/2)ε + 186ζ(5)ε^2 − (2π^6/3)ε^3 + O(ε^4)`. Write
`N(ε)Y_*(−ε) = −N(ε)ε^3/2 + N(ε)N(−ε)( −(Y/N)(−ε)/16 + (Z/N)(−ε)/32 )` and insert Prop. 8 of U2 and (6.1).
Coefficient extraction (`u3_algebra.py`, exact) gives, with `E_1 := e_1/π^3 − 2τ_1/π^4`:

    192 μ_4^* = −32 [ε^1](N(ε)Y_*(−ε)) = E_1,
    g_*(0)   = −2 [ε^3](Ξ(ε)Y_*(−ε)) = (93/4) ζ(5) E_1 + (π^2/16)(e_0 + e_2 − 4τ_0)

(the terms in `ζ(3)` cancel; also `[ε^3](N(ε)Y_*(−ε)) = 0`). By Prop. 6.1,
`e_0 + e_2 − 4τ_0 = 2π^2 Z_0 + 8 − 4π^2 S_0 = 2K`, so

    g_*(0) = (93/4) ζ(5) E_1 + π^2 K/8.                                                                   (7.1)

Insert (7.1), `μ_4^* = E_1/192`, `μ_4 = 4S_1/(3π^2)` and `g(0) = 192 μ_4 ((93/4)ζ(5) + π^6 S_0/(64 S_1))` into
Prop. 5.1: the terms with `ζ(5)` cancel and

    1/6 = μ_4 ( E_1 π^6 S_0/(64 S_1) − π^2 K/8 ) = E_1 π^4 S_0/48 − S_1 K/6,        i.e.    E_1 = 8 (1 + S_1 K)/(π^4 S_0).   (7.2)

*Signs.* `S_0 > 0`, `S_1 < 0` (U2 (1.1)), and `K < 0`: `2S_0 − Z_0 = Σ A_k (2 − 1/(k+1)) ≥ A_0 = 1 > 4/π^2`. Hence
`1 + S_1 K > 1`, `E_1 > 0`, `μ_4^* ≠ 0`. By Prop. 4.1, (7.1), (7.2):

    ℓ_3 = g_*(0)/E_1 − 4 = (93/4) ζ(5) − 4 + (π^2 K/8)/E_1 = (93/4) ζ(5) − 4 + (π^6/64) S_0 K/(1 + S_1 K),

and `1/PCF = ℓ_3` by Lemma 2.1 (which also gives `ℓ_3 > 0`). The other two forms: `S_0 = 16L(f,2)/π^2`,
`S_1 = −2L(f,1)` (`../Z/THEOREM.md`), and `1/U2 = (93/4)ζ(5) + π^6 S_0/(64 S_1)`. ∎

By-products: `e_1 = 8(1 + 4S_1 + π^2 S_1 Z_0)/(π S_0)`, i.e. the first four Taylor coefficients of `Z/N` at 0 are
explicit in `S_0, S_1, Z_0, π`; and `Σ_k A_k/(k+1) = 2 Σ_k A_k/(1−2k) − 8/π^2`.

## 8. Remarks

1. **Numerics.** `K = −7.7834942513...`, `1 + S_1 K = 6.5185080678...`, `1/PCF = 0.0437809038563...`; the formula agrees with the
   48-digit value of the fraction computed from the recurrence (`u3_numerics.py` (6)).
2. **What `Z_0` is (not proved to be new).** `Z_0 = ∫_0^1 y_0(x) dx`. With `a_n = ((1/2)_n/n!)^2`,
   `K(k) = (π/2)Σ a_n k^{2n}`, `∫_0^1 K'(k) k^{2n} dk = (π^2/4) a_n` (`../Z/THEOREM.md`, Lemma 2) and
   `∫_0^k m K(m) dm = E(k) − (1−k^2) K(k)`, one gets `Z_0 = (16/π^3) ∫_0^1 K'(k) (E(k) − (1−k^2)K(k)) k^{−2} dk`
   (numerically checked). In the language of Bönisch–Klemm–Scheidegger–Zagier (arXiv:2203.09426, §4) this is the
   kind of integral that produces quasi-periods of `f` (integrals of `f` times a modular function with a pole at
   a cusp), and (7.2) would then be the corresponding Legendre relation between periods and quasi-periods — this
   reading is an interpretation, not something proved here; PSLQ (60 digits)
   finds no relation between `Z_0`, `L(f,1)`, `L(f,2)`, `π`, `1`. No independence statement is proved.
3. **Structure.** For every solution `ρ` of `L4` holomorphic at the conifold point, `g_ρ(n) = ∫_0^1 (νθρ − ρθν)x^n dx`
   solves the U2 recurrence. `ρ = r` gives its minimal solution (U2); the gauge map `G` sends the U2 recurrence to
   the U3 recurrence and exchanges the roles: `G(g_{r_*})` is `O(1/n)` while `G(g)_n ~ 128 n`; presumably
   `G(g_{r_*})` is the minimal solution of the U3 recurrence (numerically it decays like `n^{−5}`), but minimality
   is not proved here and is not needed: Prop. 4.1 uses only `G(g_*)_n/w^{(3)}_n → 0`. (Numerically the convergents
   approach the limit like `n^{−6}`; no rate is claimed.) The Casoratian of `g_r` and `g_{r_*}` is the invariant pairing
   of `r` and `r_*`, which is how the quadratic relation (7.2) arises without any period computation.
