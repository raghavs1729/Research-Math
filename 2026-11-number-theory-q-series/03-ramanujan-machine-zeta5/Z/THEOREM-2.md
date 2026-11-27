# Zudilin's relations (6) for (r,t) = (1/2,1/4) and (1/2,1/3)

## 0. Source, statements, status

**Source.** W. Zudilin, *A hypergeometric version of the modularity of rigid Calabi–Yau manifolds*, SIGMA 14 (2018),
086, arXiv:1805.00544 (read 2026-11-26), eq. (5), (6), Table 1, Remark 2 ("None of the relations in (6) seem to be
proved"); quoted in `THEOREM.md` §0. For `r = 1/2` and `t ∈ (0, 1/2)` put, for `ε` near 0,

    F(ε) := (1/(π Γ(t)Γ(1−t))) Σ_{n≥0} Γ(1/2+ε+n)^2 Γ(t+ε+n) Γ(1−t+ε+n) / Γ(1+ε+n)^4 = F_0 + F_1 ε + F_2 ε^2 + F_3 ε^3 + O(ε^4)

(Zudilin's `F_j(1)` for `(r,t) = (1/2,t)`). Table 1 of the paper gives the weight-4 newforms
`f = η_4^16/(η_2^4 η_8^4)` (level 16) for `(1/2,1/4)` and `f = η_6^14/(η_2^3 η_18^3) − 3 η_2^3 η_6^2 η_18^3` (level 36) for
`(1/2,1/3)` (`η_m := η(mτ)`).

> **Theorem Z2 (`(r,t) = (1/2,1/4)`).** Let `f = η_4^16/(η_2^4 η_8^4)` and `f_2 = η_2^4 η_4^4` (level 8). Then
>
>     F_1 = −(32/π) L(f_2,2) = −8 L(f,1),      F_2 = 32π L(f_2,1) = 32 L(f,2),      F_3 = (7π^2/3) F_1 = −(448/3) L(f,3),
>
> i.e. `L(f,1)/F_1 = −1/8`, `L(f,2)/F_2 = 1/32`, `L(f,3)/F_3 = −3/448`: the three relations (6) hold.
> In terms of the `(1/2,1/2)` constants of `THEOREM.md`: `F_1 = −2π S_0`, `F_2 = −16π S_1`.

> **Theorem Z3 (`(r,t) = (1/2,1/3)`).** Let `f = η_6^14/(η_2^3 η_18^3) − 3 η_2^3 η_6^2 η_18^3` and
> `g := η_2^6 η_6^2 + 9 η_2^3 η_6^2 η_18^3 = a(q^2) η(2τ)^3 η(6τ)^3` (level 12; `a(q) = Σ_{m,n} q^{m^2+mn+n^2}`). Then `f = g ⊗ χ_{−3}`,
>
>     F_1 = −(12√3/π) L(g,2) = −3 L(f,1),      F_2 = 8√3 π L(g,1) = 18 L(f,2),      F_3 = (5π^2/3) F_1 = −90 L(f,3),
>
> i.e. `L(f,1)/F_1 = −1/3`, `L(f,2)/F_2 = 1/18`, `L(f,3)/F_3 = −1/90`: the three relations (6) hold.

In both cases `F_0` is *not* an L-value but an incomplete period:
`F_0^{(1/2,1/4)} = 64 ∫_{1/√8}^∞ f_2(it) dt`, `F_0^{(1/2,1/3)} = 12√3 ∫_{1/√3}^∞ g(it/2) dt` (cf. Zudilin's Remark 2).

**Status (updated 2026-11-27).** Refereed: `REFEREE-2-Z2Z3.md` (independent adversarial referee) reports "PROOF COMPLETE" for both theorems, with ratios reproduced from Zudilin's definition to 52–54 digits; its minor findings (levels in Lemma 4 to be taken termwise at level 36 — `verify/z2_symbolic.py` (S5) now checks the Gordon–Hughes–Newman and Ligozat conditions at level 36; the Wronskian sentence in Prop. B(i); the Stirling estimate in Lemma 3; closing arcs in Prop. C; "F_0 is not an L-value" is a heuristic) are addressed in the write-up R. Sharma, "Closed forms for the unidentified ζ(5)-level continued fractions of the Ramanujan Machine" (2026). Original status note:  Proof written below and checked by the author's own scripts only. **The adversarial audit by a fresh
reader was started and then stopped unfinished on the parent session's instruction (2026-11-26); no independent
referee has read this proof.** What the unfinished audit left: an independent evaluation of `F_0..F_3` for both cases
(`audit2/a1_frobenius.py`, 38 digits, agreeing with the values used here and with `F_3 = (csc^2 πt + 1/3)π^2F_1`); it
reported no error before it was stopped, and it had not yet examined the proof steps. Machine checks: `verify/z2_symbolic.py` (the three classical transformations in
operator form, a Taylor expansion, cusp orders of the eta quotients; exact), `verify/z2_numeric.py` (every step
numerically to 15–30 digits: the parametrisations, the three integral representations, the Fricke and
`γ`-transformation laws, the twist identities, the final ratios).
**Openness / novelty.** `LITERATURE-13.md` (survey of 2026-11-26 by a subagent: Zudilin; Rogers–Wan–Zucker;
Long–Tu–Yui–Zudilin; Bönisch–Klemm–Scheidegger–Zagier; Bönisch–Golyshev–Klemm; Au; Bao; Allen–Grove–Long–Tu and
others) found the relations with `F_k(1)` proved in *none* of the 13 cases other than `(1/2,1/2)`; for `(1/2,1/4)`
the L-values of `f` were known in terms of `S_0`, `S_1` (Rogers–Wan–Zucker (30), (31), reproved in §6 by a different
argument) and the relation with `F_k(1)` was numerical; for `(1/2,1/3)` everything was numerical. The ingredients
used here are classical (Clausen, Pfaff, Whipple, Jacobi and Borwein theta functions, Manin's trick for twists);
"new" means "not found in that survey". The constants agree with the survey's independent PARI computation.
Not treated: the other 11 cases (§7).

---

## 1. Three lemmas for r = 1/2

Let `c(s) := Γ(s+1/2)/(Γ(1/2)Γ(s+1))`, `u_t(s) := Γ(s+1/2)Γ(s+t)Γ(s+1−t)/(Γ(1/2)Γ(t)Γ(1−t)Γ(s+1)^3)`, so that
`F(ε) = Σ_{n≥0} c(n+ε) u_t(n+ε)`, and

    U(x,ε) := Σ_{n≥0} u_t(n+ε) x^{n+ε}  (0 < x < 1),      U_j(x) := [ε^j] U(x,ε),
    D_t := θ^3 − x(θ+1/2)(θ+t)(θ+1−t),   θ = x d/dx      (the operator of 3F2(1/2,t,1−t;1,1;x)).

**Lemma 1 (integral representation).** `F_j = (1/π) ∫_0^1 U_j(x) dx/√(x(1−x))` for `j ≥ 0`.

*Proof.* `c(s) = (1/π) ∫_0^1 x^s dx/√(x(1−x))` for `Re s > −1/2` (Beta integral). For complex `|ε| ≤ 1/8`
(inside the disc of holomorphy `|ε| < t` of `u_t`), `|u_t(n+ε)| ≤ C (n+1)^{−3/2}` uniformly, so `Σ_n |u_t(n+ε)| x^{n−1/8} ≤ C' x^{−1/8}`, which is integrable against
`dx/√(x(1−x))`. By Fubini `F(ε) = (1/π)∫_0^1 U(x,ε) dx/√(x(1−x))`, both sides are analytic in `ε`, and the coefficients
may be taken under the integral sign (Cauchy's formula). ∎

**Lemma 2 (the functions `U_0, U_1, U_2`).** `D_t U_j = 0` on `(0,1)` for `j = 0, 1, 2`; as `x → 0+`,

    U_0 = 1 + o(1),     U_1 = log(x/M_t) + o(1),     U_2 = (1/2) log^2(x/M_t) + c_t + o(1),

with `M_{1/4} = 256`, `M_{1/3} = 108`, `c_t = π^2/(2 sin^2 πt)` (`c_{1/4} = π^2`, `c_{1/3} = 2π^2/3`). A solution of `D_t v = 0` on
`(0,δ)` with `v(x) → 0` as `x → 0+` vanishes identically.

*Proof.* `u_t(s+1)(s+1)^3 = u_t(s)(s+1/2)(s+t)(s+1−t)` gives `D_t U = ε^3 u_t(ε) x^ε`, whose coefficients of `ε^0, ε^1, ε^2`
vanish. `U(x,ε) = u_t(ε) x^ε (1 + x R(x,ε))` with `R` bounded for `|ε| ≤ 1/8`, `x ≤ 1/2`, and
`log u_t(ε) = (ψ(1/2)+ψ(t)+ψ(1−t)−3ψ(1)) ε + (1/2)(ψ'(1/2)+ψ'(t)+ψ'(1−t)−3ψ'(1)) ε^2 + O(ε^3) = −ε log M_t + c_t ε^2 + O(ε^3)`
(`ψ(1/2) = ψ(1) − 2 log 2`; `ψ(t)+ψ(1−t) − 2ψ(1) = −6 log 2`, resp. `−3 log 3`; `ψ'(t)+ψ'(1−t) = π^2/sin^2 πt`,
`ψ'(1/2) = 3ψ'(1) = π^2/2`). The `U_j` are independent (leading terms `1, log, log^2`) and span the 3-dimensional
solution space, and no nontrivial combination tends to 0. ∎

**Lemma 3 (bilateral sum; cf. Zudilin's Remark 5).** For `0 < t < 1/2`: `F_3 = (csc^2 πt + 1/3) π^2 F_1`
(`= (7π^2/3)F_1` for `t = 1/4`, `(5π^2/3)F_1` for `t = 1/3`).

*Proof.* Let `T(u) := (1/(πΓ(t)Γ(1−t))) Σ_{n∈Z} Γ(1/2+u+n)^2 Γ(t+u+n)Γ(1−t+u+n)/Γ(1+u+n)^4`. For `n → +∞` the terms are
`O(n^{−2})`. For `n = −m`, the reflection formula gives `Γ(ρ+u−m)/Γ(1+u−m) = −(sin πu/sin π(u+ρ)) Γ(m−u)/Γ(m+1−ρ−u)`, so the
term is `(sin^4 πu / ∏_ρ sin π(u+ρ)) ∏_ρ Γ(m−u)/Γ(m+1−ρ−u)`, `ρ ∈ {1/2, 1/2, t, 1−t}`, which is `O(m^{−2})`. Hence `T` is
meromorphic, 1-periodic, with at most a double pole at `u ≡ 1/2` and simple poles at `u ≡ ±t` (mod 1); the terms with
`n < 0` vanish to order 4 at `u = 0`, so `T(u) = F(u) + O(u^4)`; and `T(u) = O(1/|Im u|)` as `|Im u| → ∞` uniformly in
`Re u` (each term is `O(|u+n|^{−2})`, the trigonometric factor being bounded). As in `../U2/THEOREM.md`, Prop. 8,
subtracting principal parts leaves a bounded entire periodic function, so

    T(u) = A π^2 sec^2 πu + B π cot π(u+1/2) + C_1 π cot π(u+t) + C_2 π cot π(u−t) + K,

and the limits at `±i∞` (`cot → ∓i`) give `K = 0`, `B + C_1 + C_2 = 0`. With `S_ρ := csc^2 πρ`,
`π cot π(u+ρ) = π cot πρ − π^2 S_ρ u + π^3 cot πρ · S_ρ u^2 − (π^4/3) S_ρ(3S_ρ − 2) u^3 + ...` (`z2_symbolic.py` (S4)), and
`sec^2` is even. Hence `F_1 = −π^2 (B + (C_1+C_2) S_t) = −π^2 (C_1+C_2)(S_t − 1)` and
`F_3 = −(π^4/3)(B + (C_1+C_2) S_t (3S_t−2)) = −(π^4/3)(C_1+C_2)(S_t−1)(3S_t+1)`. ∎

## 2. Classical facts used

Jacobi thetas `θ_2, θ_3, θ_4` at `σ = iy` (nome `e^{−πy}`), `y > 0`; `k := θ_2^2/θ_3^2`, `k' = θ_4^2/θ_3^2`.
(J1) `(2/π)K(k) = θ_3^2`, `K'(k)/K(k) = y`, `dk/dy = −(π/2) k θ_4^4` (`THEOREM.md` (Z2)); `K(√λ)` and `K(√(1−λ))` solve
`H_{1/2}: λ(1−λ)w'' + (1−2λ)w' − w/4 = 0`. So `θ_3^2` and `y θ_3^2` are solutions of `H_{1/2}` as functions of `λ = k^2`.
(J2) duplication: `θ_2^2(σ) = (θ_3^2 − θ_4^2)(σ/2)/2`, `θ_3^2(σ) = (θ_3^2 + θ_4^2)(σ/2)/2`, `θ_4^2(σ) = θ_3θ_4(σ/2)`;
`θ_3^4 = θ_2^4 + θ_4^4`. Consequently (Landen) `k(2y) = (1−k'(y))/(1+k'(y))`, and `k'(1/y) = k(y)`.
(J3) `f_2(τ) = η(2τ)^4η(4τ)^4 = θ_2^4θ_3^2θ_4^2(2τ)/16`, `f_2(−1/(8τ)) = 64τ^4 f_2(τ)` (`THEOREM.md`, Lemma 1).
Borweins' cubic thetas at `τ = iy`, `q = e^{−2πy}`: `a = Σ q^{m^2+mn+n^2}`, `b = η(τ)^3/η(3τ)`, `c = 3η(3τ)^3/η(τ)`.
(B1) `a^3 = b^3 + c^3`; `a(q) = b(q) + 3c(q^3)`; with `h := c^3/a^3 ∈ (0,1)`: `2F1(1/3,2/3;1;h) = a` for all `y > 0`.
(B2) `a(−1/(3τ)) = −i√3 τ a(τ)` and `b(−1/(3τ)) = −i√3 τ c(τ)` (the latter from `η(−1/τ) = (τ/i)^{1/2}η(τ)`); hence, with
(B1), `h(1/(3y)) = 1 − h(y)`.
(J. M. Borwein, P. B. Borwein, F. G. Garvan, *Some cubic modular identities of Ramanujan*, Trans. AMS 343 (1994);
B. C. Berndt, S. Bhargava, F. G. Garvan, Trans. AMS 347 (1995), §2. All are re-checked numerically in
`z2_numeric.py`.) By (B1), (B2), `a` and `2F1(1/3,2/3;1;1−h) = √3 y a` are solutions of
`H_{1/3}: h(1−h)w'' + (1−2h)w' − (2/9)w = 0` as functions of `h`.
Transformations, in operator form (`z2_symbolic.py` (S1)–(S3), exact):
(T1, Clausen) if `w_1, w_2` solve `H_t` (variable `h`), then `w_1w_2` is annihilated by `D_t` in the variable `x = 4h(1−h)`;
(T2, Pfaff) if `w` solves `H_{1/2}`, so does `(1−λ)^{−1/2} w(λ/(λ−1))`;
(T3, Whipple) if `v` solves `D_{1/4}` (variable `X`), then `(1−z)^{−1/2} v(−4z/(1−z)^2)` solves `D_{1/2}` (variable `z`).

## 3. The parametrisation for t = 1/4

For `y > 0` put `P := θ_3^4 + θ_2^4`, `z := −4θ_2^4θ_3^4/θ_4^8 = −4k^2/k'^4`, `X := −4z/(1−z)^2 = 16θ_2^4θ_3^4θ_4^8/P^4`
(`1 − z = P^2/θ_4^8`, since `θ_4^8 + 4θ_2^4θ_4^4 + 4θ_2^8 = (θ_4^4 + 2θ_2^4)^2`).

**Proposition A.** (i) `y ↦ X` is a decreasing bijection `(√2, ∞) → (0,1)`, and `log(X/256) = −πy + O(e^{−πy})`.
(ii) On `y > √2`: `U_0(X) = P`, `U_1(X) = −πy P`, `U_2(X) = (π^2y^2/2 + π^2) P`.
(iii) `dX/dy = −4π θ_2^2θ_3^2θ_4^4 √(X(1−X))/P`.

*Proof.* (i), (iii): by (J1), `d log k^2/dy = −πθ_4^4`, so `dz/dy = z · (d log k^2/dy)(1+k^2)/(1−k^2) = −πzP > 0`; `z` increases
from `−∞` to 0, and `z = −1` iff `k'^2 = 2k` iff `k = √2 − 1`. By (J2), `k(√2) = k(2/√2) = (1−k'(1/√2))/(1+k'(1/√2)) =
(1−k(√2))/(1+k(√2))`, so `k(√2) = √2 − 1`: `z(√2) = −1`. On `z ∈ (−1,0)`, `dX/dz = −4(1+z)/(1−z)^3 < 0` and
`√(X(1−X)) = 2√(−z)(1+z)/(1−z)^2`; hence `dX/dy = −2π√(−z) P √(X(1−X))/(1−z)`, which is (iii) because
`√(−z) = 2θ_2^2θ_3^2/θ_4^4`. The asymptotics follow from `θ_2^4 = 16e^{−πy}(1+O(e^{−2πy}))`.
(ii): by (J1) and (T2) with `λ = k^2`, `λ/(λ−1) = −θ_2^4/θ_4^4 =: λ_−`, the functions
`(1−λ_−)^{−1/2}θ_3^2 = θ_4^2` and `yθ_4^2` solve `H_{1/2}` in the variable `λ_−` ((T2) is an involution). By (T1), `θ_4^4·{1, y, y^2}`
are annihilated by `D_{1/2}` in the variable `z = 4λ_−(1−λ_−)`. By (T3), `v ↦ (1−z)^{−1/2} v(X(z))` maps the solution space of
`D_{1/4}` on `X ∈ (0,1)` injectively, hence (dimension 3) bijectively, onto that of `D_{1/2}` on `z ∈ (−1,0)`; so
`(1−z)^{1/2}θ_4^4·{1,y,y^2} = P·{1,y,y^2}` are solutions of `D_{1/4}` as functions of `X`. By (i), as `y → ∞`:
`P → 1`, `−πyP − log(X/256) → 0`, `(π^2y^2/2 + π^2)P − (1/2)log^2(X/256) − π^2 → 0`; Lemma 2 gives (ii). ∎

**Corollary A.** With `p_0 = 1`, `p_1 = −πy`, `p_2 = π^2y^2/2 + π^2`:

    F_j = 4 ∫_{√2}^∞ p_j(y) θ_2^2θ_3^2θ_4^4(iy) dy = 64 ∫_{1/√8}^∞ p_j(4t) f_2(it) dt        (j = 0, 1, 2).

*Proof.* Lemma 1, Proposition A, and `θ_2^2θ_3^2θ_4^4(σ) = (1/4)θ_2^4θ_3^2θ_4^2(σ/2) = 4f_2(σ/4)` by (J2), (J3); `y = 4t`. ∎

## 4. Proof of Theorem Z2, first part

Let `ρ_m := ∫_0^∞ f_2(it) t^m dt = m!(2π)^{−m−1} L(f_2, m+1)` (`THEOREM.md`, Lemma 1 (ii)). The substitution `t ↦ 1/(8t)` and
(J3) give `∫_0^{1/√8} f_2(it) t dt = ∫_{1/√8}^∞ f_2(it) t dt` and `∫_0^{1/√8} f_2(it) dt = ∫_{1/√8}^∞ 8t^2 f_2(it) dt`. Hence

    F_1 = −256π ∫_{1/√8}^∞ t f_2(it) dt = −128π ρ_1 = −(32/π) L(f_2,2),
    F_2 = 64π^2 ∫_{1/√8}^∞ (8t^2 + 1) f_2(it) dt = 64π^2 ρ_0 = 32π L(f_2,1),

and `F_3 = (7π^2/3)F_1` by Lemma 3. With `L(f_2,2) = π^2S_0/16`, `L(f_2,1) = −S_1/2` (`THEOREM.md`): `F_1 = −2πS_0`,
`F_2 = −16πS_1`. The values `L(f,k)` are computed in §6.

## 5. The case t = 1/3: parametrisation and proof of Theorem Z3, first part

For `y > 0` put `x := 4h(1−h)`, `h = c^3/a^3`, and `g_0(τ) := a(q) η(τ)^3 η(3τ)^3`.

**Proposition B.** (i) `dh/dy = −2π h(1−h) a^2`; `y ↦ x` is a decreasing bijection `(1/√3, ∞) → (0,1)`;
`log(x/108) = −2πy + O(e^{−2πy})`.
(ii) On `y > 1/√3`: `U_0(x) = a^2`, `U_1(x) = −2πy a^2`, `U_2(x) = (2π^2y^2 + 2π^2/3) a^2`.
(iii) `dx/dy = −12√3 π η(iy)^3η(3iy)^3 √(x(1−x))/a`.

*Proof.* (i) `w_1 = a`, `w_2 = √3 y a` solve `H_{1/3}` in `h`, so their Wronskian `√3 a^2 dy/dh` equals `C/(h(1−h))` (Abel).
As `y → ∞`, `h = 27q(1+O(q))` and `y = −(log q)/(2π)`, so `h dy/dh → −1/(2π)` and `C = −√3/(2π)`. Thus `h` decreases, from 1
to 0; `h(1/√3) = 1/2` by (B2); on `h ∈ (0,1/2)`, `x = 4h(1−h)` increases with `h`. (iii) `dx/dy = 4(1−2h)dh/dy`,
`√(x(1−x)) = 2√(h(1−h))(1−2h)`, `√(h(1−h)) = (bc)^{3/2}/a^3 = 3√3 η(τ)^3η(3τ)^3/a^3`. (ii) By (T1), `a^2·{1,y,y^2}` are
annihilated by `D_{1/3}` in the variable `x`; compare asymptotics as in Proposition A, using (i) and Lemma 2. ∎

So, by Lemma 1, with `p_0 = 1`, `p_1 = −2πy`, `p_2 = 2π^2(y^2 + 1/3)`: `F_j = 12√3 ∫_{1/√3}^∞ p_j(y) g_0(iy) dy`.

By (B1), `g_0 = (b(q) + 3c(q^3))η_1^3η_3^3 = η_1^6η_3^2 + 9η_1^3η_3^2η_9^3`, so `g(τ) := g_0(2τ)` is the form of the statement;
write `g = Σ α_n q^n` (`α_n = 0` for even `n`). By (B2) and `η(−1/τ) = (τ/i)^{1/2}η(τ)`: `g_0(−1/(3τ)) = 9τ^4 g_0(τ)`. Hence
`g_0(iy) = O(e^{−πy})` at `∞` and `O(y^{−4}e^{−π/(3y)})` at 0, `ρ_m := ∫_0^∞ g_0(iy) y^m dy = m! π^{−m−1} L(g, m+1)` (as in
`THEOREM.md`, Lemma 1 (ii); this defines the continuation of `L(g,s) = Σ α_n n^{−s}`), and `y ↦ 1/(3y)` gives
`∫_0^{1/√3} g_0 y dy = ∫_{1/√3}^∞ g_0 y dy`, `∫_0^{1/√3} g_0 dy = ∫_{1/√3}^∞ 3y^2 g_0 dy`, `ρ_2 = ρ_0/3`. Therefore

    F_1 = −24√3 π ∫_{1/√3}^∞ y g_0(iy) dy = −12√3 π ρ_1 = −(12√3/π) L(g,2),
    F_2 = 24√3 π^2 ∫_{1/√3}^∞ (y^2 + 1/3) g_0(iy) dy = 8√3 π^2 ρ_0 = 8√3 π L(g,1),

and `F_3 = (5π^2/3) F_1` by Lemma 3. (`g_0(iy) > 0`, so `ρ_0, ρ_1 > 0`.)

## 6. The twists: L(f,k) in terms of L(f_2,·), L(g,·)

**Lemma 4.** (a) `η_4^16/(η_2^4η_8^4) = f_2 ⊗ χ_{−4}`, and `f_2(τ/(4τ+1)) = −(4τ+1)^4 f_2(τ)`.
(b) `η_6^14/(η_2^3η_18^3) − 3η_2^3η_6^2η_18^3 = g ⊗ χ_{−3}`, and `g_0(τ/(3τ+1)) = −(3τ+1)^4 g_0(τ)`, `g_0(τ+1) = −g_0(τ)`.

*Proof.* Transformation laws: `f_2(τ+1/2) = −f_2(τ)` (only odd powers of `q`), `g_0(τ+1) = −g_0(τ)` (`η_1^3η_3^3` has the
factor `q^{1/2}`), and `τ ↦ τ/(Nτ+1)` is `W_N T^{−1} W_N^{−1}`, resp. `W_8 T^{−1/2} W_8^{−1}`, with `W_N: τ ↦ −1/(Nτ)`; conjugate
by the Fricke relations `f_2(−1/(8τ)) = 64τ^4f_2(τ)`, `g_0(−1/(3τ)) = 9τ^4g_0(τ)`. Identities of forms: both sides of (a)
lie in `M_4(Γ_0(16))` and both sides of (b) in `M_4(Γ_0(36))` (eta quotients satisfying the standard congruence
conditions, holomorphic at all cusps by Ligozat's formula — `z2_symbolic.py` (S5); twisting a form of level 8,
resp. 12, by a character of conductor 4, resp. 3, gives level 16, resp. 36); the Sturm bounds are 8 and 24, and
the `q`-expansions agree for `n < 400` (`z2_numeric.py`). ∎

**Proposition C.** (a) `L(f,1) = (4/π)L(f_2,2)`, `L(f,2) = π L(f_2,1)`, `L(f,3) = (π/2) L(f_2,2)` for `f = f_2 ⊗ χ_{−4}`.
(These are (30), (31) of Rogers–Wan–Zucker, arXiv:1303.2259, proved there by other means.)
(b) `L(f,1) = (4√3/π) L(g,2)`, `L(f,2) = (4√3π/9) L(g,1)`, `L(f,3) = (2√3π/9) L(g,2)` for `f = g ⊗ χ_{−3}`.

*Proof.* (b) `χ_{−3}(n) = (ζ^n − ζ^{−n})/(i√3)`, `ζ = e^{2πi/3}`, so `f(τ) = (g(τ+1/3) − g(τ−1/3))/(i√3)`, and with `g(τ) = g_0(2τ)`,
`g_0(w±1) = −g_0(w)` and real coefficients: `f(it) = (2/√3) Im g_0(1/3 + 2it)`. Put `γ := (1 0; 3 1)`; `γ(0) = 0`, `γ(∞) = 1/3`,
and `g_0(1/3 + 2it) = −(i/(6t))^4 g_0(γ^{−1}(1/3+2it))` with `Im γ^{−1}(1/3+2it) = 1/(18t)`: so `f(it)` decays exponentially
as `t → 0` as well, and `Λ_f(s) := ∫_0^∞ f(it)t^{s−1}dt = Γ(s)(2π)^{−s}L(f,s)` for all `s`. For `j = 0, 1, 2`,

    Λ_f(j+1) = (2/√3) 2^{−j−1} Im I_j,     I_j := ∫_0^∞ g_0(1/3 + iy) y^j dy = i^{−j−1} ∫_{1/3}^{i∞} g_0(w)(w − 1/3)^j dw.

By Cauchy's theorem on the ideal triangle with vertices `0, 1/3, i∞` (the integrand decays exponentially at the
three cusps), `∫_{1/3}^{i∞} = ∫_0^{i∞} − ∫_0^{1/3}`, the last integral along `γ(iR_+)`. With `w = γτ`: `w − 1/3 = −1/(3(3τ+1))`,
`dw = dτ/(3τ+1)^2`, `g_0(w) = −(3τ+1)^4 g_0(τ)`, so `∫_0^{1/3} g_0(w)(w−1/3)^j dw = (−1)^{j+1} 3^{−j} ∫_0^{i∞} g_0(τ)(3τ+1)^{2−j} dτ`.
With `r_m := ∫_0^{i∞} g_0(τ)τ^m dτ = i^{m+1}ρ_m` and `ρ_2 = ρ_0/3`:

    j = 0:  ∫_{1/3}^{i∞} = r_0 + (9r_2 + 6r_1 + r_0) = −iρ_0 − 6ρ_1,                 I_0 = −ρ_0 + 6iρ_1;
    j = 1:  ∫_{1/3}^{i∞} = (r_1 − r_0/3) − (1/3)(3r_1 + r_0) = −(2/3)iρ_0,            I_1 = (2/3)iρ_0;
    j = 2:  ∫_{1/3}^{i∞} = (r_2 − (2/3)r_1 + r_0/9) + r_0/9 = −iρ_0/9 + (2/3)ρ_1,     I_2 = ρ_0/9 + (2/3)iρ_1.

Hence `Λ_f(1) = 2√3ρ_1`, `Λ_f(2) = ρ_0/(3√3)`, `Λ_f(3) = ρ_1/(6√3)`; insert `ρ_0 = L(g,1)/π`, `ρ_1 = L(g,2)/π^2` and
`Λ_f(s) = Γ(s)(2π)^{−s}L(f,s)`.
(a) Likewise `f(it) = Im f_2(1/4 + it)`, `γ := (1 0; 4 1)`, `w − 1/4 = −1/(4(4τ+1))`, `f_2(γτ) = −(4τ+1)^4f_2(τ)`,
`ρ_m := ∫_0^∞ f_2(it)t^m dt`, `ρ_2 = ρ_0/8`, `Λ_f(j+1) = Im I_j`, `I_j = i^{−j−1}∫_{1/4}^{i∞} f_2(w)(w−1/4)^j dw`:

    j = 0:  2r_0 + 8r_1 + 16r_2 = −8ρ_1,   I_0 = 8iρ_1;      j = 1:  −r_0/2,   I_1 = iρ_0/2;      j = 2:  r_2 − r_1/2 + r_0/8 = ρ_1/2,   I_2 = iρ_1/2.

So `Λ_f(1) = 8ρ_1`, `Λ_f(2) = ρ_0/2`, `Λ_f(3) = ρ_1/2`, with `ρ_0 = L(f_2,1)/(2π)`, `ρ_1 = L(f_2,2)/(4π^2)`. ∎

**End of the proofs.** Z2: `L(f,1)/F_1 = (4/π)/(−32/π) = −1/8`, `L(f,2)/F_2 = π/(32π) = 1/32`,
`L(f,3)/F_3 = (π/2)/((7π^2/3)(−32/π)) = −3/448`. Z3: `L(f,1)/F_1 = (4√3/π)/(−12√3/π) = −1/3`,
`L(f,2)/F_2 = (4√3π/9)/(8√3π) = 1/18`, `L(f,3)/F_3 = (2√3π/9)/((5π^2/3)(−12√3/π)) = −1/90`. (All `L`-values involved are
nonzero: `f_2(it), g_0(iy) > 0`.) ∎

## 7. Remarks

1. **Why these two cases.** For `r = 1/2` the summand factors as `c(s) u_t(s)`, and `c` has the moment functional of
   Lemma 1, which keeps the first *three* Frobenius coefficients inside the periods of the symmetric square;
   for `t = 1/4, 1/3` the function `3F2(1/2,t,1−t;1,1;x)` is a weight-2 modular form in a Fricke-symmetric
   Hauptmodul, the path `x ∈ (0,1)` is the segment from the Fricke point to `i∞`, and `U_1, U_2` are exactly the
   Fricke-odd and Fricke-even combinations, so `F_1`, `F_2` become complete Mellin integrals. `F_0` is the incomplete one.
2. **The form that appears is a twist of Zudilin's.** The integrals produce `f_2` (level 8) and `g` (level 12),
   and Zudilin's `f` is `f_2 ⊗ χ_{−4}`, resp. `g ⊗ χ_{−3}`; the odd twist exchanges the parity of the critical values,
   which is why `F_1` pairs with `L(f,1)` although it is an `s = 2` value of the untwisted form.
3. **Not reached.** `(1/2,1/6)`: `3F2(1/2,1/6,5/6;1,1;1728/j) = E_4^{1/2}`, and the integrand `η^{12}E_4^{−1/2}` is not a
   holomorphic modular form. The cases without `r = 1/2` need a different splitting; with two-plus-two factors
   only `F_0, F_1` stay inside periods. See `LITERATURE-13.md` for the status of all 13 cases (in 8 of them
   `L(f,2) = 0`, so the `k = 2` relation is trivially true).
