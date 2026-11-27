# The "awaiting proof" families of arXiv:2412.16818 (constant a_n, quadratic b_n; and an Euler family)

## 0. Source, statements, status

**Source (read 2026-11-24).** M. Shalyt, U. Seligmann, I. Beit Halachmi, O. David, R. Elimelech, I. Kaminer,
*Unsupervised Discovery of Formulas for Mathematical Constants*, arXiv:2412.16818v1 (NeurIPS 2024).
Convention: `PCF[a,b] = a_0 + b_1/(a_1 + b_2/(a_2 + ...))`.

* §4.1, eqs (8), (9): `PCF[4, 4n^2+2n] = 2G/(4G−3)`, `PCF[4, 4n^2+2n−2] = (4G−1)/(3G−2)` (`G` Gauss's constant),
  and a formula for `PCF[4, 4n^2−2n]` with the lemniscate constant; "*Their limits were numerically validated to a
  large precision, yet formal proofs for these formula hypotheses remain an open challenge.*"
* §4.3: "(F1) `a_n = i + 2j + 1, b_n = n^2 + (i + k)n`, with integers `i, j ≥ 0`, and `k ∈ {0, 1}`. This is expected to
  be related to π if k = 1, and to ln(2) if k = 0"; "(F2) `a_n = −2n + j − 1, b_n = −n^2 + jn + k` for integer j, k
  such that b_n has real roots that are not positive integers. This is expected to converge to a root of b_n.
  These are novel experimental results and mathematical hypotheses - awaiting proof."

**What is proved here.** All of the above, with explicit values (Theorems SR, F1, F2, G below). The paper's own
Table 5 proves only the case `j = 0, k = 0` of F1. **Standing:** low — the families are corollaries of the classical
Stieltjes–Rogers continued fraction (`∫_0^∞ e^{−zt} sech^m t dt = 1/(z + 1·m/(z + 2(m+1)/(z + ...)))`, Stieltjes 1890,
Rogers 1907) and of Euler's continued fraction with Gauss's summation theorem; I did not find the statements for
these parameter families written down as theorems, but nothing here is deep. Not independently audited.
**Checks:** `verify/sr_checks.py` (bracketing of the exact values by consecutive convergents; exact rational
evaluation of the F1 integrals for `i, j ≤ 4`; the F2 convergent identity for `|j|,|k| ≤ 4`; the three
Gauss/lemniscate values).

## 1. Theorem SR (Stieltjes–Rogers, shifted form)

> Let `z > 0`, `α ≥ 0`, `β > −1`, `m := β − α + 1`, `z + m > 0`, and
> `T_c := ∫_0^∞ e^{−zt} sech^m t · tanh^c t dt` (`c > −1`). Then
> `z + K_{n≥1} (n+α)(n+β)/z` converges, to `α T_{α−1}/T_α` if `α > 0` and to `1/T_0` if `α = 0`.

*Proof.* `d/dt (sech^m t tanh^c t) = c sech^m tanh^{c−1} − (m+c) sech^m tanh^{c+1}`. Integrating against `e^{−zt}`
(boundary terms vanish for `c > 0` because `z + m > 0`; for `c = 0` the boundary term at `t = 0` is `−1`):

    z T_c = c T_{c−1} − (m+c) T_{c+1}   (c > 0),          z T_0 = 1 − m T_1.

All `T_c > 0`. Put `ρ_n := (α+n) T_{α+n−1}/T_{α+n}` for `n ≥ 0` (`n ≥ 1` if `α = 0`). The first relation at `c = α+n`
gives `ρ_n = z + (n+1+α)(n+1+β)/ρ_{n+1}` (note `m + α + n = n + 1 + β`), and for `α = 0` the second gives
`1/T_0 = z + (1+α)(1+β)/ρ_1`. Hence, with `X := α T_{α−1}/T_α` (resp. `1/T_0`) and `b_n := (n+α)(n+β) > 0`,

    X = z + b_1/(z + b_2/( ... + b_N/ρ_N)),        ρ_N ∈ (z, ∞)   for every N.

A finite continued fraction with positive entries is monotone in its last denominator, so `X` lies between the
values obtained with `ρ_N = z` and `ρ_N = ∞`, i.e. between the convergents `C_N` and `C_{N−1}`. Finally the fraction
converges by the Seidel–Stern theorem (Lorentzen–Waadeland, *Continued Fractions* vol. 1, Thm 3.14): it is
equivalent to `z + K(1/d_n)` with `d_n > 0`, `d_n d_{n+1} = z^2/b_{n+1}`, and
`d_n + d_{n+1} ≥ 2z/√(b_{n+1})`, so `Σ d_n = ∞`. ∎

## 2. Theorem F1

> For integers `i, j ≥ 0`, `k ∈ {0,1}`: `PCF[i+2j+1, n^2+(i+k)n] = 1/T`, where
> `T = 2^{i+k+1} ∫_0^1 u^{2i+2j+k+1} (1+u^2)^{−(i+k+1)} du`, and
> * `k = 0`: `T = q + (−1)^j 2^i C(i+j, i) · log 2` with `q ∈ Q`;
> * `k = 1`: `T = q + (−1)^j 2^i C(i+j+1/2, i+1) · π` with `q ∈ Q`.
>
> In both cases the coefficient is a non-zero rational, so the limit is a non-constant Möbius transform of `log 2`,
> resp. `π`, with rational coefficients. Examples: `PCF[1, n^2] = 1/log 2`, `PCF[1, n^2+n] = 2/(π−2)`,
> `PCF[2, n^2+n] = 1/(2 log 2 − 1)`, `PCF[2, n^2+2n] = 4/(3π−8)`, `PCF[5, n^2+3n] = 12/(332 − 105π)`.

*Proof.* Theorem SR with `α = 0`, `β = i+k`, `z = i+2j+1`, `m = i+k+1`, and `u = e^{−t}`
(`sech t = 2u/(1+u^2)`): `T_0 = 2^m ∫_0^1 u^{z+m−1}(1+u^2)^{−m} du`.
`k = 0`: with `w = u^2`, `T = 2^i ∫_0^1 w^{i+j}(1+w)^{−i−1} dw`; expand `w^{i+j} = Σ_l C(i+j,l)(−1)^{i+j−l}(1+w)^l`; only
`l = i` produces a logarithm, with coefficient `2^i C(i+j,i)(−1)^j`.
`k = 1`: `T = 2^{i+2} ∫_0^1 u^{2M}(1+u^2)^{−R} du`, `M = i+j+1`, `R = i+2`. With `y = 1+u^2`,
`u^{2M} y^{−R} = Σ_l C(M,l)(−1)^{M−l} y^{l−R}`, and `∫_0^1 (1+u^2)^{−r} du = q_r + γ_r π/4` with `q_r ∈ Q`,
`γ_r = (1/2)_{r−1}/(r−1)!` (reduction formula). The coefficient of `π/4` is
`Σ_{r=1}^{R} C(M, R−r)(−1)^{M−R+r} γ_r = (−1)^{M−R+1} [y^{R−1}] (1+y)^M (1+y)^{−1/2} = (−1)^j C(M − 1/2, R − 1)`,
which is non-zero (a binomial coefficient with half-integer upper argument). ∎

## 3. Theorem F2

> Let `j, k` be real with `j^2 + 4k > 0`, `r_1 > r_2` the roots of `−n^2 + jn + k`, and assume `r_2 ∉ {1, 2, ...}`.
> For `a_n = −2n + j − 1`, `b_n = −n^2 + jn + k`:
> `p_n/q_n = r_1 + (r_2 − 1)/S_n`, `S_n = Σ_{l=0}^{n} ∏_{i=1}^{l} (i − r_1)/(i + 1 − r_2)` (whenever `S_n ≠ 0`), and
> `PCF[a,b] = r_2`, the smaller root. (If `r_1` is a positive integer the fraction terminates with the same value.)

*Proof.* `b_n = −h_1(n)h_2(n)`, `a_n = h_1(n) + h_2(n+1)` with `h_1(n) = −(n − r_1)`, `h_2(n) = −(n − r_2)`: this is
the Euler family of `../FR/THEOREM.md` (Theorem A with `H_1 = h_2`, `H_2 = h_1`), whose proof gives
`p_n/q_n = h_1(0) + h_2(1)/S_n`. `S_n` is the partial sum of `2F1(1−r_1, 1; 2−r_2; 1)`, which converges because
`(2−r_2) − (1−r_1) − 1 = r_1 − r_2 > 0`, to `Γ(2−r_2)Γ(r_1−r_2)/(Γ(1−r_2)Γ(1+r_1−r_2)) = (1−r_2)/(r_1−r_2) ≠ 0` (Gauss).
So `p_n/q_n → r_1 − (r_1 − r_2) = r_2`. ∎

## 4. Theorem G (eqs (8), (9) of the paper)

Let `ϖ = 2∫_0^1 dt/√(1−t^4) = Γ(1/4)^2/(2√(2π)) = 2.62205...` and `G = ϖ/π = 0.83462...` (Gauss's constant).

> `PCF[4, 4n^2+2n] = 2G/(4G−3)`,   `PCF[4, 4n^2+2n−2] = (4G−1)/(3G−2)`,   `PCF[4, 4n^2−2n] = 6/(4−ϖ)`.

*Proof.* `PCF[4, 4c_n] = 2·(2 + K(c_n/2))`, so Theorem SR applies with `z = 2` and
`(α,β) = (0,1/2), (1,−1/2), (0,−1/2)`, i.e. `m = 3/2, −1/2, 1/2`. With `u = e^{−t} = v^2` and
`I(a) := ∫_0^1 v^{a−1}(1+v^4)^{−1/2} dv`:

* `m = 1/2`: `T_0 = 2^{3/2} I(5)`;
* `m = 3/2`: `T_0 = 2^{5/2} ∫_0^1 v^6 (1+v^4)^{−3/2} dv = −2 + 6√2 I(3)` (integration by parts);
* `m = −1/2`: `T_0 = √2 (I(3) + I(7))`, `T_1 = √2 (I(3) − I(7))`.

Differentiating `v^c (1+v^4)^{1/2}`: `√2 = I(1) + 3 I(5) = 3 I(3) + 5 I(7)`. `I(1) = ϖ/(2√2)` (classical:
`I(1) = (1/2)∫_0^∞ (1+v^4)^{−1/2} dv = B(1/4,1/4)/8`). For `I(3)`: from `d/dv [v^{−1}(1+v^4)^{1/2}] =
(v^2 − v^{−2})(1+v^4)^{−1/2}` and `v ↦ 1/v` one gets `2 I(3) = √2 − C_0`,
`C_0 := ∫_0^∞ (1 − w^2 (1+w^4)^{−1/2}) dw = −B(3/4, −1/4)/4 = Γ(3/4)^2/√π = π/(√2 ϖ)` (analytic continuation in `a` of
`∫_0^∞ w^{a−1}(1+w^4)^{−1/2} dw = B(a/4, 1/2 − a/4)/4`). Hence `I(3) = (2G−1)/(2√2 G)`, and

    m = 1/2:  T_0 = (4 − ϖ)/3,              PCF = 2/T_0 = 6/(4 − ϖ);
    m = 3/2:  T_0 = 4 − 3/G,                PCF = 2/T_0 = 2G/(4G − 3);
    m = −1/2: T_0/T_1 = (2I(3) + √2)/(8I(3) − √2) = (4G−1)/(6G−4),   PCF = 2 T_0/T_1 = (4G−1)/(3G−2). ∎
