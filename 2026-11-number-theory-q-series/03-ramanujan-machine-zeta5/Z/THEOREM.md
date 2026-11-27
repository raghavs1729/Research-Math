# Zudilin's relations (6) for r = t = 1/2: critical L-values of η(2τ)^4 η(4τ)^4 and the Frobenius values F_j(1)

## 0. Source, statement, status

**Source (read 2026-11-26).** W. Zudilin, *A hypergeometric version of the modularity of rigid Calabi–Yau
manifolds*, SIGMA 14 (2018), 086, arXiv:1805.00544. For a pair `(r,t)` from the list of 14 hypergeometric
Calabi–Yau families, `F_0(z), ..., F_3(z)` are defined by eq. (5) there (Zudilin writes the bilateral sum over
`n ∈ Z` and notes that it equals the sum over `n ≥ 0` up to `O(ε^4)`, because of the factor `Γ(1+ε+n)^{−4}`; the
unilateral form is quoted),

    (1/(Γ(r)Γ(1−r)Γ(t)Γ(1−t))) Σ_{n≥0} Γ(r+ε+n)Γ(1−r+ε+n)Γ(t+ε+n)Γ(1−t+ε+n)/Γ(1+ε+n)^4 · z^{n+ε}
        = F_0(z) + F_1(z) ε + F_2(z) ε^2 + F_3(z) ε^3 + O(ε^4),

and eq. (6) reads: "*Then numerical calculations suggest conjectural inclusions*
`L(f,1)/F_1(1) ∈ Q, L(f,2)/F_2(1) ∈ Q and L(f,3)/F_3(1) ∈ Q.`" Remark 2: "*None of the relations in (6) seem to be
proved.*" For `(r,t) = (1/2,1/2)` the eigenform is `f(τ) = η(2τ)^4 η(4τ)^4` (weight 4, level 8; LMFDB 8.4.a.a, in Zudilin's Table 1 under the old label 8.4.1.a).

**Theorem Z.** Let `f(τ) = η(2τ)^4 η(4τ)^4 = Σ_{n≥1} c_n q^n`, `A_k = C(2k,k)^4/256^k = ((1/2)_k/k!)^4`,

    S_0 = Σ_{k≥0} A_k = 4F3(1/2,1/2,1/2,1/2; 1,1,1; 1),        S_1 = Σ_{k≥0} A_k (H_{2k} − H_k − log 2).

Then

    L(f,1) = −S_1/2 = (1/(2π^2)) ∫_0^1 K'(k)^2 dk,        L(f,2) = (π^2/16) S_0 = (1/(2π)) ∫_0^1 K(k) K'(k) dk,

and, for `(r,t) = (1/2,1/2)`, all three relations in (6) hold, with

    L(f,1)/F_1(1) = −1/16,        L(f,2)/F_2(1) = 1/32,        L(f,3)/F_3(1) = −3/256.

(`F_0(1) = S_0`, `F_1(1) = 8 S_1`, `F_2(1) = 2π^2 S_0`, `F_3(1) = (32π^2/3) S_1`.)

**Status and novelty risk — read this before quoting.** *Proved* (below; elementary, self-contained modulo the
classical facts of §1 and Prop. 8 of `../U2/THEOREM.md`), audited (`AUDIT.md`: mathematics complete, numerics
re-done independently to 40 digits). **Novelty: essentially none in the ingredients** — every analytic step is
already published, and the auditor's verdict is that only the assembly is unrecorded:

* `L`-side: M. Rogers, J. G. Wan, I. J. Zucker, *Moments of elliptic integrals and critical L-values*, Ramanujan J.
  37 (2015), arXiv:1303.2259, eqs (30), (31), for `f_2 = η^4(2τ)η^4(4τ)`: `L(f_2,2) = (π^2/16)·4F3(1/2^4;1^3;1)` and
  `L(f_2,3) = (π^2/4) L(f_2,1) = (1/4)∫_0^1 K(k)^2 dk = (π^4/128)·7F6(5/4,1/2^6;1/4,1^5;1)`; and `2∫K^2 = ∫K'^2`
  (J. Wan, *Moments of products of elliptic integrals*, Adv. Appl. Math. 2012, eq. (3)).
* The "new link" `∫_0^1 K'(k)^2 dk = −π^2 S_1` is **already in** J. M. Borwein, A. Straub, J. Wan, *Three-step and
  four-step random walk integrals*, Exp. Math. 22 (2013) (preprint 2010): the display after their (29),
  `W_4(−1) = (4/π) Σ_n C(2n,n)^4/4^{4n} Σ_{k≥2n+1} (−1)^{k+1}/k` (the inner sum is `−h_n`, so this is `−(4/π)S_1`),
  and their (65), `W_4(−1) = (4/π^3) ∫_0^1 K'(k)^2 dk`. With Rogers–Wan–Zucker (31), `L(f,1) = −S_1/2` follows in
  two lines from published formulas.
* `F_2(1) = 2π^2 F_0(1)` is stated as rigorous in Zudilin's Remark 2, and `F_3(1) ∝ π^2 F_1(1)` in his Remark 5.
* K. Bönisch, A. Klemm, E. Scheidegger, D. Zagier, arXiv:2203.09426 (2022), §4, prove (via a correspondence with
  the Kuga–Sato threefold over `X_0(8)`) that the periods of the conifold fibre of this family are periods and
  quasi-periods of this `f`, "up to a multiplicative constant"; the exact entries of their transition matrix are
  numerical identifications.

What I did not find stated anywhere (searched 2026-11-26: papers citing the SIGMA paper, Rogers–Wan–Zucker
follow-ups, K. C. Au arXiv:2509.19960, D. Bao arXiv:2608.26514, BKSZ, Bönisch–Golyshev–Klemm arXiv:2510.03939) is the
conclusion itself: that Zudilin's three relations (6) hold for `(1/2,1/2)` with the constants `−1/16, 1/32, −3/256`.
So: the statement "none of the relations in (6) seem to be proved" was, for this case, already answerable from the
2010–2015 literature; this file records the assembly and a self-contained proof. It should be described as
"a known-ingredients verification of Zudilin's relations in one case", not as a new theorem.

**Consequence.** The L-value form of Theorem U2 (`../U2/THEOREM.md`, Corollary) is now unconditional:
`PCF[a^(U2), −64 n^10] = 8 L(f,1) / (186 ζ(5) L(f,1) − π^4 L(f,2))`.

**Checks:** `verify/z_checks.py` (all steps numerically to 25+ digits; the q-series identity exactly to `O(q^400)`;
the three ratios `−1/16, 1/32, −3/256` against Zudilin's definition (5) directly).

## 1. Classical facts used

Write `K(k) = ∫_0^{π/2} (1 − k^2 sin^2 θ)^{−1/2} dθ = (π/2) Σ_n a_n k^{2n}`, `a_n := ((1/2)_n/n!)^2`, `k' = √(1−k^2)`,
`K'(k) = K(k')`, and `h_n := H_{2n} − H_n − log 2`. Thetas are Jacobi Nullwerte in the nome `e^{πiσ}`.

(Z1) `θ_2(σ) = 2η(2σ)^2/η(σ)`, `θ_3(σ) = η(σ)^5/(η(σ/2)^2 η(2σ)^2)`, `θ_4(σ) = η(σ/2)^2/η(σ)` (Jacobi triple product).
(Z2) For `σ = iy`, `y > 0`: with `k := θ_2(σ)^2/θ_3(σ)^2 ∈ (0,1)` one has `k' = θ_4^2/θ_3^2`, `θ_3(σ)^2 = 2K(k)/π`,
     `y = K'(k)/K(k)`; `y ↦ k` is a decreasing bijection `(0,∞) → (0,1)` and `dy/dk = −π/(2 k k'^2 K(k)^2)`
     (Legendre's relation). (Whittaker–Watson ch. 21–22; Borwein–Borwein, *Pi and the AGM*, ch. 2.)
(Z3) `K'(k) = Σ_{n≥0} a_n k^{2n} ( log(1/k) + d(n) )`, `0 < k < 1`, `d(n) = ψ(n+1) − ψ(n+1/2)` (DLMF 19.12.1).
     Since `ψ(n+1) = −γ + H_n` and `ψ(n+1/2) = −γ − 2 log 2 + 2H_{2n} − H_n`, `d(n) = −2 h_n`.
     Moreover `h_n < 0` for all `n ≥ 0`: `H_{2n} − H_n = Σ_{j=1}^{2n} (−1)^{j+1}/j` is a partial sum, with an even number
     of terms, of the alternating series for `log 2`, hence `< log 2`.
(Z4) `η(−1/τ) = (τ/i)^{1/2} η(τ)`.
(Z5) Gauss: `2F1(a,b;c;1) = Γ(c)Γ(c−a−b)/(Γ(c−a)Γ(c−b))` for `c > a+b`.

## 2. Lemma 1 (the L-values as moments)

> `L(f,1) = (1/(2π^2)) ∫_0^1 K'(k)^2 dk`, `L(f,2) = (1/(2π)) ∫_0^1 K(k)K'(k) dk`, `L(f,3) = (π^2/4) L(f,1)`.

*Proof.* (i) By (Z1), `θ_2θ_3θ_4 = 2η(σ)^3`, so `θ_2^4θ_3^2θ_4^2 = θ_2^2 (θ_2θ_3θ_4)^2 = 16 η(2σ)^4 η(σ)^4`; with
`σ = 2τ`: `f(τ) = θ_2(2τ)^4 θ_3(2τ)^2 θ_4(2τ)^2/16`. By (Z2), for `τ = iy/2`:

    f(iy/2) = k^2 k' θ_3^8/16 = k^2 k' K(k)^4/π^4.

(ii) By (Z4), `f(−1/(8τ)) = η(−1/(4τ))^4 η(−1/(2τ))^4 = (4τ/i)^2 (2τ/i)^2 η(4τ)^4η(2τ)^4 = 64 τ^4 f(τ)`, i.e.
`f(i/(8t)) = 64 t^4 f(it)`. Hence `f(it) = O(e^{−2πt})` as `t → ∞` and `f(it) = O(t^{−4} e^{−π/(4t)})` as `t → 0+`, the
Mellin transform `Λ(s) := ∫_0^∞ f(it) t^{s−1} dt` is entire, `Λ(s) = (2π)^{−s} Γ(s) L(f,s)` (for `Re s` large by
termwise integration, `c_n = O(n^2)`; then for all `s` by analytic continuation), and the substitution `t ↦ 1/(8t)`
gives `Λ(4−s) = 8^{s−2} Λ(s)`. At `s = 3`: `Λ(1) = 8Λ(3)`, i.e. `L(f,1)/(2π) = 8·2 L(f,3)/(2π)^3`, so
`L(f,3) = (π^2/4) L(f,1)`.

(iii) `L(f,1) = 2π Λ(1) = π ∫_0^∞ f(iy/2) dy`. Substituting `y = K'/K` (Z2):

    L(f,1) = π ∫_0^1 (k^2 k' K^4/π^4) · π/(2 k k'^2 K^2) dk = (1/(2π^2)) ∫_0^1 k K(k)^2/k' dk = (1/(2π^2)) ∫_0^1 K'(κ)^2 dκ,

the last step by `κ = k'`, `dκ = −k dk/k'`, `K(k) = K'(κ)`. Likewise `L(f,2) = (2π)^2 Λ(2) = π^2 ∫_0^∞ f(iy/2) y dy
= π^2 ∫_0^1 (k^2k'K^4/π^4)(K'/K) π/(2kk'^2K^2) dk = (1/(2π)) ∫_0^1 k K(k)K'(k)/k' dk = (1/(2π)) ∫_0^1 K'(κ)K(κ) dκ`. ∎

## 3. Lemma 2 (moments of K')

> For `s > −1`: `M(s) := ∫_0^1 K'(k) k^s dk = (π/4) Γ((s+1)/2)^2 / Γ(s/2+1)^2`. In particular, for integers `n ≥ 0`,
> `M(2n) = (π^2/4) a_n` and `M'(2n) = −∫_0^1 K'(k) k^{2n} log(1/k) dk = (π^2/2) a_n h_n`.

*Proof.* `K'(k) = (π/2) Σ_m a_m (1−k^2)^m` (positive terms). With `x = k^2`, `c := (s+1)/2 > 0`, and Tonelli,

    M(s) = (π/4) Σ_m a_m ∫_0^1 x^{c−1}(1−x)^m dx = (π/4) Σ_m ((1/2)_m^2/m!^2) · m! Γ(c)/Γ(c+m+1)
         = (π/(4c)) 2F1(1/2,1/2; c+1; 1) = (π/(4c)) Γ(c+1)Γ(c)/Γ(c+1/2)^2 = (π/4) Γ(c)^2/Γ(c+1/2)^2

by (Z5). `M(2n) = (π/4)Γ(n+1/2)^2/n!^2 = (π^2/4) a_n`. `M` is differentiable under the integral sign on `s > −1`
(the integrand is dominated locally uniformly, `K'(k) = O(log(1/k))` at `0`), and
`M'(s)/M(s) = ψ((s+1)/2) − ψ(s/2+1)`, which at `s = 2n` equals `ψ(n+1/2) − ψ(n+1) = 2h_n` (Z3). ∎

## 4. Proof of Theorem Z

**`∫_0^1 K K' = (π^3/8) S_0`.** `K(k) = (π/2) Σ a_n k^{2n}` (positive terms), so by Tonelli and Lemma 2
`∫_0^1 K(k)K'(k) dk = (π/2) Σ_n a_n M(2n) = (π^3/8) Σ a_n^2 = (π^3/8) S_0` (`A_n = a_n^2`). With Lemma 1,
`L(f,2) = (π^2/16) S_0`. (This is Rogers–Wan–Zucker's (30).)

**`∫_0^1 K'^2 = −π^2 S_1`.** By (Z3), `K'(k) = Σ_n a_n k^{2n} (log(1/k) − 2h_n)` on `(0,1)`, and every term is
positive there (`h_n < 0`). By Tonelli and Lemma 2,

    ∫_0^1 K'(k)^2 dk = Σ_n a_n ( ∫_0^1 K'(k) k^{2n} log(1/k) dk − 2h_n M(2n) )
                    = Σ_n a_n ( −(π^2/2) a_n h_n − (π^2/2) a_n h_n ) = −π^2 Σ_n a_n^2 h_n = −π^2 S_1.

(The series converges absolutely: `a_n^2 ~ 1/(π^2 n^2)`, `h_n → 0`.) With Lemma 1,
`L(f,1) = (1/(2π^2)) (−π^2 S_1) = −S_1/2`.

**The values `F_j(1)`.** For `r = t = 1/2`, `Γ(r)Γ(1−r)Γ(t)Γ(1−t) = π^2`, and at `z = 1`

    F(ε) := Σ_j F_j(1) ε^j = (1/π^2) Σ_{n≥0} Γ(1/2+ε+n)^4/Γ(1+ε+n)^4 .

The terms are `O(n^{−2})` uniformly for `ε` in a compact subset of `Re ε > −1/2`, so `F` is analytic there and may be
differentiated termwise: `F_0(1) = Σ a_n^2 = S_0` and
`F_1(1) = Σ_n 4(ψ(n+1/2) − ψ(n+1)) a_n^2 = 8 Σ a_n^2 h_n = 8 S_1` (Z3). In the notation of `../U2/THEOREM.md` §7,
`π^2 F(ε) = Y(ε)/N(ε) = Σ τ_i ε^i`, and Proposition 8 there (the bilateral-series/Liouville argument; equivalently
Zudilin's Lemma 3 and Remark 5) gives `τ_2 = 2π^2 τ_0` and `τ_3 = (4π^2/3) τ_1`, i.e.

    F_2(1) = 2π^2 S_0,        F_3(1) = (4π^2/3) F_1(1) = (32π^2/3) S_1.

**The three ratios.** `L(f,1)/F_1(1) = (−S_1/2)/(8S_1) = −1/16`; `L(f,2)/F_2(1) = (π^2S_0/16)/(2π^2S_0) = 1/32`;
`L(f,3)/F_3(1) = (π^2/4)(−S_1/2)/((32π^2/3)S_1) = −3/256` (Lemma 1 for `L(f,3)`). `S_0 > 0` and `S_1 < 0` (all
`h_n < 0`), so the denominators do not vanish. ∎

## 5. Remarks

1. Numerically `S_1 = −0.709001367461929437531...`, `L(f,1) = 0.354500683730964718765...`,
   `S_0 = 1.118636387164187068349...`, `L(f,2) = 0.690031163123397525119...`.
2. The same two lemmas give `∫_0^1 K'(k)^2 k^{2m} dk` and `∫_0^1 K(k)K'(k) k^{2m} dk` as explicit series; nothing here
   is specific to `m = 0` except Lemma 1.
3. The other 13 cases `(r,t)` of Zudilin's list are not touched: there the conifold period is not a product of two
   elliptic periods, and (6) remains open as far as I know (BKSZ treat them numerically).
