# The "still unproven" rows of Table 3 of *The Ramanujan Library* (arXiv:2412.12361), and two misprints there

## 0. Source, statements, status

**Source (read 2026-11-26).** I. Beit-Halachmi, I. Kaminer, *The Ramanujan Library — Automated Discovery on the
Hypergraph of Integer Relations*, arXiv:2412.12361v2 (ICLR 2025), §5, Table 3. Notation of the paper:
`C[f_n] := 1 + f_1/(1 + f_2/(1 + f_3/(1 + ...)))`, and `R_2 := √(πe/2)·erf(1/√2)` (`= 1.4106861346424479976908...`).
Caption: "*The rows in this table represent the 8 formulas we discovered for √(πe), generalizing Ramanujan's
original formula that was so far unique. We proved the first 4 rows using transformations on Ramanujan's
original formula. The latter 4 rows are still unproven (yet equivalent, such that proving any one of them will
prove all four). The last row of each set of four is an infinite family of formulas: for any integer k ≥ 0
there exist integers α, β, γ, δ that satisfy an equality between the two sides.*" The rows:

    C[(1−2n)/(4n(n+1))]         = (R_2+1)/(2R_2)             C[n/(2n(n+1))]      = 1/(2R_2−1)
    C[(1−2n)/(4(n+1)(n+2))]     = (2R_2+1)/4                 C[(n+1)/(2n(n+1))]  = R_2
    C[(1−2n)/(4(n+2)(n+3))]     = (R_2+1)/(6R_2−1)           C[(n+2)/(2n(n+1))]  = (2R_2+1)/(R_2+1)
    C[(1−2n)/(4(n+k)(n+k+1))]   = (αR_2+β)/(γR_2+δ)          C[(n+k)/(2n(n+1))]  = (αR_2+β)/(γR_2+δ)

> **Theorem RL (the second set of four rows).** Let `k > −1` be real and `f_n = (n+k)/(2n(n+1))`. Then `C[f_n]`
> converges and
>
>     C[(n+k)/(2n(n+1))] = 1 + ((k+1)/5) · M(k+2; 7/2; 1/2) / M(k+1; 5/2; 1/2),
>
> where `M(a;b;z) = 1F1(a;b;z)` is Kummer's function. For every integer `k ≥ 0` the right side equals
> `(α_k R_2 + β_k)/(γ_k R_2 + δ_k)` with integers given by the recurrences of §3; in particular
>
>     k = 0:  1/(2R_2 − 2)        k = 1:  R_2        k = 2:  (2R_2 + 1)/(R_2 + 1)
>     k = 3:  (15R_2 + 9)/(8R_2 + 6)        k = 4:  (76R_2 + 48)/(39R_2 + 27)        k = 5:  (955R_2 + 615)/(464R_2 + 312).

> **Two printed values are false (misprints).** `C[n/(2n(n+1))] = 1.21747475218590...` is `1/(2R_2 − 2)`, not the
> printed `1/(2R_2 − 1) = 0.549...` (this first misprint is visible without computation: a C-transform of positive
> terms is `> 1`); proved, case `k = 0` of the Theorem.
> `C[(1−2n)/(4(n+2)(n+3))] = 0.97831650145726...`, whereas the printed `(R_2+1)/(6R_2−1) = 0.3229...`; numerically
> (40 digits, PSLQ) the value is `(R_2+1)/(6R_2 − 6)` — **numerical only**, this row belongs to the first set, which
> the authors state they have proved, and is not proved here.

**Status and novelty — read this first.** Theorem RL is proved below (elementary: an explicit minimal solution of
the three-term recurrence) and was audited (`AUDIT.md`: proof correct; all values re-derived and checked to 60
digits by independent code). **Novelty: nil.** The auditor identified the formula as a direct specialization of
Perron's classical continued fraction for ratios of contiguous Kummer functions,

    M(a+1; b+1; z)/M(a; b; z) = b/( b − z + (a+1)z/( b+1 − z + (a+2)z/( b+2 − z + ... ) ) )

(O. Perron, *Die Lehre von den Kettenbrüchen*, vol. II, 1957; analysed in W. Gautschi, Math. Comp. 31 (1977)
994–999, who also shows that this fraction can converge to a wrong limit for some parameters — not here, by the
minimal-solution argument of §2): with `a = k+1`, `b = 5/2`, `z = 1/2` it reads
`1 + ((k+1)/2)/(2 + ((k+2)/2)/(3 + ...))`, and rescaling the `n`-th level by `n+1` gives exactly `C[(n+k)/(2n(n+1))]`
(the convergents agree as exact fractions; auditor's script). So the contribution of this file is only:
(i) the identification of the four "still unproven" rows with Perron's fraction, with a self-contained proof and
the explicit integer Möbius coefficients; (ii) the two misprints (confirmed by the auditor in arXiv v2 and in the
ICLR 2025 proceedings PDF; there is no corrected v3); (iii) §5 and §6 below (Conjecture 2 is Gauss's logarithm
fraction; Conjecture 1 fails as printed). No published proof of the four rows as such was found (0–2 citing
works; none treats them). Checks: `verify/rl_checks.py`, `verify/rl_conj1.py`, `verify/rl_conj2.py`; auditor's
scripts in `audit/`. (§5 and §6 were added after the audit and were not seen by the auditor.)

## 1. The minimal solution

Fix `k > −1`. For integers `n ≥ −1` put

    J_n := Σ_{j≥0} Γ(n+k+2+j) / (2^j j! Γ(n+j+7/2))  = (Γ(n+k+2)/Γ(n+7/2)) · M(n+k+2; n+7/2; 1/2),
    m_n := (−1/2)^n J_n / (n+1)! .

All terms of `J_n` are positive (`n+k+2+j ≥ k+1 > 0`), and the series converges (ratio of terms `~ 1/(2j)`).

**Lemma 1.** For `n ≥ 1`: `J_n + 2(n+1) J_{n−1} = 2(n+k) J_{n−2}`, and therefore `m_n = m_{n−1} + f_n m_{n−2}`.

*Proof.* Let `A_j := Γ(n+k+j)/Γ(n+j+3/2)` (`j ≥ 0`; `n+k+j > 0`). Then `(n+k+j) A_j = (n+j+3/2) A_{j+1}`, and
`J_{n−2} = Σ A_j/(2^j j!)`, `J_{n−1} = Σ A_{j+1}/(2^j j!)`, `J_n = Σ A_{j+2}/(2^j j!)`. In
`S := Σ_j (A_{j+2} + 2(n+1)A_{j+1} − 2(n+k)A_j)/(2^j j!)` (absolutely convergent) collect the coefficient of `A_J`:

    ( 4J(J−1) + 4(n+1)J − 2(n+k) ) / (2^J J!) = ( 4J(n+J+1/2) − 2(n+k+J) ) / (2^J J!) .

Hence `S = Σ_{J≥1} 2(n+J+1/2) A_J/(2^{J−1}(J−1)!) − Σ_{J≥0} 2(n+k+J) A_J/(2^J J!)`, and by
`(n+k+J)A_J = (n+J+3/2)A_{J+1}` the second sum is the first one with `J − 1` in place of `J`. So `S = 0`.
For the recurrence: `m_n − m_{n−1} = ((−1/2)^n/(n+1)!)(J_n + 2(n+1)J_{n−1}) = ((−1/2)^n/(n+1)!) · 2(n+k) J_{n−2}`, and
`f_n m_{n−2} = ((n+k)/(2n(n+1))) (−1/2)^{n−2} J_{n−2}/(n−1)!` is the same number. ∎

(Heuristically `J_n = ∫_0^1 e^{τ/2} τ^{n+k+1}(1−τ)^{1/2−k} dτ / Γ(3/2−k)`, in the sense of analytic continuation in
the exponent, and the three-term relation is `∫ d[e^{τ/2}τ^{n+k}(1−τ)^{3/2−k}] = 0`; the proof above avoids this.)

**Lemma 2.** `|m_n| ≤ c_k^n/n!` for some `c_k`; in particular `m_n → 0`.

*Proof.* `(n+k+2+i)/(n+7/2+i) ≤ ρ := max(1, 2(k+1)/5)` for `n ≥ −1`, `i ≥ 0`, so `M(n+k+2; n+7/2; 1/2) ≤ e^{ρ/2}`, and
`Γ(n+k+2)/Γ(n+7/2) = O(n^{k−3/2})`. ∎

## 2. Proof of the formula

Let `P_n, Q_n` be the numerators and denominators of `C[f_n] = 1 + K(f_n/1)`: both satisfy `u_n = u_{n−1} + f_n u_{n−2}`
(`n ≥ 1`), `(P_{−1}, P_0) = (1, 1)`, `(Q_{−1}, Q_0) = (0, 1)`. Since `f_n > 0`, `Q_n` is nondecreasing and `Q_n ≥ 1` for
`n ≥ 0`. By Lemma 1, `m` solves the same recurrence, so (comparing at `n = −1, 0`)

    m_n = m_{−1} P_n + (m_0 − m_{−1}) Q_n ,        m_{−1} = −2 J_{−1} ≠ 0,     m_0 = J_0 .

Dividing by `m_{−1} Q_n` and using Lemma 2: `P_n/Q_n → 1 − m_0/m_{−1} = 1 + J_0/(2 J_{−1})`. Finally
`J_0/J_{−1} = (Γ(k+2)Γ(5/2)/(Γ(k+1)Γ(7/2))) · M(k+2;7/2;1/2)/M(k+1;5/2;1/2) = (2(k+1)/5) · M(k+2;7/2;1/2)/M(k+1;5/2;1/2)`. ∎

## 3. Integer k: reduction to R_2

`R_2 = M(1; 3/2; 1/2)`: indeed `M(1;3/2;z) = Σ_j z^j/(3/2)_j`, so `M(1;3/2;1/2) = Σ_{j≥0} 1/(3·5···(2j+1))`, and
`Σ_{j≥0} x^{2j+1}/(2j+1)!! = e^{x^2/2} ∫_0^x e^{−t^2/2} dt` (both sides solve `y' = 1 + xy`, `y(0) = 0`), which at
`x = 1` is `√(πe/2)·erf(1/√2)`. Put `μ_a := M(a; 3/2; 1/2)`; then `μ_0 = 1`, `μ_1 = R_2`, and by the contiguous relations
DLMF 13.3.1 (`(b−a)M(a−1,b,z) + (2a−b+z)M(a,b,z) − aM(a+1,b,z) = 0`) and 13.3.4
(`bM(a,b,z) − bM(a−1,b,z) − zM(a,b+1,z) = 0`), at `z = 1/2`:

    a μ_{a+1} = (3/2 − a) μ_{a−1} + (2a − 1) μ_a ,
    M(a; 5/2; 1/2) = 3 (μ_a − μ_{a−1}),        M(a; 7/2; 1/2) = 5 ( M(a;5/2;1/2) − M(a−1;5/2;1/2) ).

So for integers `a ≥ 0` all three are of the form (rational) + (rational)·`R_2`, and the Theorem gives
`C[(n+k)/(2n(n+1))] ∈ Q(R_2)`, a Möbius transform of `R_2` with rational, hence (after scaling) integer,
coefficients. E.g. `μ_2 = R_2 + 1/2`, `μ_3 = (5R_2+3)/4`, `μ_4 = (19R_2+12)/12`;
`M(1;5/2) = 3(R_2−1)`, `M(2;5/2) = 3/2`, `M(3;5/2) = (3/4)(R_2+1)`; `M(2;7/2) = (15/2)(3−2R_2)`, `M(3;7/2) = (15/4)(R_2−1)`,
`M(4;7/2) = 5R_2/4` (arguments `z = 1/2` omitted), whence

    k = 0:  1 + (1/5)(15/2)(3−2R_2)/(3(R_2−1)) = 1/(2R_2−2),
    k = 1:  1 + (2/5)(15/4)(R_2−1)/(3/2) = R_2,
    k = 2:  1 + (3/5)(5R_2/4)/((3/4)(R_2+1)) = (2R_2+1)/(R_2+1).

The entries for `k = 3, 4, 5` in the Theorem, and the table for `k ≤ 8`, are obtained by propagating these three
relations in exact rational arithmetic (`verify/rl_checks.py` (5); e.g. `μ_5 = (191R_2 + 123)/96`). (The denominators
`γ_kR_2 + δ_k` are nonzero: they are positive multiples of `M(k+1;5/2;1/2) > 0`.)

## 4. Remarks

1. `k = 0` and `k = 1` are tails of one another: `C[1/(2n)] = 1 + (1/2)/C[1/(2(n+1))]`, i.e. `R_2 = 1 + 1/(2C_0)`, which
   is how the misprint in the row `C[n/(2n(n+1))]` is seen at once.
2. For half-integers `k = 3/2, 5/2, ...` the value is rational, since `M(b+m; b; z)` is `e^z` times a polynomial:
   e.g. `C[(n+3/2)/(2n(n+1))] = 3/2`, `C[(n+5/2)/(2n(n+1))] = 5/3` (checked to 35 digits).
3. The first set of rows (`f_n = (1−2n)/(4(n+k)(n+k+1))`) is not treated; its third row is misprinted as stated
   in §0 (numerical).

## 5. Conjecture 2 of the same paper is the classical logarithm fraction

The paper states (§5): "**Conjecture 2** For all k ≥ 1, `C[n^2/(k^2(1−4n^2))] = (2/k)/(ln(k+1) − ln(k−1))`."

**Proposition.** This holds for every real `k > 1`, and for `k = 1` in the sense that the fraction converges to 0.

*Proof.* The classical continued fraction (Gauss; DLMF 4.9.2)
`ln((1+z)/(1−z)) = 2z/(1 − z^2/(3 − 4z^2/(5 − 9z^2/(7 − ...))))`, valid for `z ∈ C∖((−∞,−1] ∪ [1,∞))`, has the denominator
`D(z) = 1 + K_{n≥1}( −n^2z^2 / (2n+1) )`. The equivalence transformation with multipliers `r_0 = 1`, `r_n = 1/(2n+1)`
turns the partial numerators into `r_n r_{n−1}(−n^2z^2) = n^2z^2/(1−4n^2)` and the partial denominators into 1, without
changing the approximants; so `D(z) = C[n^2z^2/(1−4n^2)]`, and at `z = 1/k`:
`C[n^2/(k^2(1−4n^2))] = (2/k)/ln((k+1)/(k−1))`. For `k = 1` the `N`-th approximant equals `1/H_{N+1}` (`H` the harmonic
numbers; checked exactly for `N < 60`, not proved here in general), consistent with the value `0`. ∎

(`verify/rl_conj2.py`: the approximants of both fractions coincide exactly for `k ≤ 6`, `N ≤ 30`; 35 digits for
`k = 2..8`.) So Conjecture 2 is not an open problem; this was already noted in `../STATUS-MAP-RM.md` §1.8.

## 6. Conjecture 1 of the same paper (convergence rules) is false as printed

The paper states (§3.1; read from the PDF page, 2026-11-26): "**Conjecture 1** Let f_n be a real sequence. The
formula C[f_n] converges in the following cases, at a rate provided by the error ε_n := |1 + f_1/(1 + f_2/(··· + f_n))) −
C[f_n]|. • If f_n = O(n^k) for some k < 0, then ε_n = O((n!)^k) • If 1 + 4f_n = C + o(1) for some constant C ≠ 1, then
ε_n = O(|(1+√C)/(1−√C)|^{−n}) • If 1 + 4f_n = Cn^k + o(n^k) for some k ∈ (0, 1], then ε_n = O(e^{−4√(n/C)}) • If
1 + 4f_n = O(n^k) with k < 0 or k ∈ (1, 2], then ε_n = O(n^m) for some m < 0. Otherwise, C[f_n] does not converge."

**Proposition (counterexamples; all elementary).**
(A) *Second case, `C < 0`.* For every constant `f < −1/4`, `C[f]` diverges, although `1 + 4f = C` is a constant `≠ 1`.
(B) *Fourth case, `k = −2`.* For every `c > 1/16`, `C[−1/4 − c/(n(n+1))]` diverges, although `1 + 4f_n = −4c/(n(n+1)) = O(n^{−2})`.
(C) *First case, literally.* For `f_n = c/n` with `c > 1`, `ε_n` is not `O(1/n!)`: `ε_n + ε_{n+1} ≥ κ c^n n^{−2c}/n!` for some `κ > 0`.

*Proof.* (A) The `N`-th approximant is `T^N(1)` with `T(x) = 1 + f/x`. `T` is a Möbius map with matrix `[[1, f],[1, 0]]`,
`trace^2/det = −1/f ∈ (0, 4)`: elliptic, with non-real fixed points `(1 ± √(1+4f))/2`, and not the identity. `T` is
conjugate to a rotation of the Riemann sphere whose fixed points are not on the real circle `R ∪ {∞}`, so the orbit
of the real point 1 is either periodic with period `> 1` or dense in a circle; in neither case does it converge.
(E.g. `f = −1`: `0, ∞, 1, 0, ∞, 1, ...`; `f = −1/2`: `1/2, 0, ∞, 1, ...`; `f = −1/3`: period 6. `verify/rl_conj1.py`.)

(B) Write `c = 1/16 + β^2/4` with `β > 0`, `σ := 1/2 + iβ`, and `z_n := Γ(n+2+σ)/Γ(n+2)` (`n ≥ −1`). With `N := n+1`,
`z_{n−1} = z_n N/(N+σ)` and `z_{n−2} = z_n N(N−1)/((N+σ)(N−1+σ))`, so
`(z_n − 2z_{n−1} + z_{n−2})/z_{n−2} = ((N+σ)(N−1+σ) − 2N(N−1+σ) + N(N−1))/(N(N−1)) = σ(σ−1)/(n(n+1)) = −4c/(n(n+1))`.
Hence `u_n := 2^{−n} z_n` satisfies `u_n = u_{n−1} + f_n u_{n−2}` for `n ≥ 1`, `f_n = −1/4 − c/(n(n+1))`, and so does `ū_n`.
`u_0/u_{−1} = (1+σ)/2 ∉ R`, so `u, ū` are independent and the real solutions `P_n, Q_n` (numerators and denominators
of the approximants) are `P_n = Re(a u_n)`, `Q_n = Re(b u_n)` with complex `a, b ≠ 0`, `a/b ∉ R` (`P, Q` are
independent). With `θ_n := arg z_n`: `P_n/Q_n = h(θ_n)`, `h(θ) := Re(a e^{iθ})/Re(b e^{iθ})`, a non-constant continuous
function on the circle minus two points. Now `θ_n − θ_{n−1} = arg(1 + σ/(n+1)) = arctan(β/(n+3/2))` is positive, tends
to 0, and has divergent sum. Therefore `θ_n mod 2π` comes arbitrarily close to every point of the circle
infinitely often, and `h(θ_n)` has at least two limit points: the approximants do not converge.
(Numerically, `c = 1/4`: approximants `0.029, −1.447, 0.711, 0.196, −0.749, 0.970, ...` at `n = 10, 30, 100, ...`;
controls `c = 1/16`, `1/32` converge. `verify/rl_conj1.py`.) This is the well-known sharpness of the parabola
theorem at `a_n → −1/4` (cf. Jacobsen–Masson; Lorentzen–Waadeland, *Continued Fractions*, ch. on limit-periodic
fractions); the explicit solution makes it self-contained.

(C) `P_n/Q_n − P_{n−1}/Q_{n−1} = (−1)^{n−1} ∏_{i≤n} f_i/(Q_n Q_{n−1}) = (−1)^{n−1} c^n/(n! Q_n Q_{n−1})`, where
`Q_n − Q_{n−1} = (c/n) Q_{n−2} > 0`; so `Q_n` is increasing and `Q_n ≤ Q_{n−1}(1 + c/n)`, whence `Q_n = O(n^c)`.
By the triangle inequality `ε_n + ε_{n+1} ≥ |P_{n+1}/Q_{n+1} − P_n/Q_n| ≥ κ c^{n+1}/((n+1)!(n+1)^{2c})`, which is not
`O(1/n!)` when `c > 1`. (The intended meaning is presumably "up to exponential factors".) ∎

**What this does and does not say.** The conjecture is a set of heuristic rules used to filter a database; the
counterexamples show that, as a mathematical statement, cases 2 and 4 need the extra hypotheses `C > 0`
(respectively: not `f_n → −1/4` from below faster than allowed by the parabola theorem), and case 1 needs
"up to `c^n`". Nothing is claimed about the rules in the range where they are meant to apply (e.g. `C > 0`:
classical, Poincaré–Perron/limit-periodic theory). Novelty: none — these are textbook facts about continued
fractions applied to a printed conjecture.
