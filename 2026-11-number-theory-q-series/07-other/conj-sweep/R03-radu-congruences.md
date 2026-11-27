# R03. Eta-quotient congruence conjectures settled by Radu's finite check

All statements below are of the form `c_r(mn+t) ≡ 0 (mod u)` for the coefficients of an
eta-quotient `∏_{δ|M} f_δ^{r_δ}`, `f_k := (q^k;q^k)_∞`. They are proved by Radu's Lemma 4.5
(S. Radu, *An algorithmic approach to Ramanujan's congruences*, Ramanujan J. 20 (2009), 215–251;
in the form of Radu–Sellers, Int. J. Number Theory 7 (2011), 2249–2259, Lemma 2.4): if
`(m, M, N, t, r) ∈ Δ*`, if an integer vector `a = (a_δ)_{δ|N}` satisfies
`p_{m,r}(γ) + p*_a(γ) ≥ 0` for a complete set of representatives `γ` of
`Γ₀(N)\SL₂(ℤ)/Γ_∞`, and if the congruence holds for all `t′ ∈ P_{m,r}(t)` and `0 ≤ n ≤ ⌊ν⌋`, with

    ν = (1/24)[(Σ r_δ + Σ a_δ)·[SL₂(ℤ):Γ₀(N)] − Σ δ a_δ] − (1/(24m)) Σ δ r_δ − t_min/m,

then it holds for all `n ≥ 0`. The implementation is `code/radu.py` (definitions of `Δ*`,
`P_{m,r}(t)`, `p_{m,r}`, `p*_a` as in the lemma; all hypotheses are verified in exact rational
arithmetic; the vector `a` is found by an integer program and then re-verified exactly; cusp
representatives: all `a/c` with `c | N`, `gcd(a,c) = 1`, `1 ≤ a ≤ c`, a superset of a complete
set). Engine sanity tests on known true and known false congruences: `cert/radu_sanity.txt`.

## R03a. arXiv:2611.05623 (D. Yaqubi, M. Mirzavaziri, 3 Nov 2026), Conjectures `conj:A6`, `conj:deep`

`A_k(n)` := number of ordered `k`-tuples of partitions with designated summands of total size `n`;
`Σ A_k(n)qⁿ = (f₆/(f₁f₂f₃))^k`.

**Printed.** "Conjecture. For all `n ≥ 0`, `A₆(12n+10) ≡ 0 (mod 3⁵)`, and `3⁵` is the exact power of
3 dividing all of them." and "Conjecture. (a) `A₃(12n+7) ≡ A₃(12n+11) ≡ 0 (mod 3⁴)`.
(b) `A₄(6n+5) ≡ 0 (mod 3³)` and `A₇(12n+11) ≡ 0 (mod 3³)`. (c) `A₈(18n+15) ≡ 0 (mod 3⁵)`."

**Result: all PROVED.**

| statement | m, t | N | ⌊ν⌋ | check |
|---|---|---|---|---|
| `A₆(12n+10) ≡ 0 (mod 243)` | 12, 10 | 6 | 26 | all zero |
| `A₃(12n+7) ≡ 0 (mod 81)` | 12, 7 | 12 | 27 | all zero |
| `A₃(12n+11) ≡ 0 (mod 81)` | 12, 11 | 12 | 26 | all zero |
| `A₄(6n+5) ≡ 0 (mod 27)` | 6, 5 | 6 | 7 | all zero |
| `A₇(12n+11) ≡ 0 (mod 27)` | 12, 11 | 12 | 63 | all zero |
| `A₈(18n+15) ≡ 0 (mod 243)` | 18, 15 | 6 | 63 | all zero |

In every case `P_{m,r}(t) = {t}`. Exactness of `3⁵` in the first statement: `A₆(22) ≡ 486 (mod 729)`,
so `v₃(A₆(22)) = 5` (`cert/radu_sanity.txt`, the test "A_6(12n+10) mod 3^6").
Certificate: `cert/radu_designated_summands.txt` (program `code/radu_ds.py`).

**Conjectures 1 and 2 of the same paper** (according to `NOVELTY-1.md`, item 2, these two — numbered
6.1 and 6.2 in the arXiv listing — have since been proved by others: a Lean proof of 6.1 in the public
GitHub project `the-omega-institute/trureturing`, PR #33580, 23 Nov 2026, and 6.2 in Zenodo record
10.5281/zenodo.24552279, 26 Nov 2026; neither was read here) (2-adic: "For every `k ≥ 1`, the largest power of 2 dividing
all `A_k(4n+3)` is `2^{v₂(A_k(3))}` if `k ≢ 2 (mod 16)`, and exactly `2³` if `k ≡ 2 (mod 16)`", and the
analogue for `B_k`, exceptions `k ≡ 1, 6 (mod 8)`): NOT proved. Survived testing for all `k ≤ 2100`
and `n < 160`, and for `k ≤ 300`, `n < 300` (the authors report `k ≤ 96`). `code/ds_2adic.py`,
`cert/ds_2adic_test.txt`.

## R03b. arXiv:2503.06228 (H. Nath, A. Sarma, 8 March 2025), Conjecture 1

`ped(n)` := number of partitions of `n` with even parts distinct, `Σ ped(n)qⁿ = f₄/f₁`.

**Printed.** "Conjecture. For every `n ≥ 0`, we have `ped(225n+43) ≡ 0`, `ped(225n+88) ≡ 0`,
`ped(225n+133) ≡ 0`, `ped(225n+223) ≡ 0 (mod 192)`."

**Result: re-proved; NOT a new result.** The novelty search (`NOVELTY-1.md`, item 6) showed that this
"Conjecture" in the introduction of the source restates a conjecture of Nath (2024) and that the
source paper itself proves it (its Theorem 1.2, by Radu's algorithm). Listing it as an open candidate
was a misreading of the extracted `Conjecture` environment. The computation below only confirms the
source's theorem with an independent implementation. `m = 225`, `N = 60`; `P_{m,r}(43) = {43, 88}`, `P_{m,r}(133) = {133, 223}`;
`⌊ν⌋ = 428` in both cases; all `2·429·2` coefficients vanish mod 192.
Certificate: `cert/radu_ped_mod192.txt` (program `code/radu_ped.py`).
(Conjecture 2 of that paper, an infinite family in a prime `p`, was not attempted.)

## R03c. arXiv:2603.08510, Conjecture 1 (overpartitions mod 7)
See `R02-overpartition-congruences.md`, Theorem C.

## Status, audit, novelty
R03a (the six congruences mod powers of 3, numbered 6.3 and 6.4 in the source): no proof found as of
26 Nov 2026 (`NOVELTY-1.md`, item 2; public issues with numerical checks only); the paper is under
active automated attack by others, so the novelty risk is high. R03b: not new (see above).
See also `RESULTS.md`.
