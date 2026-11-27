# R08. The closing conjectures of arXiv:2607.16608 (2-colour overpartitions `p̄_ℓ(n)`): all PROVED

**Source.** H. S. Sumanth Bharadwaj, N. Sujatha, S. Chandankumar, *Arithmetic properties of the
2-color overpartition function `p̄_ℓ(n)`*, arXiv:2607.16608 (18 July 2026), Section "Closing Remarks".
`p̄_ℓ(n)` counts 2-colour overpartitions in which one colour appears only in multiples of `ℓ`:

    Σ p̄_ℓ(n) qⁿ = f₂ f_{2ℓ} / (f₁² f_ℓ²) = 1/(φ(−q) φ(−q^ℓ)),   φ(−q) = 1 + 2S(q),  S(q) := Σ_{n≥1} (−1)ⁿ q^{n²}.

**Printed.**

> Conjecture. For all `ℓ ≡ 0, 1 (mod 4)` and `n ≥ 0`, `p̄_ℓ(4n+3) ≡ 0 (mod 8)`. The case `ℓ ≡ 0 (mod 8)`
> follows from Corollary [c1]. The remaining cases are supported by numerical evidence and left open.
>
> Conjecture. For all `n ≥ 0`: `p̄₈(32n+28) ≡ 0 (mod 32)`, `p̄₆(24n+17) ≡ 0 (mod 64)`,
> `p̄₆(24n+23) ≡ 0 (mod 64)`.
>
> Conjecture. For all `n ≥ 0` and `α ≥ 1`: `p̄₉(2^α(8n+7)) ≡ 0 (mod 16)`, `p̄₉(36n+33) ≡ 0 (mod 128)`.

(and, before these, "`p̄₄(36n+21) ≡ 0 (mod 9)`, whose proofs we leave to the interested reader".)

## Theorem 1 (first conjecture)
For every `ℓ ≡ 0 or 1 (mod 4)` and every `n ≥ 0`: `p̄_ℓ(4n+3) ≡ 0 (mod 8)`.

*Proof.* Put `S = S(q)`, `T = S(q^ℓ)`. Modulo 8,
`1/((1+2S)(1+2T)) ≡ (1 − 2S + 4S²)(1 − 2T + 4T²) ≡ 1 − 2S − 2T + 4S² + 4ST + 4T²`.
The exponents occurring in `S, T, S², ST, T²` are `a²`, `ℓa²`, `a²+b²`, `a²+ℓb²`, `ℓ(a²+b²)`
(`a, b ≥ 1`). For `ℓ ≡ 0 (mod 4)` these are `≡ a²` or `0 (mod 4)`, or `a²+b²`, hence in `{0,1,2}` mod 4.
For `ℓ ≡ 1 (mod 4)` they are `≡ a²`, `a²+b²` or `ℓ(a²+b²) ≡ a²+b² (mod 4)`, again in `{0,1,2}`.
So no exponent is `≡ 3 (mod 4)`. ∎

## Theorem 2 (third conjecture, first part; also true for `α = 0`)
For all `α ≥ 0`, `n ≥ 0`: `p̄₉(2^α(8n+7)) ≡ 0 (mod 16)`.

*Proof.* With `T = S(q⁹)`, modulo 16,
`1/((1+2S)(1+2T)) ≡ 1 − 2(S+T) + 4(S²+ST+T²) − 8(S³+S²T+ST²+T³)`.
Let `m = 2^α u`, `u ≡ 7 (mod 8)`.
* `S`, `T`: `m` or `m/9` would be a square; impossible, since the odd part of a square is `≡ 1 (mod 8)`
  (and `9^{−1}·7 ≡ 7 (mod 8)`).
* `S²`, `ST`, `T²`: `m` would be a sum of two positive squares (`a²+b²`, `a²+(3b)²`, `(3a)²+(3b)²`);
  impossible, because `u ≡ 3 (mod 4)` contains a prime `≡ 3 (mod 4)` to an odd power.
* Cubic terms: it suffices that the total number `N` of solutions in positive integers of
  `a²+b²+c² = m`, `a²+b²+9c² = m`, `a²+9b²+9c² = m`, `9(a²+b²+c²) = m` is even (signs do not matter
  mod 2). If `α` is even, `m = 4^{α/2}u` is not a sum of three squares at all, so `N = 0`. If `α` is
  odd: pairing solutions under the permutations of the variables with equal coefficients
  (orbits of size 6 or 2 are even; orbits with exactly two equal variables have size 3 resp. 1,
  and those with three equal variables size 1), `N ≡ #{2a²+c² = m} + #{2a²+9c² = m} + #{a²+18b² = m}
  + #{18a²+9c² = m} (mod 2)`. Each of these counts solutions of `x² + 2y² = m`. But an integer
  represented by `x²+2y²` has every prime `≡ 5, 7 (mod 8)` to an even power, which forces its odd
  part to be `≡ 1 or 3 (mod 8)`; here it is `≡ 7`. So all four sets are empty and `N` is even. ∎

## Theorem 3 (second conjecture, and the second part of the third)
For all `n ≥ 0`: `p̄₈(32n+28) ≡ 0 (mod 32)`, `p̄₆(24n+17) ≡ p̄₆(24n+23) ≡ 0 (mod 64)`,
`p̄₉(36n+33) ≡ 0 (mod 128)`; also `p̄₄(36n+21) ≡ 0 (mod 9)`.

*Proof.* Radu's finite check (see `R03-radu-congruences.md` for the lemma) applied to the
eta-quotients `f₂f_{2ℓ}/(f₁²f_ℓ²)`; in each case `P_{m,r}(t) = {t}` and

| statement | `⌊ν⌋` | verified for `0 ≤ n ≤ ⌊ν⌋` |
|---|---|---|
| `p̄₈(32n+28)` mod 32 | 16 | all zero |
| `p̄₆(24n+17)` mod 64 | 41 | all zero |
| `p̄₆(24n+23)` mod 64 | 41 | all zero |
| `p̄₉(36n+33)` mod 128 | 28 | all zero |
| `p̄₄(36n+21)` mod 9 | 31 | all zero |

Certificate: `cert/radu_2color_overpartitions.txt` (program `code/radu_2color.py`). ∎

## Numerical sanity checks (not part of the proofs)
`code/radu_2color.py`: Theorem 1 for all `ℓ < 80`, `ℓ ≡ 0, 1 (mod 4)`, arguments `< 60000`; Theorem 2
for `α < 12`, arguments `< 60000`; controls `ℓ = 2, 3, 6, 7` fail as expected.

## Novelty (26 Nov 2026, `NOVELTY-2.md`, item 4)
No proof found: the source has only v1, no citing paper, and no public issue or record mentions these
conjectures (GitHub issue search for the arXiv id: 0 hits).
