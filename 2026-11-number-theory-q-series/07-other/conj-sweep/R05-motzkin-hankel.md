# R05. Hankel determinants of `M(x)^r` (Motzkin): the conjectures are false as printed (sign errors); repaired forms survive

**Source.** Ying Wang, Yingrui Zhang, *Hankel determinants for convolution powers of Motzkin numbers*,
arXiv:2502.21050v1 (28 Feb 2025; only version). `M(x)` = Motzkin generating function,
`F(x,r) = M(x)^r`, `H_n(A) = det(a_{i+j})_{0≤i,j≤n−1}`, `H_0 = 1`.

**Printed.**

> Conjecture. For `r ≡ 0 (mod 3)`, we have `H_{rn}(F(x,r)) = H_{rn+1}(F(x,r)) = (−1)^n (n+1)^{r−1}`,
> `H_{rn+2}(F(x,r)) + H_{rn−1}(F(x,r)) = 3²α(n+1)^{r−1}`, where `|α| = (1/18) r(r−3)`.
>
> Conjecture. For `r ≡ 1 or 2 (mod 3)`, we have
> `H_{3rn} = H_{3rn+1} = H_{3rn+r} = H_{3rn+r+1} = α`; `H_{3rn+2r} = H_{3rn+2r+1} = 0`;
> `H_{3rn+2r−1} + H_{3rn+2r+2} = β((2r)(n+1))^{r−2}`; `H_{3rn+2} + H_{3rn−1} = ½γ r(r−3)`,
> where `|α| = |β| = |γ| = 1`.

followed by "Theorem. The above Conjectures hold for `r ≤ 27`."

## Counterexamples to the printed statements (exact; two independent programs)

1. **First conjecture, first line, `r = 9`, `n = 1`:** `H_9(M^9) = H_{10}(M^9) = +256`, while the
   printed formula gives `(−1)^1·2^8 = −256`. The same failure occurs for every odd `n` when
   `r ∈ {9, 12, 21, 24}` (e.g. `H_{12}(M^{12}) = +2048`). This is inside the range `r ≤ 27` for which
   the paper states the conjectures as a theorem.
2. **Second conjecture, third line, `r = 4`, `n = 0`:** `H_7(M^4) = 64`, `H_{10}(M^4) = −64`, so the
   sum is `0`, not `±(2·4)^2 = ±64`. Likewise `r = 5`: `H_9 = −1000`, `H_{12} = +1000`. The sum is `0` for
   every `4 ≤ r ≤ 22` with `3 ∤ r` and `n = 0, 1`.
   (The paper's own Theorem for `r = 4` prints `H_{12n+7} = H_{12n+10} = −64(n+1)²`; the true value of
   `H_7(M^4)` is `+64`. First values: `H_0,…,H_{13}(M^4) = 1, 1, −2, −18, −1, −1, 20, 64, 0, 0, −64, −44, 1, 1`.)
3. Second conjecture, first line read literally (four equal quantities): `r = 4`, `n = 0` gives
   `H_0 = H_1 = 1`, `H_4 = H_5 = −1`. (Probably only `|·| = 1` was meant; the paper's own `r = 4`
   theorem has `H_{12n+4} = −1`.)

Programs: `code/motzkin_hankel.py` (python-flint, exact integer determinants; output
`cert/motzkin_hankel_printed.txt`) and `code/motzkin_hankel.gp` (PARI/GP, series from the closed form
of `M`, `matdet`; output `cert/motzkin_hankel_gp.txt`). They agree.

## Repaired statements (tested, NOT proved)

* `3 | r`: `H_{rn}(M^r) = H_{rn+1}(M^r) = (−1)^{n·r(r−1)/2}(n+1)^{r−1}` and
  `|H_{rn+2} + H_{rn−1}| = ½ r(r−3)(n+1)^{r−1}`. Tested for `r = 3, 6, …, 21` (`n ≤ 5`) and
  `r = 24, 27, 30, 33` (`n ≤ 2`).
* `3 ∤ r`: `H_{3rn} = H_{3rn+1} = ±1`, `H_{3rn+r} = H_{3rn+r+1} = ±1`, `H_{3rn+2r} = H_{3rn+2r+1} = 0`,
  `H_{3rn+2r−1} = −H_{3rn+2r+2} = ±(2r(n+1))^{r−2}`, `|H_{3rn+2} + H_{3rn−1}| = ½ r(r−3)`.
  Tested for `4 ≤ r ≤ 22`, `n = 0, 1`.

`cert/motzkin_hankel_repaired.txt`. The errors are sign errors; the structural content of the
conjectures (which determinants are `±1`, `0`, pure powers) is not affected.

## Novelty (26 Nov 2026, `NOVELTY-2.md`, item 3)
Partly known: `the-omega-institute/trureturing` PR #6686 (9 Sep 2026) refutes the printed form of the
second conjecture at `r = 4` (the identically vanishing sum). The sign error of the first conjecture
at `r = 9, 12, 21, 24` was not found anywhere. No corrected or published version of the source exists.
