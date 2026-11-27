# R01. Unimodality of `[a_1]_q … [a_k]_q [b]_{q^r}`: the necessity clause is false

**Source.** B. B. Connelly, E. Ito, T. C. Martinez, O. Shevchenko, K. Yang,
*Unimodality of q-Fibonomial coefficients for small cases*, arXiv:2605.12822 (12 May 2026),
Conjecture labelled `conj:unimodalitylemma` (stated after Proposition `prop:unimodal`).

**Statement as printed** (their `\pow` is written `r` here; `[a]_q = 1+q+…+q^{a-1}`):

> Let `r ≥ 2, k ≥ 1` and let `a_1,…,a_k,b` be positive integers. Then if `r | a_i` for any
> `1 ≤ i ≤ k` or `b ≤ 1 + Σ_{i=1}^k ⌊a_i/r⌋`, the polynomial `[a_1]_q…[a_k]_q[b]_{q^r}` is
> unimodal. Moreover, if `k ≤ 3` or `r ≤ 3`, this condition is also necessary.

The authors add: "We have computationally verified [the conjecture] for `k ≤ 5`, `r ≤ 6` and
`max{a_i,b} ≤ 15`. In general, the condition is not necessary: `([3]_q)^4[2]_{q^4}` is unimodal".

## Result: the second sentence is FALSE for `r = 3`, `k = 6`

**Counterexample.** `r = 3, k = 6, a_1 = … = a_6 = 2, b = 2`. Then `3 ∤ 2` and
`b = 2 > 1 + 6·⌊2/3⌋ = 1`, so the condition fails; `r ≤ 3`, so the conjecture asserts that the
polynomial is not unimodal. But

    [2]_q^6 [2]_{q^3} = (1+q)^6 (1+q^3)
                      = 1 + 6q + 15q^2 + 21q^3 + 21q^4 + 21q^5 + 21q^6 + 15q^7 + 6q^8 + q^9,

which is unimodal (`(1+q)^6` has coefficients 1, 6, 15, 20, 15, 6, 1; adding the copy shifted by 3
gives the list above). This is a hand check. It lies just outside the authors' tested range
(`k ≤ 5`).

More counterexamples with `r = 3`: all 167 parameter sets found for `k = 6`, `a_i ≤ 13`, `b ≤ 20`
(e.g. `a = (2,2,2,2,2,5), b ∈ {3,4}`; `a = (2,2,2,2,2,2), b = 3`), and 184 more for `k = 7`,
`a_i ∈ {2,4,5,7,8}`. Each was confirmed by two independent implementations
(numpy sliding-window products in `code/qint_unimodal.py`; plain integer convolution in the
session log, reproduced by `code/qint_check2.py`).

## What survives (tested, NOT proved)

* The first sentence (sufficiency): no failure for
  `k = 1` (`r ≤ 8`, `a,b ≤ 18`), `k = 2` (`r ≤ 30`, `a ≤ 70`, `b ≤ 12`), `k = 3` (`r ≤ 16`, `a ≤ 40`,
  `b ≤ 10`), `k = 4` (`r ≤ 8`, `a ≤ 20`, `b ≤ 10`; `r ≤ 3`, `a,b ≤ 30`), `k = 5` (`r ≤ 3`, `a ≤ 22`,
  `b ≤ 24`), `k = 6` (`r ≤ 3`, `a ≤ 13`, `b ≤ 20`), `k = 7, 8` (`r = 2`, `a ≤ 9`, `b ≤ 14`).
* Necessity for `k ≤ 3`: no failure in the ranges above (`k ≤ 3`).
* Obvious repair "necessary if `k ≤ 5` or `r = 2`": no failure for `r = 3`, `k ≤ 5`
  (`a ≤ 22`, `b ≤ 24`), nor for `r = 2`, `k ≤ 14`, `a_i ∈ {3,5,7,9}`, `b` up to 7 above the bound.
  (`k = 4` already fails for `r = 4`, by the authors' own example.)

Certificates: `cert/qint_*.txt`.

## Novelty check
NOT NEW. `NOVELTY-1.md`, item 3: the same counterexample `(1+q)^6(1+q^3)` is recorded, with a Lean
formalisation, in the public GitHub project `the-omega-institute/trureturing` (PR #6805, merged
10 Sep 2026). It is not in any arXiv version of the source (v1 is the only one) nor in a citing
paper. The larger families of counterexamples (`k = 6, 7`) and the tests of the repaired statement
were not found elsewhere. A Lean check of the counterexample is also in `lean/ConjSweep.lean`.
