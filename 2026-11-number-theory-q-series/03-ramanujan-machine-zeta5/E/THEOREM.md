# The three "new and unproven" continued fractions for e in Table 3 of the Ramanujan Machine paper

## 0. Source, statements, status

**Source (read 2026-11-26).** G. Raayoni et al., *Generating conjectures on fundamental constants with the Ramanujan
Machine*, arXiv:1907.00205v4 (Nature 590 (2021) 67–73), Appendix A, Table 3 ("Conjectures of e"). Three rows carry
the label "new and unproven" (defined there: "we have not found this result in the literature, we consider it as
a new conjecture, until proven or until an equivalent form that is unknown to us is found"):

    6e/(−3+2e)  = 7 − 4/(14 − 20/(23 − 54/(34 − 112/...))),      a_n = n(n+6)+7,    b_n = −(n+3)n^2
    3/(−10+4e)  = 3 + 4/(8 + 20/(15 + 54/(24 + 112/...))),       a_n = (n+1)(n+3),  b_n = −(n+3)n^2   [as printed]
    4e/(−1+2e)  = 3 − 3/(7 − 16/(13 − 45/(21 − 96/...))),        a_n = n(n+3)+3,    b_n = −(n+2)n^2

Convention: `PCF[a,b] = a_0 + b_1/(a_1 + b_2/(a_2 + ...)) = lim p_n/q_n`, `u_n = a_n u_{n−1} + b_n u_{n−2}`,
`(p_{−1},p_0) = (1,a_0)`, `(q_{−1},q_0) = (0,1)`. In the second row the displayed fraction has plus signs, i.e.
`b_n = +(n+3)n^2`; the printed `b_n = −(n+3)n^2` is a misprint (with the minus sign the value would be below 3,
whereas `3/(4e−10) = 3.4359...`).

> **Theorem E1.** `PCF[n^2+6n+7, −n^2(n+3)] = 6e/(2e−3)`.
> **Theorem E2.** `PCF[(n+1)(n+3), +n^2(n+3)] = 3/(4e−10)`.
> **Theorem E3.** `PCF[n^2+3n+3, −n^2(n+2)] = 4e/(2e−1)`.
> In each case `p_n/q_n = a_0 − U_0/U_{−1} + 1/(U_{−1}^2 T_n)` with the explicit `U_n`, `T_n` of §2.

**Status and novelty — read this first.** The proofs below are correct and were audited (`AUDIT.md`: "the
mathematics is correct"; every step re-derived and re-computed independently to 60 digits). **They do not settle
open problems: novelty is nil.** The auditor found:

* All three are instances of a published general theorem with proof: O. David, *On Euler polynomial continued
  fractions*, arXiv:2308.02567 (2023), Theorem 1 (`b = −h_1 h_2`, `f(x)a(x) = f(x−1)h_1(x) + f(x+1)h_2(x+1)` with a
  polynomial `f`). The triples are `(h_1, h_2, f) = (x+3, x^2, x^3+3x^2+2x+2)` for E1, `(−x−3, x^2, (x+2)(x^2+5x+5))`
  for E2, `(x+2, x^2, 1)` for E3; the `Π(n)` below is David's `f(n+1)` and `U_n` is his `F_n`, so the method is the
  same as well. (My survey had recognised only the case `f = 1`.)
* E2 is a rescaling of `PCF[n+2, n] = 1/(2e−5)`: `3/(4e−10) = (3/2)·PCF[n+2, n]`, which is a case of
  Dougherty-Bliss–Zeilberger, arXiv:2004.00090 (March 2020), Theorem 1, and is listed as "new and proven" in the
  same Table 3. So E2 was proved in substance before v4 appeared.
* E3 (`f = 1`) is a direct case of Euler's classical series-to-continued-fraction identity.
* H. Cohen's dictionary of continued fractions (arXiv:2607.06581; database described in arXiv:2409.06086) states
  E2 as entry 5.10.5 (`e = 5/2 + (3/4)/E2`) and a tail-equivalent of E3 as entry 5.10.8 (`k = 1`), both without proof.
* For **E1** no explicit statement or proof was found anywhere; it follows from David's theorem by a short
  computation. The Nature supplementary table (Feb 2021) still labels all three "new and unproven".

So this file is an explicit, self-contained write-up of three instances of a known theorem; the only item without
a located explicit prior statement is E1. Checks: `verify/e_checks.py` (polynomial identities with sympy; the
convergent formula with rationals for `n ≤ 40`; 50+ digit agreement); auditor's scripts in `audit/`.
(v4's own title is "The Ramanujan Machine: Automatically Generated Conjectures on Fundamental Constants".)

## 1. Lemma (reduction of order)

> Let `a_n, b_n` be sequences and `U_n` (`n ≥ −1`) a solution of `u_n = a_n u_{n−1} + b_n u_{n−2}` (`n ≥ 1`) with
> `U_n ≠ 0` for all `n ≥ −1`. Put `T_n := Σ_{j=0}^{n} ∏_{i=1}^{j}(−b_i) / (U_j U_{j−1})`. Then
> `q_n = U_n U_{−1} T_n` and `p_n = U_n (1/U_{−1} + (a_0 U_{−1} − U_0) T_n)`; hence, if `T_n → T ∉ {0, ∞}`,
> `PCF[a,b] = a_0 − U_0/U_{−1} + 1/(U_{−1}^2 T)`.
> (If `T_n = 0` for some `n` then `q_n = 0` and that convergent is undefined; in the three applications all terms of
> `T_n` are positive (E1, E3) or `T_n` is an alternating sum with decreasing terms and positive first term (E2), so
> `T_n ≠ 0`.)

*Proof.* For a solution `u`, `W_n := u_n U_{n−1} − u_{n−1} U_n` satisfies `W_n = −b_n W_{n−1}`, so with `s_n := u_n/U_n`:
`s_n − s_{n−1} = W_0 ∏_{i=1}^n(−b_i)/(U_n U_{n−1})`, i.e. `s_n = s_{−1} + W_0 T_n` (the `j = 0` term of `T_n` is `1/(U_0U_{−1})`
and `s_0 − s_{−1} = W_0/(U_0U_{−1})`). For `q`: `s_{−1} = 0`, `W_0 = U_{−1}`. For `p`: `s_{−1} = 1/U_{−1}`, `W_0 = a_0U_{−1} − U_0`. ∎

## 2. Proofs

In all three cases `U_n = (n+1)!^2 Π(n)` with a polynomial `Π`; that `U` solves the recurrence is the polynomial
identity `(n+1)^2 n^2 Π(n) = a_n n^2 Π(n−1) + b_n Π(n−2)` (checked in `e_checks.py` (1)), and `U_n > 0` for `n ≥ −1`.

**E3.** `Π = 1`, `U_{−1} = U_0 = 1`, `∏_{i≤j}(−b_i) = j!^2 (j+2)!/2`, so
`T = Σ_{j≥0} (j+2)/(2 (j+1)!) = (1/2) Σ_{m≥1} (m+1)/m! = (1/2)(e + e − 1) = e − 1/2`, and
`PCF = 3 − 1 + 1/(e − 1/2) = 4e/(2e−1)`.

**E1.** `Π(n) = n^3 + 6n^2 + 11n + 8`, `U_{−1} = Π(−1) = 2`, `U_0 = Π(0) = 8`, `∏_{i≤j}(−b_i) = j!^2 (j+3)!/6`, so

    T = (1/6) Σ_{j≥0} (j+2)(j+3) / ( (j+1)! Π(j) Π(j−1) ).

With `N(j) := −(j^2 + 3j + 1)/3` one has the polynomial identity
`(j+2)(j+3) = (j+1) N(j) Π(j) − N(j+1) Π(j−1) + (1/3) Π(j)Π(j−1)`, i.e.

    (j+2)(j+3)/((j+1)! Π(j)Π(j−1)) = N(j)/(j! Π(j−1)) − N(j+1)/((j+1)! Π(j)) + (1/3)/(j+1)! .

The first two terms telescope (`N(j)/(j!Π(j−1)) → 0`), so `6T = N(0)/Π(−1) + (e−1)/3 = −1/6 + (e−1)/3 = (2e−3)/6`,
`T = (2e−3)/36`, and `PCF = 7 − 8/2 + 1/(4T) = 3 + 9/(2e−3) = 6e/(2e−3)`.

**E2.** `Π(n) = (n+3)(n^2+7n+11)`, `U_{−1} = 2·5 = 10`, `U_0 = 33`, `∏_{i≤j}(−b_i) = (−1)^j j!^2 (j+3)!/6`. With
`c_m := m^2 + 3m + 1` (so `Π(j) = (j+3) c_{j+2}`, `Π(j−1) = (j+2) c_{j+1}`),

    T = (1/6) Σ_{j≥0} (−1)^j / ( (j+1)! c_{j+1} c_{j+2} ).

With `N(j) := −(j+3)/2`: `1 = (j+1) N(j) c_{j+2} + N(j+1) c_{j+1} + (1/2) c_{j+1}c_{j+2}`, i.e.

    (−1)^j/((j+1)! c_{j+1}c_{j+2}) = (−1)^j N(j)/(j! c_{j+1}) − (−1)^{j+1} N(j+1)/((j+1)! c_{j+2}) + (1/2)(−1)^j/(j+1)! ,

so `6T = N(0)/c_1 + (1/2)(1 − 1/e) = −3/10 + 1/2 − 1/(2e) = (2e−5)/(10e)`, `T = (2e−5)/(60e)`, and
`PCF = 3 − 33/10 + 60e/(100(2e−5)) = (−3(2e−5) + 6e)/(10(2e−5)) = 3/(4e−10)`. ∎

(`T ≠ 0` in each case, and the series converge absolutely.) The series for E2 is, after the shift `m = j+1`, the
one recorded in Cohen's database.
