# ESMA tables (arXiv:2212.09470v2, ICML 2023): numerical status of the 26 "new and unproven" rows — work stopped here

Stopped on instruction of the parent session (the rows are Hurwitz-type signed continued fractions of Möbius
transforms of `e`, `e^2`, `√e`, `tanh(1/4)`, `tan 1`, `J_0(1)/J_1(1)`; expected novelty low). **Nothing is proved
here.** What exists: `verify/rows.py` (transcription of Tables 1–3: sign period, initial values, printed linear
recurrence, printed value as a Möbius matrix in the base constant) and `verify/esma_numeric.py` (60-digit
evaluation of each fraction, depth 600 and 900, compared with the printed value; PSLQ on mismatch).

* 24 of the 26 rows agree with their printed values to 60 digits (numerical only).
* **Table 3, row "(−12047 J_1(1) + 6928 J_0(1)) / (8 (1777 J_0(1) + 3090 J_1(1)))"** (signs `{−1,−1,−1,1,1}` repeated,
  `a = [1,112,2,143,1,2,2,45,1,1,3,2,53,1,1, 3,240,4,271,1,4,2,77,1,1,5,2,85,1,1]`, period 15): **false as printed.**
  The fraction equals `0.991031249607386591...`; the printed expression equals `7.63·10^{−9}`. PSLQ (60 digits):
  the value is `(6928 J_0(1) − 12047 J_1(1)) / (8 (3090 J_1(1) − 1777 J_0(1)))`, i.e. the sign of `1777 J_0(1)` in the
  denominator is misprinted (`3090 J_1(1) − 1777 J_0(1) = J_6(1)`). Numerical only. Read from the rendered PDF page,
  not only from the text extraction.
* **Table 1, row "(−3+5e)/(−6+6e) = 1 + 1/(36 + 1/(2 − 1/(4 − ...)))"** (signs `{1,1,−1,−1}` repeated): **internally
  inconsistent.** The printed coefficient data (`[2,2,18,3,3,4,3,5,2,2,6,2, 2,6,54,7,3,8,3,9,2,2,10,2]` and the
  matching piecewise formula for `a_j`) give `2.48620940761916...`, which is not a Möbius transform of `e` with small
  integer coefficients for any of the 4096 sign patterns of period 12 (PSLQ, coefficients ≤ 2000). The displayed
  fraction is consistent with the value: expanding `(5e−3)/(6e−6) = 1.0273255689...` with the printed sign pattern
  gives `a = [1,36,2,4,1,8,28,2,1,12,2,2, 5,180,6,4,1,24,64,2,1,28,2,2, 9,324,10,4,1,40,100,2,1,44,2,2, ...]`
  (period 12, increments `4,144,4,0,0,16,36,0,0,16,0,0`; computed from 400-digit arithmetic, first 60 terms). So
  the coefficient list printed in that row belongs to some other fraction. Numerical only.

A uniform proof strategy (not carried out): for each row find a constant matrix `U` with
`U·B_k = ±(E_{dk+r} ··· E_{dk+r+d−1})·U`, where `B_k` is the product of the row's elementary matrices over one period
and `E_j = [[αj+β, ±1],[1, 0]]` are those of Lambert's fractions for `coth(1/m)`, `cot 1`, or of the Bessel-ratio
fraction; this is a polynomial identity in `k`, and gives the value as a Möbius transform of the base constant.
