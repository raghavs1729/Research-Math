# R06. The 3-sumfree sequence `S_{1,d−1,2d−1}`: printed conjecture false (modulus; `d = 4`); corrected statement PROVED

**Source.** *On t-sumfree sequences*, arXiv:2609.16843 (16 Sep 2026; citing D. van Berkel and
W. Bosma), Conjecture `conj:except`. The greedy (strict) 3-sumfree sequence `S_{f,g,h}` starts with
`f < g < h` and each further entry is the least integer exceeding the current entries that is not
the sum of 3 distinct previous entries. "Periodic modulo `m` with preperiod length `k` and period
length `p`" means (source's definition): the difference sequence `D_i = S[i+1] − S[i]` is ultimately
periodic with minimal preperiod of length `k` and minimal period of length `p`, the period summing
to `m`; equivalently `S` consists of its first `k` entries `a_1 < … < a_k` followed by all integers
`b_i + jm` (`1 ≤ i ≤ p`, `j ≥ 0`) with `a_k < b_1`.

**Printed.**

> Conjecture. For every `d ≥ 3` the 3-sumfree sequence `S_{1,d−1,2d−1}` is ultimately periodic, and,
> with the exception of `d = 3` or `d = 5`, the preperiod length is `d+3` and the period length is `d`
> modulo `8d−1`.

(preceded by: "although our proof methods apply here and we are convinced the following result is
correct, we have not written down a full proof and so it is correct to state this case as a
conjecture.")

## Theorem (corrected statement)
For every `d ≥ 6`, with `m := 8d − 3`,

    S_{1,d−1,2d−1} = {1, d−1} ∪ [2d−1, 3d−2] ∪ {4d−1} ∪ ⋃_{j≥0} [10d−5+jm, 11d−6+jm].

Hence `S_{1,d−1,2d−1}` is ultimately periodic **modulo `8d − 3`** with preperiod length `d+3` and
period length `d`; its difference sequence is `(d−2, d, 1^{d−1}, d+1, 6d−4, \overline{1^{d−1}, 7d−2})`.
These values are minimal: the block `(1^{d−1}, 7d−2)` contains exactly one entry `≠ 1`, so it is
primitive (least period `d`), and the entry preceding the periodic part is `6d−4 ≠ 7d−2`, so the
preperiod cannot be shortened.

**Consequences for the printed conjecture.**
* For every `d ≥ 6` the modulus is `8d−3`, not `8d−1`: the printed statement is false as written
  (e.g. `d = 6`: `S = 1, 5, 11,…,16, 23, 55,…,60, 100,…,105, …`, modulus `45`, not `47`).
  The paper's own table of small cases is consistent with `8d−3` (e.g. `(2,4,10)`: `m = 43`,
  `(2,8,18)`: `m = 75`, i.e. `m = 8d−1−2f` for `f+g = d`).
* `d = 4` is a further exception, not listed: `S_{1,3,7} = 1, 3, 7, 8, 9, 10, 15, 29, 35, 36, 56, 62, 63, 83, …`,
  difference sequence `(2, 4, 1, 1, 1, 5, 14, \overline{6, 1, 20})`: preperiod length 7 `= d+3`, but
  period length `3 ≠ d`, modulus `27` (by simulation of 60 terms; ultimate periodicity for `d = 4` is
  observed, not proved here; the failure of the printed statement does not depend on it, since
already `D_8 = 6 ≠ D_{12} = 1` contradicts period length 4 after 7 terms). `d = 3`: differences `(1, 3, 1, 1, 4, \overline{5})`; `d = 5`: preperiod 13,
  period 10, modulus 104 (simulation).

## Proof of the Theorem
Let `T` be the right-hand side; `C := [2d−1, 3d−2]`, `E := 4d−1`, `B_j := [10d−5+jm, 11d−6+jm]`.
The greedy rule decides membership of `x > 2d−1` from the entries `< x`. So `S = T` follows by
induction on `x` from: (I) no `x ∈ T`, `x > 2d−1`, is a sum of three distinct smaller elements of `T`;
(II) every `x ∉ T`, `x > 2d−1`, is such a sum. (Seeds: `1 < d−1 < 2d−1` since `d ≥ 3`.)

**(I).**
* `x ∈ [2d, 3d−2]`: the least sum of three distinct elements is `1 + (d−1) + (2d−1) = 3d−1 > x`.
* `x = 4d−1`: sums of three distinct elements of `{1, d−1} ∪ C`: `1 + (d−1) + c = d + c ≤ 4d−2`;
  any sum with two elements of `C` is `≥ 1 + (2d−1) + 2d = 4d`. Neither equals `4d−1`.
* `x ∈ B_j`. Reduce modulo `m`: `B_j ≡ [2d−2, 3d−3]` (as `10d−5 = (2d−2) + m`). The residues (in
  `[0,m)`) of the elements of `T` are `a = 1`, `b = d−1`, `e = 4d−1` (one element each) and
  `U := [2d−2, 3d−2]` (elements of `C` and of the blocks; the residue `3d−2` is attained only by the
  single element `3d−2`). If `x = y_1+y_2+y_3`, the sum `σ` of the three residues satisfies
  `σ ≡ x (mod m)` and `0 < σ < 2m`, so `σ ∈ [2d−2, 3d−3] ∪ [10d−5, 11d−6]`. But, by type:
  `UUU`: `σ ∈ [6d−6, 9d−6]`; `aUU`: `[4d−3, 6d−3]`; `bUU`: `[5d−5, 7d−5]`;
  `eUU`: `σ ≤ (4d−1) + (3d−2) + (3d−3) = 10d−6` (two distinct elements cannot both have residue
  `3d−2`) and `σ ≥ 8d−5`; `abU`: `[3d−2, 4d−2]`; `aeU`: `[6d−2, 7d−2]`; `beU`: `[7d−4, 8d−4]`;
  `abe`: `5d−1`. None of these meets `[2d−2, 3d−3] ∪ [10d−5, 11d−6]` (for `d ≥ 2`). Contradiction.

**(II).** Sums of two distinct elements of `C` fill `[4d−1, 6d−5]`, of three fill `[6d, 9d−9]`.
* `x ∈ [3d−1, 4d−2]`: `x = 1 + (d−1) + c` with `c = x − d ∈ C`.
* `x ∈ [4d, 10d−6]`: `1 + c_1 + c_2` gives `[4d, 6d−4]`; `(d−1) + c_1 + c_2` gives `[5d−2, 7d−6]`;
  `c_1+c_2+c_3` gives `[6d, 9d−9]`; `E + c_1 + c_2` gives `[8d−2, 10d−6]`. These four intervals cover
  `[4d, 10d−6]` because `5d−2 ≤ 6d−3`, `6d ≤ 7d−5` (`d ≥ 5`) and `8d−2 ≤ 9d−8` (**`d ≥ 6`**). All
  summands are `≤ 4d−1 < x`, and `E` is used only for `x ≥ 8d−2`.
* `x` in the gap after `B_j`: `x = y + (j+1)m` with `y ∈ [3d−2, 10d−6]`. Take `β ∈ B_j`, i.e.
  `β = ρ + (j+1)m`, `ρ ∈ [2d−2, 3d−3]`, and two distinct elements `s_1, s_2` of `{1, d−1} ∪ C ∪ {E}`:
  `s_1+s_2 = d` gives `y ∈ [3d−2, 4d−3]`; `1 + c` gives `[4d−2, 6d−4]`; `(d−1) + c` gives
  `[5d−4, 7d−6]`; `c_1 + c_2` gives `[6d−3, 9d−8]`; `E + c` gives `[8d−4, 10d−6]`. The union is
  `[3d−2, 10d−6]` (`d ≥ 3`). `β < x`, and `β` exceeds `s_1, s_2`.

These cases exhaust the complement of `T` above `2d−1`: `C`, `[3d−1,4d−2]`, `E`, `[4d,10d−6]`, `B_0`,
gap, `B_1`, … (`11d−5+jm = (3d−2) + (j+1)m`, `10d−6+(j+1)m = 10d−5+(j+1)m − 1`). ∎

For `d = 5` the last covering fails at `x = 37 = 8d−3`, and indeed `37 ∈ S_{1,4,9}`. For `d = 4` the four
families of sums used above leave `x = 23` uncovered, but `23 = 1+7+15` is still a 3-sum; the first
actual deviation from the pattern for `d = 4` is `x = 29 ∈ S_{1,3,7}` (found by the referee,
`AUDIT-2.md`). This is why `d ≥ 6` is needed.

## Theorem 2 (the whole boundary family `f + g = d`)
The source's main theorem (`thm:gend>1`) covers `S_{f,g,g+d}` for `f+g ≥ d+1`; the conjecture above is
the case `f = 1` of the boundary `f + g = d`, and the other boundary cases appear in the source's
tables of unexplained small cases (e.g. `(2,4,10)`, `(2,8,18)`, `(3,7,17)`, `(4,6,16)`). The same
argument gives all of them:

**Theorem 2.** Let `f ≥ 1`, `d ≥ max(6, 2f+1)`, `g = d − f` and `m := 8d − 2f − 1`. Then

    S_{f, d−f, 2d−f} = {f, d−f} ∪ [2d−f, 3d−f−1] ∪ {4d−f} ∪ ⋃_{j≥0} [10d−3f−2+jm, 11d−3f−3+jm].

Hence `S_{f,d−f,2d−f}` is ultimately periodic modulo `8d−2f−1` with (minimal) preperiod length `d+3` and
period length `d`; its difference sequence is
`(d−2f, d, 1^{d−1}, d+1, 6d−2f−2, \overline{1^{d−1}, 7d−2f})`. The pattern fails for
`(f,d) = (1,3), (1,4), (1,5), (2,5)`, the only pairs with `2f+1 ≤ d < 6`.

*Proof.* As for the Theorem, with `C := [2d−f, 3d−f−1]`, `E := 4d−f`,
`B_j := [10d−3f−2+jm, 11d−3f−3+jm]`; `g = d−f > f` since `d ≥ 2f+1`.

**(I).** `x ∈ C`, `x > 2d−f`: the least sum of three distinct elements is `f + g + (2d−f) = 3d−f > x`.
`x = E`: `f+g+c = d+c ∈ [3d−f, 4d−f−1]`, and any sum with two elements of `C` is
`≥ f + (2d−f) + (2d−f+1) = 4d−f+1`. `x ∈ B_j`: modulo `m` the blocks are `≡ [2d−f−1, 3d−f−2]`, the
elements of `T` have residues `a = f`, `b = d−f`, `e = 4d−f` and `U = [2d−f−1, 3d−f−1]` (the residue
`3d−f−1` only from one element); a representation would give a residue sum
`σ ∈ [2d−f−1, 3d−f−2] ∪ [10d−3f−2, 11d−3f−3]` (`σ < 2m`). By type:
`UUU`: `[6d−3f−3, 9d−3f−3]`; `aUU`: `[4d−f−2, 6d−f−2]`; `bUU`: `[5d−3f−2, 7d−3f−2]`;
`eUU`: `[8d−3f−2, 10d−3f−3]` (two distinct elements cannot both have residue `3d−f−1`);
`abU`: `[3d−f−1, 4d−f−1]`; `aeU`: `[6d−f−1, 7d−f−1]`; `beU`: `[7d−3f−1, 8d−3f−1]`; `abe`: `5d−f`.
All lie strictly between `3d−f−2` and `10d−3f−2` (using `d ≥ 2f+1`). Contradiction.

**(II).** `x ∈ [3d−f, 4d−f−1]`: `x = f + g + c`, `c = x−d ∈ C`.
`x ∈ [4d−f+1, 10d−3f−3]`: `f+c_1+c_2` gives `[4d−f+1, 6d−f−3]`; `g+c_1+c_2` gives `[5d−3f+1, 7d−3f−3]`;
`c_1+c_2+c_3` gives `[6d−3f+3, 9d−3f−6]`; `E+c_1+c_2` gives `[8d−3f+1, 10d−3f−3]`. They cover the
interval because `5d−3f+1 ≤ 6d−f−2`, `6d−3f+3 ≤ 7d−3f−2` (`d ≥ 5`) and `8d−3f+1 ≤ 9d−3f−5` (`d ≥ 6`).
Gap after `B_j`, `x = y + (j+1)m`, `y ∈ [3d−f−1, 10d−3f−3]`: with `β ∈ B_j` (residue
`ρ ∈ [2d−f−1, 3d−f−2]`) and two distinct small elements: `f+g` gives `y ∈ [3d−f−1, 4d−f−2]`; `f+c` gives
`[4d−f−1, 6d−f−3]`; `g+c` gives `[5d−3f−1, 7d−3f−3]`; `c_1+c_2` gives `[6d−3f, 9d−3f−5]`; `E+c` gives
`[8d−3f−1, 10d−3f−3]`; these cover `[3d−f−1, 10d−3f−3]` (`d ≥ 3`). All summands are smaller than `x`.
Minimality as before (`6d−2f−2 ≠ 7d−2f`). ∎

For `f = 1` this is the Theorem above. Checks: `code/sumfree3_family_verify.py` — simulation of all
1186 sequences with `f ≤ 14`, `max(6,2f+1) ≤ d ≤ 100` (first `12d+40` terms), and a brute-force
check of the statements (I), (II) for `f ≤ 8`, `d ≤ 40` (`cert/sumfree3_family_verify.txt`). The
moduli agree with the source's table: `(2,4,10) ↦ 43`, `(2,8,18) ↦ 75`, `(3,7,17) ↦ 73`, `(4,6,16) ↦ 71`.
Theorem 2 was refereed separately (`AUDIT-3.md`): **proof complete**; the referee checked all 46
linear inequalities on the whole region `f ≥ 1`, `d ≥ max(6, 2f+1)` exactly, and simulated all 22 346
pairs with `f ≤ 149`, `d ≤ 300` (own C program) with no mismatch. Remarks from that report: the upper
ends printed for the types `UUU`, `aUU`, `bUU` are enclosures, not attained values; the covering in
(II) uses that the right ends increase (`d ≥ 2f`); (I) for the blocks only needs `d ≥ f+1`; `d ≥ 6` is
used exactly once (`8d−3f+1 ≤ 9d−3f−5`), and for `d = 5` exactly the integer `8d−3f` is uncovered.

**Novelty caveat (from `AUDIT-3.md`).** The source's tables list all 120 boundary triples with `f ≥ 2`,
`6 ≤ d ≤ 25`, each with `(k,p,m) = (d+3, d, 8d−2f−1)`, and the source states that its tabulated values
are proved by its lemma plus computation. So those 120 instances are already established there. New
here: the general formula `8d−2f−1`, the case `f = 1` (the printed conjecture, corrected), all `d > 25`,
and a computation-free proof. The formula appears nowhere in the source.

## Checks
`code/sumfree3_verify.py`: for all `6 ≤ d ≤ 400` the greedy simulation (bitset implementation,
which reproduces the paper's example `D_{1,2,3} = (1⁴, 8, 1, \overline{1, 10, 1², 10})`) agrees with
`T` on the first `25d+50` terms (`cert/sumfree3_verify.txt`); `code/sumfree3b.py` (period finder).

## Novelty (26 Nov 2026, `NOVELTY-2.md`, item 2)
No correction or proof found. The conjecture appears in v2 of the source (11 Nov 2026), which still
prints "modulo `8d−1`" and the exceptions `d = 3, 5`. A public issue (`ndesmo/leaning` #24, opened
25 Nov 2026) targets the conjecture as printed. Time-sensitive.
