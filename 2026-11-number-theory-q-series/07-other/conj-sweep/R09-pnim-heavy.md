# R09. P-NIM on partitions: the "chopped heavy rectangle" conjecture is FALSE (`a = 8`, `b = 7`)

**Source.** E. Gottlieb, M. Krnc, P. Muršič, *Nim on Integer Partitions and Hyperrectangles*,
arXiv:2506.04991v1 (5 June 2025; only version), Conjecture `conj:heavychoppedrect`.
P-NIM on one partition: "Players alternate choosing a nonempty partition and removing either a
positive number of rows or a positive number of columns from its Young diagram" (the removed rows
or columns need not be adjacent; the rest is merged); normal play. A position is **heavy** if its
Sprague–Grundy value equals the length of the longest play, which for `λ = [λ_1,…,λ_r]` is
`λ_1 + r − 1`. `≤` is the order of Young's lattice (`ℓ ≤ m` and `λ_j ≤ μ_j`).

**Printed.**

> Conjecture. Let `a, b ∈ ℤ` be such that `[(a+1)^{b+1}]` is heavy. Then any `λ` satisfying
> `[a+1, a, …, a−b+1] ≤ λ ≤ [(a+1)^{b+1}]` is heavy.

(The authors' data: "In [the appendix] we list all heavy partitions of `n ≤ 8`".)
By the source's theorem on rectangles, `[c^r]` is heavy iff `C(c+r−2, r−1)` is odd.

## Counterexample
`a = 8`, `b = 7`. The rectangle `[9^8]` is heavy: `C(15,7) = 6435` is odd, and directly
`SG([9^8]) = 16 = 9 + 8 − 1`. The partition

    λ = [9, 9, 7, 7, 6, 5, 4, 2]     (size 49)

satisfies `[9,8,7,6,5,4,3,2] ≤ λ ≤ [9^8]`, its longest play has length 16, but **`SG(λ) = 3`**. So `λ`
is not heavy.

In total 40 of the 1430 partitions in this interval are not heavy (Grundy values 3 or 4), e.g.
`[9,9,7,7,6,5,5,2]`, `[9,9,7,7,6,6,4,2]`, `[9,8,8,8,8,6,5,2]` (`SG = 4`), `[9,9,9,9,9,7,6,3]` (`SG = 4`).

## Verification
Two independent programs compute the Grundy values:
* `code/pnim.py`: positions as boundary words, a move deletes a nonempty set of `N`'s or of `E`'s;
* `code/pnim_check2.py`: positions as tuples, a nonempty sub-multiset of the parts of `λ` or of
  its conjugate is removed.

Both give the values above (`code/pnim_cex.py`, output `cert/pnim_counterexample.txt`). Sanity
checks against the source's own theorems: the characterisation of Grundy value 1
(`SG(λ) = 1` iff `λ = [1]` or `[r,r,r−1,…,2] ≤ λ ≤ [r^r]`) is reproduced for all 11 731 partitions of
`n ≤ 26`, the two programs agree on all partitions of `n ≤ 18`, and the criterion for heavy
rectangles is reproduced for `r, c ≤ 9`.

## What survives (tested only)
* The same conjecture for all other pairs `(a,b)` with `b ≤ a ≤ 9` and `C(a+b,b)` odd, `a+b ≤ 16`:
  `(a,b) ∈ {(2,1), (4,1), (4,2), (4,3), (5,2), (6,1), (8,1), …, (8,6), (9,2), (9,4), (9,6)}` and `b = 0`:
  no failure. So `(8,7)` is the smallest failure with `a ≥ b`; no repaired statement is proposed.
* The companion conjecture `conj:shallowstair` (`[i+(k−1)s, …, i+s, i]` is heavy): no failure among the
  420 such partitions with `rows + columns ≤ 22`.

## Novelty
arXiv:2506.04991 has only v1; its single citing paper on Semantic Scholar (26 Nov 2026) is
arXiv:2606.16828 (CRIM), which does not discuss this conjecture. No further search was made
(see `NOVELTY-2.md` if it lists this item).

## Novelty update: NOT NEW
A search of public GitHub issues (26 Nov 2026) shows that the same failure at `a = 8`, `b = 7` (with
`λ = [9,9,9,9,9,7,5,3]`) was recorded in `the-omega-institute/trureturing`, issue #8467 and pull requests
#8543 / #8565 ("PnimHeavyIntervalRefutation", 17–18 Sep 2026). The counterexample above is an
independent rediscovery; `[9,9,7,7,6,5,4,2]` is the smallest (by size) of the 40 in that interval.
