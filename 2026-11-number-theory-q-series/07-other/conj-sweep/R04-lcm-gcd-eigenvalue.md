# R04. `−1` as a generalized eigenvalue of the LCM matrix to the GCD matrix: "if" direction PROVED

**Source.** arXiv:2601.16626 (23 Jan 2026), *On generalized eigenvalues of MAX matrices to MIN
matrices and of LCM matrices to GCD matrices*, Conjecture `conj`:

> Let `T = {1,…,n}`, `n > 3`. Then `−1` is a g-eigenvalue of `L_T` to `G_T` if and only if the binary
> representation of `n` begins with `10`.

Here `L_T = [lcm(i,j)]`, `G_T = [gcd(i,j)]` (`1 ≤ i,j ≤ n`), and `−1` is a g-eigenvalue iff
`det(L_T + G_T) = 0`. The authors report computer experiments for `1 ≤ n ≤ 1000`.

## Theorem (the "if" direction)
If `2^k ≤ n < 3·2^{k−1}` for some `k ≥ 2` (i.e. `n ≥ 4` and the binary expansion of `n` begins with
`10`), then the rows `R_i` of `M_n := [lcm(i,j) + gcd(i,j)]` satisfy

    R_{2^k} − 3 R_{2^{k−1}} + 2 R_{2^{k−2}} = 0,

so `det M_n = 0` and `−1` is a g-eigenvalue of `L_T` to `G_T`.

*Proof.* Put `a = 2^{k−1}`, `b = 2^k`, `c = 2^{k−2}`; all three are `≤ n`. Let `1 ≤ j ≤ n`, `j = 2^e v`
with `v` odd. Since `n < 3·2^{k−1}`, the only `j ≤ n` with `e ≥ k−1` are `j = a` and `j = b`.

* `e ≤ k−2`: `gcd(b,j) = gcd(a,j) = gcd(c,j) = 2^e =: g` and `lcm(x,j) = xj/g` for `x ∈ {a,b,c}`. The
  combination is `g(1 − 3 + 2) + (j/g)(b − 3a + 2c) = 0 + (j/g)·2^{k−2}(4 − 6 + 2) = 0`.
* `j = a`: entries `lcm+gcd` are `b + a = 3a` (row `b`), `2a` (row `a`), `a + c = 3a/2` (row `c`):
  `3a − 6a + 3a = 0`.
* `j = b`: entries are `2b = 4a` (row `b`), `b + a = 3a` (row `a`), `b + c = 5a/2` (row `c`):
  `4a − 9a + 5a = 0`. ∎

(For `n = 3` one also has `det M_3 = 0`; `n = 3` is excluded in the conjecture and is not covered by this relation, which needs `k ≥ 2`.)

## The "only if" direction: NOT proved
For `3·2^{k−1} ≤ n < 2^{k+1}` the conjecture asserts `det M_n ≠ 0`. A structural reduction: with
`i = 2^a u`, `u` odd, the row operations `R_{2^a u} − 3R_{2^{a−1}u} + 2R_{2^{a−2}u}` (`a ≥ 2`) leave in
row `(a,u)` only the entries in the columns `2^{a−1}v` (`v` odd), equal to
`2^{a−1}(lcm(u,v) − gcd(u,v))`. So `M_n` is equivalent to a block matrix whose rows of 2-adic level
`a ≥ 2` are supported on the columns of level `a−1`; the determinant vanishes iff this structured
matrix is singular. The determinants for `n = 6, 7, 12,…,15, 24,…,30` have large prime factors
(e.g. `7407624403 | det M_26`), so no product formula exists; we found no argument for
non-singularity. Verified (`det ≠ 0`) only in the authors' range.

Numerical check of the theorem and of the small determinants: `code/lcmgcd.gp`.

## Novelty
NOT NEW. `NOVELTY-1.md`, item 4: the same relation (kernel vector `(−2, 3, −1)` on the rows
`r, 2r, 4r`) is recorded in a public GitHub issue of the project `the-omega-institute/trureturing`
(#9180, 20 Sep 2026). The "only if" direction is open there too; no proof was found anywhere. The
source is published as J. Integer Seq. 28 (2025), Art. 25.7.1, with the conjecture unchanged.
