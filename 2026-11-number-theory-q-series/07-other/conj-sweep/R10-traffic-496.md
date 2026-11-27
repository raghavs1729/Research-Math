# R10. Points of maximal traffic with an obstruction on the antidiagonal: `ρ(n) < 1` for all `n ≥ 496` — PROVED

**Source.** *Points of maximal traffic on a grid with obstruction*, arXiv:2609.01562 (1 Sep 2026),
Conjecture `conj:496`. Lattice paths from `(0,0)` to `(n,n)` with unit north/east steps; `f_B(A)` =
number of paths through `A` avoiding `B`.

**Printed.**

> Conjecture. For every `n ≥ 496` and every obstruction `B` on the antidiagonal `x+y = n`, the maximum
> of `f_B` is attained at `(1,1)` and `(n−1,n−1)`. Equivalently, `ρ(n) < 1` for all `n ≥ 496`.

Here (source, Section on the antidiagonal) `ρ(n) = max_{1≤a<n/2} R_n(a)`, `R_n(a) = G(n,a)/D(n)`,

    G(n,a) = [2 − (n−1)/(n−a)]·C(n,a)·C(n−2,a−1),      D(n) = 2C(2n−2,n−1) − C(2n−1,n−1).

The source proves that some antidiagonal obstruction has its maximum at boundary points precisely
when `ρ(n) > 1`, reports `ρ(n) > 1` for `8 ≤ n ≤ 375` and for 82 values in `376..495`, and "For
`496 ≤ n ≤ 2000` we have verified that `ρ(n) < 1`."

## Theorem
`ρ(n) < 1` for every `n ≥ 496`. (By the equivalence stated in the source, this is `conj:496`. We prove
the second, explicit form; the equivalence itself is the source's.)

Moreover `ρ(n) ≤ (4/√(2πe))·e^{4.128/√n}` for all `n ≥ 12`; `4/√(2πe) = 0.96788…`.

## Proof

**Step 0 (a simpler formula).** `D(n) = C(2n−2,n−1)/n`, `C(n−2,a−1) = C(n,a)·a(n−a)/(n(n−1))` and
`C(2n−2,n−1) = C(2n,n)·n/(2(2n−1))`. Hence, with `u := n − 2a ≥ 1`,

    R_n(a) = (u+1) · a/(n−1) · 2(2n−1)/n · C(n,a)² / C(2n,n).                                   (1)

(Checked against the definition for all `3 ≤ n < 60` in `code/traffic496.py`.)

**Step 1 (binomial bounds).**
(B1) For `m ≥ 1`: `4^m/√(π(m+½)) < C(2m,m) < 4^m/√(π(m+¼))`. [With `w_m = C(2m,m)/4^m`,
`w_{m+1}/w_m = (2m+1)/(2m+2)`; then `w_m²(m+¼)` is increasing and `w_m²(m+½)` is decreasing
(compare `(2m+1)²(m+5/4)` with `(2m+2)²(m+¼)`, resp. `(2m+1)²(m+3/2)` with `(2m+2)²(m+½)`), and both
tend to `1/π` by Wallis' product.]
(B2) `n = 2m` even, `a = m−t`, `u = 2t`: `C(n,a)/C(2m,m) = ∏_{i=1}^{t}(1 − t/(m+i)) ≤ exp(−t²/(m+t))`.
(B3) `n = 2m+1` odd, `a = m−j`, `u = 2j+1`: `C(n,a)/C(n,m) = ∏_{i=1}^{j}(1 − (j+1)/(m+1+i)) ≤
exp(−j(j+1)/(m+1+j))`, and `C(2m+1,m) = ½C(2m+2,m+1)`.

From these, for both parities of `n` and all `1 ≤ a < n/2`:

    C(n,a)²/C(2n,n) < (2/√(π(n+½)))·exp(−(u²−1)/(n+u)).                                        (2)

[Even: `C(2m,m)²/C(4m,2m) < 16^m/(π(m+¼)) · √(π(2m+½))/16^m = 2/√(π(n+½))`, exponent
`2t²/(m+t) = u²/(n+u)`. Odd: `C(2m+1,m)²/C(4m+2,2m+1) < (4·16^m/(π(m+5/4)))·√(π(2m+3/2))/(4·16^m)
= 2√(n+½)/(√π(n+3/2)) < 2/√(π(n+½))`, exponent `2j(j+1)/(m+1+j) = (u²−1)/(n+u)`.]

**Step 2 (uniform bound).** In (1), `a/(n−1) · 2(2n−1)/n = (n−u)(2n−1)/(n(n−1)) < 2` since `u ≥ 1`. So by (2)

    R_n(a) < (4/√(π(n+½))) · φ(u),      φ(u) := (u+1)·exp(−(u²−1)/(n+u)).                      (3)

Put `v = u+1 ≥ 2`. Since `(u²−1)/(n+u) = (v²−2v)/(n+v−1) ≥ v²/(n+v) − 2v/n ≥ v²/n − v³/n² − 2v/n`,

    φ(u) ≤ v e^{−v²/n} · exp(v³/n² + 2v/n).

* If `v ≤ 1.2√n`: `v e^{−v²/n} ≤ √(n/(2e))` (maximum over all real `v`) and
  `v³/n² + 2v/n ≤ (1.728 + 2.4)/√n = 4.128/√n`. So `φ(u) ≤ √(n/(2e))·e^{4.128/√n}`.
* If `u ≥ 1.2√n − 1`: `(log φ)′(u) = 1/(u+1) − (u²+2un+1)/(n+u)² < 0` as soon as `(n+u)² ≤ 2u²n`, i.e.
  `u ≥ n/(√(2n) − 1)`, which holds for `u ≥ 1.2√n − 1` when `n ≥ 12`
  (`(1.2√n − 1)(√(2n) − 1) ≥ n ⇔ 0.697n − 2.615√n + 1 ≥ 0`). So `φ` is decreasing there and bounded
  by its value at `v = 1.2√n`, covered by the first case.

Hence for `n ≥ 12`:  `ρ(n) < (4/√(π(n+½)))·√(n/(2e))·e^{4.128/√n} < (4/√(2πe))·e^{4.128/√n}`.

**Step 3 (large `n`).** `(4/√(2πe))·e^{4.128/√n} < 1` iff `4.128/√n < ln(√(2πe)/4) = 0.03264417…`,
i.e. iff `n > 15990.72…`. So `ρ(n) < 1` for all `n ≥ 15991`.

**Step 4 (`496 ≤ n ≤ 16100`, exact).** By (1), `R_n(a) < 1` is the integer inequality
`(u+1)·a·2(2n−1)·C(n,a)² < (n−1)·n·C(2n,n)`. Also by (1),
`R_n(a−1) < R_n(a) ⇔ a(a−1)(n−2a+3) < (n−a+1)²(n−2a+1)` (this is the source's Lemma `lem:monotone`);
so if `a_min` is the least `a ≥ 2` violating this inequality, the maximum of `R_n` is attained on
`[a_min − 1, n/2)`. The program `code/traffic496.py` checks, in exact integer arithmetic (FLINT
`fmpz`), the first inequality for every `a` in that window, for every `496 ≤ n ≤ 16100`: it holds
strictly in all cases. ∎

## Consistency with the source's data
The same program reproduces the source's statements exactly: `ρ(n) > 1` for all `8 ≤ n ≤ 375`; in
`376..495` the values with `ρ(n) ≤ 1` are `376, 378, 423, 425, …, 463` (odd), `466, 468, …, 494` (even);
the last `n` with `ρ(n) > 1` is `495`; `ρ(n) = 1` never occurs for `8 ≤ n ≤ 16100`
(`cert/traffic496_cert.txt`).

## Status
**Audit** (`AUDIT-4.md`, independent referee with own programs): PROOF COMPLETE, no error. The referee
re-verified Step 4 in exact integers for **every** `1 ≤ a < n/2` and every `3 ≤ n ≤ 16100`
(64 794 450 pairs, no monotonicity used), checked all displayed inequalities in ball arithmetic for
`n ≤ 3000`, and confirmed by brute-force path counting that the verbal form holds at
`n = 496, 497, 498` and fails at `n = 9, 20, 375, 377, 495`. Closest call above 495: `ρ(497) = 0.9999553`.
Remarks from that report: `G`, `D` are quoted here in the form used inside the source's proof of its
criterion proposition (same values as its displayed `G = (n−2a+1)/(n−a)·C(n,a)C(n−2,a−1)`,
`D = C(2n−2,n−1)/n`); the source reports `ρ(n) > 1` already from `n = 5`; the chain in `v` of Step 2 is
valid for all real `v ≥ 1`, which is what the monotonicity case uses; the verbal form of the
conjecture rests on the source's own reductions (classification of maximisers for `n ≥ 9`, the
diagonal case `2a = n`, the corner cases `a ∈ {0, n}`), which are not re-proved here.

**Novelty: NOT NEW.** The same inequality (`G(n,a)/D(n) < 1` for all `n ≥ 496`, `1 ≤ a`, `2a < n`) was
formalised in Lean 4 in the public project `the-omega-institute/trureturing` (issue #8249, PR #8934,
merged 20 Sep 2026; 652 lines, a Stirling/entropy envelope plus local certificates; not compiled or
read here). The earlier search by arXiv id missed it because the issue writes `2609.01562v1`.
The proof above is an independent, short human-readable proof by a different route (uniform
Gaussian bound + exact finite check).
