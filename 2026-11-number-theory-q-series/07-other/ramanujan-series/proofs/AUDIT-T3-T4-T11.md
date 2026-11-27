# AUDIT of RESULTS.md — Theorems 11, 3, 4 (and section 0)

Referee report, 2026-11-24. Nothing in the repository was modified. All scratch scripts and their outputs are in
`proofs/audit1/` (added to the repository afterwards; outputs carry the suffix `.out.txt`):

| file | what it does |
|---|---|
| `t11_partA.py` | Part A of T11 in exact `Fraction` arithmetic, independent of the author's sympy code |
| `t11_B1.py` | (B1) exactly, with `r(k) = sigma(N,k)` from the theorem, all N = 1..1000 |
| `t11_B2345.py` (+ `.out.txt`) | (B2)–(B5): symbolic re-derivations, O-constants for N up to 10^5–10^6 |
| `t3_t4.py` (+ `.out.txt`), `extra.py` | Theorems 3 and 4 |
| `rerun_verify_*.out.txt/.diff` | re-run of the author's three scripts: all pass, output identical to the stored `.out.txt` |

**Overall verdict: no error and no real gap found in any of the three proofs.** Every item is either "correct as
written" or "correct, text needs a small fix". The fixes are listed per item; the most useful ones are collected at the end.

---

## 1. Theorem 11, Part B

### (B1) shape — correct, text needs three one-line additions

* Positivity of `a(N,k)`, `0 <= k < N`: **correct as written.** `4k-2N+9/2 = 1/2 + 2(2k-N+2)`, so when negative it is
  `1/2 - j` with `j` even and `Gamma > 0`; `C_N > 0`. Exact check: `a(n,k) > 0` for all `1 <= n <= 60`, `k < n`.
* Log-concavity on `I = {N/2-1 <= k <= N-1}`: **correct as written.** On `I`, `4k-2N+9/2 >= 1/2`; the displayed product form of
  `r(k)` is an identity (checked exactly), each factor is nonnegative and decreasing. Exact check: `r(k+1) < r(k)` for every
  consecutive pair in `I`, all `N <= 1000`, no failure.
* `r(k) >= 1` below `I`: **correct, but three assertions are stated without proof.**
  - `r(k) >= ((N-k-1)(N+2k)^2/((N+k+1)(N-2k)^2))^2` for `0 <= k <= N/2-2`: true. Reason (not in the text): the four
    numerators `4k+2N+1/2+i >= 2N+4k`, the four denominators satisfy `0 < 2N-4k-9/2-i <= 2N-4k`, and `(2k+3)/(2k+1) >= 1`.
  - `6N^2 k - 8k^3 - 2N^2 - 8k^2 >= 0` for `1 <= k <= N/2`: true but only asserted. Proof: `8k^3 <= 2N^2 k` and
    `8k^2 <= 2N^2`, so the polynomial is `>= 4N^2 (k-1) >= 0`. (The polynomial identity itself is correct.)
  - `k = 0`: `r(0) >= 3((N-1)/(N+1))^2 >= 1` needs `N >= 2+sqrt3`, i.e. `N >= 4` — which is exactly when `k = 0 <= N/2-2`. Fine.
  - Odd `N`, `k = (N-3)/2`: the four numbers are indeed `-3/2,-1/2,1/2,3/2` (product `9/16 > 0`), and this is the only index
    covered by neither case. "`r(k) >= 1` too" is asserted without proof; proof:
    `r = (N/(N-2)) ((N+1)/(3N-1))^2 (16/9) prod_{i<4}(4N-11/2+i) >= (16/81)(4N-11/2)^4 >= 1` for `N >= 3`.
  - Exact check, all `N <= 1000`: `r(k) >= 1` for **every** `k < min I`; the stated lower bound, the polynomial inequality, the
    `k=0` case and the odd-`N` index all hold with no failure. (Smallest margin `r(k)-1` is 0.416 at N=400, 0.245 at N=1000;
    it tends to 0 like `N^{-1/2}`, so numerics alone would not prove it, but the argument above does.)

### (B2) ratio expansion — correct, text needs one fix

* `L(x*) = 0`, `L'(x) = 12/(4x^4-5x^2+1)`, `L'(x*) = -24`, `phi' = L`: confirmed symbolically. `sup |L''| ~ 388` on `[0.8,0.93]`.
* Uniformity `log r(k) = L(k/N) + O(1/N)` on `[0.8,0.93]`: `N |log r(k) - L(k/N)| <= 43.2, 41.3, 41.11, 41.09, 41.09` for
  `N = 10^2 ... 10^6`. Constant in `log r(k0+j) = -24j/N + O(1/N + j^2/N^2)`: `<= 65.2` (N=10^3), `47.0`, `42.3`, `41.5` (N=10^6).
* Summed form: the constant in `O(|j|/N + |j|^3/N^2)` is `<= 20.3` and tends to `|L''(x*)|/6 = 13.86`; the worst error times
  `N^{1/5}` is `26.3, 18.0, 15.1, 14.7` for `N = 10^3 ... 10^6`. So `O(N^{-1/5})` is right — but note it is an honest
  `~14 N^{-1/5}`, i.e. about 1.5 at `N = 10^5` at the edge `|j| = J`; the statement is purely asymptotic. No error.
* **Fix:** (B4) uses the expansion at `j = -J-1` (`r(k0-J-1)`), outside the stated range `|j| <= J`. State (B2) for
  `|j| <= J+1` (or `2J`); the proof is unchanged.

### (B3) centre — correct as written
`|sum_Z f - int f| <= max f = 1` for unimodal `f`; the discarded part is `<= 2 int_J^oo`, fine. Numerically
`centre/(a(k0) sqrt(pi N/12)) - 1 = 6.6e-3, 3.3e-4, 9.1e-6` for `N = 10^3, 10^4, 10^5`.

### (B4) tails — correct as written (for N large), modulo the range remark above

* Directions of the inequalities are right: `r` is decreasing on `I`, so `a(k) = a(k0+J) prod_{i=k0+J}^{k-1} r(i) <= a(k0+J) r(k0+J)^{k-k0-J}`,
  and `a(k) = a(k0-J)/prod_{i=k}^{k0-J-1} r(i) <= a(k0-J) r(k0-J-1)^{-(k0-J-k)}` because every `r(i)` there is `>= r(k0-J-1) > 1`.
* `r(k0+J) <= e^{-12J/N}` and `r(k0-J-1) >= e^{12J/N}`: hold for every `N` from 165 to 20000 and at `10^5, 10^6, 10^7`
  (the only exceptions below 20000 are small `N`, the last being 164, where `k0+J > N-2` and the statement is void; irrelevant).
* `sum_m e^{-12Jm/N} <= N/(6J)` uses `1-e^{-u} >= u/2`, valid for `u = 12J/N <= 1.59`, i.e. `N >= 151`. Fine.
* `a(min I) <= a(k0) e^{-4 N^{3/5}}`: follows from the chain as claimed (`12 (x*-1/2) = 4.39 > 4`); numerically true for
  all `14 <= N <= 3000` and at `10^4, 10^5` (it fails for a few small `N`, the last being `N = 13`; harmless). In fact `log(a(min I)/a(k0)) = N log(4/27) + O(1)`.
* Measured tails over `a(k0) sqrt N`: upper `2e-26, 1e-36, 2e-55`, lower `7e-22, 2e-34, 4e-54` for `N = 10^3, 10^4, 10^5`,
  all below the stated `N^{2/5} e^{-11 N^{1/5}}`.

### (B5) Stirling — correct as written

* Re-derived symbolically (`t11_B2345.py`): inserting `log Gamma(z+a) = (z+a-1/2) log z - z + (1/2) log 2pi` for the eight
  Gamma functions (including the four in `C_N`: `log C_N = -10N log2 + (15/2) log2 + (5/2) log N - (1/2) log pi + O(1/N)`)
  gives exactly `N phi(x) - (1/2) log N + log psi0(x)`; difference simplifies to 0. Also done by hand.
* `phi(x*) = 0` exactly: with `2 log(s+1) = log 2 + log(2+s)`, `2 log(s-1) = log 2 + log(2-s)` (`s = sqrt3`) the coefficient of
  `log 2` is `-6+(s+1)+(2-s)+(2+s)-(s-1) = 0` and what remains is `-log(2+s) - log(2-s) = -log 1 = 0`. Correct.
* `psi0(x*) = sqrt3/(sqrt2 pi^{5/2})` and `psi0(x*) sqrt(pi/12) = sqrt2/(4 pi^2)`: confirmed (`x(1-x)/((1+x)(2x-1)^4) = sqrt3/8` at `x*`).
* O(1/N) uniformly near `x*`: `N (log a(N,k) - [N phi - (1/2)log N + log psi0])` at `x = 0.8, 0.83, x*, 0.9, 0.93` is
  `-7.75, -7.32, -7.01, -6.99, -7.38` at `N = 10^5` (stable from `N = 100` on).
* `a(k0) sqrt N / psi0(x*) - 1 = -6.5e-3, -3.3e-4, -9.1e-6` for `N = 10^3, 10^4, 10^5`. This is not monotone in `N` because the
  `1/N` coefficient depends on `theta` (`~ -7 + 17.8 theta - 12 theta^2`); consistent with the claim.

### Combination (B3)+(B4)+(B5) — correct as written, no circularity
`T_N = a(k0) sqrt N [ sqrt(pi/12)(1+o(1)) + o(1) ]` and `a(k0) sqrt N -> psi0(x*)`. (B4) uses only (B1), (B2); (B5) is independent;
the value of the series is never used. `T_N/(sqrt2/(4pi^2)) - 1` is `< 1e-34` for `N = 10^3, 10^4, 10^5` (as it must be, since
by (A) `T_N` differs from the limit by `~64^{-N}`).

### Script vs. text for Part B — the text overstates what the script checks
* `verify_T11...py`, block B3: log-concavity is tested only for indices `k >= N//2+2` and monotonicity only for `k < N//2-2`
  (and the docstring says `I = [N/2+1, N-1]`). Four indices around `N/2` — including the first elements of `I`, the last index
  below `I`, and the odd-`N` exceptional index `(N-3)/2` — are **not** tested, and only for 5 values of `N`, in floating point.
  RESULTS.md says "log-concavity on `I`, monotonicity below". My exact test (`t11_B1.py`, all `N <= 1000`) closes this.
* RESULTS.md lists `phi(x*) = 0` among the items checked "exactly"; the script checks it numerically (`< 1e-55`). The hand
  proof above is exact.

---

## 2. Theorem 11, Part A — correct as written

* Telescoping over `1 <= n <= N-1`, `0 <= k <= N`: left side `sum_k (a(N,k) - a(1,k)) = T_N - a(1,0)`, right side
  `sum_n (gam(n,N+1) - gam(n,0)) = -sum_{n=1}^{N-1} gam(n,0)` since `N+1 > n`. With `a(1,0) = t(0)/256` (`C_1 = 9/256`, by hand)
  and `gam(n,0) = -t(n)/256` this is (A). Correct, including `N = 1`.
* Four cases of the WZ equation: `k <= n-2` — division by `a(n,k) > 0`, no denominator of `sigma, rho, R(n,k), R(n,k+1)` vanishes;
  `k = n-1` — needs `rho(n,n-1) - 1 = -rho(n,n-1) sigma(n+1,n-1) - R(n,n-1)`, which is exactly what the script checks;
  `k = n` — `a(n+1,n) - 0 = 0 - gam(n,n)`; `k > n` — all zero. Correct.
* The script does check what the text claims (both rational identities symbolically; `gam(n,0)` by `n = 1..8` plus a ratio
  identity, which is a valid induction). Only weakness: `sigma`/`rho` vs. the Gamma closed form is tested for `n <= 4` only.
* Independent exact recomputation (rational `C_n = 2^{6-8n} (3/4)_n^2 (1/4)_n / (-4 (-1/4)_n (2n-2)!)`): ratio formulas for all
  `n <= 40`; WZ equation at all `1 <= n <= 40`, `0 <= k <= 45` (all four cases), no failure; `gam(n,0) = -t(n)/256` for `n <= 40`;
  **(A) for all `N = 1..40`**. Values `sum_{n<N} t(n) = 256 T_N`:

  | N | exact | decimal |
  |---|---|---|
  | 1 | 9 | 9.000000000000 |
  | 2 | 1201725/131072 | 9.168434143066 |
  | 3 | 5041545945975/549755813888 | 9.170518653218 |
  | 4 | 660807502747917525/72057594037927936 | 9.170546304947 |
  | 5 | 88692085354689833196255075/9671406556917033397649408 | 9.170546686537 |
  | 6 | 11625049018450350423933547694139/1267650600228229401496703205376 | 9.170546691933 |

  (`64 sqrt2/pi^2 = 9.170546692012`.)

---

## 3. Theorem 3 — correct; three small text fixes

* Step 1: correct (zeros of `F(n,.)` are isolated, continuity suffices).
* Step 2: correct. `F' = O(x^{2n-1/2} (3/4)^x)`, `G' = O(x^{2n+3/2} (3/4)^x)`, so the stated `O(x^{2n+3}(3/4)^x)` is true
  (generous); the ratios are numerically bounded for `n = 0,1,3`, `k <= 2000`.
* Step 3: correct. **Fix:** the "exact ratio identity" covers `n -> n+1` for `n >= 1` only, so its base case is `n = 1`,
  which the script checks only numerically (`mp.diff`, 40 digits). By hand: `G'(1,0) = -7/(160 pi)`, `t(2) = 7/120`. OK.
* Step 4: correct as written. Beta integral, `c_0 = 3/(8 sqrt pi)`, the substitutions, `3cos^2 - sin^2 = 4cos(th+pi/6)cos(th-pi/6)`,
  the three pieces `-(1/2)Cl_2(2pi/3)`, `Cl_2(2pi/3)`, `Cl_2(pi/3)`, duplication, `Cl_2(2pi/3) = (sqrt3/2)K` all re-derived by
  hand and numerically (`I + 4K ~ 1e-39`; `D_0 + 3K/(2pi) ~ 1e-41`). The interchange is monotone convergence (all terms
  have the sign of `log(3t/4) < 0`).
* Step 5: correct. Closed form of `Psi_N` re-derived (Leibniz with `M = 2N`; the sign `(-1)^M = +1` is where evenness is used) and
  confirmed symbolically for `N <= 4` and at `y = 1/20, 3/10, 3/4, 99/100`; bound `(2N)!/((2N+1)y)` holds.
  **Fix (justification):** "dominated convergence" needs a dominating function; one is
  `sum_k |y^{2N} d^{2N}/dy^{2N}(y^k log y)| <= C_N (1 + |log y|)` on `0 < y <= 3/4` (for `k < 2N` the term is a constant times
  `y^k`; for `k >= 2N` it is `<= k^{2N} y^k (|log y| + 2N)`).
  **Fix (size):** `|D_N| = O(N 9^{-N})` is true but two powers of `N` too weak: the displayed bound equals
  `(sqrt3/6 + o(1)) 9^{-N}/N` (bound `* 9^N N = 0.2871, 0.2885, 0.28866` at `N = 50, 500, 5000`), and `|D_N| 9^N N ~ 0.13 ... 0.20`
  for `N = 1..8`. Telescoped identity `sum_{n<N} G'(n,0) = D_0 - D_N` verified to 1e-41 for `N = 1..8`.

---

## 4. Theorem 4 — correct; one line should be added

* `1 - sigma - (3/4)(n+1)/(k+3n+1) = (4nk+4k+16n+1)/(16(k+n+1)(k+3n+1)) > 0` (exact). Then
  `h(n,k)/h(n,0) <= prod_{j<k} exp(-c/(j+m)) <= (m/(k+m))^c` (`c = 3(n+1)/4`, `m = 3n+1`) by `sum_{j<k} 1/(j+m) >= log((k+m)/m)`.
  Correct; numerically the ratio never exceeds 1 (`n < 12`, `k < 4000`).
* `sum_k h(n,k) <= 7 h(n,0)` for `n >= 2`: true, asserted without proof. **Add:** `sum_{k>=0} (m/(k+m))^c <= 1 + m/(c-1) = 1 + 4(3n+1)/(3n-1) <= 33/5`
  for `n >= 2` (it is 9 at `n = 1`, which is why `n >= 2`). True values of `sum_k h(n,k)/h(n,0)`: 2.95, 3.06, 3.84, 3.96 for `n = 1, 2, 40, 200`.
* `K (1/4)_K (3/4)_K/K!^2 -> 1/(Gamma(1/4)Gamma(3/4)) = 1/(sqrt2 pi)`: correct (`(a)_K/K! ~ K^{a-1}/Gamma(a)`); `Q(0,K)/K -> 1`.
* Telescoping with `u = n h`, `u(0,k) = 0`, `K -> oo` (needs `Q(n,K)h(n,K) -> 0` for `1 <= n < N`, exponent `>= 3/2 > 1`, and convergence
  of `sum_k h(N,k)` for `N >= 1`): correct. Checked: `sum_{n<N} Q(n,0)h(n,0) + N sum_k h(N,k) = 1/(sqrt2 pi)` to 1e-41 for `N = 1..5`.
  The bound `7N(27/32)^N` relies on the "7", proved for `N >= 2` only (the sentence reads as if for `N >= 1`; it is numerically true there too); irrelevant for the limit.

---

## Section 0
Telescoping lemma and the equivalence with (WZ) (for `f != 0`): correct as written.

## Suggested edits to RESULTS.md, most useful first

1. T11 "Check" paragraph: say that log-concavity/monotonicity are tested by the script only away from `k ~ N/2` (or fix the
   script's index ranges `kI = N//2+2`, `range(0, N//2-2)` and the docstring `I = [N/2+1, N-1]`), and that `phi(x*) = 0` is
   checked numerically there (the exact argument is the one in the text).
2. T11 (B2): state the expansion for `|j| <= J+1`, since (B4) uses `r(k0-J-1)`.
3. T11 (B1): add the three one-line proofs given above (lower bound for the Gamma quotient, the polynomial inequality, the odd-`N` index).
4. T4: add `sum_k (m/(k+m))^c <= 1 + m/(c-1) <= 33/5`.
5. T3 step 5: name the dominating function; replace `O(N 9^{-N})` by `O(9^{-N}/N)` (optional); step 3: note the exact base case `n = 1`.
