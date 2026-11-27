# Referee report on `monotonicity/hkmn23/HKMN23.md`

Date: 2026-11-24. Scope: every claim of §§0–6 of HKMN23.md as it stood at 08:29 (the file
still contains two placeholders, see R1). Every identity, residue, sign and inequality was
re-derived by hand before reading the note's proof. All numerics below were done with my own
scripts in
`/tmp/claude-0/-home-user-weighted-prime-races/603595b8-b464-579a-a8f3-65c6b0cc4976/scratchpad/referee/`
(`lemA.py`, `fin.py`, `chars.py`, `phi_check.py`, `zeros.py`, `uchk.py`, `ef.py`, `lb.py`,
`t1.py`, `wit.py`). The author's scripts were read line by line but not run; only their output
files were read. One core, each run under 3 minutes. No existing file was edited.

**Verdict in one line.** I could not break it: no FATAL error. Lemmas A, E, I, S, N and
Theorems 1–5 are correct as far as I can re-derive and re-compute them. There are three
REPAIRABLE items (the note is unfinished in two places and Theorem 6(b) is not yet a fixed
statement; two small gaps in §4.4; unverifiable priority/comparison claims) and a list of
MINOR points, including one wrong number in a remark (the low zero mod 19).

---------------------------------------------------------------------------

## 1. Lemma A (§2.1) — correct

Re-derived from scratch. With ν = φ Σ log p·δ_p − dt on [2,x] and W(t) = ν([2,t]):
`∫_{[2,x]} F dν = F(x)W(x) − ∫_2^x F'W` (the possible atom at t = 2 is inside the closed
interval, W(2^-) = 0), and `−∫F'W = −F'(x)W_1(x) + ∫W_1F''` because W_1(2) = 0. For
h_x = (x−t)/log t I get h_x(x) = 0, h_x'(x) = −1/log x and
`h_x'' = 2/(t log²t) + (x−t)(1/(t² log²t) + 2/(t² log³t)) = (x+t)/(t² log²t) + 2(x−t)/(t² log³t)`,
positive on (1,x]; this matches. The left side is `φΣ_{p≤x}(x−p) − ∫_2^x (li(t) − li(2))dt`
(Fubini), so the constant is `−li(2)(x−2)` with the helpful sign. The reduction of q = 1 to
(2,1) is correct: `∫_2^x(π − li) = A(x;2,1) + (x−2)` and li(2) − 1 = 0.04516 > 0.

Numerics (30 digits, piecewise quadrature, `lemA.py`), left − right:

    (12,1), x = 300.7 :  −4029.94660382879172681687597…  both sides, diff 1.6e-27
    (20,9), x = 123.4 :  diff 2.0e-28
    (3,2),  x = 77.7  :  diff 1.8e-28     (class containing p = 2: tests the atom at the base point)
    (1,·),  x = 50.5  :  diff 0;   and ∫_2^x(π−li) = A(x;2,1) + (x−2) to all 30 digits

The trivial reductions (6,1) = (3,1), (10,1) = (5,1), (10,9) = (5,4), (30,1) = (15,1),
(30,19) = (15,4) are correct (in none of them is 2 ≡ a mod m).

## 2. Lemma E, the functional-equation constants, Lemma I (§2.2) — correct

Residues of `(−L'/L)(s,χ) x^{s+1}/(s(s+1))`, re-derived:

- at ρ: `−x^{ρ+1}/(ρ(ρ+1))`;
- odd χ: at 0, `−bx`; at −1 (double pole, w = s+1:
  `(1/w²)(1 + b_0 w)(1 + w(log x + 1))`), `log x + 1 + b_0`; at −(2n+1), `−x^{−2n}/(2n(2n+1))`;
- even χ: at 0 (double pole, `−(x/s²)(1 + bs)(1 + s(log x − 1))`), `−x log x + x − bx`;
  at −1 (pole of the kernel only, residue of 1/(s(s+1)) is −1), `+(L'/L)(−1,χ)`;
  at −2r, `−x^{1−2r}/(2r(2r−1))`;
- ζ: `x²/2`, `−x log 2π`, `+(ζ'/ζ)(−1)`, `−T_e`.

All match the note. Closed forms of T_e, T_o: correct (series checked term by term);
T_e(1) = log 2, T_o(1) = 1 − log 2; x·T_e(x) = 1/2 + u²/12 + … ≤ 0.51 and
x²·T_o(x) = 1/6 + u²/20 + … ≤ 0.173 at x = 3, so the bounds 0.6/x and 0.2/x² hold.

Functional equation: from `(L'/L)(s,χ) = −(L'/L)(1−s,χ̄) − log(f/π) − ½ψ((s+a)/2) − ½ψ((1−s+a)/2)`
I get, for both parities, `b = −(L'/L)(1,χ̄) − log(f/π) + γ_E + log 2` and
`b_0 = c_{−1} = −(L'/L)(2,χ̄) − log(f/π) + γ_E + log 2 − 1`. Matches.

Lemma I: correct, including the sign (`φψ_1(x;q,a) = Σχ̄(a)ψ_1(x,χ*) − imp`, imp ≥ 0) and
the orthogonality step (χ*(p) ≠ 0 ⟺ f_χ | m_p; the characters mod q of conductor dividing m_p
are exactly those induced from mod m_p). `imp ≤ I_q x log x` is correct.
Formula (E), the integers n_e = (φ/2)·1_{a≡±1} − 1, n_o = (φ/2)(1_{a≡1} − 1_{a≡−1}), and B, C:
correct. All twelve (n_e, n_o) pairs of the table in §2.4 agree with these formulas.

Numerics (`ef.py`: zeros with |γ| < 150 located by my own scan of Hardy's Z; constants b and
b_0/c_{−1} computed by me from the functional equation with mpmath, i.e. by a route different
from the note's Taylor coefficients at 0 and −1). "direct prime-power sum − formula":

    χ = 5.4 (even, real)      x = 2.5: 3.9e-5   x = 7.3: 5.6e-4   x = 20.6: 5.9e-5   x = 1: 1.23e-2 (= omitted tail 1.23e-2)
    χ = 7.2 (even, cubic)     x = 7.3: 6.9e-4   x = 20.6: 9.5e-3                     x = 1: 1.30e-2 (= tail)
    χ = 5.2 (odd, quartic)    x = 7.3: 1.8e-3   x = 20.6: 2.6e-3                     x = 1: 1.23e-2 (= tail)

A wrong residue or constant would show as O(1) or O(x). The even-character double pole is
therefore confirmed numerically as well.

## 3. Lemma S (§2.3) — correct

`|ρ|² = A`, `|ρ+1|² = A + 2` for ρ = 1/2 + iγ; `∫_0^π dθ/(c − cos θ) = π(c²−1)^{−1/2}` with
c = A + 1 gives `(A(A+2))^{−1/2}`; `A + 1 − cos θ = α(θ)² + γ²`. V_χ(α) is Davenport (17)–(18)
with Re B(χ) = −Σ Re 1/ρ; for ζ, `½ψ(s/2+1) = ½ψ(s/2) + 1/s` gives the stated form.
Trapezoid: the Poisson-kernel expansion gives `a_0 = S`, `0 ≤ a_k ≤ 2^{1−k}S` (r ≤ 1/2 because
c ≥ 5/4), the n-interval rule on [0,π] aliases exactly k ∈ 2nℤ, so
`S ≤ T_n ≤ S(1 + 2/(4^n − 1))`, i.e. `T_n(4^n−1)/(4^n+1) ≤ S ≤ T_n`. Correct. The sandwich
`1/(A+1) ≤ (A(A+2))^{−1/2} ≤ ½(1/A + 1/(A+2))` is correct.

Numerics. (a) Same formula, different software (`chars.py`, mpmath, my own character
construction): S(ζ) = 0.046117644421509, S(4.3) = 0.15379470538051,
S(5.2) = 0.19897259073964, S(5.4) = 0.15517198067265, L_1(20) = 2.28170575077 — all digits
agree with `l1const.out`. (b) Different route (`zeros.py`: actual zeros with |γ| < 120 plus a
smooth-density tail): S(4.3) = 0.139664 + 0.014153 = 0.153817 (closed form 0.153795),
S(5.4) = 0.140447 + 0.014746 = 0.155193 (closed form 0.155172). The differences 2·10⁻⁵ are
the size of the neglected fluctuation of N(T) at T = 120.

## 4. Theorem 2 (§2.4) and the scripts — correct

Every step re-derived.

- Step 1: `ψ(u;q,a) ≥ θ(u;q,a) + Σ_{b∈R}θ(u^{1/2};q,b)`; correct, terms k ≥ 3 dropped with the
  right sign (they make ψ_1 − θ_1 larger).
- Step 2: imp ≥ 0 dropped with the right sign; `−(1+n_e)T_e ≤ 0` dropped; `−n_oT_o ≤ max(0,−n_o)·0.2/t²`. Correct.
- Step 3: `∫_0^{√t}(Ny − c_1√y)2y dy = (2/3)Nt^{3/2} − (4/5)c_1t^{5/4}`. Correct. The reduction of
  the supremum to left limits at primes of R plus the endpoint Y_1 is correct (Θ ≥ 0).
- Step 4: `Sq = 2YΘ_1(Y) − 2∫_0^YΘ_1`. Upper bound for Θ_1: correct, and E_R is exactly what
  one gets term by term for y ≥ 3. Lower bound: per class b, the powers k ≥ 2 contribute at
  most `φ log 4 Σ_{k=2}^{⌊log_2 Y⌋} Y^{1+1/k}/(1+1/k) ≤ φ log 4[(2/3)Y^{3/2} + (3/4)Y^{4/3} + (log Y/log 2)Y^{5/4}]`
  (at most log_2Y − 3 terms with k ≥ 4, each ≤ Y^{5/4}); this gives λ(Y) as printed (the
  factor N comes from bounding each class separately; crude but valid). Assembling, I get
  `K_B(t) = 2λ(√t) + (4/5)NL_1 + E_R t^{−1/4}(log t)/2 + φ log 2·t^{−5/4}`; the note has
  `2φ log 2·t^{−5/4}`, i.e. a harmless factor 2 of slack in the last term. K_B is decreasing
  for t ≥ e⁸.
- Step 5: `−(t−2)²/2 = −t²/2 + 2t − 2`; r_+ is what one gets; every term decreases for t ≥ e² < 8.

`w1_analytic.py` implements exactly these formulas (read line by line: `pair_constants`,
`rplus`, `KB`, E_R, I_q, the regime switch at Y_1², the Arb assertions `UA(X0) < 0`,
`UB(Y1²) < 0`). `w1_cert.py`: M_i = (φϑ_i)²/2 + 2φϑ_i − φΛ_i is the maximum over all real t
of the concave quadratic piece (at t* = φϑ_i + 2); the test is `M < 0` on Arb balls (true only
if certain); the last piece runs to the first prime ≡ a beyond 10⁹, which is what
`valid_up_to` records; part B takes left limits and the endpoint, and inflates the float by
1 + 10⁻¹⁵. All as §2.5 says.

Independent re-evaluation (`uchk.py`, mpmath; B, C, E_R from my functional-equation constants,
not from the note's JSON):

    (20,1): L_1 = 2.28170575077, m = 0.3849609159, n_e = 3, n_o = 4, B = 1.1433438, C = 0.3024841, E_R = 9.280959, I_q = 6
            U(8.16e5) = +3.7e-6,  U(8.17e5) = −1.15e-4,  U(1e6) = −0.019241,  U(1e12−) = −0.373492,
            K_B(1e12) = 148.32694,  U(1e12+) = −0.236631,  K_B(∞) = 84.703664
    (2,1):  B = −1.8378771, C = 1.9850537, E_R = 2.3358817;  U(58.4) = +5.3e-5, U(58.5) = −2.2e-4, U(1e12+) = −0.616491
    (5,4):  B = 0.4797282, C = −0.0240460, E_R = 8.9457776;  U(689) = −4.8e-5,  U(1e12+) = −0.699883
    (15,4): B = −3.8324662, C = 3.7149050, E_R = 17.891555;  U(1.78e4) = −1.3e-4, U(1e12+) = −0.596209

All agree with Table 1, the table of §2.4 and `w1_analytic.out` to the digits printed. X_0 is
sharp to three digits in the two cases where I tested the neighbour.

Finite range (`fin.py`, numpy longdouble, primes ≤ 10⁶; not ball arithmetic, margins are
≥ 0.49 in absolute terms):

    (20,1): 9758 primes ≡ 1 (20) up to 1e6, all M_i < 0, max M_i = −717.33 (p = 41); decade sups identical to the JSON to 10 digits;
            c_1 = 14.332099223989 attained before 379; value at Y_1: 9.366793911          (JSON: 14.332099223989076, 379, 9.36679391143)
    (2,1):  78497 primes, max M_i = −0.49514 (p = 3); c_1 = 2.071192618039 before 1423      (JSON: 2.071192618038794, 1423)
    (20,9), (15,1), (5,1): same agreement; c_1 = 9.103522434555 (37), 14.332099223989 (379), 7.166049611995 (379)

c_1(15,1) = c_1(20,1) = 2c_1(5,1) is not a coincidence: in all three cases R is "p ≡ ±1 mod 5".
Since X_0 ≤ 8.17·10⁵ < 10⁶ in every case, my recomputation to 10⁶ already covers the range the
proof of Theorem 1 needs for these five pairs; the certificate to 10⁹ is 1200 times more than
needed in the worst case.

## 5. Corollary 1.3 / §3 — correct

`d/dx[x li(x) − li(x²)] = li(x)`; both sides → 0 at 0+; continuity at 1 (limit −log 2).
g' = −li, so g increases on (0,μ) and decreases after; unique zero x_0. For q ≥ 2 the first
prime ≡ a is ≥ 3 > x_0 in all 17 cases, and `c_0 − li(2)(3−2) = −0.168 < 0`. For q = 1,
π(p) < li(p) at the primes ≤ 31 implies A_0 decreasing on [2,31] (li is increasing, π constant
between primes), and 21.42 < 31. Numerics (`lemA.py`, 30 digits): c_0 = 0.87725753480406530832,
x_0 = 2.61113944566016890519, x_1 = 3.56517637359934422071, g(3) = −0.76952803, A_0(3) = 0.2304720,
A_0(5⁻) = −1.6599851, li(p) − π(p) > 0.045 at all primes ≤ 31. All as stated. In particular
the remark on the problem list is right: with lower limit 0, `∫_0^x(π − li) > 0` on (2, 3.5652).

## 6. §4

**Theorem 3.** Correct. F is real on all of ℝ (zeros of χ̄* are the conjugates), Σ|c_j| = L_1,
the Dirichlet argument gives arbitrarily large almost periods. The asymptotic formula: P_{≥3}
bound correct; `h_x'' ≤ (2 + 2/log 2)x/(t² log²t) = 4.89 x/(t² log²t)`; the upper bound for Sq
"with roles exchanged" works. One point suffices, and y_0 may be negative.

**Lemma N.** Correct (second derivative 1 + (3/4)c x^{−1/2} + ℓ''(x) → 1 against a piecewise
linear left side).

**Theorem 4.** Correct. The Bohr mean is the coefficient at frequency 0,
`−(4/3)Σχ̄(a)m_{χ*}(1/2)` (1/(ρ(ρ+1)) = 4/3 at ρ = 1/2); a real almost periodic function of
mean 0 that is not identically 0 takes both signs. NV(q) for imprimitive χ is equivalent to
NV for χ* (the Euler factors do not vanish at 1/2).

**Theorem 5(a).** Correct, both parities and ζ: `c(ζ) = log 2π + log 2 − 1/2 − (ζ'/ζ)(−1) = 0.0459705226`;
I also get it as `Σ1/(ρ+1) − Σ1/ρ = 0.069066 − 0.023096`. Zero-side check (`zeros.py`):
c(4.3) = 0.136178 + 0.014153 = 0.150331 against the closed form 0.150309; c(5.4) = 0.152459
against 0.152439. Independent values of Φ (`phi_check.py`, mpmath):
Φ(3) = 0.1571905, Φ(4) = 0.1962790, Φ(7) = 1.2542288, Φ(16) = 2.3257355, Φ(20) = 2.0541722,
**Φ(21) = 3.6294975, Φ(60) = 5.4858331** (thresholds 2.6667 and 5.3333). They agree with
`t2_summary.out` to the ~9 digits my method gives (mpmath's `dirichlet` is unreliable at s = 1
for complex characters; I evaluate at s = 1 + 10⁻¹⁰). Also, for q = 16, 27, 32, 48, 60, 96,
120 the prime-side value `F_at_1` in `raw/t2_witness_120_*.jsonl` (Taylor constants at 0 and
−1) equals the closed form assembled from `raw/t2_closed_400_*.jsonl` (L'/L at 1 and 2) to
15 digits; this tests the corrected identification of inducing characters at exactly the
moduli where the Conrey-label bug of §6.4 occurred.

**Theorem 5(b), q > 400 (§4.4).** Re-derived step by step.

- (i) Correct. `Σ_{j≥1}log(j+1)/j² = 1.80074…`; `p^{2o_p} − 1 ≥ m_p(m_p+2)`; each p | q costs ≤ (log p)/3.
- (ii) Correct: `s∫_1^∞φψ x^{−s−1} = s(s+1)/(2(s−1)) + s(s+1)∫E_1x^{−s−2}`,
  `s(s+1)/(2(s−1)) = 1/(s−1) + 3/2 + O(s−1)`, `ζ'/ζ = −1/(s−1) + γ + …`; hence
  `−Σ_{χ≠χ_0}(L'/L)(1,χ) = 1 + γ + Σ log p/(p−1) + I`. The passage to primitive characters and
  the identity `T_1 = 1 + γ + I + Σφ(m_p)log p/(p^{o_p} − 1)` are correct. Numerics (`t1.py`,
  my own code, exact truncated integral at X = 4·10⁶): q = 7: −1.510378 vs −1.510292;
  q = 12: −1.090641 vs −1.090582; q = 20: −1.685065 vs −1.685157 (truncation error O(L_1X^{−1/2}) ≈ 5·10⁻³).
  The formula for T_2 in (i) agrees to 6 digits for q = 7, 12, 20.
- (iii) R(q) is exactly what the six integrals give. `V(1/2) ≤ 3V(3/2)` is right under GRH.
  The bounds |Re(L'/L)(1,χ*)| ≤ log q + 1.62, L_1 ≤ φ(3 log q + 0.11), |B| ≤ φ(2 log q + 4.4),
  |C| ≤ φ(log q + 2.4) are all true for q ≥ 3 (I get the sharper |Re b| ≤ 2 log q + 2.94,
  |b_0|, |c_{−1}| ≤ log q + 0.9), but see R2 for the ζ term.
- (iv) The constant is exactly 1: `log 2π − (ζ'/ζ)(−1) − (ζ'/ζ)(2) + γ = 1` follows from the
  functional equation at s = 2 (`ζ'/ζ(2) + ζ'/ζ(−1) = log 2π + γ − 1`); numerically
  1.837877 − 1.985054 + 0.569961 + 0.577216 = 1.000000.
- 400 < q ≤ 10⁶ (`lb.py`, my own sieve): min of LB − (2/3)N = 18.6056 at q = 420, ratio 2.744,
  max R = 3.94·10⁻⁴, max of the two small terms 0.0291. Identical to `t2_summary.out`.
- `N(q,1) ≤ φ(q)/4` for q > 168: correct (the exceptions are exactly the 30 moduli
  1,…,10, 12, 14, 15, 16, 18, 20, 21, 24, 28, 30, 36, 40, 42, 48, 56, 60, 72, 84, 120, 168);
  φ(q) ≥ (q/2)^{1/2} is classical; g(707) = 208.67, g' > 0.17 for φ ≥ 20. The final inequality
  LB(q) > φ/6 ≥ (2/3)N holds.

**Theorem 6.** Counts re-derived: 160 classes, 44 squares, 116 non-squares; 18 + 121 + 21 = 160;
the list of 21 undecided squares is the complement of the right sets. Spot checks (`wit.py`,
my own evaluation of (F)): F_{4,3}(log 7381) = 0.0831465174503 and F_{21,13}(log 30.8149…) =
1.2247452 agree with the certified values. The zero-side method (`t2_zeros.py`, read only) is
sound: missed zeros go into the remainder, bracket widths are carried as balls, and the rows
of the table in §4.5 match `raw/t2_zeros_{21,40,48,60}.jsonl`; `Φ − F_Z(≈0)` = 0.0664, 0.0902
is just below the remainders 0.0699, 0.0968, as it must be.

## 7. Statement hygiene

Hypotheses are stated where they are used; GRH(q) includes RH and this is said. Nothing in
§1.1 is stronger than what §§2–3 prove. In §1.2, Theorem 6(b) is not a finished statement
(R1). The comparison claims are R3.

---------------------------------------------------------------------------

## Issues

### FATAL

None found.

### REPAIRABLE

**R1. The note is unfinished, and Theorem 6(b) is not a fixed statement.**
Line 9 is `SUMMARY_PLACEHOLDER`; line 654 (§4.5, "zero side at the other moduli") ends in
`T6Q_PLACEHOLDER`. Theorem 6(b) defines its exceptional set as "the square classes listed in
`t2_summary.out` under undecided (§4.5 gives the list at the time of writing)", but §4.5 gives
no list, and `t2_summary.out` (08:20) is stale: it reports "zero-side Z = 0" and lists
(11,3), (11,4), (11,5) as undecided although `raw/t2_zeros_11.jsonl` (08:29) contains
certified witnesses for them; `t2_zeros.py` jobs for q = 19, 39, 44 were still running while I
read, and `raw/t2_zeros_{19,28,32}.jsonl` are empty. §7 lists a `SHA256SUMS` file that does
not exist. *Correction:* finish the runs, re-run `t2_assemble.py`, print the final undecided
list in §4.5 and state 6(b) with that explicit list (or drop 6(b)); write the summary; create
or remove `SHA256SUMS`.

**R2. Two small gaps in §4.4 (Theorem 5(b), q > 400). The conclusion stands.**
(a) "`0 ≤ V_χ(1/2) ≤ 3V_χ(3/2) ≤ 3[0.59 + ½log(f/π)]`" is false for f = 1:
V_ζ(3/2) = 0.0691 > 0.59 − ½ log π = 0.0176. The consequence `L_1(q) ≤ φ(3 log q + 0.11)` is
still true, because the ζ term is 2V_ζ(1/2) = 0.0462 and the other φ − 1 terms are each
≤ 3 log q + 0.106; say so. (b) For q > 10⁶, "the last two terms of LB are < 0.04" is asserted,
not shown; R(q) contains ω(q). One line suffices (ω(q) ≤ log_2 q and φ ≥ (q/2)^{1/2} give
R(q) + (φ/q²)(…) < 10⁻⁴; my grid evaluation gives sup ≈ 7.9·10⁻⁵). Also state that the bounds
for |B|, |C| are derived from the functional-equation formulas of §2.2 with
|(L'/L)(2,χ)| ≤ −(ζ'/ζ)(2) = 0.570; the derivation is currently one clause.
(c) "certified for q ≤ 400" in Theorem 5(b): the Arb part is the blocks c_f; the final
assembly `Phi(q)` in `t2_assemble.py` is floating point (ball midpoints and radii added as
floats, c(ζ) hard-coded as a float). The smallest margin is 0.15, so nothing is at risk, but
the label should say this, or the sum should be done in Arb.

**R3. Claims about the sources that this referee could not verify, and that carry the novelty
statement.** §0 and §1.3 assert: Johnston needs x ≥ 10¹⁹ / Büthe's range; "no threshold is
given or follows from" LWS; LWS's numbers of excluded classes (67, 49, 75); "the remark of LWS
§8.5 is too pessimistic"; "a correction of detail to LWS Proposition 8.3"; "every value they
print is reproduced". I did not fetch the papers. LWS is an unrefereed Zenodo preprint; the
note says so in §6.1 but §1.3 is written as if its content were settled. *Correction:* give
page/equation references for each of these, and phrase the two critical remarks on LWS as
observations ("their proof as written uses A(q,a) ≠ 0; Lemma N removes this"). The
mathematical statements of the note do not depend on any of this.

### MINOR

M1. §4.2, Remark. The low zero mod 19 is at **γ = −0.018956** for the Conrey character 19.3
(equivalently +0.018956 for 19.13 = its conjugate), not "γ ≈ 0.015". |L(1/2,χ)| = 0.019140 is
right (the note seems to have mixed the two numbers). "Essentially constant for |y| ≤ 100" is
then wrong as stated (the phase reaches 1.9 rad); it is true on the scan range y ≤ log 10⁷ =
16.1 (phase ≤ 0.31). A sign flip needs |γy| ≈ π, i.e. |y| ≈ 166, not 200.

M2. GRH(q) (§0) is defined as "every nontrivial zero of every L(s,χ), χ mod q, has real part
1/2". For imprimitive χ the Euler factors give zeros on Re s = 0. Say "zeros in 0 < Re s < 1"
or define it by the primitive characters (the "Equivalently" sentence is the correct one).

M3. Corollary 1.4: "below that W_1(t;2,1) < 0 is a check at 15 primes". W_1 is concave between
primes, so values at primes do not suffice; what is checked (correctly, by `w1_cert.py`) is
M_i < 0 on the 15 quadratic pieces starting at 3, …, 53. Also "it uses the one computed
constant c_1" understates the inputs: L_1(1) = 0.04612 and the explicit B, C are used too
(for q = 2 these are classical constants).

M4. "Rigor class" under Table 1: besides GRH(q), Chebyshev and c_1, the analytic range uses
the Arb values of L_1(q), B, C, E_R (L'/L on [1,2], Taylor coefficients at 0 and −1). List them.

M5. Table 1, column "sup G [1e8,1e9]", and the "consistency test" of §2.4: the per-decade
number is max of M_i/p_{i+1}^{3/2} over pieces *starting* at primes p_i ∈ [10^k, 10^{k+1}).
It is an upper bound for sup W_1/t^{3/2} on [first prime ≥ 10^k, first prime ≥ 10^{k+1}], not
on [10^k, 10^{k+1}); and in the low decades it is far from "what the truth looks like"
(e.g. (20,1), decade 2: bound −0.74, true sup −1.76). Harmless; reword.

M6. `w1_analytic.py`, function `B2`: balls are rebuilt from a 40-digit decimal midpoint and a
radius printed to 5 significant digits (rounded to nearest, not up). This is not strictly
outward; the possible loss is ~10⁻⁵ of a radius of ~10⁻²⁶, against margins ≥ 10⁻⁴. Print the
radius rounded up, or pass Arb strings.

M7. §2.4 Step 4: the last term of K_B can be φ log 2·t^{−5/4} (the note has twice that; safe).

M8. `t2_zeros.py`: brackets are constrained to [t − 10⁻¹², t₂ + 10⁻¹²], so two brackets from
adjacent scan cells may overlap by 2·10⁻¹²; "disjoint" is not enforced. With the certified end
signs the two brackets still contain two distinct zeros, so no double counting can occur, but
the docstring claim should match the code (clip to the cell).

M9. Stale cross-references in docstrings: `t2_assemble.py` ("Theorems 4–7"), `t2_witness.py`
and `t2_zeros.py` ("Theorem 4" for the witness principle, which is Theorem 3), `w1_analytic.py`
("H(q)").

M10. §4.3, last paragraph: "Φ(q)/L_1(q) is 0.98–0.99 for q ≤ 4": it is 0.9968 for q = 1, 2.

M11. The Hardy–Wright theorem number for θ(x) < x log 4 was not checked by me; it is stated
there for integers n, which suffices.

---------------------------------------------------------------------------

## Statement-by-statement verdict

    Lemma A                 proved as written (re-derived; 4 numerical checks to 27–30 digits)
    Lemma E                 proved as written (residues re-derived; even/odd/complex checked against zeros)
    functional-eq. formulas proved as written
    Lemma I                 proved as written
    Lemma S                 proved as written (incl. the trapezoid enclosure)
    Theorem 2               proved as written, up to the certified constants L_1, B, C, c_1, which I reproduced
                            independently for (20,1), (2,1), (5,4), (15,4) (and c_1 also for (20,9), (15,1), (5,1))
    Theorem 1               proved as written, conditional on GRH(q) and on the Arb computations; for (20,1), (2,1),
                            (20,9), (15,1), (5,1) I reproduced the finite range up to 10⁶ ≥ X_0 in floating point
    Corollary 1.3           proved as written (constants reproduced to 20 digits)
    Corollary 1.4           proved up to M3 (wording); the comparison with Johnston's proof is R3
    Theorem 3               proved as written
    Lemma N                 proved as written
    Theorem 4               proved as written; NV(q ≤ 400) not re-verified (one character checked)
    Theorem 5(a)            proved as written
    Theorem 5(b)            q ≤ 400: proved up to R2(c) (7 values of Φ reproduced);
                            q > 400: proved up to R2(a),(b) (both repairable in a line; scan 400 < q ≤ 10⁶ reproduced)
    Theorem 5(c)            follows
    Theorem 6(a)            supported by the raw files; counts and lists consistent; 2 witnesses re-evaluated; not re-run
    Theorem 6(b)            not established as a statement (R1): the exceptional list is not fixed
    §1.3 novelty claims     not established by this report (R3)

## What I did NOT check

- The certificate `W_1 < 0` between 10⁶ and 10⁹ (not re-run; not needed for the proof since
  X_0 ≤ 8.17·10⁵), and the finite range for the seven cases (3,1), (4,1), (5,4), (8,1), (12,1),
  (15,4), (24,1) (their JSON files report 0 failures; c_1 values are consistent with the
  obvious relations c_1(4,1) = 2c_1(2,1), c_1(8,1) = 4c_1(2,1), c_1(24,1) = 2c_1(12,1)).
- B, C, E_R, X_0 for the eight pairs other than (20,1), (2,1), (5,4), (15,4).
- Central non-vanishing for the 29 468 characters (only 19.3: |L(1/2)| = 0.0191397, confirmed).
- Φ(q) for all q ≤ 400 (seven moduli recomputed; seven more cross-checked between the two raw
  data sets).
- The 4266 prime-side witnesses (two re-evaluated) and the zero-side certificates (code read,
  not run).
- Any quotation from HKMN, Johnston or Liu–Wang–Shi; `EXPLICIT.md`; "triage §6";
  `verify_formula.py` and `verify_identity.py` (outputs read only).
- Correctness of Arb / python-flint itself. My cross-checks used mpmath and numpy, except that
  the zero ordinates in `zeros.py`/`ef.py` come from sign changes of flint's `hardy_z`.
