# tau(p^4) is never 461 or -599

Date: 2026-11-26. Files: `tau_p4.gp` (PARI/GP generator), `certificate_tau_p4.txt` (its output),
`check_certificate.py` (independent Python checker of the certificate), `sieve_tau_p4.py` (the
first, Python, implementation), `thue_check.gp`, `thue_results.txt` (Thue equations, Section 4 only), `pair_table.py` (table in Section 4).

## 1. Statements

Source. A. Goto, "On some values which do not belong to the image of Ramanujan's tau-function",
arXiv:2405.16723v1 (26 May 2024). With `L_1^+ := {461}, L_1^- := {599}`:

> "Let l < 1000 be an odd prime and eps = +-1. [...] Theorem 0.3. If tau(n) = eps l, we must have
> n = p^4, where p is prime and l in L_1^eps."

and Section 2.3: "When d = 5, the equation (2.1) for a = 281, 461, -599, -919 cannot be ruled out by
the congruences. In this subsection, we will show that Dembner-Jain's method makes it possible ruling
out 281, -919." So the paper leaves open exactly: is tau(p^4) = 461, is tau(p^4) = -599, for a prime p.

**Theorem 1.** Neither of the equations

    x^22 - 3 x^11 y^2 + y^4 = 461,        x^22 - 3 x^11 y^2 + y^4 = -599

has a solution in integers x, y. Consequently tau(p^4) is not 461 and not -599 for every prime p.

**Corollary 2.** tau(n) is not in {461, -461, 599, -599} for every n >= 1. (Section 4; this uses the
literature for the reduction to prime powers, but not Goto's paper.)

## 2. Proof of Theorem 1

**Step 0 (Hecke).** For a prime p put t = tau(p). From tau(p^{j+1}) = t tau(p^j) - p^11 tau(p^{j-1}):
tau(p^2) = t^2 - p^11, tau(p^3) = t^3 - 2 p^11 t, tau(p^4) = t^4 - 3 p^11 t^2 + p^22. So
(U, V) = (p^11, t^2) is an integer solution of

    U^2 - 3UV + V^2 = alpha,   alpha = tau(p^4).                                   (*)

Theorem 1 is the statement that (*) with alpha in {461, -599} has no integer solution in which U is
an 11th power and V is a square. (The papers of Balakrishnan-Craig-Ono-Tsai, Dembner-Jain, Goto and
arXiv:2603.29970 pass from (*) to u^2 = 5x^22 + 4 alpha, u = 2V - 3U, which forgets that V is a
square. Keeping it is what settles 461; for -599 it is not needed, see the remark after Step 3.)

**Step 1 (all solutions lie on finitely many recurrence sequences).** Fix alpha != 0 and let S be the
set of integer solutions of (*). The map T(U, V) = (V, 3V - U) is a bijection of S (inverse
(U, V) -> (3U - V, U); V^2 - 3V(3V - U) + (3V - U)^2 = U^2 - 3UV + V^2). For (U, V) in S define
Z_0 = U, Z_1 = V and Z_{n+1} = 3 Z_n - Z_{n-1} for all n in Z; then T^n(U, V) = (Z_n, Z_{n+1}).

*Lemma A (seed bound).* Every T-orbit in S contains a pair (U, V) with U^2 <= alpha/2 if alpha > 0
and U^2 <= |alpha| if alpha < 0.

Proof. Let (Z_n) be the sequence of an orbit and choose m with |Z_m| minimal. Since
Z_{m+1} + Z_{m-1} = 3 Z_m and (Z_m, Z_{m+1}) is in S,
Z_{m+1} Z_{m-1} = Z_{m+1}(3 Z_m - Z_{m+1}) = Z_m^2 - alpha. By minimality |Z_{m+1} Z_{m-1}| >= Z_m^2.
If alpha > 0: Z_m^2 >= alpha would give 0 <= Z_m^2 - alpha < Z_m^2, a contradiction; so
Z_m^2 < alpha and alpha - Z_m^2 >= Z_m^2. If alpha < 0: Z_{m+1} Z_{m-1} = Z_m^2 + |alpha| > 0, so
Z_m != 0 (the sum of the two neighbours would be 0) and both neighbours have the sign of Z_m; after
changing the sign of the whole sequence, a = Z_{m-1} and b = Z_{m+1} satisfy a, b >= z := Z_m > 0 and
a + b = 3z, hence a, b in [z, 2z] and ab = a(3z - a) >= 2 z^2 (a concave function, minimal at the
endpoints). So z^2 + |alpha| >= 2 z^2. The pair (Z_m, Z_{m+1}) is in the orbit. QED

So S is the union of the orbits of the finitely many "seeds" (U, V) in S with |U| <= B_0,
B_0 = floor(sqrt(alpha/2)) or floor(sqrt|alpha|). For such U the V are the integer roots of
V^2 - 3UV + U^2 - alpha, i.e. V = (3U +- s)/2 with s^2 = 5U^2 + 4 alpha.

- alpha = 461: B_0 = 15; 5U^2 + 1844 is a square only for U = +-1 (43^2). Seeds
  (1, 23), (1, -20), (-1, 20), (-1, -23).
- alpha = -599: B_0 = 24; 5U^2 - 2396 is a square only for U = +-24 (22^2). Seeds
  (24, 47), (24, 25), (-24, -25), (-24, -47).

(The seed lists are recomputed by `tau_p4.gp`, by `sieve_tau_p4.py`, and, by brute force over U and
V, by `check_certificate.py`.) Hence: if x^22 - 3x^11 y^2 + y^4 = alpha then for one of the four
seeds and some n in Z,

    Z_n = x^11   and   Z_{n+1} = y^2.                                            (**)

**Step 2 (periodicity).** Let q be a prime and A = [0, 1; -1, 3], so that
(Z_n, Z_{n+1})^t = A^n (Z_0, Z_1)^t for all n in Z (det A = 1). If A^M = I mod q then
n -> (Z_n, Z_{n+1}) mod q has period dividing M, for every seed and in both directions. The orders
T_q of A mod q used below are

    q    3   7   43   89   199   967
    T_q  4   8   44   22    11    88

(computed in `tau_p4.gp` as the matrix order; `check_certificate.py` verifies A^M = I mod q by an
independent exponentiation; `tau_p4.gp` in addition compares exact integers Z_n with the periodic
table for -2M <= n <= 2M + 1). All divide 88; T_89 = 22.

**Step 3 (the sieve).** If (**) holds then for every prime q: Z_n mod q lies in
P_q = {x^11 mod q} and Z_{n+1} mod q lies in Q_q = {y^2 mod q} (0 included in both). A residue class
n mod M is *killed* by q (with T_q | M) if Z_n mod q is not in P_q (code P) or Z_{n+1} mod q is not
in Q_q (code S); this depends only on n mod M by Step 2.

- alpha = -599, M = 22, q = 89 alone: for each of the 4 seeds all 22 classes are killed
  (88 lines of the certificate). P_89 has 9 elements (89 = 8*11 + 1).
- alpha = 461, M = 88, q in {3, 7, 43, 89, 199, 967}: for each of the 4 seeds all 88 classes are
  killed (352 lines).

`certificate_tau_p4.txt` lists, for every seed and every class n mod M, a witness (q, code,
Z_n mod q, Z_{n+1} mod q). No class survives, so (**) is impossible and Theorem 1 follows. QED

**Remark (which condition does the work).**
- alpha = -599: the 11th-power condition alone suffices. With M = 22 and q in {89, 199}, no class
  has Z_n in P_89 and in P_199 simultaneously (checked by me, by the literature scout, and by the
  referee, each with separate code). So u^2 = 5x^22 - 2396 has no integer point at all, and this is
  Goto's own Lemma 2.15 method (Dembner-Jain) applied to the right sequence; why Goto's v1 lists
  -599 as not excluded I do not know. The case -599 is therefore an oversight in v1 rather than a
  hard case.
- alpha = 461: here the square condition is essential. The sequence contains Z_0 = +-1, an 11th
  power, so no congruence on Z_n alone can work; with the 11th-power condition only, the classes
  n = 0, 11 mod 22 survive on every seed. The neighbours of +-1 are 23, -20 (resp. -23, 20), none
  a square, and the square condition kills those classes. (The scout's reading of Dembner-Jain is
  that such cases needed linear forms in logarithms there; I did not read that part myself.)
- The six witness primes for 461 are not a minimal set (the referee notes that 43 or 89 can be
  dropped); they are simply what the generator picks.

## 3. Verification

1. `sieve_tau_p4.py` (Python/numpy; first implementation): modulus 55440, primes up to 100000 with
   period dividing 55440; zero survivors for 461 and -599 (`cert_461.json`, `cert_m599.json`).
   Planted solutions (U, V) = (2^11, 9) and (3^11, 25) (alpha = 4139089, 31367774209) survive in
   exactly one class, as they must; brute force over |n| < 400 with exact integers finds no n with
   Z_n an 11th power and Z_{n+1} a square.
2. `tau_p4.gp` (PARI/GP 2.15.4; written separately; residue tests by Euler's criterion and the
   Kronecker symbol instead of enumeration): regenerates seeds, periods and the full witness table
   with M = 88 (461) and M = 22 (-599); zero survivors. Run: `gp -q tau_p4.gp` (12 s).
3. `check_certificate.py` (Python, exact integers, brute-force residue sets, brute-force seeds):
   re-verifies every line of the certificate: `CERTIFICATE OK`, 1376 class lines over 10 blocks.

## 4. Proof of Corollary 2, and what is now known for odd primes l < 1000

**Reduction (literature).** Let l be an odd prime with l not dividing tau(l) (for l < 1000 this
excludes only 3, 5, 7; checked in `tau_p4.gp`) and tau(n) = +-l. Two papers are cited; to avoid
confusion (their lemma numbers collide):
- BCOT = Balakrishnan-Craig-Ono-Tsai, Adv. Math. 428 (2023), arXiv:2005.10354. Theorem 1.1 (proved
  in Section 3 together with Theorem 3.2): n = p^{d-1} with p and d odd primes, d | l(l^2 - 1).
  Section 5 (the paragraph after (5.1)): for primes p not dividing tau(p), {tau(p^m)} is a Lucas
  sequence with no defective term equal to +-1 or +-l, "d + 1 must be an odd prime, and tau(p^d) is
  the very first term that is divisible by l". (In BCOT, Lemma 2.1 only classifies the defective
  terms; the no-defective conclusion for tau is drawn in Section 5.)
- BCO = Balakrishnan-Craig-Ono, J. Number Theory 237 (2022), arXiv:2005.10345. Its Lemma 2.1 states
  the same no-defective fact directly ("there are no defective Lucas numbers u_n(alpha, beta) in
  {+-1, +-l}" for X^2 - AX + p^11), and its Section 3 gives the same reduction.

So u_d = tau(p^{d-1}) = +-l is not defective: l is a primitive prime divisor of u_d, hence l does
not divide (alpha - beta)^2 p^11, and therefore d | l - ((alpha-beta)^2 / l), in particular
d | l^2 - 1 (this step is written out in Goto, Section 3.1). This excludes d = l. I rely on these
published statements; I did not re-derive the Bilu-Hanrot-Voutier/Abouzaid case analysis behind
them.

**Congruences (recomputed, `tau_p4.gp` Part B).** For n = p^{d-1} with p an odd prime not in
{3, 5, 7, 23, 691} (so n = 1 mod 8, n = 1 mod 3):
tau(n) = sigma_11(n) mod 691 and mod 2^11; tau(n) = n^{-610} sigma_1231(n) mod 3^6;
tau(n) = n^{-30} sigma_71(n) mod 5^3; tau(n) = n sigma_3(n) mod 7; tau(p) = 0 mod 23 if (p/23) = -1
and tau(p) in {2, -1} mod 23 if (p/23) = 1 (then the Hecke recurrence mod 23). These are the
classical congruences (Ramanujan, Wilton, Bambah-Chowla-Gupta-Lahiri-Lehmer-Kolberg; see
Swinnerton-Dyer, LNM 350 (1973)); the script re-checks each of them for all n <= 3000 in the stated
class, and the resulting residue sets against exact values tau(p^{d-1}) for p <= 113 and
d in {3,5,7,11,13,23}: 0 violations. For p in {3, 5, 7, 23, 691} the script checks
tau(p^{d-1}) != +-l directly.

Result: among all (eps l, d), l < 1000 an odd prime not in {3,5,7}, d an odd prime dividing l^2 - 1,
the congruences exclude every pair except these 17:

    d = 3 :  -137, -227, 463, 601, 691, -827
    d = 5 :  -139, -179, 281, 461, -599, -919
    d = 23:  -139, -367        d = 59: 827        d = 83: -829        d = 107: -643

(The same list was obtained independently by the literature scout with a separate Python script.)
For 461 the only pair is (461, 5); for -599 only (-599, 5); for -461 and 599 there is none.
In detail (`pair_table.py`, a separate Python computation; 461^2 - 1 = 2^3 3 5 7 11 23,
599^2 - 1 = 2^4 3 5^2 13 23):

| value | d | excluded modulo | not excluded modulo |
|---|---|---|---|
| 461 | 3 | 3^6 | 691, 2^11, 5^3, 7, 23 |
| 461 | 5 | **none** | 691, 2^11, 3^6, 5^3, 7, 23 |
| 461 | 7 | 691, 3^6, 7 | 2^11, 5^3, 23 |
| 461 | 11 | 7 | 691, 2^11, 3^6, 5^3, 23 |
| 461 | 23 | 7 | 691, 2^11, 3^6, 5^3, 23 |
| -461 | 3 | 691, 5^3, 23 | 2^11, 3^6, 7 |
| -461 | 5 | 691, 2^11, 5^3 | 3^6, 7, 23 |
| -461 | 7 | 3^6, 5^3, 23 | 691, 2^11, 7 |
| -461 | 11 | 5^3 | 691, 2^11, 3^6, 7, 23 |
| -461 | 23 | 5^3 | 691, 2^11, 3^6, 7, 23 |
| 599 | 3 | 3^6, 5^3 | 691, 2^11, 7, 23 |
| 599 | 5 | 2^11, 5^3 | 691, 3^6, 7, 23 |
| 599 | 13 | 691, 2^11, 3^6, 5^3, 7 | 23 |
| 599 | 23 | 5^3 | 691, 2^11, 3^6, 7, 23 |
| -599 | 3 | 691, 23 | 2^11, 3^6, 5^3, 7 |
| -599 | 5 | **none** | 691, 2^11, 3^6, 5^3, 7, 23 |
| -599 | 13 | 3^6, 7, 23 | 691, 2^11, 5^3 |
| -599 | 23 | 7 | 691, 2^11, 3^6, 5^3, 23 |

With Theorem 1 (the two rows marked "none") this proves Corollary 2. QED

**Status of each surviving pair.**

| pair | status | how |
|---|---|---|
| (461, 5), (-599, 5) | **excluded here**; open in Goto's arXiv v1 (published version not seen) | Theorem 1 |
| (281, 5), (-919, 5) | excluded | Goto Lemma 2.15 (11th powers mod 89, 199); re-proved by the same sieve, M = 22 |
| (-139, 5), (-179, 5) | excluded | l < 100..252: Dembner-Jain / Lin-Ma; re-proved by the same sieve, M = 22 |
| (-827, 3) | **excluded here** (not in Goto v1) | t^2 + 827 = p^11: h(Q(sqrt(-827))) = 7 is prime to 11, so t + sqrt(-827) = ((a + b sqrt(-827))/2)^11; comparing coefficients, b F(a,b) = 2^11 with b = +-2^j; no integer solution (`tau_p4.gp` Part C). Coprimality of (t +- sqrt(-827)): p is odd and p != 827 (827^11 - 827 is divisible by 827 exactly once). |
| (-137, 3), (-227, 3) | excluded | same argument: D = 227 as for 827 (disc -227, h = 5); D = 137: ring Z[sqrt(-137)], disc -548, h = 8, t + sqrt(-137) = (a + b sqrt(-137))^11 and b F(a,b) = 1; p != D since D^11 - D is divisible by D exactly once. Also inside Lin-Ma's range |l| <= 252. Class numbers 7, 5, 8 confirmed by `bnfcertify` (and by the referee by enumerating reduced forms). |
| (691, 3) | excluded | BCOT Theorem 1.2(1): tau(n) != +-691 |
| (-139, 23), (-367, 23) | excluded, certified | F_22(X,Y) = -139, -367 has no integer solution: `thue_check.gp`, `thueinit(P, 1)` (no GRH) |
| (827, 59), (-829, 83), (-643, 107) | **excluded assuming the uncertified PARI computation** | Goto Lemma 2.14; re-run here: F_58 = 827, F_82 = -829, F_106 = -643 have no integer solution, PARI `thue` with uncertified `thueinit`, tentative class number 1 in all three fields (`thue_results.txt`; see "On GRH" below) |
| **(463, 3), (601, 3)** | **OPEN** | tau(p^2) = 463 or 601, i.e. t^2 = p^11 + 463, + 601. No prime p < 10^5. Not treated in Goto's arXiv v1. |

**Conclusion.** For every odd prime l < 1000 and every n: tau(n) != +-l, **except possibly
tau(p^2) = 463 and tau(p^2) = 601**, and subject to the cited results: for l = 3, 5, 7 (not
ordinary, not covered by Part B) BCO Theorem 1.1, tau(n) not in {+-1, +-3, +-5, +-7, +-691}; the
BCO/BCOT reduction to n = p^{d-1}; the absence of defective terms (BCO Lemma 2.1 = BCOT Section 5),
which is the only thing excluding d = l (the congruences do not exclude (+l, d = l) for any l); BCOT Theorem 1.2(1) for
(691, 3). Lin-Ma
(Ramanujan J. 63 (2024), Cor. 1.3, |tau(n)| > 252 for odd tau(n), n > 3; known to me only through
Goto's citation) is *not* needed, because Part B was run from l = 11 up and every surviving pair
with l <= 251 is handled in the table; Bennett-Gherga-Patel-Siksek is likewise not needed for
this table. The three largest Thue equations (d = 59, 83, 107) rest on PARI's uncertified `thueinit`: for
l = 827, 829 (minus sign), 643 (minus sign) the exclusion is conditional on that computation.

**On GRH.**
- Bennett-Gherga-Patel-Siksek (Math. Ann. 382 (2022), arXiv:2101.02933) is unconditional and says
  so: the remark at the end of their Section 11 notes that the Thue-Mahler solver needs the class
  and unit groups of Q(zeta_m)^+, that h_m^+ = 1 is known for all primes m <= 67, and that the
  units are then the cyclotomic units.
- Goto's Lemma 2.14 says only "Calculating by PARI/GP"; whether `thueinit(P, 1)` was used is not
  stated. PARI's documentation of `thue` says an uncertified result "is conditional on the GRH if
  a != +-1 and P has a single irreducible rational factor, whose attached tentative class number h
  and regulator R (as computed assuming the GRH) satisfy h > 1, R/0.2 > 1.5". The fields here are
  Q(zeta_d)^+; I found tentative h = 1 for d = 59, 83 and 107, so by PARI's stated criterion those
  outputs are unconditional. This is PARI's criterion, not something I proved. For d = 23 I ran
  the certified version. (Independently, h^+ = 1 is a theorem for these conductors: van der Linden
  for d <= 67 as cited by BGPS; for 83 and 107 the scout cites J. C. Miller, which I did not read.)
- Nothing in Theorem 1, Part B or Part C uses GRH (the class numbers in Part C are certified with
  `bnfcertify`; `qfbclassno` alone carries a correctness caveat in PARI's documentation).

**Goto's arXiv v1 versus the published version.** In v1 the proof of Theorem 0.3 lists the pairs
(l, d) not excluded by "Lemma A.1"; the list contains no pair with d = 3. But (601, 3) passes
every congruence of Lemma A.1 as printed (601 = 1 mod 12, 1 mod 5, 6 mod 7, 3 mod 23), and the
mod 12 line of Lemma A.1 ("tau(p^{d-1}) = 0, 1, d mod 12") is false for d = 3: tau(5^2) = -25499225
= 7 mod 12 (the true set for d = 3 is {1, 3, 7, 9}; the referee finds it false for d = 11 too),
and 463 = 7 mod 12. So v1 does not exclude tau(p^2) = 463, 601 (nor -827, which Part C
excludes). The published version (Arch. Math. 125 (2025) 157-172, doi 10.1007/s00013-025-02139-5)
is behind a paywall; I could read only its abstract on the Springer page, which no longer says
"except 14 cases" but: "We prove that tau(n) != +-l, +-2l, +-4l, +-8l for l in L, where L is an
explicit finite subset of odd primes less than 1000." I do not know the set L, so I cannot say
whether 461, 599, 463, 601, 827 are in it.

## 5. Other cases by the same sieve

All in `certificate_tau_p4.txt`, checked by `check_certificate.py`.

- alpha = 281, -919 (Goto Lemma 2.15), -139, -179: all classes killed with M = 22 (q = 89, 199).
- Hanada-Madhukara (Acta Arith. 200 (2021), arXiv:2007.08683), Remark 6.2, "were not able to
  determine" the integer points of H^+_{8,33}, H^+_{8,55}, H^+_{13,33}, H^+_{13,55}, where
  H^+_{k,c}: Y^2 = 5 X^{2(2k-1)} + 4c comes from tau_{2k}(p^4) = c (weights 16 and 26).
  - c = 33: U^2 - 3UV + V^2 = 33 has **no integer solution at all** (mod 3 it forces 3 | U and
    3 | V, then 9 | 33). So H^+_{8,33} and H^+_{13,33} have no integer points; no sieve needed.
  - c = 55: seeds (+-1, ...); with e = 15 (M = 40) and e = 25 (M = 84) every class is killed except
    n = 0 mod M for the seed (Z_0, Z_1) = (1, 9), which is the genuine solution 1 - 27 + 81 = 55
    with U = 1. Suppose Z_n = x^e and Z_{n+1} = y^2 with n != 0, n = 0 mod 4. Write n = 2Kk,
    K = 2^a, a >= 1, k odd. With V_K = lambda^K + lambda^{-K}, lambda = (3 + sqrt 5)/2 (V_1 = 3,
    V_{2K} = V_K^2 - 2), one has Z_{m+2K} = V_K Z_{m+K} - Z_m, hence Z_{1+2Kk} = -Z_1 = -9 mod V_K.
    V_2 = 7 and V_{2K} = V_K^2 - 2 give V_K = 7 mod 8 and V_K = 1 or 2 mod 3 for all K = 2^a,
    a >= 1. So the Jacobi symbol (-9 / V_K) = (-1 / V_K) = -1, while y^2 = -9 mod V_K with
    gcd(y, V_K) = 1 (a common prime would divide 9) would make it +1. Contradiction (Cohn's
    argument). So the only solutions have x^e = 1, and x is not a prime: tau_16(p^4) != 55 and
    tau_26(p^4) != 55.
  (This determines only the points with V a square, i.e. of the form needed for tau_{2k}(p^4);
  it does not list all integer points of H^+_{k,55}, which has the point X = 1.)
  What this gives for Hanada-Madhukara's Theorem 1.2: by their proof, +33 and +55 were excluded
  from the lists for weights 16 and 26 only because of these four curves, so they can now be
  included, *provided* the other steps of their proof for 33 and 55 (curves C^+_{k,c}, Thue
  equations, Lemma 3.1) are as stated there. I did not re-verify those steps and their text does
  not say whether they used GRH for these particular c.

## 6. Novelty check (2026-11-26)

- arXiv API queries (titles/abstracts on Ramanujan tau values, Lehmer, prime values; sorted by
  date) on 2026-11-24 and 2026-11-26: the only later papers are arXiv:2603.29970 (Angdinata et
  al., "ABC implies that Ramanujan's tau function misses almost all primes", v3 25 Apr 2026),
  which states tau(p^4) = Y^2 - 3XY + X^2 with (X, Y) = (p^11, tau(p)^2) (its (2.5)) but then
  passes to u^2 = 5x^22 +- 4l and only counts; arXiv:2311.12073v4 (Xiong); arXiv:2512.02345;
  arXiv:2604.13365. None mentions 461 or 599 or cites Goto (text search of the first three).
- Crossref lists 0 citations of the published Goto paper. Semantic Scholar's API was rate-limited.
- Not checked: MathSciNet/zbMATH reviews, the published version of Goto, non-arXiv preprints.

## 7. Adversarial review (fresh subagent, 2026-11-26)

A separate agent was asked to break this document, with its own scripts (none of the repository
code reused). Verdict: no fatal error; Theorem 1 and Corollary 2 hold up. It confirmed
independently: the Hecke identity (p <= 59); Lemma A by hand (zero terms, ties, sign flip) and by
a stress test (for every alpha with 0 < |alpha| <= 700, brute-force seeds equal the formula seeds
and all 14,548 solutions with |U|, |V| <= 4000 lie on seed orbits); the seed lists for 461 and
-599 (all solutions with |U| <= 10^6 lie on the four sequences); the periods; zero survivors with
its own sieve for (461, M = 88, the six primes) and (-599, M = 22, q = 89), and again with
M = 55440 and 41 other primes; non-vacuity on seven planted solutions; the six congruences for
n <= 4000; the same 17 surviving pairs; the Part C argument (with positive controls for other
exponents); the Section 5 arguments; the quotations from Goto v1, BCO, BCOT, BGPS,
Hanada-Madhukara and the PARI manual; and that `tau_p4.gp` reproduces the certificate byte for
byte.

Its objections, all addressed above: (1) the d = 107 row was stale (fixed; `thue_results.txt`);
(2) the Conclusion did not cite anything for l = 3, 5, 7, nor BCO Lemma 2.1 and BCOT 1.2(1)
(added); (3) -599 does not need the square condition (now stated in the remark after Step 3);
(4) "new" was claimed although the published version of Goto was not read (reworded);
(5) `qfbclassno` caveat (class numbers now certified); (6)-(8) wording.

**What remains unverified by anyone here:** the published version of Goto's paper (the set L);
the Bilu-Hanrot-Voutier analysis behind BCO Lemma 2.1 / BCOT Section 5; the correctness of PARI's "tentative
h = 1 implies unconditional" criterion for the degree 29, 41, 53 Thue equations; the other steps
of Hanada-Madhukara's proof for c = 33, 55.

## 8. Second referee (parent session, 2026-11-26)

An independent referee commissioned by the parent session (report `REFEREE-2-TAU4.md` on the parent
branch) found Theorem 1 complete and Corollary 2 complete modulo the cited literature, and asked
for three changes, made above: the BCO/BCOT citation is now unambiguous (Section 4, "Reduction");
the per-pair table for (+-461, d) and (+-599, d) was added; the rows for d = 59, 83, 107 are
relabelled as resting on an uncertified PARI computation.
