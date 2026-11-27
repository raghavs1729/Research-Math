# Referee report 2: tau(p^4) != 461, -599 (PROOF.md)

Date: 2026-11-26. Independent adversarial check. I reused none of the repository code. My scripts are in
`/tmp/claude-0/-home-user-weighted-prime-races/7b9224b4-6348-502d-a364-b106f526c180/scratchpad/tau4ref/`
(`seeds_sieve.py`, `sieve.py`, `brute.py`, `congr.py`, `hecke.py`).

## Verdicts

- **(T) Theorem 1: PROOF COMPLETE.** I found no mathematical gap. Every step was recomputed with my own code.
- **(C) Corollary 2: PROOF COMPLETE, modulo the cited literature.** It relies on Balakrishnan-Craig-Ono-Tsai
  (BCOT) Theorem 1.1 and its Lucas-sequence or primitive-divisor input, and on the classical congruences for
  tau. The reduction is quoted correctly in substance and is unconditional. One citation pointer is imprecise
  (Finding 1); this is not a gap.
- **(S) Side claims:** the status labels are honest. tau(p^2) = 463 and 601 are clearly marked OPEN. One label
  could be sharper (Finding 3).

## What I recomputed

1. **Hecke identities.** I computed tau from the q-expansion of Delta and checked tau(p^4) = t^4 - 3p^11 t^2 + p^22
   and tau(p^2) = t^2 - p^11 for p = 2, 3, 5.
2. **Lemma A (seed bound), checked by hand.**
   - Case alpha > 0: Z_{m+1} Z_{m-1} = Z_m^2 - alpha and minimality give Z_m^2 <= alpha/2. This is correct,
     including the case Z_m = 0.
   - Case alpha < 0: the product is positive, so Z_m != 0 and both neighbours have the sign of Z_m. Then
     a + b = 3z with a, b >= z gives ab >= 2z^2, hence z^2 <= |alpha|. This is correct.
   - A minimum of |Z_n| exists because the values are integers.
   - Non-primitive solutions: gcd(U, V)^2 divides alpha, and 461 and 599 are prime, so every solution is
     primitive. The bound does not use primitivity anyway.
   - Orderings and signs: S contains both (U, V) and (V, U), and the seed list consists of all pairs with
     |U| <= B_0, so both orderings and both signs are covered automatically.
3. **Seeds, by brute force.** I took every U with |U| <= 2*10^6 and solved for V exactly.
   - alpha = 461: 96 solutions with |U|, |V| <= 2*10^6.
   - alpha = -599: 100 solutions in the same range.
   - Every solution lies on the orbit of one of the four listed seeds: (+-1, 23), (+-1, -20) for 461, and
     (+-24, +-25), (+-24, +-47) for -599. The four orbits are pairwise distinct.
   - The seed lists in PROOF.md and the certificate are exact.
4. **Sieve, rebuilt from scratch.**
   - I computed the matrix order of A = [0, 1; -1, 3] mod q: T_3 = 4, T_7 = 8, T_43 = 44, T_89 = 22,
     T_199 = 11, T_967 = 88. These are as stated.
   - The primes q < 1000 with T_q | 88 are 3, 7, 43, 89, 199, 263, 307, 881 and 967. This matches the
     certificate header.
   - I built the residue sets by brute-force enumeration (11th powers and squares, 0 included).
   - alpha = 461, M = 88, the six primes: 0 survivors over 4 seeds x 88 classes.
   - alpha = -599, M = 22, q = 89 alone: 0 survivors.
   - Each of the 440 witness lines in `certificate_tau_p4.txt` for 461 and -599 is a valid kill in my
     recomputation: correct residues Z_n, Z_{n+1} mod q, T_q | M, and the stated code (P or S) holds.
     There were 0 errors.
   - **Remark after Step 3 confirmed.**
     - -599: the 11th-power condition with q in {89, 199} alone kills everything. With q = 89 alone, 16
       classes survive without the square condition.
     - 461: with the 11th-power condition only, the classes n = 0 and n = 11 mod 22 survive on all four
       seeds, as stated.
     - Prime-dropping test: dropping 43 or 89 still leaves 0 survivors, as stated. Dropping 3, 7, 199 or 967
       leaves 1 survivor each.
5. **All integers x, y.** The theorem is proved in full generality, not only for x = p > 0. The set P_q of
   11th-power residues contains the residues of negative 11th powers, since (-x)^11 = -x^11. The sign of y is
   irrelevant. The seeds include U < 0.
   - Trivial cases y = 0, x = 0 and x = +-1 have no solutions; I checked them directly.
   - Brute force over |x| <= 50 (solving exactly for y^2 = (3x^11 +- sqrt(5x^22 + 4 alpha))/2) found no
     solution for either alpha.
6. **Corollary reduction.**
   - tau(461) = 101 mod 461 and tau(599) = 280 mod 599, so both primes are ordinary and BCOT Theorem 1.1
     applies.
   - The odd primes d dividing l(l^2 - 1) are:
     - l = 461: d in {3, 5, 7, 11, 23, 461}.
     - l = 599: d in {3, 5, 13, 23, 599}.
   - My own congruence sets (mod 691, 2^11, 3^6, 5^3, 7 and 23) exclude these pairs, by these moduli:

     | pair | excluded by |
     |---|---|
     | (+461, 3) | 3^6 only |
     | (-461, 3) | 691, 5^3, 23 |
     | (-461, 5) | 691, 2^11, 5^3 |
     | (+-461, 7) | several moduli |
     | (+461, 11) | 7 |
     | (-461, 11) | 5^3 |
     | (+461, 23) | 7 |
     | (-461, 23) | 5^3 |
     | (+599, 3) | 3^6, 5^3 |
     | (-599, 3) | 691, 23 |
     | (+599, 5) | 2^11, 5^3 |
     | (+-599, 13) | several moduli |
     | (+599, 23) | 5^3 |
     | (-599, 23) | 7 |
     | (-461, 461) | 691, 2^11, 5^3 |
     | (-599, 599) | 7 |

   - I checked the special primes p in {3, 5, 7, 23, 691} exactly.
   - Survivors of the congruences:
     - (461, 5) and (-599, 5): these are Theorem 1.
     - (461, 461) and (599, 599): these are excluded only by d | l^2 - 1, i.e. by the primitive-divisor
       argument. That argument is correct: u_d = +-l is prime and non-defective, so l is a primitive divisor,
       so l does not divide (alpha - beta)^2 by definition, and the rank of apparition d divides
       l - (D/l), which is not 0 mod l. For d = l > 30, non-defectiveness is just Bilu-Hanrot-Voutier.
   - **d = 3** is therefore handled by congruences: 461 by 3^6 alone, -599 by 691 and 23. I verified the
     congruences numerically for n < 760 (0 violations).
   - Composite n: BCOT reduces this via multiplicativity and |tau(m)| != 1 for m > 1 (BCOT Proposition 3.1,
     Lemma 2.1(1)). PROOF.md quotes this correctly.
   - Running the full computation for odd primes 11 <= l < 1000, I reproduced exactly the 17 surviving pairs
     with d != l listed in Section 4.
7. **Sources read.**
   - BCOT arXiv:2005.10354v4: Theorem 1.1 is quoted correctly ("Suppose that l is an odd prime for which
     l does not divide tau(l). If tau(n) = +-l^m ... then n = p^{d-1}, where p and d | l(l^2-1) are odd
     primes"). It is unconditional; only Theorem 1.2(2) assumes GRH. v4 added the "ordinary" hypothesis in
     several statements; it is satisfied here.
   - Goto arXiv:2405.16723: there is only v1 (26 May 2024), and it has no journal reference on arXiv.
     - Theorem 0.3, L_1^+ = {461}, L_1^- = {599}, and the Section 2.3 sentence are quoted correctly.
     - The Lemma A.1 mod 12 line "tau(p^{d-1}) = 0, 1, d (12)" is false for d = 3:
       tau(25) = -25499225 = 7 mod 12, and 463 = 7 mod 12.
     - So the referee point in PROOF.md is right: v1's proof of Theorem 0.3 silently drops (463, 3), (601, 3)
       and (-827, 3).
   - Goto, published version: Crossref gives Arch. Math. 125, 157-172, 2025-06-21, with no abstract
     deposited and 0 citations. Springer redirects to a login page, so I could not read the published
     version.
   - BGPS arXiv:2101.02933: the abstract mentions no conditional hypothesis. PROOF.md does not need it.
   - Two web searches found no 2025-2026 result on tau(n) = +-461, +-599, or on tau(p^2) = 463 or 601. They
     did find 2603.29970, 2512.02345 and 2311.12073, which PROOF.md already lists.

## Findings (by severity; none fatal)

1. **(Minor, citation; Section 4, "Reduction" and "Conclusion").**
   - "BCO Lemma 2.1 (no defective terms)" is an imprecise pointer. In BCOT (arXiv:2005.10354v4), Lemma 2.1 is
     the classification of defective Lucas numbers with |u_n| = 1 or l. The conclusion that the sequences
     {tau(p^j)} have no defective term with |u_n| = l odd is drawn in BCOT Section 5, from Lemma 2.1(2) and
     Lemma 4.2 (no points on Y^2 = X^11 +- 3), together with BHV Theorem 2.1 for n > 30.
   - Also, "(Section 3)" is attributed to the JNT exposition arXiv:2005.10345. The reduction to
     n = p^{d-1} is in BCOT Section 3 (Theorem 3.2 and the proof of Theorems 1.1 and 3.2).
   - Recommend citing BCOT v4 Theorem 1.1, Proposition 3.1, Lemma 2.1 and Section 5 directly.
   - The mathematics is correct. I re-derived the d | l^2 - 1 step myself (item 6).
2. **(Minor, presentation; Section 4).** The corollary gives the surviving list but does not say which
   congruence kills each (+-461, d) and (+-599, d). In particular, (461, 3) rests on a single modulus, the
   3^6 congruence tau(n) = n^{-610} sigma_1231(n) for n = 1 mod 3.
   - Recommend a small table, like the one in item 6, so that a reader can check (C) without running Part B.
   - This congruence is classical (Swinnerton-Dyer, LNM 350), and I confirmed it numerically.
3. **(Minor, labelling; (S) table).** The rows (827, 59), (-829, 83) and (-643, 107) are labelled "excluded",
   but rest on uncertified `thueinit` and PARI's stated GRH criterion. The prose and Conclusion say so.
   - Suggest the status read "excluded (uncertified PARI / Goto Lemma 2.14)".
   - This does not affect (T) or (C).
4. **(Minor, Lemma A wording).** For alpha > 0 the proof could say explicitly that Z_m = 0 is allowed; then
   alpha = Z_{m+1}^2, which is impossible for 461. The inequality argument already covers this case.
5. **(Informational).** -599 indeed needs no square condition, and 461 does. Both claims were confirmed
   (item 4).

## Not checked

- The published version of Goto (paywalled): whether its set L contains 461, 599, 463 or 601.
- The BHV / Abouzaid tables behind BCOT Lemma 2.1.
- Hanada-Madhukara's other steps for c = 33 and c = 55.
- Certification of the d = 59, 83, 107 Thue equations.
- MathSciNet and zbMATH.
