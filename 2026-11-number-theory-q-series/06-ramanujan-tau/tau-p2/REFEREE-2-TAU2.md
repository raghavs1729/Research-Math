# Referee report 2 on TAU2.md (tau(p^2) = 601 and 463)

Date: 2026-11-26. Independent adversarial check. No existing file was edited. My scripts are in
`/tmp/claude-0/-home-user-weighted-prime-races/7b9224b4-6348-502d-a364-b106f526c180/scratchpad/tau601/`
(`ref601.py`, `ref_tau23_and_search.py`, `ref463_form.py`). They share no code with `local_k.gp` or
`local_k_check.py`. No `thueinit(P, 1)` was run, and I did not rerun any Thue computation.

## Verdict

- **Claim 1 (tau(p^2) = 601 is impossible, unconditionally): PROOF COMPLETE.** The argument rests
  on two classical facts: the Hecke relation tau(p^2) = tau(p)^2 - p^11 and Wilton's congruence for
  tau mod 23. Everything else is a finite computation, and I recomputed each one independently. I
  found no gap. The few inaccuracies are cosmetic (findings 8 and 9).
- **Claim 2 (463), status only:** the document labels it honestly. It says the case rests on an
  uncertified PARI computation and that the certified run did not finish. I add one caveat about the
  manual's "right-hand side +-1" exemption (finding 10). The positive control makes sense for what it
  tests, but it says nothing about certification (finding 11).

## Findings: Claim 1

1. **Descent.** The setting is K = Q(sqrt 601) with t^2 - p^11 = 601.
   - **Coprimality.** Suppose a prime ideal divides both t + sqrt 601 and t - sqrt 601. Then its
     rational prime divides p and N(2 sqrt 601) = 4*601, so p = 2 or p = 601. Both are excluded:
     tau(4) = -1472, and 601^11 + 601 has 601-adic valuation exactly 1. This also rules out 601 | t,
     which would force p = 601. So (t + sqrt 601) is the 11th power of an ideal. Since h = 1 it is
     principal, and the sign is absorbed because -1 = (-1)^11. Equation (D) is correct and k mod 11
     is exhaustive. (Parity of t plays no role: t is even, which is harmless.)
   - **Half-integral g.** O_K = Z[(1 + sqrt 601)/2], so g may have half-integral coordinates. For odd
     q, O_K/q is isomorphic to F_q[y]/(y^2 - 601), so working with pairs mod q covers every g. My own
     code works in the basis (1, omega), omega = (1 + sqrt 601)/2, and gets identical results.
   - **Units.** The mod-8 argument holds: for X, Y odd, X^2 - 601 Y^2 = 0 mod 8, never +-4. So every
     unit lies in Z[sqrt 601].
2. **Fundamental unit.** I recomputed it from the continued fraction of omega = (1 + sqrt 601)/2, as
   the first convergent h/k with N(h - k*omega) = +-1. The result is
   eps = (278936607359064 + 11378061539690 sqrt 601)/2 = 139468303679532 + 5689030769845 sqrt 601,
   with norm -1. This equals eps_0 in TAU2.md. Because it is the least unit > 1 of O_K, it is
   fundamental, not a power.
3. **Class number, independent of PARI.** I used the finite formula
   h log eps = -1/2 sum_{a=1}^{600} (a/601) log sin(pi a/601), evaluated with mpmath at 50 digits. It
   gives h = 1.0000... This agrees with PARI's `bnfcertify`, and 11 does not divide h.
4. **Congruence mod 23.**
   - I computed tau(n) mod 23 from Delta = q (prod (1 - q^n)^3)^8, using Jacobi's identity, and
     checked it against known values of tau(2), ..., tau(23).
   - For all primes p < 6000 with p != 23, Wilton's classification holds with no exceptions:
     - tau(p) = 0 if (p/23) = -1;
     - tau(p) = 2 if p = x^2 + xy + 6y^2;
     - tau(p) = -1 otherwise.
   - The deduction is correct. 601 = 3 mod 23 (also 463 = 3 mod 23), and p^11 = (p/23) mod 23.
     t^2 - p^11 = 3 forces (p/23) = 1 and t^2 = 4, so t = +-2. Since -2 is not in {2, -1},
     t = 2 mod 23.
   - The case p = 23 is excluded: 23^11 + 601 is not a square (I checked).
5. **Local test, recomputed.** I used the power-residue criterion: (z1^((q-1)/11), z2^((q-1)/11)) =
   (1, 1) for split q, and z^((q^2-1)/11) = 1 for inert q. This differs from the brute-force
   enumeration of 11th powers in `local_k_check.py`. With eps_0 as orientation:
   - A(23, k), with only t = 2 allowed, is empty for every k except k = 7.
   - A(67, k) is empty exactly for k in {0, 1, 4, 7, 10}.
   - Both rows match the table in TAU2.md. So every k in Z/11 is killed.
   - 23 and 67 both split. Neither 23^11 + 601 nor 67^11 + 601 is a square, so p != q.
   - I also ran the document's own `local_k_check.py` and `local_k.gp`; they reproduce the stated
     output.
6. **Robustness, and where the proof is fragile.**
   - The class k = 7 is also killed by q = 89. So the class that escapes q = 23 is killed twice.
   - Without the restriction from step 4, A(23, k) is nonempty for k = 3, 4, 7, 8 (t = 5, 21, 2, 18).
     The classes k = 3 and k = 8 then survive every prime q < 700 (the only other killers are
     q = 67, 89, 199).
   - So the proof genuinely needs Wilton's congruence, not just local data. That is fine, since the
     congruence is classical and verified in step 4, but it is the load-bearing step.
   - A consistency check: N(eps_0) = -1 makes t -> -t swap class k with class -k. This matches
     A(23, 3) = {5} against A(23, 8) = {18}, and A(23, 4) = {21} against A(23, 7) = {2}.
7. **Direct search.** I checked every odd x < 10^7 (the document went to 2*10^6). No x^11 + 601 and
   no x^11 + 463 is a perfect square. I used a quadratic-residue prefilter over 14 moduli followed by
   an exact `isqrt`. The prefilter passes the known c = 94 solution x = 3.
8. **Cosmetic.** TAU2.md Section 2 calls k = 4 (PARI orientation; k = 7 for eps_0) "the class killed
   only by 67". It is also killed by 89. "Only by 67 among 23 and 67" is what is meant.
9. **Cosmetic.** Section 3 says the positive control finds the solution "for exactly one unit class"
   and then lists two classes, k = -4 and k = 4. These are the classes of t = 421 and t = -421, which
   are conjugate. The wording should say one class per sign of t.

## Findings: Claim 2 (status only)

10. **Honesty of labelling: adequate.** The document never calls 463 unconditional:
    - the table says "excluded, given one PARI Thue computation";
    - Sections 3 and 5 say the runs are uncertified, the certified run stopped after 3000 s, and "I
      have not proved that criterion myself";
    - Section 4 lists 463 among the PARI dependencies.

    I checked the manual quote against the `??thue` entry of the installed PARI 2.15.4. It says
    "conditional on the GRH if a != +-1 and ..." and also "all results computed using this
    particular tnf are likewise conditional, except for a right-hand side of +-1". So the document
    reports the manual correctly.

    Caveat: the manual gives no reason for the +-1 exemption. For right-hand side +-1 the class group
    is irrelevant, but the answer depends entirely on having a full system of fundamental units. An
    uncertified `bnfinit` does not prove that, and the manual's own criterion R/0.2 > 1.5 is a
    unit-index condition. Here R is about 7.6*10^9. I could not get the PARI source to see what
    `thue` actually checks.

    Recommendation: keep the label "uncertified". Also note that only the unit part of the
    certification is needed. The class-group part of `bnfcertify` is irrelevant for +-1 and is
    plausibly the infeasible step: the Minkowski bound for disc -2^10 11^11 463^5 is about 3.7*10^10.
    So a units-only certification, for example with a regulator lower bound and an index check,
    could be a cheaper route.

    I also verified the form exactly in `ref463_form.py`. Q(a, b) as printed equals the coefficient
    of sqrt 463 in eps_0^(-2) (a w1 + b w2)^11, and not in eps_0^(+2). The basis is unimodular
    (det = -1), and k = -2 = 9 is the class that survives the local test (I reproduced that A(23, k)
    is empty for all k != 9). PARI's `K.fu` is -eps_0, which does not matter. Also, because
    N(eps_0) = +1, Q_k = 1 and Q_{-k} = 1 are equivalent under conjugation, so the orientation cannot
    affect the Thue step.
11. **The positive control (c = 94) makes sense, with limited reach.**
    - The setup is right: 421^2 - 94 = 177147 = 3^11, and h(Q(sqrt 94)) = 1.
    - My local test, using eps = 2143295 + 221064 sqrt 94, independently leaves exactly the classes
      {4, 7}: t = 421 in class 4 and t = -421 in class 7. This matches the document's k = -+4 in
      PARI's orientation.
    - It confirms that `setup()` builds the correct form and orientation, and that the local test
      keeps a true class alive.
    - It does not test the completeness of PARI's `thue`, which is the only open point. The solution
      (a, b) = (1, 0) is trivially small. The control also does not exercise the mod-23 filter,
      which has no analogue for c = 94; that filter is verified separately in finding 4.

## What I recomputed

- The fundamental units of Q(sqrt 601), Q(sqrt 463) and Q(sqrt 94), from continued fractions.
- h(Q(sqrt 601)) = 1 from the analytic class number formula.
- Wilton's congruence for p < 6000, and the deduction t = 2 mod 23.
- A(q, k) for every k in Z/11 and every q < 700 with 11 | q^2 - 1, for c = 601, 463 and 94.
- That q^11 + c is not a square for q = 23 and q = 67.
- The direct search for odd x < 10^7.
- The exact binary form Q for 463.
- I also ran the document's `local_k_check.py` and `local_k.gp`, which agree with all of the above.
