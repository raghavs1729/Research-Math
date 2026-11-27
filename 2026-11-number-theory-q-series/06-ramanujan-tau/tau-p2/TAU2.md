# tau(p^2) = 601 and tau(p^2) = 463

Date: 2026-11-26. These are the two pairs left open in `../tau/PROOF.md`, Section 4.
Files: `local_k.gp` (PARI/GP), `local_k_check.py` (independent Python check), `thue463.gp`,
`thue463_runs.gp`, `thue463_results.txt`, `control94.gp`.

## Result

| case | status | method |
|---|---|---|
| tau(p^2) = 601 | **excluded, unconditional** | descent in Q(sqrt 601) + local obstruction at 23 and 67; no Thue equation needed |
| tau(p^2) = 463 | **excluded, given one PARI Thue computation** (see Section 3 for its certification status) | descent in Q(sqrt 463) + local obstruction at 23 for ten of the eleven unit classes + one Thue equation of degree 11 |

## 1. Set-up (both cases)

Let c in {463, 601} and suppose tau(p^2) = c for a prime p. With t = tau(p): t^2 - p^11 = c.

Facts used about p and t:
- p is odd: tau(4) = -1472.
- p != c: c^11 + c is divisible by c exactly once, so it is not a square.
- p != 23 and **t = 2 mod 23**: c = 3 mod 23 for both values. For p != 23, Ramanujan's congruence
  gives tau(p) = 0 mod 23 if (p/23) = -1, and tau(p) in {2, -1} mod 23 if (p/23) = 1, while
  p^11 = (p/23) mod 23. Then t^2 - p^11 mod 23 is 1, 3 or 0 respectively; only t = 2 gives 3.
  For p = 23: 23^11 + c is not a square (checked).

Let K = Q(sqrt c), O_K its ring of integers. PARI (`bnfinit` + `bnfcertify`): h(K) = 1 for both.
Fundamental units (also recomputed in Python from the continued fraction of sqrt c):

    c = 463:  eps_0 = 247512720456368 + 11502891625161 sqrt(463),  norm +1   (O_K = Z[sqrt 463])
    c = 601:  eps_0 = 139468303679532 +  5689030769845 sqrt(601),  norm -1

(For c = 601 = 1 mod 8 every unit lies in Z[sqrt 601]: X^2 - 601 Y^2 = +-4 with X, Y odd is
impossible mod 8.)

Descent. (t + sqrt c)(t - sqrt c) = p^11. A prime ideal dividing both factors divides 2 sqrt c and
p, which is impossible (p odd, p != c). So (t + sqrt c) = a^11 for an ideal a; h = 1 gives a = (g),
and t + sqrt c = +- eps_0^j g^11. Absorbing -1 = (-1)^11 and eps_0^11 into g:

    t + sqrt c = eps_0^k g^11,   g in O_K,   k in {0, 1, ..., 10},   N(g) = +-p.          (D)

Local test. Let q be a prime, q != p. Then g is a unit mod q, so (t + sqrt c) eps_0^(-k) is an
11th power in (O_K/q)^*. For given k and q let A(q, k) be the set of t mod q with this property
(at q = 23 only t = 2 is allowed). If A(q, k) is empty, (D) is impossible for that k. Only primes
with 11 | q^2 - 1 can give information. For every q used, q^11 + c is not a square, so p != q.

## 2. c = 601: every k is impossible

    k (exponent of eps_0)   0   1   2   3   4   5   6   7   8   9   10
    A(23, k) empty          yes yes yes yes yes yes yes no  yes yes yes
    A(67, k) empty          yes yes no  no  yes no  no  yes no  no  yes

(`local_k_check.py`; `local_k.gp` gives the same with PARI's generator eps = 1/eps_0, i.e. with
k replaced by -k: there k = 4 is the class killed only by 67.) Details for the two primes:
23 splits (601 = 3 = 7^2 mod 23), the 11th powers in F_23^* are +-1, and t = 2 maps to
(2 + 7, 2 - 7) = (9, 18); 67 splits (601 = -2 mod 67, a square), 67 = 1 mod 11, and the 11th powers
in F_67^* form the subgroup of order 6.

Every k in Z/11 is killed by 23 or by 67. Hence **there is no prime p with tau(p^2) = 601.** This
uses: h(Q(sqrt 601)) = 1 (PARI, certified; only 11 not dividing h is needed), the fundamental
unit, Ramanujan's congruence mod 23, and finite computations mod 23 and 67. No GRH, no Thue
equation.

## 3. c = 463: one Thue equation

Local tests (`local_k.gp`, `local_k_check.py`): A(23, k) is empty for every k except k = 9.
The class k = 9 survives everything I tried: all q < 700 with 11 | q^2 - 1, all split q = 1 mod 11
up to 20000, and the 11-adic test (11 splits; one admissible t mod 121). So k = 9 looks locally
soluble everywhere and a global argument is needed.

For k = 9 = -2 mod 11 write t + sqrt 463 = eps_0^(-2) g^11 and g = a w1 + b w2 with the Z-basis

    w1 = 624 + 29 sqrt 463,   w2 = 1829 + 85 sqrt 463        (624*85 - 1829*29 = -1)

chosen by LLL so that the coefficients are small (the choice of basis does not affect the
statement; the form is computed exactly). Comparing coefficients of sqrt 463:

    Q(a, b) = 1,
    Q = 2069 a^11 + 143 a^10 b + 1072720 a^9 b^2 - 899250 a^8 b^3 + 61199160 a^7 b^4
        - 48219864 a^6 b^5 + 821606016 a^5 b^6 - 492420720 a^4 b^7 + 2836978320 a^3 b^8
        - 1043992400 a^2 b^9 + 1842900224 a b^10 - 226837664 b^11.

Q(x, 1) is irreducible with one real root; its field has signature (1, 5) (unit rank 5) and
discriminant -2^10 11^11 463^5.

PARI's `thue` finds **no integer solution of Q(a, b) = 1** (nor of = -1):
- run A: P(x) = Q(x, 1); run B: P(x) = Q(1, x); run C: Q(a + b, b) -- three different inputs,
  uncertified `thueinit`, a few seconds each (`thue463_results.txt`);
- brute force: no solution with |a|, |b| <= 500;
- tentative class number of the degree-11 field: 1.

Certification status. PARI's manual says an uncertified `thue` result "is conditional on the GRH
if a != +-1 and [...]"; here the right-hand side is 1, so by the manual the result is
unconditional. I have not proved that criterion myself. A certified run (`thueinit(P, 1)`, which
calls `bnfcertify` on the degree-11 field) was started; its outcome is recorded in Section 5.

Positive control (`control94.gp`): the same code applied to c = 94, where 421^2 - 94 = 3^11,
finds the solution for exactly one unit class (k = -4: (a, b) = (1, 0), t = 421; and its conjugate
at k = 4 with t = -421), and the local tests leave alive exactly the classes k = 4 and 7 = -4.
So the pipeline finds a solution when there is one and the local test does not kill a true class.

Hence **there is no prime p with tau(p^2) = 463, given the Thue computation.**

## 4. Consequence

Together with `../tau/PROOF.md`, Section 4: **tau(n) is never +-l for any odd prime l < 1000**,
subject to
- the literature used for the reduction (Balakrishnan-Craig-Ono, Balakrishnan-Craig-Ono-Tsai) and
  the classical congruences for tau;
- PARI Thue computations: for l = 463 the degree-11 equation above (right-hand side 1); for
  +827, -829, -643 the equations of degree 29, 41, 53 over real cyclotomic fields (uncertified
  `thueinit`, tentative class number 1); the degree-11 cyclotomic cases are certified.
Everything else (the sieve for the six tau(p^4) cases, the three tau(p^2) = -D cases, and
tau(p^2) = 601) is elementary and unconditional.

## 5. Verification log

- `local_k.gp` (PARI) and `local_k_check.py` (Python, own continued-fraction unit, brute-force
  11th powers in F_q[y]/(y^2 - c)) agree: 601: no k survives; 463: only k = 9 survives.
- Thue equation: three inputs, see Section 3. Certified run: see below.
- Small search: no odd x < 2*10^6 with x^11 + 463 or x^11 + 601 a square (referee of
  `../tau/PROOF.md`).

**Certified run: unfinished after 50 minutes.** `thueinit(P, 1)` for P(x) = Q(x, 1) (PARI/GP
2.15.4, one core, parisizemax 10 GB) printed nothing and was stopped by its 3000 s time limit
(`thue463_certified_run.txt`). So the 463 case rests on the uncertified runs A, B, C above (which
PARI's manual declares unconditional for a right-hand side of +-1). To resume, from this
directory, without a time limit:

    (echo '\r thue463.gp'; echo "S = setup(463, -2); P = subst(subst(S[1], 'a, 'x), 'b, 1); tnf = thueinit(P, 1); print(thue(tnf, 1)); print(thue(tnf, -1));") | gp -q > thue463_certified_run.txt 2>&1

The slow step is `bnfcertify` on the degree-11 field of discriminant -2^10 11^11 463^5 (tentative
regulator 7.6*10^9); `bnfcertify(bnfinit(P, 1))` alone is the same computation.

**Second certified attempt (coordinating session, 26 Nov 2026, 19:30–21:24 UTC):** the same
`thueinit(P, 1)` command, `parisizemax` 10 GB, one core, stopped by a 7000 s time limit without
output (total certified time now about 2 h 47 min, no result). The independent referee
(`REFEREE-2-TAU2.md`) estimates the full certification as plausibly infeasible here (Minkowski
bound ≈ 3.7·10¹⁰) and suggests certifying only the fundamental units, which is all a right-hand
side of ±1 needs. Status unchanged: the 463 case rests on the uncertified runs A, B, C.
