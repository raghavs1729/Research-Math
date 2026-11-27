\\ tau_p4.gp -- independent PARI/GP implementation (shares no code with sieve_tau_p4.py).
\\
\\ Part A.  For alpha and an exponent e, all integer solutions of U^2 - 3UV + V^2 = alpha are
\\          (Z_n, Z_{n+1}), n in Z, for finitely many "seed" sequences Z_{n+1} = 3 Z_n - Z_{n-1}.
\\          For every seed and every residue class n mod M the script finds a prime q, with the
\\          period of Z mod q dividing M, such that Z_n is not an e-th power residue mod q
\\          (code P) or Z_{n+1} is not a square mod q (code S).  It writes the witness table.
\\ Part B.  Classical congruences for tau: which (eps*l, d) with l < 1000 an odd prime and d an odd
\\          prime dividing l(l^2-1) are NOT excluded by congruences for tau(p^(d-1)).
\\ Part C.  tau(p^2) = -827: the Lebesgue-Nagell argument in Q(sqrt(-827)) (class number 7).
\\
\\ Run:  gp -q tau_p4.gp      (writes certificate_tau_p4.txt in the current directory)

OUT = "certificate_tau_p4.txt";
wr(s) = write(OUT, s);

\\ ---------- Part A ----------
seeds(al) = {
  my(B0 = if (al > 0, sqrtint(al\2), sqrtint(-al)), L = List(), D, s);
  for (U = -B0, B0,
    D = 5*U^2 + 4*al;
    if (D >= 0 && issquare(D, &s),
      if ((3*U + s) % 2 == 0, listput(L, [U, (3*U + s)/2]));
      if (s != 0 && (3*U - s) % 2 == 0, listput(L, [U, (3*U - s)/2]))));
  L = Vec(L);
  for (i = 1, #L, if (L[i][1]^2 - 3*L[i][1]*L[i][2] + L[i][2]^2 != al, error("bad seed")));
  [B0, L]
};

\\ period of n -> (Z_n, Z_{n+1}) mod q for every seed: the order of the matrix [0,1;-1,3] mod q
matperiod(q) = {
  my(A = Mod([0, 1; -1, 3], q), P = A, k = 1);
  while (P != matid(2), P = P*A; k++);
  k
};

isepow(z, e, q) = { z = z % q; z == 0 || Mod(z, q)^((q - 1)/gcd(e, q - 1)) == 1 };
issq(z, q) = { z = z % q; z == 0 || kronecker(z, q) == 1 };

\\ returns the list of surviving [seed index, n]; writes the witness table
sieve(al, e, M, qmax) = {
  my(sd = seeds(al), B0 = sd[1], S = sd[2], qs = List(), T, surv = List(), z, w, found, zz, nn);
  forprime (q = 3, qmax, if (q != 5, T = matperiod(q); if (M % T == 0, listput(qs, [q, T]))));
  qs = Vec(qs);
  wr(Str("# alpha = ", al, "  exponent e = ", e, "  modulus M = ", M, "  seed bound |U| <= ", B0));
  wr(Str("# seeds (Z_0, Z_1): ", S));
  wr(Str("# primes q <= ", qmax, " with period T_q | M, as [q, T_q]: ", qs));
  for (i = 1, #S,
    \\ exact integers Z_n for -2M <= n <= 2M+1, to confirm periodicity in both directions
    zz = vector(4*M + 2); nn = 2*M + 1;  \\ zz[n + nn] = Z_n
    zz[nn] = S[i][1]; zz[nn + 1] = S[i][2];
    for (n = 2, 2*M + 1, zz[n + nn] = 3*zz[n - 1 + nn] - zz[n - 2 + nn]);
    forstep (n = -1, -2*M, -1, zz[n + nn] = 3*zz[n + 1 + nn] - zz[n + 2 + nn]);
    for (j = 1, #qs,
      for (n = -2*M, 2*M + 1,
        if ((zz[n + nn] - zz[(n % qs[j][2]) + nn]) % qs[j][1], error("period check failed"))));
    wr(Str("seed ", S[i], ":  n mod M, witness q, code (P: Z_n not an e-th power residue mod q; S: Z_{n+1} not a square mod q), Z_n mod q, Z_{n+1} mod q"));
    for (n = 0, M - 1,
      found = 0;
      for (j = 1, #qs,
        if (!found,
          z = zz[n + nn] % qs[j][1]; w = zz[n + 1 + nn] % qs[j][1];
          if (!isepow(z, e, qs[j][1]),
            found = 1; wr(Str(n, " ", qs[j][1], " P ", z, " ", w)),
            if (!issq(w, qs[j][1]),
              found = 1; wr(Str(n, " ", qs[j][1], " S ", z, " ", w))))));
      if (!found,
        listput(surv, [S[i], n]);
        wr(Str(n, " SURVIVES  Z_n = ", zz[n + nn], "  Z_{n+1} = ", zz[n + 1 + nn])))));
  wr(Str("# survivors: ", Vec(surv)));
  wr("");
  Vec(surv)
};

\\ ---------- Part B ----------
\\ residues of tau(p^(d-1)) mod m allowed by the classical congruences, p prime not dividing 2m
allowed(d, m) = {
  my(R = List(), a, t, u0, u1, u2);
  if (m == 691, for (a = 1, 690, listput(R, lift(sum(i = 0, d - 1, Mod(a, 691)^(11*i))))));
  if (m == 2048, forstep (a = 1, 2047, 2, listput(R, lift(sum(i = 0, d - 1, Mod(a, 2048)^(11*i))))));
  if (m == 729, for (a = 1, 728, if (a % 3, listput(R, lift(Mod(a, 729)^(-610*(d - 1))*sum(i = 0, d - 1, Mod(a, 729)^(1231*i)))))));
  if (m == 125, for (a = 1, 124, if (a % 5, listput(R, lift(Mod(a, 125)^(-30*(d - 1))*sum(i = 0, d - 1, Mod(a, 125)^(71*i)))))));
  if (m == 7, for (a = 1, 6, listput(R, lift(Mod(a, 7)^(d - 1)*sum(i = 0, d - 1, Mod(a, 7)^(3*i))))));
  if (m == 23,
    for (a = 1, 22,
      foreach (if (kronecker(a, 23) == -1, [0], [2, -1]), t,
        u0 = Mod(1, 23); u1 = Mod(t, 23);
        for (j = 2, d - 1, u2 = t*u1 - Mod(a, 23)^11*u0; u0 = u1; u1 = u2);
        listput(R, lift(u1)))));
  Set(Vec(R))
};

\\ exact tau(p^(d-1)) by the Hecke recurrence
taupow(p, d) = { my(t = ramanujantau(p), u0 = 1, u1 = t, u2); for (j = 2, d - 1, u2 = t*u1 - p^11*u0; u0 = u1; u1 = u2); if (d == 1, 1, u1) };

MODS = [691, 2048, 729, 125, 7, 23];

validate() = {
  my(v, bad = 0);
  \\ (i) the congruences themselves, on all n <= 3000 in the stated classes
  for (n = 1, 3000,
    v = ramanujantau(n);
    if ((v - sigma(n, 11)) % 691, bad++);
    if (n % 8 == 1 && (v - sigma(n, 11)) % 2048, bad++);
    if (n % 3 == 1 && Mod(v, 729) != Mod(n, 729)^(-610)*sigma(n, 1231), bad++);
    if (n % 5 && Mod(v, 125) != Mod(n, 125)^(-30)*sigma(n, 71), bad++);
    if ((v - n*sigma(n, 3)) % 7, bad++));
  forprime (p = 2, 3000, if (p != 23,
    v = ramanujantau(p) % 23;
    if (kronecker(p, 23) == -1, if (v != 0, bad++), if (v != 2 && v != 22, bad++))));
  \\ (ii) the residue sets against exact values tau(p^(d-1))
  forprime (p = 3, 120, if (p != 3 && p != 5 && p != 7 && p != 23,
    foreach ([3, 5, 7, 11, 13, 23], d,
      v = taupow(p, d);
      for (k = 1, #MODS, if (!setsearch(allowed(d, MODS[k]), v % MODS[k]), bad++)))));
  bad
};

partB() = {
  my(Q = nextprime(10^30), ds, al, alive, killers, surv = List(), nord = 0);
  wr("# Part B: pairs (eps*l, d), l < 1000 an odd prime with l not dividing tau(l), d an odd prime dividing l^2-1,");
  wr("# NOT excluded by the congruences mod 691, 2^11, 3^6, 5^3, 7, 23 for tau(p^(d-1)), p not in {3,5,7,23,691};");
  wr("# for p in {3,5,7,23,691} the inequality tau(p^(d-1)) != eps*l is checked directly (mod a 31-digit prime).");
  forprime (l = 3, 999,
    if (ramanujantau(l) % l == 0, nord++; wr(Str("# non-ordinary prime (l | tau(l)): ", l)); next);
    ds = select(d -> d > 2, factor(l^2 - 1)[, 1]~);
    foreach ([1, -1], eps,
      al = eps*l;
      foreach (ds, d,
        foreach ([3, 5, 7, 23, 691], p, if (Mod(taupow(p, d), Q) == Mod(al, Q), error("tau(p^(d-1)) = alpha ?!")));
        alive = 1;
        for (k = 1, #MODS, if (alive && !setsearch(allowed(d, MODS[k]), al % MODS[k]), alive = 0));
        if (alive, listput(surv, [al, d])))));
  wr(Str("# survivors [eps*l, d]: ", Vec(surv)));
  wr("");
  Vec(surv)
};

\\ ---------- Part C ----------
\\ t^2 + D = p^11 (D prime, p an odd prime != D), i.e. tau(p^2) = -D.  If 11 does not divide the class
\\ number of Q(sqrt(-D)), then t + sqrt(-D) = g^11 with g an integer of the field.
\\ D = 3 mod 4: g = (a + b sqrt(-D))/2, a = b mod 2, and comparing coefficients of sqrt(-D):
\\      2^11 = b * F(a, b),   F(a,b) = sum_{k=0}^{5} binomial(11, 2k+1) a^(10-2k) (-D)^k b^(2k);
\\ D = 1 mod 4: g = a + b sqrt(-D) and 1 = b * F(a, b).
\\ In both cases b = +-2^j; all integer roots a are found exactly (nfroots over Q).
partC(D) = {
  my(disc = if (D % 4 == 3, -D, -4*D), h = qfbclassno(disc), sols = List(), b, r, rhs = if (D % 4 == 3, 2^11, 1), pol, K = bnfinit('y^2 + D, 1));
  if (bnfcertify(K) != 1 || K.no != h || K.disc != disc, error("class number not certified"));
  for (j = 0, if (D % 4 == 3, 11, 0), foreach ([1, -1], s,
    b = s*2^j;
    pol = b*sum(k = 0, 5, binomial(11, 2*k + 1)*'x^(10 - 2*k)*(-D)^k*b^(2*k)) - rhs;
    r = nfroots(, pol);
    for (i = 1, #r, if (denominator(r[i]) == 1, listput(sols, [r[i], b])))));
  wr(Str("# Part C, D = ", D, ": field discriminant ", disc, ", class number ", h, " (bnfcertify), 11 | h: ", h % 11 == 0,
         "; integer solutions (a, b) of b*F(a,b) = ", rhs, ": ", Vec(sols)));
  [D, h, Vec(sols)]
};

\\ ---------- main ----------
system(Str("rm -f ", OUT));
wr("# certificate_tau_p4.txt -- generated by tau_p4.gp (PARI/GP); see PROOF.md");
wr("");
s1 = sieve(461, 11, 88, 1000);
s2 = sieve(-599, 11, 22, 1000);
s3 = sieve(281, 11, 22, 1000);
s4 = sieve(-919, 11, 22, 1000);
s4a = sieve(-139, 11, 22, 1000);
s4b = sieve(-179, 11, 22, 1000);
s5 = sieve(33, 15, 4, 100);
s6 = sieve(33, 25, 4, 100);
s7 = sieve(55, 15, 40, 1000);
s8 = sieve(55, 25, 84, 1000);
bad = validate();
wr(Str("# validation of the congruences and residue sets: ", bad, " violations"));
wr("");
sB = partB();
sC = [partC(827), partC(227), partC(137)];
print("survivors: 461: ", s1, "; -599: ", s2, "; 281: ", s3, "; -919: ", s4, "; -139: ", s4a, "; -179: ", s4b);
print("33 (e=15): ", s5, "; 33 (e=25): ", s6);
print("55 (e=15): ", s7, "; 55 (e=25): ", s8);
print("validation violations: ", bad);
print("Part B survivors: ", sB);
print("Part C [D, h, solutions]: ", sC);
quit
