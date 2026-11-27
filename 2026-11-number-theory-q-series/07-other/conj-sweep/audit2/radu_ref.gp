\\ Referee (audit2): generic Radu-Sellers finite check, written from the statement of the lemma in arXiv:2302.11830,
\\ Section "An algorithmic approach by Radu and Sellers" (Definition of Delta^*, p_{m,r}, p^*_a, Lemma "finite check").
\\ Independent of code/radu.py.  The vectors a are taken from the session's certificates and RE-VERIFIED here.
default(parisizemax, 2000000000);
\\ p_{m,r}(gamma) and p^*_a(gamma) for gamma with first column (a, c)
PMR(rd, rv, m, kap, a, c) = { my(best = oo, v); for (l = 0, m-1, v = sum(i = 1, #rd, rv[i]*gcd(rd[i]*a + rd[i]*kap*l*c, m*c)^2/(rd[i]*m))/24; if (v < best, best = v)); best; }
PST(ad, av, c) = sum(i = 1, #ad, av[i]*gcd(ad[i], c)^2/ad[i])/24;
\\ rd, rv : divisors delta of M with r_delta (others 0);  ad = divisors(N), av = a-vector in that order.
radu(name, M, rd, rv, m, N, t, u, av) =
{
  my(ad = divisors(N), kap = gcd(m^2-1, 24), sdr, sr, ca, cb, cc, cd, ce, cf, pr, s2, jj, P, idx, nu, B, D, ser, bad, cnt, w1, w2, w3, n1, ncusp, At);
  print("## ", name);
  if (#av != #ad, error("a-vector length"));
  for (i = 1, #rd, if (M % rd[i], error("delta does not divide M")));
  if (t < 0 || t >= m, error("t range"));
  sdr = sum(i = 1, #rd, rd[i]*rv[i]); sr = vecsum(rv);
  \\ Delta^*
  ca = 1; fordiv(m, q, if (isprime(q) && N % q, ca = 0));
  cb = 1; for (i = 1, #rd, if (rv[i] != 0 && (m*N) % rd[i], cb = 0));
  cc = (kap*N*sum(i = 1, #rd, rv[i]*m*N/rd[i])) % 24 == 0;
  cd = (kap*N*sr) % 8 == 0;
  At = 24*m/gcd(-24*kap*t - kap*sdr, 24*m);
  ce = (N % At == 0);
  pr = prod(i = 1, #rd, rd[i]^abs(rv[i])); s2 = valuation(pr, 2); jj = pr/2^s2;
  cf = if (m % 2, 1, ((kap*N) % 4 == 0 && (s2*N) % 8 == 0) || (s2 % 2 == 0 && ((1-jj)*N) % 8 == 0));
  print("   M = ", M, ", r = ", vector(#rd, i, [rd[i], rv[i]]), ", m = ", m, ", t = ", t, ", N = ", N, ", kappa = ", kap, ", u = ", u);
  print("   Delta^* (a)..(f): ", [ca, cb, cc, cd, ce, cf], "   [24m/gcd = ", At, ", s = ", s2, ", j = ", jj, "]");
  \\ P_{m,r}(t): s runs over the squares of units mod 24m
  P = Set();
  for (y = 1, 24*m, if (gcd(y, 24*m) == 1, my(s = y^2 % (24*m)); if ((s-1) % 24, error("s != 1 mod 24")); P = setunion(P, [(t*s + (s-1)/24*sdr) % m])));
  print("   P_{m,r}(t) = ", P);
  \\ cusp inequality
  w1 = oo; ncusp = 0;
  fordiv(N, c, my(g = gcd(c, N/c), seen = Set()); for (a = 1, c, if (gcd(a, c) == 1 && !setsearch(seen, a % g), seen = setunion(seen, [a % g]); ncusp++; w1 = min(w1, PMR(rd, rv, m, kap, a, c) + PST(ad, av, c)))));
  w2 = oo; n1 = 0; fordiv(N, c, for (a = 1, c, if (gcd(a, c) == 1, n1++; w2 = min(w2, PMR(rd, rv, m, kap, a, c) + PST(ad, av, c)))));
  w3 = PMR(rd, rv, m, kap, 1, 0) + PST(ad, av, 0); for (c = 1, 3*N, for (a = 0, c-1, if (gcd(a, c) == 1, w3 = min(w3, PMR(rd, rv, m, kap, a, c) + PST(ad, av, c)))));
  print("   a = ", vector(#ad, i, [ad[i], av[i]]));
  print("   min of p_mr + p*_a: over ", ncusp, " cusps of Gamma_0(N) [expected ", sumdiv(N, c, eulerphi(gcd(c, N/c))), "]: ", w1, ";  over the ", n1, " pairs a/c, c | N: ", w2, ";  over all first columns with c <= 3N: ", w3);
  idx = N*prod(i = 1, #factor(N)~, 1 + 1/factor(N)[i,1]);
  nu = ((sr + vecsum(av))*idx - sum(i = 1, #ad, ad[i]*av[i]))/24 - sdr/(24*m) - vecmin(P)/m;
  B = floor(nu);
  print("   [Gamma:Gamma_0(N)] = ", idx, ", nu = ", nu, ", floor(nu) = ", B);
  D = m*(B+3) + 5;
  ser = prod(i = 1, #rd, eta(x^rd[i] + O(x^D))^rv[i]);       \\ exact integers
  bad = 0; cnt = 0; for (i = 1, #P, for (n = 0, B, cnt++; if (polcoeff(ser, m*n + P[i]) % u != 0, bad++)));
  print("   c_r(m n + t') mod u for t' in P, 0 <= n <= ", B, ": ", cnt, " coefficients, nonzero: ", bad, "   [first coefficients of the series: ", Vec(ser + O(x^6)), "]");
  my(ok = ca && cb && cc && cd && ce && cf && w1 >= 0 && w2 >= 0 && w3 >= 0 && bad == 0 && B >= 0);
  print("   VERDICT: ", if (ok, "all hypotheses verified, congruence holds for all n", "FAILED"));
  print();
  [ok, ser];
}
print("=== R08, Theorem 3 (arXiv:2607.16608), pbar_l = f_2 f_{2l}/(f_1^2 f_l^2) ===");
\\ divisors(8)=[1,2,4,8]; divisors(24)=[1,2,3,4,6,8,12,24]; divisors(12)=[1,2,3,4,6,12]
radu("pbar_8(32n+28) = 0 mod 32", 16, [1,2,8,16], [-2,1,-2,1], 32, 8, 28, 32, [72,-36,-30,60]);
radu("pbar_6(24n+17) = 0 mod 64", 12, [1,2,6,12], [-2,1,-2,1], 24, 24, 17, 64, [61,-30,-16,9,8,-18,-27,54]);
radu("pbar_6(24n+23) = 0 mod 64", 12, [1,2,6,12], [-2,1,-2,1], 24, 24, 23, 64, [61,-30,-16,9,8,-18,-27,54]);
radu("pbar_9(36n+33) = 0 mod 128", 18, [1,2,9,18], [-2,1,-2,1], 36, 12, 33, 128, [89,-41,-27,-7,3,22]);
radu("pbar_4(36n+21) = 0 mod 9", 8, [1,2,4,8], [-2,1,-2,1], 36, 12, 21, 9, [101,-50,-33,-1,15,3]);
print("=== R03a (arXiv:2611.05623), A_k = (f_6/(f_1 f_2 f_3))^k ===");
R = radu("A_6(12n+10) = 0 mod 3^5", 6, [1,2,3,6], [-6,-6,-6,6], 12, 6, 10, 243, [181,-98,-75,78]);
a22 = polcoeff(R[2], 22); print("   A_6(22) = ", a22, ",  v_3 = ", valuation(a22, 3), ",  A_6(22) mod 729 = ", a22 % 729);
print("   v_3(A_6(12n+10)), n = 0..10: ", vector(11, n, valuation(polcoeff(R[2], 12*(n-1)+10), 3))); print();
radu("A_3(12n+7) = 0 mod 3^4", 6, [1,2,3,6], [-3,-3,-3,3], 12, 12, 7, 81, [85,-36,-23,-4,5,10]);
radu("A_3(12n+11) = 0 mod 3^4", 6, [1,2,3,6], [-3,-3,-3,3], 12, 12, 11, 81, [85,-36,-23,-4,5,10]);
radu("A_4(6n+5) = 0 mod 3^3", 6, [1,2,3,6], [-4,-4,-4,4], 6, 6, 5, 27, [62,-36,-42,60]);
radu("A_7(12n+11) = 0 mod 3^3", 6, [1,2,3,6], [-7,-7,-7,7], 12, 12, 11, 27, [200,-91,-54,-2,21,7]);
radu("A_8(18n+15) = 0 mod 3^5", 6, [1,2,3,6], [-8,-8,-8,8], 18, 6, 15, 243, [329,-119,-123,70]);
print("=== negative control: a false congruence must fail the coefficient check ===");
radu("pbar_6(24n+17) = 0 mod 128 (FALSE?)", 12, [1,2,6,12], [-2,1,-2,1], 24, 24, 17, 128, [61,-30,-16,9,8,-18,-27,54]);
quit
