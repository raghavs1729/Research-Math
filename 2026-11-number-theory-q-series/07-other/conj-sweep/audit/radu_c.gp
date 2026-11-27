\\ AUDIT (independent of code/radu.py): Theorem C of R02 via the Radu-Sellers lemma, written from the statement
\\ in arXiv:2302.11830 (Section "An algorithmic approach by Radu and Sellers").
default(parisizemax, 4000000000);
m = 896; MM = 2; N = 56; t = 176; u = 7;
rd = [1, 2]; rv = [-2, 1];                       \\ r_1 = -2, r_2 = 1  (f_2/f_1^2)
ad = divisors(N); av = [1829, -915, 7143, 0, -256, -14284, 128, -49993, 99986];
\\ divisors(56) = [1,2,4,7,8,14,28,56]: re-map the certificate's vector  {1:1829, 2:-915, 4:7143, 7:-256, 8:-14284, 14:128, 28:-49993, 56:99986}
av = [1829, -915, 7143, -256, -14284, 128, -49993, 99986];
if (ad != [1,2,4,7,8,14,28,56], error("divisors"));
kap = gcd(m^2 - 1, 24);
print("kappa = ", kap);
sdr = sum(i = 1, #rd, rd[i]*rv[i]);  sr = vecsum(rv);
\\ ---- Delta^* conditions (a)-(f)
ca = 1; fordiv(m, q, if (isprime(q) && N % q, ca = 0));
cb = 1; for (i = 1, #rd, if (rv[i] != 0 && (m*N) % rd[i], cb = 0));
cc = (kap*N*sum(i = 1, #rd, rv[i]*m*N/rd[i])) % 24 == 0;
cd = (kap*N*sr) % 8 == 0;
ce = N % (24*m/gcd(-24*kap*t - kap*sdr, 24*m)) == 0;
pr = prod(i = 1, #rd, rd[i]^abs(rv[i])); s2 = valuation(pr, 2); jj = pr/2^s2;
cf = if (m % 2, 1, ((kap*N) % 4 == 0 && (s2*N) % 8 == 0) || (s2 % 2 == 0 && ((1-jj)*N) % 8 == 0));
print("Delta^* conditions (a)..(f): ", [ca, cb, cc, cd, ce, cf], "   [24m/gcd(...) = ", 24*m/gcd(-24*kap*t - kap*sdr, 24*m), ", s = ", s2, ", j = ", jj, "]");
\\ ---- P_{m,r}(t)
P = Set(); for (y = 1, 24*m, if (gcd(y, 24*m) == 1, my(s = y^2 % (24*m)); if ((s-1) % 24, error("s != 1 mod 24")); P = setunion(P, [(t*s + (s-1)/24*sdr) % m])));
print("P_{m,r}(t) = ", P, ";  {16k : k in 11,43,51} = ", [176, 688, 816]);
\\ ---- cusp inequality
pmr(a, c) = { my(best = oo, v); for (l = 0, m-1, v = sum(i = 1, #rd, rv[i]*gcd(rd[i]*a + rd[i]*kap*l*c, m*c)^2/(rd[i]*m))/24; if (v < best, best = v)); best; }
pst(c) = sum(i = 1, #ad, av[i]*gcd(ad[i], c)^2/ad[i])/24;
\\ (1) a genuine complete set of representatives of the cusps of Gamma_0(N): a/c, c | N, a mod gcd(c, N/c), gcd(a,c)=1
worst = oo; ncusp = 0;
fordiv(N, c, my(g = gcd(c, N/c), seen = Set()); for (a = 1, c, if (gcd(a, c) == 1 && !setsearch(seen, a % g), seen = setunion(seen, [a % g]); ncusp++; my(v = pmr(a, c) + pst(c)); print("  cusp ", a, "/", c, ": p_mr = ", pmr(a, c), ", p* = ", pst(c), ", sum = ", v); worst = min(worst, v))));
print("number of cusps of Gamma_0(56) used: ", ncusp, " (expected sum_{c|N} phi(gcd(c,N/c)) = ", sumdiv(N, c, eulerphi(gcd(c, N/c))), "); min of p_mr + p* = ", worst);
\\ (2) the session's superset: all (a,c), c | N, 1 <= a <= c, gcd(a,c) = 1
worst2 = oo; n2 = 0; fordiv(N, c, for (a = 1, c, if (gcd(a, c) == 1, n2++; worst2 = min(worst2, pmr(a, c) + pst(c)))));
print("superset of ", n2, " pairs (a,c): min = ", worst2);
\\ (3) brute force over ALL first columns (a,c), 1 <= c <= 120, 0 <= a < c (or a = 1 for c = 1), and c = 0, a = 1:  tests the double-coset invariance too
worst3 = oo; n3 = 0; for (c = 1, 120, for (a = 0, c-1, if (gcd(a, c) == 1, n3++; worst3 = min(worst3, pmr(a, c) + pst(c)))));
worst3 = min(worst3, pmr(1, 0) + pst(0));
print("brute force over ", n3+1, " first columns (c <= 120 and the identity): min = ", worst3);
\\ ---- nu
idx = N*prod(i = 1, #factor(N)~, 1 + 1/factor(N)[i,1]);
nu = ((sr + vecsum(av))*idx - sum(i = 1, #ad, ad[i]*av[i]))/24 - sdr/(24*m) - vecmin(P)/m;
print("[Gamma:Gamma_0(56)] = ", idx, ", nu = ", nu, ", floor(nu) = ", floor(nu));
\\ ---- the finite check, pbar from f_2/f_1^2 mod 7 (PARI eta)
B = floor(nu); D = m*(B+1) + 2;
pb = (eta(x^2 + O(x^D))*Mod(1,u))/(eta(x + O(x^D))*Mod(1,u))^2;
if (lift(Vec(pb + O(x^8))) != [1,2,4,8,14,24,40,64] % u, error("pbar"));
bad = 0; cnt = 0; for (i = 1, #P, for (n = 0, B, cnt++; if (polcoeff(pb, m*n + P[i]) != 0, bad++)));
print("pbar(896 n + t') mod 7, t' in P, 0 <= n <= ", B, ": ", cnt, " values, nonzero: ", bad);
\\ beyond the bound, as far as the series goes; and a negative control (another residue)
print("negative control: nonzero among pbar(896 n + 16*3), n <= 200: ", sum(n = 0, 200, polcoeff(pb, m*n + 48) != 0));
quit
