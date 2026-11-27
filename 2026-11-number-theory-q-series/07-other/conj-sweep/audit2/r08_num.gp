\\ Referee numerical tests of R08 Theorems 1, 2 (PARI; independent of code/radu_2color.py which uses flint).
default(parisizemax, 3000000000);
D = 30000; U = 2^12;
\\ pbar_l via theta series: 1/(phi(-q) phi(-q^l)), coefficients mod U
phim(k) = { my(v = vector(D), c0 = 1); for (n = 1, sqrtint((D-1)\k), v[k*n^2] = 2*(-1)^n); (Mod(1,U) + Ser(concat([0], v), q)*Mod(1,U)) + O(q^D); }
pb(l) = 1/(phim(1)*phim(l));
\\ sanity against the eta-quotient definition (exact integers), and pbar_3(3) = 10 (source)
{ for (l = 1, 12, my(e = eta(q^2 + O(q^400))*eta(q^(2*l) + O(q^400))/(eta(q + O(q^400))^2*eta(q^l + O(q^400))^2), s = pb(l));
    for (n = 0, 399, if (Mod(polcoeff(e, n), U) != polcoeff(s, n), error("theta vs eta ", l, " ", n)));
    if (l == 3 && polcoeff(e, 3) != 10, error("pbar_3(3)")));
  print("1/(phi(-q)phi(-q^l)) = f_2 f_{2l}/(f_1^2 f_l^2) to q^399 for l = 1..12; pbar_3(3) = 10: OK"); }
\\ Theorem 1
{ my(bad = 0, cnt = 0, ls = List());
  for (l = 1, 100, if (l % 4 <= 1, listput(ls, l); my(s = pb(l)); forstep (n = 3, D-1, 4, cnt++; if (lift(polcoeff(s, n)) % 8, bad++))));
  print("Theorem 1: pbar_l(4n+3) mod 8, l = 0,1 mod 4, l <= 100 (", #ls, " values of l), 4n+3 < ", D, ": ", cnt, " coefficients, nonzero: ", bad); }
{ for (l = 2, 11, if (l % 4 >= 2, my(s = pb(l), f = 0); forstep (n = 3, 2000, 4, if (lift(polcoeff(s, n)) % 8, f++)); print("   control l = ", l, ": nonzero pbar_l(4n+3) mod 8 for 4n+3 <= 2000: ", f))); }
\\ Theorem 2
{ my(s = pb(9), bad = 0, cnt = 0, minv = oo);
  for (al = 0, 14, for (n = 0, D, my(x = 2^al*(8*n+7)); if (x >= D, break); cnt++; my(c = lift(polcoeff(s, x))); if (c % 16, bad++); minv = min(minv, valuation(c + U*(c == 0), 2))));
  print("Theorem 2: pbar_9(2^a (8n+7)) mod 16, a = 0..14, argument < ", D, ": ", cnt, " coefficients, nonzero: ", bad, "; minimal 2-adic valuation seen: ", minv);
  my(f = 0); for (n = 0, 500, if (lift(polcoeff(s, 8*n+3)) % 16, f++)); print("   control: nonzero pbar_9(8n+3) mod 16, n <= 500: ", f);
  my(g = 0, s5 = pb(5)); for (n = 0, 500, if (lift(polcoeff(s5, 2*(8*n+7))) % 16, g++)); print("   control: nonzero pbar_5(2(8n+7)) mod 16, n <= 500: ", g); }
\\ the cubic-term parity claim of Theorem 2, tested directly: N(m) = number of positive solutions of the four ternary equations
{ my(bad = 0, cnt = 0);
  for (al = 0, 6, for (n = 0, 200, my(mm = 2^al*(8*n+7), Nn = 0); if (mm > 3000, break);
     for (a = 1, sqrtint(mm), for (b = 1, sqrtint(mm), my(r1 = mm - a^2 - b^2);
        if (r1 > 0 && issquare(r1), Nn++);                                  \\ a^2+b^2+c^2
        if (r1 > 0 && r1 % 9 == 0 && issquare(r1/9), Nn++);                 \\ a^2+b^2+9c^2
        my(r2 = mm - a^2 - 9*b^2); if (r2 > 0 && r2 % 9 == 0 && issquare(r2/9), Nn++);   \\ a^2+9b^2+9c^2
        my(r3 = mm - 9*a^2 - 9*b^2); if (r3 > 0 && r3 % 9 == 0 && issquare(r3/9), Nn++))); \\ 9(a^2+b^2+c^2)
     cnt++; if (Nn % 2, bad++)));
  print("Theorem 2, cubic terms: N(m) even for all ", cnt, " m = 2^a(8n+7) <= 3000: odd cases ", bad); }
quit
