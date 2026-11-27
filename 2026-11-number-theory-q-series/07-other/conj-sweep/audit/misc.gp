\\ AUDIT: small checks. (1) periodicity (c/d') = (c/d) for d' = d mod c, 4 | c, used in the direct proof of (F4);
\\ (2) the Remark after Theorem A; (3) arithmetic constants; (4) residue classes; (5) index of Gamma_1(256);
\\ (6) the eta-quotient/E_2 identities of (F8) and E_2 = E_{l+1}, E_{l-1} = 1 mod l.
main() =
{ my(bad = 0, D, ph, f, phm, R, E);
  forstep (c = -400, 400, 4, if (c, forstep (d = -801, 801, 2, if (gcd(c, d) == 1, for (k = -3, 3, if (kronecker(c, d) != kronecker(c, d + k*c), bad++))))));
  print("(1) violations of kronecker(c,d) = kronecker(c,d+kc), 4|c, |c| <= 400, |d| <= 801 odd: ", bad);
  D = 2000; ph = sum(n = -50, 50, x^(n^2)) + O(x^D); phm = subst(ph, x, -x); f = sum(n = 1, D-1, if (n % 2, sigma(n)*x^n)) + O(x^D);
  print("(2) phi(q)^4 - phi(-q)^4 = 16 F: ", ph^4 - phm^4 == 16*f);
  R = f*ph^7*phm^4; print("    coefficients of F phi(q)^7 phi(-q)^4 at n = 3 mod 8, n < ", D, " all zero (exactly): ", sum(n = 0, D-1, (n % 8 == 3) && polcoeff(R, n) != 0) == 0);
  print("    12 F phi^7 (phi^4 + F) = 12 F phi^7 phi(-q)^4 mod 17: ", (12*f*ph^7*(ph^4 + f) - 12*R)*Mod(1,17) == 0);
  print("(3) 11^6 mod 17 = ", 11^6 % 17, "; 13^9 mod 23 = ", 13^9 % 23, "; 17^358 mod 19 = ", lift(Mod(17,19)^358), "; (-1/11)^7 = ", kronecker(-1,11)^7, "; (-1/13)^10 = ", kronecker(-1,13)^10, "; (-1/17)^359 = ", kronecker(-1,17)^359);
  print("    Sturm bounds: ", [8, 11, 360]*49152/12, "; 9+18*20 = ", 9+18*20, ", 9+9*20 = ", 9+9*20, ", 369-189 = ", 369-189, " = 10*18; dim Lambda_369 = ", 738\4+1, ", dim Lambda_{719/2} = ", 719\4+1);
  print("(4) A: {m mod 88 : m = 3 (8), (m/11) = -1} = ", select(m -> m % 8 == 3 && kronecker(m, 11) == -1, [0..87]), "  (source: 19,35,43,51,83)");
  print("    B: {l mod 104 : l = 5 (8), (l/13) = +1} = ", select(m -> m % 8 == 5 && kronecker(m, 13) == 1, [0..103]), "  (source: 29,53,61,69,77,101)");
  print("    C: {k mod 56 : k = 3 (8), (k/7) = +1} = ", select(m -> m % 8 == 3 && kronecker(m, 7) == 1, [0..55]), "  (source: 11,43,51)");
  print("    D: number of classes mod 2584 = 8*17*19: ", #select(m -> m % 8 == 3 && kronecker(m, 17) == -1 && kronecker(m, 19) == 1, [0..2583]), "  (source: 72)");
  print("(5) [SL2(Z):Gamma_1(256)] = 256^2 * prod(1 - p^-2) = ", 256^2*(1 - 1/4));
  D = 600; E = ((k, d) -> 1 - 2*k/bernfrac(k)*sum(n = 1, (D-1)\d, sigma(n, k-1)*x^(d*n)) + O(x^D));
  ph = sum(n = -30, 30, x^(n^2)) + O(x^D); f = sum(n = 1, D-1, if (n % 2, sigma(n)*x^n)) + O(x^D);
  print("(6) Theta F / F = (4 E2(q^4) - E2(q^2))/3: ", x*f'/f == (4*E(2,4) - E(2,2))/3, ";  Theta phi / phi = (10 E2(q^2) - 2 E2(q) - 8 E2(q^4))/24: ", x*ph'/ph == (10*E(2,2) - 2*E(2,1) - 8*E(2,4))/24);
  print("    E_20 = E_2 mod 19: ", (E(20,1) - E(2,1))*Mod(1,19) == 0, ";  E_18 = 1 mod 19: ", (E(18,1) - 1)*Mod(1,19) == 0, " (both 19-integral)");
}
main();
quit
