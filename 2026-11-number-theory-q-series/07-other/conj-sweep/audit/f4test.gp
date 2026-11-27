\\ AUDIT: numerical test of (F4) [U_4 preserves M_{K/2}(Gamma~_0(4))] and of (F3) on the basis elements, in exact
\\ rational arithmetic: for every K <= KMAX and every basis element B_j = F^j phi^(K-4j), the series U_4 B_j
\\ (resp. T(p^2) B_j for odd K) is expanded in the basis using the first dim coefficients and the expansion is then
\\ compared on EXTRA further coefficients.  Negative control: U_2 does NOT preserve the space.
default(parisizemax, 4000000000);
et(d, D) = eta(x^d + O(x^D));
PHI(D) = et(2,D)^5/(et(1,D)^2*et(4,D)^2);
FF(D)  = x*et(4,D)^8/et(2,D)^4;
inspan(v, K, D) =
{ my(dim = K\4+1, ph = PHI(D), f = FF(D), L = #v, M, c);
  M = matrix(L, dim, i, j, polcoeff(f^(j-1)*ph^(K-4*(j-1)) + O(x^D), i-1));
  c = matsolve(matrix(dim, dim, i, j, M[i,j]), vector(dim, i, v[i])~);
  [M*c == v~, denominator(c) == 1];
}
main(KMAX, EXTRA) =
{ my(ok4 = 0, tot4 = 0, ok2 = [0,0], tot2 = [0,0], okT = 0, totT = 0, int4 = 1, dim, L, D, ph, f, B, r, lam, pp);
  for (K = 1, KMAX, dim = K\4 + 1; L = dim + EXTRA;
    D = 9*L + 2; ph = PHI(D); f = FF(D);
    for (j = 0, dim-1, B = f^j*ph^(K-4*j) + O(x^D);
      r = inspan(vector(L, i, polcoeff(B, 4*(i-1))), K, D); tot4++; ok4 += r[1]; int4 = int4 && r[2];
      r = inspan(vector(L, i, polcoeff(B, 2*(i-1))), K, D); tot2[K % 2 + 1]++; ok2[K % 2 + 1] += r[1];
      if (K % 2 && K >= 3, lam = (K-1)/2; pp = 3;
        r = inspan(vector(L, i, my(n = i-1); polcoeff(B, 9*n) + kronecker(-1,pp)^lam*kronecker(n,pp)*pp^(lam-1)*polcoeff(B, n) + if(n % 9 == 0, pp^(2*lam-1)*polcoeff(B, n/9), 0)), K, D);
        totT++; okT += r[1])));
  print("K <= ", KMAX, ", ", EXTRA, " extra coefficients each:");
  print("  U_4 B_j in span:    ", ok4, " of ", tot4, "  (integral coordinates: ", int4, ")");
  print("  T(3^2) B_j in span: ", okT, " of ", totT, "  (odd K >= 3, Shimura's formula as in F3)");
  print("  U_2 B_j in span (negative control):  even K (integer weight, where U_2 does act): ", ok2[1], " of ", tot2[1], ";  odd K (half-integral weight, U_2 changes the character by (2/.)): ", ok2[2], " of ", tot2[2]);
}
main(48, 40);
quit
