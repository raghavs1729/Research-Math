\\ AUDIT (independent of the session's code): Theorem D of R02 (mod 19).
\\ Differences from overp19_proof.py: PARI instead of FLINT; phi, F, pbar from Euler pentagonal series / eta-quotients;
\\ W is taken directly from overpartition numbers (pbar = f_2/f_1^2), not from the basis expansion;
\\ coordinates in the basis F^j phi^(K-4j) are obtained as "Taylor coefficients" in the local parameter
\\ t = F/phi^4 = q + ... :      f = sum_j c_j F^j phi^(K-4j)   <=>   f/phi^K = sum_j c_j t^j ,
\\ by repeated division by t, and always to MORE than dim terms, so that the claim "f lies in the span"
\\ is itself tested (the coefficients of t^j for j >= dim must vanish).
\\ Everything is inside main(), so any error aborts the run (no RESULT line is printed after an error).
default(parisizemax, 9000000000);
ell = 19; p = 17;
pent(d, D) = { my(v = vector(D), k = 0, e); while (1, e = d*k*(3*k-1)/2; if (e >= D, break); v[e+1] += (-1)^k; e = d*k*(3*k+1)/2; if (k > 0 && e < D, v[e+1] += (-1)^k); k++); Ser(v*Mod(1,ell), x) + O(x^D); }
PHI(D) = pent(2,D)^5/(pent(1,D)^2*pent(4,D)^2);
FF(D)  = x*pent(4,D)^8/pent(2,D)^4;
PBAR(D) = pent(2,D)/pent(1,D)^2;
\\ first L coefficients c_0..c_{L-1} of the t-expansion of f/phi^K; f must be known to O(x^L)
tcoords(f, K, L) =
{ my(ph = PHI(L), tinv = ph^4*x/FF(L), G = (f + O(x^L))/ph^K, c = vector(L), c0);
  for (j = 1, L, c0 = polcoeff(G, 0); c[j] = lift(c0); G = ((G - c0)/x)*tinv);
  c;
}
\\ evaluate sum_j c_j F^j phi^(K-4j) to O(x^D)
evalform(c, K, D) = { my(ph = PHI(D), t = FF(D)/ph^4, acc = c[#c] + O(x^D)); forstep (j = #c-1, 1, -1, acc = acc*t + c[j]); acc*ph^K; }
kr(a,q) = kronecker(a,q);
sieve19(s, D) = Ser(vector(D, i, if (kr(i-1, ell) == 1, polcoeff(s, i-1), 0))*Mod(1, ell), x) + O(x^D);

main() =
{ my(D, ph, f, L1, A, Ap, u, T18, T9, w0, L3, D3, pb, W, w0s, lam, eps, c1, c2, bser, bet, etav, L5, D5, hs, gser, gco, Bd, D6, gs, bad, cnt, hd);
  \\ ---- sanity of phi, F against their definitions (squares, odd divisor sums)
  D = 3000; ph = PHI(D); f = FF(D);
  for (n = 0, D-1, if (polcoeff(ph,n) != if(n==0, 1, issquare(n), 2, 0), error("phi")); if (polcoeff(f,n) != if(n % 2, sigma(n), 0), error("F")));
  print("phi, F (eta-quotients) agree with sum q^(n^2), sum_{n odd} sigma(n) q^n mod 19 for n < ", D);
  \\ ---- steps 1-2
  L1 = 1200;
  A = PHI(L1)^18;
  Ap = sieve19(A, L1);
  u = tcoords(Ap, 738, L1);
  print("step 2: u_180..u_184 = ", u[181..185], ";  coefficients of t^j in A_+/phi^738 for 185 <= j < ", L1, " all zero: ", vecmax(u[186..L1]) == 0);
  if (vecmax(u[181..L1]), error("A_+ not in phi^19 * Lambda_{719/2}"));
  \\ Theta-lemma (F8) numerically: Theta^18 A in Lambda_369, Theta^9 A in Lambda_189 (mod 19)
  T18 = Ser(vector(L1, i, Mod(i-1, ell)^18*polcoeff(A, i-1)), x) + O(x^L1);
  T9  = Ser(vector(L1, i, Mod(i-1, ell)^9*polcoeff(A, i-1)), x) + O(x^L1);
  print("        (F8 test) t-coefficients j >= 185 of Theta^18 A/phi^738 vanish: ", vecmax(tcoords(T18, 738, L1)[186..L1]) == 0,
        ";  t-coefficients j >= 95 of Theta^9 A/phi^378 vanish: ", vecmax(tcoords(T9, 378, L1)[96..L1]) == 0);
  w0 = u[1..180];
  \\ ---- step 3
  L3 = 220; D3 = p^2*L3 + 10;
  pb = PBAR(D3);
  if (lift(Vec(pb + O(x^8))) != [1,2,4,8,14,24,40,64] % ell, error("pbar"));
  W = sieve19(subst(pb, x, -x), D3);
  w0s = evalform(w0, 719, D3);
  print("        W(n) = [(n/19)=1] (-1)^n pbar(n) equals the basis expansion w0 for all n < ", D3, ": ", W == w0s);
  if (W != w0s, error("W != w0"));
  lam = 359; eps = kr(-1, p)^lam; c1 = Mod(p, ell)^(lam-1); c2 = Mod(p, ell)^(2*lam-1);
  print("step 3: eps = ", eps, ", 17^358 = ", lift(c1), ", 17^717 = ", lift(c2), " mod 19");
  bser = Ser(vector(L3, i, my(n = i-1); polcoeff(W, p^2*n) + eps*kr(n,p)*c1*polcoeff(W, n) + if(n % p^2 == 0, c2*polcoeff(W, n/p^2), 0)), x) + O(x^L3);
  bet = tcoords(bser, 719, L3);
  print("        t-coefficients of b/phi^719 for 180 <= j < ", L3, " all zero (test of F3): ", vecmax(bet[181..L3]) == 0);
  if (vecmax(bet[181..L3]), error("b not in Lambda"));
  bet = bet[1..180];
  \\ ---- step 4
  etav = lift(Mod(bet, ell) + c1*w0);
  \\ ---- step 5
  L5 = 700; D5 = 4*L5 + 4;
  hs = evalform(etav, 719, D5);
  \\ h against its direct definition from W: h(n) = W(289 n) + 5((n/17)+1) W(n) + 17^717 W(n/289), where W is available
  hd = 0; for (n = 0, L3-1, if (polcoeff(hs, n) != polcoeff(W, p^2*n) + c1*(eps*kr(n,p) + 1)*polcoeff(W, n) + if(n % p^2 == 0, c2*polcoeff(W, n/p^2), 0), hd++));
  print("step 4: h (basis expansion) agrees with (T(17^2) + 5) W coefficientwise for n < ", L3, ": ", hd == 0);
  if (hd, error("h"));
  gser = Ser(vector(L5, i, polcoeff(hs, 4*(i-1))), x) + O(x^L5);
  gco = tcoords(gser, 719, L5);
  print("step 5: t-coefficients of (U_4 h)/phi^719 for 180 <= j < ", L5, " all zero (numerical test of F4): ", vecmax(gco[181..L5]) == 0);
  if (vecmax(gco[181..L5]), error("U_4 h not in Lambda"));
  gco = gco[1..180];
  write("d_coords.txt", "u=", w0); write("d_coords.txt", "beta=", bet); write("d_coords.txt", "eta=", etav); write("d_coords.txt", "g=", gco);
  \\ ---- step 6
  Bd = 360*49152/12; D6 = Bd + 10;
  gettime();
  gs = evalform(gco, 719, D6);
  print("step 6: g expanded to O(q^", D6, ") in ", gettime()\1000, " s");
  if (type(gs) != "t_SER" || serprec(gs, x) < D6, error("precision"));
  bad = 0; cnt = 0; forstep (n = 3, Bd, 8, cnt++; if (polcoeff(gs, n) != 0, bad++));
  print("        Sturm bound ", Bd, "; n = 3 mod 8, n <= bound: ", cnt, " values; nonzero g(n): ", bad);
  print("        nonzero g(n), n < 4000, by n mod 8: ", vector(8, i, sum(n = 0, 3999, (n % 8 == i-1) && polcoeff(gs, n) != 0)));
  print("        g(n) = h(4n) for n < ", L5, ": ", sum(n = 0, L5-1, polcoeff(gs, n) != polcoeff(hs, 4*n)) == 0);
  write("d_g.txt", lift(vector(200001, i, polcoeff(gs, i-1))));
  print(if (bad, "RESULT: FAILED", "RESULT: g(n) = 0 mod 19 for all n = 3 mod 8 up to the Sturm bound"));
}
main();
quit
