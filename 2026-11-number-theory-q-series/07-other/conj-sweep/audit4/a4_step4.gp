\\ AUDIT-4 part 4 (own code, PARI/GP exact integers, direct binomial() calls -- no recurrences).
\\ Re-does Step 4 of R10 with the window given by the monotonicity criterion, which is re-derived in AUDIT-4.md:
\\   R_n(a-1) < R_n(a)  <=>  g(n,a) := (n-a+1)^2 (n-2a+1) - a(a-1)(n-2a+3) > 0     (2 <= a < n/2).
\\ For each n: V = { a in [2,amax] : g(n,a) <= 0 } is computed completely (every a), amin = min V;
\\ R_n is strictly increasing on [1, amin-1], so max_{1<=a<n/2} R_n = max over [amin-1, amax].
\\ On that window the exact integer inequality  (n-2a+1) a 2(2n-1) C(n,a)^2 < n(n-1) C(2n,n)  is tested.
\\ Also records whether V is an interval ending at amax (unimodality of a -> R_n(a)).
default(parisize, 256000000);
g(n,a) = (n-a+1)^2*(n-2*a+1) - a*(a-1)*(n-2*a+3);
{
N0 = 3; N1 = 16100;
gt = List(); eq = List(); nonuni = List(); wmax = 0; totw = 0; nearest = [0, 0, 0.];
for(n = N0, N1,
  amax = (n-1)\2;
  amin = 0; uni = 1;
  for(a = 2, amax, if(g(n,a) <= 0, if(amin == 0, amin = a), if(amin, uni = 0)));
  if(!uni, listput(nonuni, n));
  lo = if(amin, amin-1, amax); if(lo < 1, lo = 1);
  rhs = n*(n-1)*binomial(2*n, n);
  isgt = 0; iseq = 0; best = 0;
  for(a = lo, amax,
    lhs = (n-2*a+1)*a*2*(2*n-1)*binomial(n,a)^2;
    if(lhs > rhs, isgt = 1); if(lhs == rhs, iseq = 1); if(lhs > best, best = lhs);
  );
  totw += amax-lo+1; if(amax-lo+1 > wmax, wmax = amax-lo+1);
  if(isgt, listput(gt, n)); if(iseq, listput(eq, n));
  if(n >= 496, r = best*1./rhs; if(r > nearest[3], nearest = [n, 0, r]));
);
gt = Vec(gt); eq = Vec(eq);
print("n range ", N0, "..", N1, " ; window pairs tested exactly: ", totw, " ; largest window: ", wmax);
print("n for which {a: R_n(a-1) >= R_n(a)} is NOT an interval [amin, amax]: ", Vec(nonuni));
print("#n with rho(n) > 1: ", #gt, " ; largest: ", vecmax(gt));
print("n with rho(n) = 1: ", eq);
print("rho(n) > 1 for all 5 <= n <= 375: ", sum(i = 1, #gt, gt[i] >= 5 && gt[i] <= 375) == 371);
S = Set(gt); print("n in 376..495 with rho(n) <= 1: ", select(n -> !setsearch(S, n), [376..495]));
print("rho(495) > 1: ", setsearch(S, 495) > 0, " ; any n >= 496 with rho(n) >= 1: ", #select(n -> n >= 496, concat(gt, eq)) > 0);
print("largest rho(n) for 496 <= n <= ", N1, ": n = ", nearest[1], ", rho = ", nearest[3]);
}
\\ a few exact values of rho near the threshold
rho(n) = vecmax(vector((n-1)\2, a, (n-2*a+1)/(n-a)*binomial(n,a)*binomial(n-2,a-1)/(binomial(2*n-2,n-1)/n)));
print("rho(494), rho(495), rho(496), rho(497), rho(498) = ", vector(5, i, rho(493+i)*1.));
print("rho(2000), rho(16100) = ", rho(2000)*1., " ", rho(16100)*1.);
quit
