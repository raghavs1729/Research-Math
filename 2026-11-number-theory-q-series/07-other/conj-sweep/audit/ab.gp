\\ AUDIT (independent of the session's code): Theorems A and B of R02.
\\ phi and F are built from eta-quotients (not from squares / divisor sums);
\\ coordinates are obtained by solving an OVERDETERMINED linear system (matsolve on a square part,
\\ then checked on many further coefficients), all in exact integer arithmetic; h is expanded mod ell.
default(parisizemax, 8000000000);
default(seriesprecision, 16);
et(d, D) = eta(x^d + O(x^D));     \\ prod (1-x^(dn)) to O(x^D)
PHI(D) = et(2,D)^5/(et(1,D)^2*et(4,D)^2);
FF(D)  = x*et(4,D)^8/et(2,D)^4;
kr(a,p) = kronecker(a,p);
\\ coefficient vector (length n) of series s
cf(s, n) = vector(n, i, polcoeff(s, i-1));
\\ coordinates of the form with coefficient vector v (length L >= dim) in basis F^j phi^(K-4j), by matsolve
\\ on the first dim rows; then verify on all L rows.  Returns coordinates (exact rationals/integers).
coords(v, K, D) =
{ my(dim = K\4+1, ph = PHI(D), f = FF(D), L = #v, B, M, c);
  B = vector(dim, j, cf(f^(j-1)*ph^(K-4*(j-1)) + O(x^D), L));
  M = matrix(L, dim, i, j, B[j][i]);
  c = matsolve(matrix(dim, dim, i, j, M[i,j]), vector(dim, i, v[i])~);
  if (M*c != v~, error("form is NOT in the span on the extended range"));
  c~;
}
chkc0(ell, ga, K2, dim2) =
{ my(Dp = 3000*ell, pb = (et(2,Dp)*Mod(1,ell))/(et(1,Dp)*Mod(1,ell))^2, c0l = sum(j=1, dim2, ga[j]*FF(3000)^(j-1)*PHI(3000)^(K2-4*(j-1))) + O(x^3000));
  for (n = 0, 2999, if (polcoeff(c0l, n) - (-1)^n*polcoeff(pb, ell*n) != 0, error("c0 != U_ell(1/phi) at ", n)));
  print("check: c0(n) = (-1)^n pbar(ell n) mod ell for n < 3000 (pbar from eta(2tau)/eta(tau)^2): OK");
}
doit(ell, p, r, s) =
{ my(lam = (ell-3)/2, K1 = ell-1, dim1 = K1\4+1, L1, D1, ph, f, w1, chi, T, al, ga, K2, dim2, L2, D2, c0, eps, b, be, mu, et_, Bd, D, h, bad, hs);
  print("## ell=", ell, " p=", p, " r=", r, " s=", s, " lam=", lam);
  \\ step 1: T_ell(phi^(ell-1)), checked on L1 = dim1+40 coefficients
  L1 = dim1 + 40; D1 = ell*L1 + 2; ph = PHI(D1); f = ph^(ell-1);
  w1 = (ell-1)/2; chi = kr(-4, ell)^w1;
  T = vector(L1, i, my(n=i-1); polcoeff(f, ell*n) + if(n % ell == 0, chi*ell^(w1-1)*polcoeff(f, n/ell), 0));
  al = coords(T, K1, D1);
  print("alpha mod ell = ", al % ell, "   (integral: ", denominator(al) == 1, "; expansion verified on ", L1, " coefficients)");
  \\ direct check: U_ell(1/phi) = sum pbar(ell n)(-1)^n q^n  vs  sum_j alpha_j F^j phi^(ell-2-4j)  mod ell
  if (K1 % 4 == 0, if (al[dim1] % ell != 0, error("top alpha nonzero")); ga = vector(dim1-1, j, al[j] % ell), ga = al % ell);
  print("gamma = ", ga);
  K2 = ell-2; dim2 = K2\4+1; if (#ga != dim2, error("dim"));
  L2 = dim2 + 25; D2 = p^2*L2 + 2; ph = PHI(D2); f = FF(D2);
  c0 = sum(j=1, dim2, ga[j]*f^(j-1)*ph^(K2-4*(j-1))) + O(x^D2);
  \\ independent check of step 1: c0(n) = (-1)^n pbar(ell n) mod ell, pbar from eta(2)/eta(1)^2
  chkc0(ell, ga, K2, dim2);
  eps = kr(-1, p)^lam;
  b = vector(L2, i, my(n=i-1); polcoeff(c0, p^2*n) + eps*kr(n,p)*p^(lam-1)*polcoeff(c0, n) + if(n % (p^2) == 0, p^(2*lam-1)*polcoeff(c0, n/p^2), 0));
  be = coords(b, K2, D2);
  print("beta mod ell = ", be % ell, "  eps = (-1/p)^lam = ", eps, "   (integral: ", denominator(be) == 1, "; expansion verified on ", L2, " coefficients)");
  mu = s*eps*p^(lam-1) % ell;
  print("s*eps = ", s*eps, ", p^(lam-1) mod ell = ", lift(Mod(p, ell)^(lam-1)));
  et_ = (be - mu*ga) % ell;
  print("eta = ", et_);
  Bd = (lam+1)*49152/12; D = Bd + 10;
  ph = PHI(D)*Mod(1, ell); f = FF(D)*Mod(1, ell);
  h = sum(j=1, dim2, et_[j]*f^(j-1)*ph^(K2-4*(j-1))) + O(x^D);
  bad = 0; forstep(n = r, Bd, 8, if (polcoeff(h, n) != 0, bad++));
  print("Sturm bound ", Bd, "; number of n = ", r, " mod 8, n <= bound: ", #[n | n <- [0..Bd], n % 8 == r], "; nonzero h(n): ", bad);
  print("nonzero counts n<2000 by class mod 8: ", vector(8, i, sum(n=0, 1999, (n % 8 == i-1) && polcoeff(h,n) != 0)));
  \\ dump h for the cross-check against direct pbar values (python)
  write(Str("h_", ell, ".txt"), lift(vector(Bd+1, i, polcoeff(h, i-1))));
  if (bad, print("RESULT ell=", ell, ": FAILED"), print("RESULT ell=", ell, ": all zero up to the Sturm bound"));
}
doit(17, 11, 3, -1);
doit(23, 13, 5, 1);
quit
