\\ thue463.gp -- the one unit exponent that survives the local tests for c = 463 (k = 9 = -2 mod 11):
\\      t + sqrt(463) = eta^(-2) g^11,   eta = |eps| the fundamental unit > 1,  g in Z[sqrt 463].
\\ Write g = a w1 + b w2 with (w1, w2) a Z-basis of Z[sqrt 463] chosen (by LLL on the lattice
\\ x -> (eta^(-2/11) x, eta^(2/11) x')) so that the binary form
\\      Q(a, b) := coefficient of sqrt(463) in eta^(-2) (a w1 + b w2)^11
\\ has small coefficients.  The Thue equation is Q(a, b) = 1.  The basis only affects the size of the
\\ coefficients; Q is computed exactly.
default(parisizemax, 10*10^9);
setup(c, k) = {
  my(K = bnfinit('y^2 - c, 1), eps, eta, A, s, G, M, w1, w2, g, Q, e);
  if (bnfcertify(K) != 1 || K.no != 1, error("class group"));
  eps = K.fu[1];
  default(realprecision, 200);
  s = sqrt(c); e = subst(lift(eps), 'y, s); eta = abs(e);
  \\ embedding of eps^k: (e^k, e'^k) with e' = norm/e ; use absolute values (signs are absorbed in g)
  A = abs(e)^(k/11);
  G = [A, A*s; 1/A, -s/A];                 \\ columns = images of 1 and sqrt(c)
  M = qflll(G);                             \\ unimodular 2x2 integer matrix
  if (abs(matdet(M)) != 1, error("not unimodular"));
  w1 = M[1,1] + M[2,1]*'y; w2 = M[1,2] + M[2,2]*'y;
  g = Mod(lift(eps), 'y^2 - c)^k * Mod('a*w1 + 'b*w2, 'y^2 - c)^11;
  Q = polcoef(lift(g), 1, 'y);
  [Q, polcoef(lift(g), 0, 'y), w1, w2, M]
};
