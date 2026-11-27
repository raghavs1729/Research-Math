\\ local_k.gp -- tau(p^2) = c, c in {463, 601}:  t^2 - c = p^11,  t + sqrt c = eps^k g^11  (k mod 11).
\\ For each k and each prime q with 11 | q^2 - 1 (q != p), count the t mod q for which
\\ (t + sqrt c) eps^(-k) is an 11th power in (O_K/q)^*.  If there is none, k is impossible.
\\ At q = 23 the Ramanujan congruence forces t = 2 mod 23 (see TAU2.md), so only t = 2 is tested there.

\\ number of admissible t for (c, k, q); eps = u + v sqrt(c)
cnt_split(c, u, v, k, q, T) = {
  my(s = lift(sqrt(Mod(c, q))), e1 = Mod(u + v*s, q)^k, e2 = Mod(u - v*s, q)^k, n = 0, z, zz);
  foreach (T, t,
    z = Mod(t + s, q)/e1; zz = Mod(t - s, q)/e2;
    if (z != 0 && zz != 0 && z^((q - 1)/11) == 1 && zz^((q - 1)/11) == 1, n++));
  n
};
cnt_inert(c, u, v, k, q, T) = {
  my(P = Mod(1, q)*('y^2 - c), ee = Mod(Mod(1, q)*(u + v*'y), P)^k, n = 0, z);
  foreach (T, t,
    z = Mod(Mod(1, q)*(t + 'y), P)/ee;
    if (z^((q^2 - 1)/11) == 1, n++));
  n
};
dok(c, qmax) = {
  my(K = bnfinit('y^2 - c, 1), eps, u, v, alive = vector(11, i, 1), killer = vector(11, i, 0), T, n, spl);
  if (bnfcertify(K) != 1 || K.no != 1, error("class group"));
  eps = lift(K.fu[1]); u = polcoef(eps, 0); v = polcoef(eps, 1);
  print("c = ", c, ": h = 1 (bnfcertify), eps = ", eps, ", norm ", u^2 - c*v^2);
  forprime (q = 3, qmax,
    if ((q^2 - 1) % 11 == 0 && q != c,
      if (issquare(q^11 + c), error("p = q is a solution"));
      spl = (kronecker(c, q) == 1);
      if (!spl || q % 11 == 1,
        T = if (q == 23, [2], vector(q, i, i - 1));
        for (k = 0, 10,
          n = if (spl, cnt_split(c, u, v, k, q, T), cnt_inert(c, u, v, k, q, T));
          if (n == 0 && alive[k + 1], alive[k + 1] = 0; killer[k + 1] = q)))));
  print("  k = 0..10 alive: ", alive, "   first killing prime: ", killer);
  [alive, killer]
};
