\\ control94.gp -- positive control for the pipeline: 421^2 - 94 = 3^11 in Q(sqrt 94) (h = 1).
\\ Run:  (echo '\\r thue463.gp'; echo '\\r local_k.gp'; cat control94.gp) | gp -q
\\ The solution must appear for exactly one unit exponent k, and that k must survive the local tests.
{
my(c = 94, K = bnfinit('y^2 - c, 1), S, P, sol, t, eps, u, v, good = List(), al);
print("control c = 94: h = ", K.no, ", bnfcertify = ", bnfcertify(K));
for (k = -5, 5,
  S = setup(c, k); P = subst(subst(S[1], 'a, 'x), 'b, 1);
  sol = thue(thueinit(P), 1);
  for (i = 1, #sol, t = subst(subst(S[2], 'a, sol[i][1]), 'b, sol[i][2]); listput(good, [k, sol[i], t, factor(t^2 - c)])));
print("control solutions [k, (a,b), t, factor(t^2 - c)]: ", Vec(good));
eps = lift(K.fu[1]); u = polcoef(eps, 0); v = polcoef(eps, 1); al = vector(11, i, 1);
forprime (q = 5, 700, if ((q^2 - 1) % 11 == 0 && c % q,
  for (k = 0, 10, if (if (kronecker(c, q) == 1, if (q % 11 == 1, cnt_split(c, u, v, k, q, vector(q, i, i - 1)), 1), cnt_inert(c, u, v, k, q, vector(q, i, i - 1))) == 0, al[k + 1] = 0))));
print("control: k = 0..10 alive after the local tests (no t-restriction at 23): ", al);
}
