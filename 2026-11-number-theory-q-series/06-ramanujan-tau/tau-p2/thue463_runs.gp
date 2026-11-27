\\ thue463_runs.gp -- run as:  (echo '\r thue463.gp'; cat thue463_runs.gp) | gp -q
\\ Solves Q(a, b) = +-1 for the form of thue463.gp in several independent ways and writes
\\ thue463_results.txt.  FLAG = 0 (default here): uncertified thueinit; set FLAG = 1 for certification.
{
my(OUT = "thue463_results.txt", S, Q, P, R, Q2, P2, bnf, L, t0, FLAG = 0);
system(Str("rm -f ", OUT));
S = setup(463, -2); Q = S[1];
write(OUT, "# thue463_results.txt -- PARI/GP ", version(), ", written by thue463_runs.gp");
write(OUT, "basis: w1 = ", S[3], ", w2 = ", S[4], "   (y = sqrt 463)");
write(OUT, "Q(a,b) = ", Q);
P = subst(subst(Q, 'a, 'x), 'b, 1);
bnf = bnfinit(P, 1);
write(OUT, "degree-11 field: irreducible ", polisirreducible(P), ", signature ", bnf.sign, ", disc ", bnf.disc, " = ", factor(bnf.disc), ", tentative h = ", bnf.no, ", tentative R = ", precision(bnf.reg, 19));
t0 = getabstime();
write(OUT, "run A (P(x) = Q(x,1), flag ", FLAG, "):   Q = 1: ", thue(thueinit(P, FLAG), 1), "   Q = -1: ", thue(thueinit(P, FLAG), -1), "   [", getabstime() - t0, " ms]");
R = subst(subst(Q, 'a, 1), 'b, 'x); t0 = getabstime();
write(OUT, "run B (P(x) = Q(1,x), flag ", FLAG, "):   Q = 1: ", thue(thueinit(R, FLAG), 1), "   Q = -1: ", thue(thueinit(R, FLAG), -1), "   [", getabstime() - t0, " ms]");
Q2 = subst(Q, 'a, 'a + 'b);      \\ unimodular change of variables (a, b) -> (a + b, b)
P2 = subst(subst(Q2, 'a, 'x), 'b, 1); t0 = getabstime();
write(OUT, "run C (Q(a+b, b), flag ", FLAG, "):   = 1: ", thue(thueinit(P2, FLAG), 1), "   = -1: ", thue(thueinit(P2, FLAG), -1), "   [", getabstime() - t0, " ms]");
L = List();
for (b = -500, 500, my(Qb = subst(Q, 'b, b)); for (a = -500, 500, if (abs(subst(Qb, 'a, a)) == 1, listput(L, [a, b]))));
write(OUT, "brute force |a|, |b| <= 500 with Q = +-1: ", Vec(L));
}
