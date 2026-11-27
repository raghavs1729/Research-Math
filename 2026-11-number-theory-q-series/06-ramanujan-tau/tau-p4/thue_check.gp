\\ thue_check.gp -- certified (thueinit flag 1: no GRH) solution of F_{d-1}(X, Y) = c,
\\ F_{d-1}(X, Y) = prod_{k=1}^{(d-1)/2} (Y - 4 cos^2(pi k/d) X),  tau(p^(d-1)) = F_{d-1}(p^11, tau(p)^2).
\\ PARI convention: thue(thueinit(P), c) lists (x, y) with y^deg(P) * P(x/y) = c; here x = Y, y = X.
\\ usage: gp -q, then  \r thue_check.gp ; run(d, c)
Fpol(d) = subst(polsubcyclo(d, (d - 1)/2), 'x, 'x - 2);
\\ sanity: F_{d-1}(X, t^2) against the Hecke/Lucas recurrence u_1 = 1, u_2 = t, u_{j+1} = t u_j - X u_{j-1}
chk(d) = { my(P = Fpol(d), t = 7, X = 3, u0 = 1, u1 = t, u2, n = (d - 1)/2);
  for (j = 2, d - 1, u2 = t*u1 - X*u0; u0 = u1; u1 = u2);
  X^n*subst(P, 'x, t^2/X) == u1 };
run(d, c) = { my(P = Fpol(d), tnf, S, T = getabstime());
  if (!chk(d), error("form mismatch"));
  tnf = thueinit(P, 1);
  S = thue(tnf, c);
  print("d = ", d, "  c = ", c, "  certified solutions (Y, X): ", S, "   [", getabstime() - T, " ms]");
  S };
\\ uncertified variant (thueinit flag 0).  By the PARI documentation of thue(), the output is then
\\ unconditional when the tentative class number is 1; the tentative class number and regulator are printed.
run0(d, c) = { my(P = Fpol(d), bnf, tnf, S, T = getabstime());
  if (!chk(d), error("form mismatch"));
  bnf = bnfinit(P, 1);
  tnf = thueinit(P);
  S = thue(tnf, c);
  print("d = ", d, "  c = ", c, "  flag 0: tentative h = ", bnf.no, ", R = ", bnf.reg, ";  solutions (Y, X): ", S, "   [", getabstime() - T, " ms]");
  S };
