WZ certificates (plain text, PARI/GP and sympy readable; ^ = power).
For each theorem: sg = f(n,k+1)/f(n,k)  ("sigma" in RESULTS.md; renamed because sigma is a GP built-in), rho = f(n+1,k)/f(n,k), R = certificate, g = R*f.
WZ equation  f(n+1,k) - f(n,k) = g(n,k+1) - g(n,k)  <=>  rho - 1 == subst(R,k,k+1)*sg - R   (rational identity).
Check in GP, e.g.:   read("T1.gp"); print(rho - 1 - (subst(R,k,k+1)*sg - R) == 0)
For T4, T8, T10 the pair is (n*h, Q*h):  (n+1)*rho - n == subst(Q,k,k+1)*sg - Q.
For T7 and T11 the pair lives on the triangle k < n; besides the generic identity there is the boundary row k = n-1:
   subst(rho,k,n-1) - 1 == -subst(rho,k,n-1)*subst(subst(sigma,k,n-1),n,n+1)... (see proofs/verify_T7*, verify_T11*).
T12 (Guillera's shifted series, z = -1/48): standard form; T12.gp also checks the compatibility of the two ratios.
The proofs (boundary terms, limits) are in ../RESULTS.md; exact checks in ../proofs/verify_T*.py.
