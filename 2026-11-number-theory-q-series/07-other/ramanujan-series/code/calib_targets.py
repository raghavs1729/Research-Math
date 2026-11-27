"""classical published identities used for calibration (same format as survey/verify_*.py)"""
from mpmath import mp, mpf, binomial, pi, zeta, log, sqrt, catalan, nsum, inf
E = []
def add(id_, start, summand, rhs, src): E.append(dict(id=id_, start=start, summand=summand, rhs=rhs, src=src))
def C(n, k): return binomial(n, k)
add("C1", 0, lambda k: (6*k+1)*C(2*k,k)**3/mpf(256)**k, lambda: 4/pi, "Ramanujan 1914; WZ proof: Guillera 2006")
add("C2", 0, lambda k: (42*k+5)*C(2*k,k)**3/mpf(4096)**k, lambda: 16/pi, "Ramanujan 1914; WZ proof: Guillera 2002/2006")
add("C3", 0, lambda k: (20*k**2+8*k+1)*C(2*k,k)**5/mpf(-4096)**k, lambda: 8/pi**2, "Guillera 2002 (WZ)")
add("C4", 0, lambda k: (4*k+1)*C(2*k,k)**3/mpf(-64)**k, lambda: 2/pi, "Bauer 1859; Ekhad-Zeilberger 1994 (WZ)")
ENTRIES = E
