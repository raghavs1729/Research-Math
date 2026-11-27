"""AUDIT a5. (i) The convention of the source (eq. (3) of arXiv:2308.11829v2) tested on the IDENTIFIED rows of its
Table 1, to make sure that 'PCF(a,b)' in THEOREM.md is the quantity tabulated in the paper.
(ii) The value printed in the paper for our row is '1.20426...', the true limit is 1.2042534...: find the depth n at
which the convergent p_n/q_n shows 1.20426 (the paper's value is a finite-depth convergent of a 1/n-convergent PCF)."""
from mpmath import mp, mpf, zeta, pi, nstr, floor
mp.dps = 40
def pcf(a, b, N, marks=()):
    pm, p, qm, q = mpf(1), mpf(a(0)), mpf(0), mpf(1); out = {}
    for n in range(1, N + 1):
        pm, p = p, a(n) * p + b(n) * pm; qm, q = q, a(n) * q + b(n) * qm
        s = abs(q); pm /= s; p /= s; qm /= s; q /= s
        if n in marks: out[n] = p / q
    return p / q, out
z3, z5 = zeta(3), zeta(5)
rows = [("n^5+(n+1)^5+6(n^3+(n+1)^3), -n^10", lambda n: n**5 + (n + 1)**5 + 6 * (n**3 + (n + 1)**3), lambda n: -n**10,
         [("2/(2z5+6z3-9)", 2 / (2 * z5 + 6 * z3 - 9))]),
        ("... - 4(2n+1), -n^10", lambda n: n**5 + (n + 1)**5 + 6 * (n**3 + (n + 1)**3) - 4 * (2 * n + 1), lambda n: -n**10,
         [("2/(2z5-2z3+1) [Table 1]", 2 / (2 * z5 - 2 * z3 + 1)), ("2/(2z5-2z3-1) [Table 4]", 2 / (2 * z5 - 2 * z3 - 1))]),
        ("n^5+(n+1)^5+16(n^3+(n+1)^3)-4(2n+1), -n^10", lambda n: n**5 + (n + 1)**5 + 16 * (n**3 + (n + 1)**3) - 4 * (2 * n + 1), lambda n: -n**10,
         [("64/(64z5+176z3-273)", 64 / (64 * z5 + 176 * z3 - 273))])]
for name, a, b, cands in rows:
    val, _ = pcf(a, b, 3000)
    print("%-48s depth 3000: %s" % (name, nstr(val, 15)))
    for cn, cv in cands: print("        %-28s = %s   diff %s" % (cn, nstr(cv, 15), nstr(val - cv, 3)))
a = lambda n: 8 * (n**5 + (n + 1)**5) - 15 * (n**3 + (n + 1)**3) + 9 * (2 * n + 1)
b = lambda n: -64 * n**10
L = 96 / (372 * z5 - pi**5)
marks = (30, 1000, 10000, 20000, 43000, 50000, 62000, 72000, 100000, 110000, 120000)
_, out = pcf(a, b, 120000, marks)
print("U1 limit = %s ; paper prints 1.20426..." % nstr(L, 12))
for n in marks: print("   depth %6d: %s" % (n, nstr(out[n], 10)))
