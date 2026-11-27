"""sw_yoshida.py -- tests of the numerically stated identities of Y. Yoshida, arXiv:2608.08780, for the (E8,T1) system
(Nahm sums with matrix C(E8)), and of Sun-Wang (3.21).

Claims tested (Yoshida (2.68)-(2.73)):  for h in {0, 5/11, 8/11, 10/11, 12/11}
     q^{gamma_h} sum_{m in N^8} q^{m^T C(E8) m/2 - Q_h.m}/prod (q)_{m_i}  =  chi_h of T_10 M_eff(11,2),
(2.72): the displayed expansions;  and the bilinear relation (1.9)/(1.10)
     sum_h chi_h chi^vee_{sigma(h)} = E_4/eta^8.

The characters of T_10 M_eff(11,2) are computed from the Harvey-Wu Hecke formula for conductor N = 33:
     b^{(p)}_i(M) = sum_{a | gcd(p,M)} (p/a) eps_a(i) b_{pi_a(i)}(p M / a^2),      p = 10,
where b_i(n) are the coefficients of the M_eff(11,2) characters in q^{n/33}, pi_a is forced by the support
(n_{pi_a(i)} = a^{-2} n_i mod 33) and eps_a(i) = +-1 are the signs of the Galois action on the S-matrix,
S_{s,s'} = (-1)^{s+s'} (zeta^{ss'} - zeta^{-ss'})/sqrt(-11):  sigma_a(S)_{s,s'} = (a/11) (-1)^{s+[as]} sgn(as) S_{[as],s'}.
The script tries eps from sigma_a and from sigma_{a^{-1}} and reports which one (if any) matches.

usage: python3 sw_yoshida.py N
"""
import sys, time
from fractions import Fraction as Fr
from sw_lib import Ctx, Ser, compare, report
from fusion import inv, cartanE

N = int(sys.argv[1]) if len(sys.argv) > 1 else 12
C8 = cartanE(8)          # chain 1-...-7, node 8 attached to node 5  (= Yoshida's labelling, his (2.11))
Q = {Fr(0): (0,) * 8, Fr(5, 11): (2, -2, 2, -2, 2, -1, 1, -2), Fr(8, 11): (1, -1, 1, -1, 2, -2, 2, -2),
     Fr(10, 11): (0, 0, 1, -1, 1, -1, 1, -1), Fr(12, 11): (1, -1, 1, -1, 1, -1, 1, -1)}
gamma = {Fr(0): Fr(-10, 33), Fr(5, 11): Fr(104, 33), Fr(8, 11): Fr(80, 33), Fr(10, 11): Fr(20, 33), Fr(12, 11): Fr(26, 33)}
printed = {Fr(5, 11): (-3, [8, 245, 2360, 15100, 75600, 321552, 1210000]),
           Fr(8, 11): (-2, [35, 680, 5628, 33440, 159020, 652240, 2386140]),
           Fr(10, 11): (0, [84, 1240, 9570, 54160, 249880, 1001296, 3599710]),
           Fr(12, 11): (0, [50, 560, 4180, 22352, 100695, 393960, 1393820]),
           Fr(0): (0, [1, 120, 1660, 12320])}
Ci = inv(C8)
print("Q^vee = C^-1 Q:", {str(h): [str(sum(Ci[i][j] * Q[h][j] for j in range(8))) for i in range(8)] for h in Q})

# ---- M_eff(11,2) characters: chi_s = q^{n_s/33} prod_{n != 0,+-s mod 11} 1/(1-q^n),  n_s = 33 h_eff - 1
LB = 10 * (N + 12) + 40          # number of q-coefficients needed of the M(11,2) characters
cM = Ctx(1, LB)
heff = {5: Fr(0), 4: Fr(1, 11), 3: Fr(3, 11), 2: Fr(6, 11), 1: Fr(10, 11)}
nres = {s: int(33 * heff[s]) - 1 for s in heff}
bser = {s: cM.prod(lambda n, s=s: -1 if n % 11 not in (0, s, 11 - s) else 0).coeffs() for s in heff}


def b(s, n):
    """coefficient of q^{n/33} in chi_s of M_eff(11,2)"""
    if (n - nres[s]) % 33:
        return 0
    k = (n - nres[s]) // 33
    return bser[s].get(k, 0) if k >= 0 else 0


def leg(a):
    a %= 11
    return 1 if a in (1, 3, 4, 5, 9) else -1


def fold(x):
    x %= 11
    return (x, 1) if x <= 5 else (11 - x, -1)


def galois(a):
    """eps_a(s), pi_a(s) with sigma_a(S)_{s,.} = eps_a(s) S_{pi_a(s),.}"""
    out = {}
    for s in range(1, 6):
        t, sg = fold(a * s)
        out[s] = (leg(a) * (-1) ** (s + t) * sg, t)
    return out


def hecke10(variant):
    """returns dict s -> Ser (in q^{1/33}) : component of T_10 chi with exponents = 10 n_s mod 33"""
    c = Ctx(33, N + 4)
    res = {}
    for s in range(1, 6):
        cls = (10 * nres[s]) % 33
        d = {}
        M = cls - 33 * 2
        while M < 33 * (N + 4):
            val = 0
            for a in (1, 2, 5, 10):
                if M % a:
                    continue
                if a == 1:
                    val += 10 * b(s, 10 * M)
                    continue
                abar = pow(a, -1, 33)
                g = galois(a if variant == 0 else pow(a, -1, 11))
                eps, t = g[s]
                # support check: n_t = abar^2 n_s mod 33
                if (nres[t] - abar * abar * nres[s]) % 33:
                    # permutation from the support
                    cands = [u for u in range(1, 6) if (nres[u] - abar * abar * nres[s]) % 33 == 0]
                    assert len(cands) == 1 and cands[0] == t, ("support/galois mismatch", a, s, t, cands)
                val += (10 // a) * eps * b(t, 10 * M // (a * a))
            if val:
                d[Fr(M, 33)] = val
            M += 33
        res[s] = c.from_dict(d)
    return res


t0 = time.time()
c1 = Ctx(1, N)
nahm = {}
for h in Q:
    s = c1.nahm(C8, [-x for x in Q[h]])
    nahm[h] = s
    cs = s.coeffs()
    v0 = min(cs)
    ser = [cs.get(v0 + i, 0) for i in range(min(N - v0, 12))]
    pv, pc = printed[h]
    okp = (v0 == pv and ser[:len(pc)] == pc[:len(ser)])
    print("h=%-6s Nahm sum (C(E8), B=-Q_h): q^%d * %s ...   printed expansion (Yoshida (2.4)/(2.72)): %s   [%.0fs]" % (
        h, v0, ser[:9], "AGREES (%d coefficients)" % min(len(pc), len(ser)) if okp else "DISAGREES: printed q^%d %s" % (pv, pc), time.time() - t0), flush=True)

# sector h of T_10 <-> s with 10 n_s = 33(h - 10/33) mod 33
c33 = Ctx(33, N)
for variant in (0, 1):
    try:
        T = hecke10(variant)
    except AssertionError as e:
        print("variant %d: %s" % (variant, e)); continue
    print("--- Hecke images T_10 M_eff(11,2), signs from sigma_a%s :" % ("" if variant == 0 else "^{-1}"))
    allok = True
    for h in Q:
        s = [u for u in range(1, 6) if (10 * nres[u] - (33 * h - 10)) % 33 == 0][0]
        lhs_c = nahm[h].coeffs()
        lhs = c33.from_dict({Fr(k) + gamma[h]: v for k, v in lhs_c.items()}).cut(33 * N + int(33 * gamma[h]) if gamma[h] < 0 else 33 * N)
        ok = report("  chi_{%s} = q^{%s} * Nahm  [M(11,2) label s=%d]" % (h, gamma[h], s), lhs, T[s])
        allok = allok and ok
    print("   variant %d: %s" % (variant, "ALL FIVE AGREE" if allok else "mismatch"))

# ---- bilinear relation with the (T1,E8) characters: pairing (Yoshida (1.9)):
#   chi_0 chi^v_0 + chi_{5/11} chi^v_{6/11} + chi_{8/11} chi^v_{3/11} + chi_{10/11} chi^v_{1/11} + chi_{12/11} chi^v_{10/11}
pairing = {Fr(0): Fr(0), Fr(5, 11): Fr(6, 11), Fr(8, 11): Fr(3, 11), Fr(10, 11): Fr(1, 11), Fr(12, 11): Fr(10, 11)}
inv_heff = {v: k for k, v in heff.items()}
tot = 0
for h in Q:
    lhs_c = nahm[h].coeffs()
    chiA = c33.from_dict({Fr(k) + gamma[h]: v for k, v in lhs_c.items()})
    s = inv_heff[pairing[h]]
    chiB = c33.prod(lambda n, s=s: -1 if n % 11 not in (0, s, 11 - s) else 0).shift(Fr(nres[s], 33))
    tot = tot + chiA * chiB
# E_4/eta^8
sig3 = lambda n: sum(d ** 3 for d in range(1, n + 1) if n % d == 0)
E4 = c33.ps([1] + [240 * sig3(n) for n in range(1, N + 2)])
e = c33.euler_inv()
e8 = e * e; e8 = e8 * e8; e8 = e8 * e8
rhs = (E4 * e8).shift(Fr(-1, 3))
P = min(tot.P, rhs.P) - 33
report("bilinear: sum_h chi_h(Nahm) chi^vee = E_4/eta^8", tot.cut(P), rhs.cut(P))
