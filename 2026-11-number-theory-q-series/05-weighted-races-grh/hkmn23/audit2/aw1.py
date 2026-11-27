"""Second audit, Theorem 1/2 for selected pairs.
(1) constants n_e, n_o, B, C, E_R, I_q from my own characters and functional-equation constants (al1.Char);
(2) c_1 from my own sieve (numpy integers), checked at the maximising primes in Arb;
(3) U(t), X_0, U(1e12 +-) re-implemented from the formulas of section 2.4 as I re-derived them;
(4) W_1 < 0 certified in Arb (python-flint) on every piece between consecutive primes = a mod q up to T,
    by the concave-quadratic maximum M_i, plus the exact location of the true maximum of each piece;
(5) prime-side F(log x) on a grid (formula (F)) to test |F| <= L_1 (a wrong constant shows at small x).
Usage: python3 aw1.py q a T
"""
import sys, math, json
import numpy as np
import mpmath as mp
import flint
from flint import arb
from sympy import divisors, primefactors, totient
from achars import primitive_characters
from al1 import Char

flint.ctx.prec = 128
mp.mp.dps = 30


def sieve(n):
    s = np.ones(n + 1, dtype=bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return np.nonzero(s)[0]


def all_chars(q):
    """primitive characters inducing the characters mod q (one per character mod q)."""
    out = []
    for f in divisors(q):
        if f == 1:
            z = Char([0], 1); z.a = 0; out.append(z)
        else:
            out += [Char(v, f) for v in primitive_characters(f)]
    return out


def chibar(ch, a):
    if ch.f == 1:
        return mp.mpf(1)
    return mp.conj(ch.c[a % ch.f])


def constants(q, a, chars):
    ne = mp.mpf(0); no = mp.mpf(0)
    B = -mp.log(2 * mp.pi); C = mp.zeta(-1, 1, 1) / mp.zeta(-1)
    for ch in chars:
        if ch.f == 1:
            continue
        cb = chibar(ch, a)
        b, b2 = ch.fe_constants()
        if ch.a == 0:
            ne += mp.re(cb); B += mp.re(cb * (1 - b)); C += mp.re(cb * b2)
        else:
            no += mp.re(cb); B -= mp.re(cb * b); C += mp.re(cb * (1 + b2))
    return int(mp.nint(ne)), int(mp.nint(no)), B, C


def main(q, a, T):
    from math import gcd
    phi = int(totient(q))
    R = [b for b in range(1, q + 1) if gcd(b, q) == 1 and (b * b - a) % q == 0] if q > 1 else [0]
    N = len(R)
    chars = all_chars(q)
    L1 = sum(ch.S() for ch in chars)
    m = mp.mpf(2) / 3 * N - L1
    ne, no, B, C = constants(q, a, chars)
    # R-sums
    neR = noR = 0; BR = CR = mp.mpf(0)
    for b in R:
        x = constants(q, b, chars); neR += x[0]; noR += x[1]; BR += x[2]; CR += x[3]
    l3 = mp.log(3)
    ER = abs(neR) + abs(BR) / l3 + abs(noR) / 3 + abs(CR) / (3 * l3) + mp.mpf('0.6') * (N + abs(neR)) / (9 * l3) + mp.mpf('0.2') * abs(noR) / (27 * l3)
    Iq = sum(int(totient(q // p ** (q.bit_length() and _v(q, p)))) for p in primefactors(q))
    print(f'({q},{a}) phi={phi} N={N} R={R} L1={mp.nstr(L1,15)} m={mp.nstr(m,12)}')
    print(f'  n_e={ne} n_o={no} B={mp.nstr(B,10)} C={mp.nstr(C,10)}  n_e^R={neR} n_o^R={noR} B^R={mp.nstr(BR,10)} C^R={mp.nstr(CR,10)} E_R={mp.nstr(ER,10)} I_q={Iq}')

    P = sieve(int(T * 1.01) + 1000)
    # ---- c_1 : sup_{0<y<=1e6} (N y - Theta(y))/sqrt(y); Theta(y) = phi * sum_{p<=y, p^2=a} log p
    Y1 = 10 ** 6
    PR = P[P <= Y1]
    PR = PR[np.isin(PR % q, R)] if q > 1 else PR
    Th_before = phi * np.concatenate([[0.0], np.cumsum(np.log(PR.astype(float)))])  # Theta just before PR[i] is Th_before[i]
    cand_y = np.concatenate([PR.astype(float), [float(Y1)]])
    vals = (N * cand_y - Th_before) / np.sqrt(cand_y)
    i = int(np.argmax(vals))
    # certify that point in Arb
    yA = arb(int(cand_y[i]))
    ThA = sum((arb(int(p)).log() for p in PR[:i]), arb(0)) * phi
    c1A = (N * yA - ThA) / yA.sqrt()
    c1 = float(vals[i])
    print(f'  c_1 = {c1:.12f} at y -> {int(cand_y[i])}^-  (Arb {c1A.str(15)})   K_A = (4/5)c_1 = {0.8*c1:.6f}')

    # ---- U(t)
    l4 = math.log(4); l2 = math.log(2)
    L1f = float(L1); mf = float(m)
    Bf, Cf, ERf = float(B), float(C), float(ER)

    def rplus(t):
        return (max(0, -ne) * t ** -0.5 * math.log(t) + max(0, Bf + 2) * t ** -0.5 + max(0, no) * t ** -1.5 * math.log(t)
                + max(0, Cf - 2) * t ** -1.5 + max(0, -no) * 0.2 * t ** -3.5)

    def lam(Y):
        return (N * (L1f + 2 / 3 * phi * l4) + N * phi * l4 * (0.75 * Y ** (-1 / 6) + math.log(Y) / l2 * Y ** -0.25)
                + (N * Iq + ERf) * Y ** -0.5 * math.log(Y))

    def KB(t):
        return 2 * lam(math.sqrt(t)) + 0.8 * N * L1f + ERf * t ** -0.25 * math.log(t) / 2 + 2 * phi * l2 * t ** -1.25

    def U(t):
        K = 0.8 * c1 if t <= 1e12 else KB(t)
        return -mf + K * t ** -0.25 + rplus(t)

    lo, hi = 8.0, 1e12
    if U(lo) < 0:
        X0 = lo
    else:
        for _ in range(200):
            mid = math.sqrt(lo * hi)
            if U(mid) < 0: hi = mid
            else: lo = mid
        X0 = hi
    print(f'  X_0 (exact root of U on the K_A branch) = {X0:.6g};  U(X0*1.001)={U(X0*1.001):.3e}  U(1e6)={U(1e6):.5f}  U(1e12-)={U(1e12):.5f}  U(1e12+)={U(1e12*1.0000001):.5f}  K_B(1e12)={KB(1e12):.4f}  K_B(1e30)={KB(1e30):.4f}  -m={-mf:.6f}')

    # ---- certificate of W_1 < 0 on pieces up to T
    Pa = P[(P % q == a % q)] if q > 1 else P
    Pa = Pa[Pa >= 2]
    nxt_idx = np.searchsorted(Pa, T, side='right')  # pieces starting at primes <= T
    vt = arb(0); lam_ = arb(0)
    worst = None; worst_rel = -1e300; fail = 0; interior = 0; nint = 0
    true_sup_rel = -1e300; true_sup_at = None
    decade = {}
    tf_th = 0.0; tf_La = 0.0
    for idx in range(nxt_idx):
        p = int(Pa[idx]); pn = int(Pa[idx + 1])
        lp = arb(p).log()
        vt += lp; lam_ += p * lp
        M = (phi * vt) ** 2 / 2 + 2 * phi * vt - phi * lam_
        if not (M < 0):
            fail += 1
            print('   FAIL/undecided at', p, M)
        # float exploration of the true maximum on [p, pn]
        tf_th += math.log(p); tf_La += p * math.log(p)
        tstar = 2 + phi * tf_th
        def W(t):
            return phi * (tf_th * t - tf_La) - (t - 2) ** 2 / 2
        if p <= tstar <= pn:
            nint += 1
            tmax = tstar
        else:
            tmax = p if W(p) > W(pn) else pn
        rel = W(tmax) / tmax ** 1.5
        if rel > true_sup_rel:
            true_sup_rel, true_sup_at = rel, (p, tmax, p <= tstar <= pn)
        Mrel = float(M.mid()) / pn ** 1.5
        dk = int(math.log10(p))
        if dk not in decade or Mrel > decade[dk][0]:
            decade[dk] = (Mrel, p)
    print(f'  certificate: {nxt_idx} pieces (primes = {a} mod {q} up to {int(Pa[nxt_idx-1])}), W_1<0 valid up to {int(Pa[nxt_idx])}, failures {fail}')
    print(f'  pieces whose concave maximum t*=2+phi*theta lies inside the piece: {nint} of {nxt_idx}')
    print(f'  true sup of W_1/t^1.5 over (p_1, {int(Pa[nxt_idx])}] = {true_sup_rel:.5f} at piece starting {true_sup_at[0]}, t = {true_sup_at[1]:.2f}, interior={true_sup_at[2]}')
    print('  per decade max M_i/p_{i+1}^1.5 :', {k: (round(v[0], 4), int(v[1])) for k, v in sorted(decade.items())})

    # ---- (F) on a grid: F(log x) = x^{-3/2}[phi psi_1(x;q,a) - x^2/2 - l(x) + imp(x)]
    # prime powers n = a mod q; imp: p | q, p^k = a mod m_p
    Xmax = min(T, 2e6)
    ppw = []
    for p in P[P <= Xmax]:
        p = int(p); pk = p; k = 1
        while pk <= Xmax:
            ppw.append((pk, p)); pk *= p
    ppw.sort()
    grid = np.unique(np.concatenate([np.geomspace(1, Xmax, 4000), np.array([n for n, _ in ppw if n <= 2000], float)]))
    maxabs = 0; arg = None
    ps1_terms = [(n, math.log(p)) for n, p in ppw if q == 1 or n % q == a % q]
    imp_terms = []
    for p in primefactors(q):
        mpp = q
        while mpp % p == 0: mpp //= p
        for n, pp in ppw:
            if pp == p and (n - a) % mpp == 0:
                imp_terms.append((n, int(totient(mpp)) * math.log(p)))
    ns = np.array([n for n, _ in ps1_terms], float); ws = np.array([w for _, w in ps1_terms])
    ni = np.array([n for n, _ in imp_terms], float); wi = np.array([w for _, w in imp_terms])
    for x in grid:
        ps1 = phi * np.sum(ws[ns <= x] * (x - ns[ns <= x]))
        imp = np.sum(wi[ni <= x] * (x - ni[ni <= x])) if len(ni) else 0.0
        u = 1 / x
        Te = (x / 2) * ((1 + u) * math.log(1 + u) + ((1 - u) * math.log(1 - u) if u < 1 else 0))
        To = (1 - 0.5 * math.log(1 - u * u) - math.atanh(u) / u) if u < 1 else 1 - math.log(2)
        ell = -ne * x * math.log(x) + Bf * x + no * math.log(x) + Cf - (1 + ne) * Te - no * To
        F = (ps1 - x * x / 2 - ell + imp) / x ** 1.5
        if abs(F) > maxabs:
            maxabs, arg = abs(F), x
    print(f'  prime-side F on grid 1<=x<={Xmax:.0e}: max|F| = {maxabs:.6f} at x={arg:.4g}  (L_1 = {L1f:.6f}; must be <=)')
    return dict(q=q, a=a, L1=str(L1), m=str(m), ne=ne, no=no, B=str(B), C=str(C), ER=str(ER), Iq=Iq, c1=c1, X0=X0,
                cert_pieces=int(nxt_idx), valid_to=int(Pa[nxt_idx]), failures=fail, true_sup=true_sup_rel, maxabsF=maxabs)


def _v(q, p):
    v = 0
    while q % p == 0:
        q //= p; v += 1
    return v


if __name__ == '__main__':
    q, a, T = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3])
    res = main(q, a, T)
    json.dump(res, open(f'aw1_{q}_{a}.json', 'w'), indent=1, default=str)
