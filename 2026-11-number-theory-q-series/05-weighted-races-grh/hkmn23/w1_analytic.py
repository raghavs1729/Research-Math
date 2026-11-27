#!/usr/bin/env python3
"""
w1_analytic.py -- the analytic bound of HKMN23.md (Theorem 2) made explicit, in ball arithmetic.

Under GRH(q) (RH and GRH for every Dirichlet L-function of conductor dividing q), for t >= 8:

    W_1(t) = phi(q) theta_1(t;q,a) - (t-2)^2/2  <=  t^(3/2) * U(t),
    U(t)   = -m + K(t) t^(-1/4) + r_plus(t),

    m       = (2/3) N(q,a) - L_1(q)                      [l1const.json]
    K(t)    = K_A = (4/5) c_1                            for t <= Y1^2   (c_1 from w1_cert.py, part B)
            = K_B(t) (see below)                         for t >  Y1^2
    r_plus  = max(0,-n_e) t^(-1/2) log t + max(0,B+2) t^(-1/2) + max(0,n_o) t^(-3/2) log t
              + max(0,C-2) t^(-3/2) + max(0,-n_o) * 0.2 * t^(-7/2)

with, for chi mod q induced by the primitive chi* (constants b, b0, cm1 of chi* from l1const.json):
    n_e = sum_{chi even, chi != chi_0} conj(chi(a)),     n_o = sum_{chi odd} conj(chi(a)),
    B   = -log(2 pi) + sum_{even, != chi_0} Re[conj(chi(a)) (1 - b)] - sum_{odd} Re[conj(chi(a)) b],
    C   = (zeta'/zeta)(-1) + sum_{odd} Re[conj(chi(a)) (1 + b0)] + sum_{even, != chi_0} Re[conj(chi(a)) cm1].

K_B (prime squares beyond the computed range, from H(q) itself; c_theta = log 4 (Chebyshev), Y = sqrt t):
    lambda(Y) = N [L_1 + (2/3) phi c_theta] + N phi c_theta [(3/4) Y^(-1/6) + (log Y/log 2) Y^(-1/4)]
                + (N I_q + E_R) Y^(-1/2) log Y
    K_B(t)    = 2 lambda(Y) + (4/5) N L_1 + E_R t^(-1/4) (log t)/2 + 2 phi log 2 * t^(-5/4)
    I_q = sum_{p | q} phi(q / p^{v_p(q)}),
    E_R = |n_e^R| + |B^R|/log 3 + |n_o^R|/3 + |C^R|/(3 log 3) + 0.6 (N + |n_e^R|)/(9 log 3) + 0.2 |n_o^R|/(27 log 3),
    where X^R = sum over b in R = {b : b^2 = a mod q} of the constant X for the pair (q, b).

U is decreasing on [8, Y1^2] and on (Y1^2, oo).  The script finds X0 with U(X0) < 0 on the first
range (certified in Arb at the printed X0), checks U(Y1^2 +) < 0, and compares with the finite
certificate (W_1 < 0 on (2, P'] from w1_cert.py).  Output: table + w1_analytic.json.
"""
import json, math, sys
from flint import arb, acb, dirichlet_char, ctx

ctx.prec = 128
D = json.load(open("l1const.json"))
CASES = [(2, 1), (3, 1), (4, 1), (5, 1), (5, 4), (8, 1), (12, 1), (15, 1), (15, 4), (20, 1), (20, 9), (24, 1)]
CTHETA = arb(4).log()


def B2(x):
    """Ball from [midpoint string (40 digits), radius string (5 digits, rounded to nearest)].
    Rounded OUTWARD: the radius is inflated by 1% and by the truncation error of the midpoint."""
    mid = arb(x[0])
    rad = arb(x[1]) * arb("1.01") + abs(mid) * arb(10) ** (-38) + arb(10) ** (-60)
    return mid + arb(0, 1) * 0 + (rad * arb(-1)).union(rad)


def C2(e):
    return acb(B2(e["re"]), B2(e["im"]))


def phi_of(n):
    return sum(1 for k in range(1, n + 1) if math.gcd(k, n) == 1)


def chars(q):
    """list of (conj(chi(.)) evaluator, primitive key 'f.n', parity) for all chi mod q."""
    out = []
    for n in range(1, q + 1):
        if math.gcd(n, q) != 1:
            continue
        if q == 1:
            out.append((None, "1.1", 0, True))
            continue
        chi = dirichlet_char(q, n)
        f = chi.conductor()
        if f == 1:
            out.append((chi, "1.1", 0, True))
            continue
        nf = n % f
        prim = dirichlet_char(f, nf)
        assert prim.is_primitive()
        # check that chi is induced by prim (Conrey labels are compatible with induction)
        for k in range(1, 2 * q):
            if math.gcd(k, q) == 1:
                assert abs(complex(chi(k)) - complex(prim(k))) < 1e-25
        out.append((chi, f"{f}.{nf}", prim.parity(), False))
    return out


def chibar(chi, a):
    if chi is None:
        return acb(1)
    return acb(chi(a)).conjugate()


def pair_constants(q, a):
    ne = arb(0); no = arb(0)
    Bc = -(2 * arb.pi()).log()
    Cc = C2(D["chars"]["1.1"]["zeta_logder_m1"]).real
    L1 = arb(0)
    for chi, key, par, principal in chars(q):
        e = D["chars"][key]
        L1 += B2(e["S"])
        if principal:
            continue
        w = chibar(chi, a)
        b = C2(e["b"])
        if par == 0:
            ne += w.real
            Bc += (w * (1 - b)).real
            Cc += (w * C2(e["cm1"])).real
        else:
            no += w.real
            Bc -= (w * b).real
            Cc += (w * (1 + C2(e["b0"]))).real
    return {"n_e": ne, "n_o": no, "B": Bc, "C": Cc, "L1": L1}


def pos(x):
    return x if x > 0 else (arb(0) if x < 0 else x.union(arb(0)).upper() * 1 if False else abs(x))


def up(x):
    """float upper bound of an arb"""
    return float(arb(x).upper())


def analyse(q, a):
    phi = phi_of(q)
    R = [b for b in range(q) if math.gcd(b, q) == 1 and (b * b - a) % q == 0]
    N = len(R)
    pc = pair_constants(q, a)
    L1 = pc["L1"]
    m = arb(2) * N / 3 - L1
    assert m > 0
    # integer checks for n_e, n_o
    ne_i, no_i = round(float(pc["n_e"].mid())), round(float(pc["n_o"].mid()))
    assert pc["n_e"].contains(ne_i) and pc["n_o"].contains(no_i) and pc["n_e"].rad() < 1e-20
    Bc, Cc = pc["B"], pc["C"]
    # set constants for R
    neR = arb(0); noR = arb(0); BR = arb(0); CR = arb(0)
    for b in R:
        pb = pair_constants(q, b)
        neR += pb["n_e"]; noR += pb["n_o"]; BR += pb["B"]; CR += pb["C"]
    l3 = arb(3).log()
    ER = abs(neR) + abs(BR) / l3 + abs(noR) / 3 + abs(CR) / (3 * l3) + arb("0.6") * (N + abs(neR)) / (9 * l3) + arb("0.2") * abs(noR) / (27 * l3)
    Iq = 0
    qq = q
    p = 2
    while qq > 1:
        if qq % p == 0:
            pe = 1
            while qq % p == 0:
                qq //= p; pe *= p
            Iq += phi_of(q // pe)
        p += 1
    cert = json.load(open(f"raw/w1_cert_{q}_{a}.json"))
    assert cert["A"]["failures"] == 0 and cert["N"] == N and cert["phi"] == phi
    c1 = arb(cert["B"]["c1_upper"])
    Y1 = cert["Y1"]
    KA = arb(4) / 5 * c1

    def rplus(t):
        t = arb(t)
        lt = t.log()
        r = arb(0)
        if ne_i < 0:
            r += (-ne_i) * lt / t.sqrt()
        Bp = Bc + 2
        if not (Bp < 0):
            r += abs(Bp) / t.sqrt()
        if no_i > 0:
            r += no_i * lt / t ** 1.5
        Cp = Cc - 2
        if not (Cp < 0):
            r += abs(Cp) / t ** 1.5
        if no_i < 0:
            r += (-no_i) * arb("0.2") / t ** 3.5
        return r

    def KB(t):
        t = arb(t)
        Y = t.sqrt()
        lam = N * (L1 + arb(2) / 3 * phi * CTHETA) + N * phi * CTHETA * (arb(3) / 4 * Y ** (-arb(1) / 6) + (Y.log() / arb(2).log()) * Y ** (-arb(1) / 4)) \
            + (N * Iq + ER) * Y.log() / Y.sqrt()
        return 2 * lam + arb(4) / 5 * N * L1 + ER * t.log() / 2 / t ** 0.25 + 2 * phi * arb(2).log() / t ** 1.25

    def UA(t):
        t = arb(t)
        return -m + KA / t ** 0.25 + rplus(t)

    def UB(t):
        t = arb(t)
        return -m + KB(t) / t ** 0.25 + rplus(t)

    # X0: bisection in floats on [8, Y1^2], then certify at a rounded-up point
    lo, hi = 8.0, float(Y1) ** 2
    assert UA(hi) < 0, "regime A bound not negative at Y1^2"
    if UA(lo) < 0:
        X0 = 8.0
    else:
        for _ in range(200):
            mid = math.sqrt(lo * hi)
            if UA(mid) < 0:
                hi = mid
            else:
                lo = mid
        # round up to 3 significant digits
        e = math.floor(math.log10(hi)) - 2
        X0 = math.ceil(hi / 10 ** e * (1 + 1e-9)) * 10 ** e
    assert UA(X0) < 0
    UBv = UB(arb(Y1) ** 2)
    assert UBv < 0, f"regime B fails at Y1^2: {UBv}"
    lim = -m      # U(t) -> -m as t -> oo
    PA_valid = cert["A"]["valid_up_to"]
    ok = X0 < cert["PA"]
    # consistency: certified sup of W_1/t^{3/2} per decade must lie below U at the left end of the decade (for decades >= X0)
    cons = []
    for k, v in cert["A"]["per_decade_sup"].items():
        k = int(k)
        left = 10.0 ** k
        if left >= X0:
            cons.append((k, v["sup"], up(UA(left))))
            # U(t) >= W_1/t^{3/2} pointwise and U is decreasing; so sup over decade <= U(left)
            assert v["sup"] <= up(UA(left)) + 1e-12, (q, a, k, v, up(UA(left)))
    out = {
        "q": q, "a": a, "phi": phi, "N": N, "R": R,
        "L1": L1.str(20, radius=False), "margin_m": m.str(15, radius=False),
        "n_e": ne_i, "n_o": no_i, "B": Bc.str(15, radius=False), "C": Cc.str(15, radius=False),
        "n_e_R": neR.str(6, radius=False), "n_o_R": noR.str(6, radius=False), "B_R": BR.str(12, radius=False), "C_R": CR.str(12, radius=False),
        "E_R": ER.str(12, radius=False), "I_q": Iq,
        "c1": cert["B"]["c1_upper"], "c1_attained_before": cert["B"]["attained_before"], "Y1": Y1, "K_A": up(KA),
        "K_B_at_Y1sq": up(KB(arb(Y1) ** 2)), "K_B_limit": up(2 * N * (L1 + arb(2) / 3 * phi * CTHETA) + arb(4) / 5 * N * L1),
        "X0": X0, "U_A_at_X0": up(UA(X0)), "U_A_at_Y1sq": up(UA(arb(Y1) ** 2)), "U_B_at_Y1sq": up(UBv),
        "U_limit": up(lim),
        "finite_cert_valid_up_to": PA_valid, "finite_cert_num_primes": cert["A"]["num_primes"],
        "finite_sup_W1_over_t32": cert["A"]["sup_W1_over_t32_upper"], "finite_sup_at_p": cert["A"]["sup_at_p"],
        "overlap_ok": ok, "overlap_factor": PA_valid / X0,
        "U_table": {f"1e{k}": up(UA(10.0 ** k) if 10.0 ** k <= float(Y1) ** 2 else UB(10.0 ** k)) for k in (4, 5, 6, 7, 8, 9, 10, 12, 13, 16, 20, 30)},
        "consistency_decades": cons,
    }
    return out


def main():
    res = []
    print("# w1_analytic.py : W_1(t) <= t^{3/2} U(t) under GRH(q);  U(t) = -m + K(t) t^{-1/4} + r_plus(t)")
    print("# (q,a) phi N  L_1(q)      m=(2/3)N-L_1  n_e n_o    B        C       c_1(Y1=1e6)  X0(analytic)  finite cert up to   overlap   U(1e6)   U(1e12-)  U(1e12+)  -m")
    for (q, a) in CASES:
        o = analyse(q, a)
        res.append(o)
        print(f"({q:2d},{a:2d}) {o['phi']:2d} {o['N']:1d}  {float(o['L1']):.7f}  {float(o['margin_m']):.7f}   {o['n_e']:2d} {o['n_o']:2d}  {float(o['B']):+8.4f} {float(o['C']):+8.4f}  {o['c1']:8.4f}   {o['X0']:10.4g}   {o['finite_cert_valid_up_to']:12d}  {o['overlap_factor']:9.3g}  {o['U_table']['1e6']:+.4f}  {o['U_A_at_Y1sq']:+.4f}  {o['U_B_at_Y1sq']:+.4f}  {o['U_limit']:+.4f}")
        assert o["overlap_ok"]
    json.dump(res, open("w1_analytic.json", "w"), indent=1)
    print("# all cases: analytic range [X0, oo) and certified finite range (2, P'] overlap: OK")


if __name__ == "__main__":
    main()
