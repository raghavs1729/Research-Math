#!/usr/bin/env python3
"""Rigorous verification of the Riemann hypothesis for L(s, Sym^2 Delta) up to height H.

Inputs: the table Phi_k = Phi(e^{kh}) of sym2_phi.py (h = 2 pi/M), exact tau(p) for the Euler product.

Step 1 (critical line).  Z(t) = Lambda_S(1/2 + it) is real, and for |Im s| <= TMAX, -1 <= Re s <= 2,
    Lambda_S(s) = h [ Phi_0 + sum_{k=1}^{kmax} Phi_k (e^{s k h} + e^{(1-s) k h}) ] + E,   |E| <= E_tail + E_alias
  (Poisson summation; E_alias = sum_{m != 0} |Lambda_S(s + i M m)|; see TAU-RACE.md, Lemma V1).
  Sign changes of Z certified in ball arithmetic give zeros gamma_j in isolating intervals [l_j, r_j].
Step 2 (completeness, no contour integral).  For real T,
    Re Lambda_S'/Lambda_S (2 + iT) = sum_rho Re 1/(2 + iT - rho)      (all zeros, each term > 0),
  the left side is computed from the Euler product (Re s = 2) and the gamma factor.  Put
    Rem(T) = Re Lambda_S'/Lambda_S(2+iT) - sum_j [ w(T - gamma_j) + w(T + gamma_j) ],  w(u) = (3/2)/(9/4 + u^2).
  Rem >= 0 is the contribution of all zeros not among the found ones, |Rem'| <= Rem, and a zero (or a
  symmetric pair of zeros off the line) at height T0 that is not among the found ones forces Rem(T0) >= 2/3.
  So  max_grid Rem * e^{step/2} < 2/3  on [0, H] proves: the zeros with |Im| <= H are exactly +-gamma_j,
  they are simple and lie on the critical line.
Usage: sym2_verify.py [M=640] [H=60] [H2=72] [P=200000]
"""
import sys, json, hashlib, time, os, math
from fractions import Fraction
from flint import arb, acb, ctx, fmpz
from arbutil import dec_down, dec_up, encl, primes_upto
DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data')

def load_table(M):
    fn = DATA + "/sym2_phi_M%d.json" % M
    raw = open(fn).read()
    d = json.loads(raw)
    ctx.prec = d["prec"]
    tab = []
    for k in range(d["kmax"]+1):
        mm, me, rm, re_ = d["table"][str(k)]
        mid = arb(fmpz(int(mm))) * arb(2)**me
        rad = arb(fmpz(int(rm))) * arb(2)**re_
        tab.append(mid + arb(0, rad.upper()))
    return d, tab, hashlib.sha256(raw.encode()).hexdigest()

class LambdaS:
    def __init__(self, M, TMAX):
        self.d, self.tab, self.sha = load_table(M)
        self.M = M; self.TMAX = TMAX
        self.h = 2*arb.pi()/M
        self.kmax = self.d["kmax"]
        pi = arb.pi()
        # ---- k-tail: sum_{k>kmax} h |Phi(e^{kh})| 2 e^{2kh} <= 2 C (1+1.5 c0) (3/2) Gamma(9/2, b y0^{2/3}) / b^{9/2}
        eta = arb(1)/10
        C = 4*arb(362880)/(2*pi*eta)**10; b = 3*pi*(1-eta)**(arb(2)/3)
        c0 = b.gamma_upper(arb(6))*b.exp()/b**6
        y0 = (self.h*self.kmax).exp()
        self.E_tail = 2*C*(1+arb(3)/2*c0)*arb(3)/2*(b*y0**(arb(2)/3)).gamma_upper(arb(9)/2)/b**(arb(9)/2)
        # ---- alias: sum_{m != 0} |Lambda_S(s + iMm)| <= 2 M_theta e^{theta TMAX} sum_{m>=1} e^{-theta M m}
        th = pi/8; eta2 = arb(1)/5
        C2 = 4*arb(362880)/(2*pi*eta2)**10
        b2 = 3*pi*((2*th).cos()*(1-eta2)**2)**(arb(1)/3)
        assert b2 > arb("4.5")
        c02 = b2.gamma_upper(arb(6))*b2.exp()/b2**6
        Mth = 2*C2*(1+arb(3)/2*c02)*arb(3)/2*b2.gamma_upper(arb(9)/2)/b2**(arb(9)/2)
        q = (-th*M).exp()
        self.E_alias = 2*Mth*(th*TMAX).exp()*q/(1-q)
        self.M_theta = Mth
        self.E = arb(0, (self.E_tail + self.E_alias).upper())
        self.hk = [self.h*k for k in range(self.kmax+1)]

    def __call__(self, s):
        """s: acb with |Im s| <= TMAX, -1 <= Re s <= 2"""
        s = acb(s)
        S = acb(self.tab[0])
        for k in range(1, self.kmax+1):
            v = self.hk[k]
            S += self.tab[k]*((s*v).exp() + ((1-s)*v).exp())
        S = S*self.h
        return S + acb(self.E, self.E)

    def Z(self, t):
        """Z(t) = Lambda_S(1/2+it), real: = h[Phi_0 + sum 2 Phi_k e^{kh/2} cos(t k h)] + E"""
        t = arb(t)
        S = arb(self.tab[0])
        for k in range(1, self.kmax+1):
            v = self.hk[k]
            S += self.tab[k]*2*(v/2).exp()*(t*v).cos()
        return S*self.h + self.E

def gammaS(s):
    s = acb(s); pi = acb.pi()
    return pi**(-(s+1)/2)*((s+1)/2).gamma()*2*(2*pi)**(-(s+11))*(s+11).gamma()

def gammaS_logderiv(s):
    s = acb(s); pi = acb.pi()
    return -pi.log()/2 + ((s+1)/2).digamma()/2 - (2*pi).log() + (s+11).digamma()

def euler_logderiv(T, plist, tail):
    """rigorous enclosure of (L_S'/L_S)(2 + iT) from primes in plist = [(p, log p, c2_p)]"""
    S = acb(0)
    s = acb(2, T)
    for (p, lp, c2) in plist:
        x = (-s*lp).exp()
        S += lp*((c2*x - 2*x*x)/(1 - c2*x + x*x) + x/(1-x))
    return -S + acb(arb(0, tail.upper()), arb(0, tail.upper()))

def main(M, H, H2, P):
    t0 = time.time()
    Lam = LambdaS(M, H2 + 1)
    ctx.prec = Lam.d["prec"]
    out = dict(M=M, H=H, H2=H2, P=P, phi_table_sha256=Lam.sha,
               E_tail=dec_up(Lam.E_tail*arb(10)**80, 3) + "e-80", E_alias=dec_up(Lam.E_alias*arb(10)**80, 3) + "e-80")
    print("error bounds: tail", Lam.E_tail.str(5), " alias", Lam.E_alias.str(5), " M_theta", Lam.M_theta.str(5))
    maxrad = max(float(t.rad().upper()) for t in Lam.tab)
    zr = [float(Lam.Z(arb(t)).rad().upper()) for t in (0, 30, 60, 72)]
    print("largest radius of a table entry: %.3e;  radius of Z(t) at t = 0, 30, 60, 72: %s" % (maxrad, ["%.2e" % r for r in zr]))
    out.update(max_table_radius="%.3e" % maxrad, Z_radius_at_0_30_60_72=["%.2e" % r for r in zr])
    # ---- sanity: Lambda_S(2) against gamma_S(2) * Euler product; Lambda_S(1/2), Lambda_S(1)
    from tau_table import tau_list
    tau = tau_list(P+1)
    primes = primes_upto(P)
    plist = []
    for p in primes:
        lp = arb(p).log(); c2 = arb(tau[p])**2/arb(p)**11 - 2
        plist.append((p, lp, c2))
    tailP = arb("3.03")*(arb(P).log()+1)/P
    logL2 = arb(0)
    for (p, lp, c2) in plist:
        x = 1/arb(p)**2
        logL2 += -((1 - c2*x + x*x).log()) - (1-x).log()
    tl = 3*arb("1.01")/P     # |log tail| <= 3 sum_{p>P} p^-2/(1-p^-2) <= 3.03/P
    logL2 += arb(0, tl.upper())
    L2_euler = logL2.exp()
    L2_poisson = (Lam(acb(2))/gammaS(2)).real
    print("L_S(2): Euler", L2_euler.str(12), " Poisson", L2_poisson.str(40))
    assert L2_euler.overlaps(L2_poisson)
    LS_half = (Lam(acb(1)/2)/gammaS(acb(1)/2)).real; LS_one = (Lam(acb(1))/gammaS(1)).real
    print("L_S(1/2) =", LS_half.str(30), "  L_S(1) =", LS_one.str(30))
    out.update(LS_2=encl(L2_poisson, 30), LS_half=encl(LS_half, 30), LS_1=encl(LS_one, 30))
    # ---- Step 1: scan
    step = Fraction(1, 16)
    pts = []; t = Fraction(0)
    while t <= H2:
        pts.append(t); t += step
    def sgn(tq):
        z = Lam.Z(arb(tq.numerator)/tq.denominator)
        if z > 0: return 1
        if z < 0: return -1
        return 0
    signs = []
    for tq in pts:
        sg = sgn(tq); j = 0
        while sg == 0:           # move slightly (never needed in practice)
            j += 1; tq = tq + Fraction(1, 1024*j); sg = sgn(tq)
            assert j < 20, "cannot decide the sign of Z"
        signs.append((tq, sg))
    zeros = []
    for (a, sa), (b, sb) in zip(signs, signs[1:]):
        if sa != sb:
            lo, hi, slo = a, b, sa
            for _ in range(34):          # bisection to width 2^-38
                mid = (lo+hi)/2; sm = sgn(mid)
                if sm == 0: break
                if sm == slo: lo = mid
                else: hi = mid
            zeros.append((lo, hi))
    print(len(zeros), "sign changes on [0, %d]; first:" % H2, [float(z[0]) for z in zeros[:6]])
    print("widest isolating interval: %.3e" % max(float(hi-lo) for lo, hi in zeros))
    out["max_interval_width"] = "%.3e" % max(float(hi-lo) for lo, hi in zeros)
    out["zeros"] = [[str(lo), str(hi)] for lo, hi in zeros]
    out["n_zeros_upto_H"] = sum(1 for lo, hi in zeros if hi <= H)
    out["n_zeros_upto_H2"] = len(zeros)
    assert all(hi <= H or lo >= H for lo, hi in zeros)
    # ---- Step 2: completeness on [0, H]
    gstep = Fraction(1, 8)
    def w_lower(T, lo, hi):
        """lower bound of (3/2)/(9/4 + (T-g)^2) + (3/2)/(9/4 + (T+g)^2) for g in [lo, hi]"""
        g = arb(lo.numerator)/lo.denominator; g = g.union(arb(hi.numerator)/hi.denominator)
        T = arb(T.numerator)/T.denominator
        return arb(3)/2/(arb(9)/4 + (T-g)**2) + arb(3)/2/(arb(9)/4 + (T+g)**2)
    worst = None; T = Fraction(0); rem_rows = []
    ctx.prec = 128
    plist = [(p, arb(p).log(), arb(tau[p])**2/arb(p)**11 - 2) for p in primes]
    while T <= H:
        s = acb(2, arb(T.numerator)/T.denominator)
        tot = (gammaS_logderiv(s) + euler_logderiv(arb(T.numerator)/T.denominator, plist, tailP)).real
        found = arb(0)
        for lo, hi in zeros: found += w_lower(T, lo, hi)
        rem = tot - found
        if worst is None or not (rem < worst[0]): 
            if worst is None or rem > worst[0]: worst = (rem, T)
            else: worst = (worst[0].union(rem), T)
        if T.denominator == 1 and T % 5 == 0: rem_rows.append((int(T), dec_down(tot, 6), dec_up(rem, 6)))
        T += gstep
    bound = worst[0]*(arb(gstep.numerator)/gstep.denominator/2).exp()
    print("max Rem on grid <= %s at T = %s;  with Lipschitz factor: %s  (must be < 2/3)" % (dec_up(worst[0], 6), worst[1], dec_up(bound, 6)))
    for r in rem_rows: print("   T = %3d  Re Lambda'/Lambda(2+iT) >= %s   Rem <= %s" % r)
    ok = bool(bound < arb(2)/3)
    out.update(rem_max=dec_up(worst[0], 6), rem_bound=dec_up(bound, 6), complete=ok)
    js = json.dumps(out, indent=1, sort_keys=True) + "\n"
    open(DATA + "/sym2_rh_H%d.json" % H, "w").write(js)
    print("RESULT:", "PASS: all zeros of L(s,Sym^2 Delta) with |Im s| <= %d are simple and on the critical line; there are %d with 0 < Im <= %d" % (H, out["n_zeros_upto_H"], H) if ok else "FAIL")
    print("sha256 =", hashlib.sha256(js.encode()).hexdigest(), " time %.1f s" % (time.time()-t0))

if __name__ == "__main__":
    a = sys.argv[1:]
    main(int(a[0]) if len(a) > 0 else 640, int(a[1]) if len(a) > 1 else 60, int(a[2]) if len(a) > 2 else 72, int(a[3]) if len(a) > 3 else 200000)
