#!/usr/bin/env python3
"""Rigorous table of  Phi(e^{kh}),  k = 0..kmax,  for the symmetric square of Delta.

  L_S(s)      = L(s, Sym^2 Delta) = sum a_S(n) n^-s   (analytic normalisation, centre 1/2, conductor 1)
  gamma_S(s)  = Gamma_R(s+1) Gamma_C(s+11) = pi^{-(s+1)/2} Gamma((s+1)/2) * 2 (2 pi)^{-(s+11)} Gamma(s+11)
  Lambda_S(s) = gamma_S(s) L_S(s) = Lambda_S(1-s)
  K(x)        = 4 x int_0^oo u^9 exp(-pi x^2/u^2 - 2 pi u) du      (inverse Mellin transform of gamma_S)
  Phi(y)      = sum_n a_S(n) K(n y),      Lambda_S(s) = int_0^oo Phi(y) y^{s-1} dy,   Phi(1/y) = y Phi(y).

a_S(n) is rational: multiplicative, a_S(p^k) = c a_S(p^{k-1}) - c a_S(p^{k-2}) + a_S(p^{k-3}),
c = tau(p)^2/p^11 - 1.

Each Phi_k = Phi(e^{kh}) is enclosed in an Arb ball containing:
   * the Arb integral over [u_lo, u_hi] for every n <= XMAX/y,
   * the tails of the u-integral   (int_0^{u_lo} <= u_lo^10 e^{-pi x^2/u_lo^2},  int_{u_hi}^oo = Gamma(10, 2 pi u_hi)/(2 pi)^10),
   * the tail n > N = max(1, [XMAX/y]):  sum_{n>N} n^2 K(n y) <= C_eta y [ (N+1)^3 e^{-beta (N+1)^{2/3}} + (3/2) Gamma(6, beta (N+1)^{2/3})/beta^6 ],
     beta = b_eta y^{2/3},
     from K(x) <= C_eta x exp(-b_eta x^{2/3}),  C_eta = 4*9!/(2 pi eta)^10, b_eta = 3 pi (1-eta)^{2/3}.
Usage: sym2_phi.py [M=640] [XMAX=160] [nproc=4]      (h = 2 pi / M)
"""
import sys, json, hashlib, time, os
from fractions import Fraction
from flint import arb, acb, ctx, fmpq
DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data')
PREC = 330
ETA = Fraction(1, 10)

def aS_table(N):
    """exact a_S(n), 1 <= n <= N, as Fractions"""
    from tau_table import tau_list
    tau = tau_list(N+2)
    a = [Fraction(0)]*(N+1); a[1] = Fraction(1)
    spf = list(range(N+1))
    for i in range(2, int(N**0.5)+1):
        if spf[i] == i:
            for j in range(i*i, N+1, i):
                if spf[j] == j: spf[j] = i
    pp = {}
    for n in range(2, N+1):
        p = spf[n]; m = n; k = 0
        while m % p == 0: m //= p; k += 1
        if m == 1:
            c = Fraction(tau[p]**2, p**11) - 1
            seq = pp.setdefault(p, [Fraction(1)])
            while len(seq) <= k:
                j = len(seq)
                v = c*seq[j-1] - (c*seq[j-2] if j >= 2 else 0) + (seq[j-3] if j >= 3 else 0)
                seq.append(v)
            a[n] = seq[k]
        else:
            a[n] = a[m]*a[n//m]
    return a

def K_arb(x):
    """rigorous enclosure of K(x), x a positive arb (exact or thin)"""
    pi = arb.pi()
    # exact (dyadic) endpoints, so that the integration path is a fixed interval
    u_lo = arb(float((x*(pi/230).sqrt()).mid()))          # pi x^2/u_lo^2 ~ 230
    u_hi = arb(float(arb(46).max(2*x**(arb(2)/3) + 30).mid()))
    assert u_lo.rad() == 0 and u_hi.rad() == 0 and u_lo > 0 and u_hi > u_lo
    px2 = pi*x*x
    def f(u, _):
        return u**9 * (-(px2/(u*u)) - 2*pi*u).exp()
    I = acb.integral(f, acb(u_lo), acb(u_hi), rel_tol=arb(2)**(-PREC+20), abs_tol=arb(10)**(-95))
    I = I.real
    t_lo = u_lo**10 * (-(px2/(u_lo*u_lo))).exp()       # int_0^{u_lo} u^9 e^{-pi x^2/u^2 - 2 pi u} du <= u_lo^10 e^{-pi x^2/u_lo^2}
    t_hi = (2*pi*u_hi).gamma_upper(arb(10)) / (2*pi)**10
    I = I + arb(0, (t_lo + t_hi).upper())
    return 4*x*I

def tail_n(y, N):
    """bound for sum_{n>N} n^2 K(n y)  (y >= 1)"""
    pi = arb.pi(); eta = arb(ETA.numerator)/ETA.denominator
    Ceta = 4*arb(362880)/(2*pi*eta)**10
    b = 3*pi*(1-eta)**(arb(2)/3)
    beta = b*y**(arb(2)/3)
    # n -> n^3 exp(-beta n^{2/3}) decreases for n >= 1 (beta >= 9/2), so
    # sum_{n>N} <= (N+1)^3 exp(-beta (N+1)^{2/3}) + int_{N+1}^oo = first + (3/2) Gamma(6, beta (N+1)^{2/3})/beta^6
    assert beta > arb("4.5")
    w = arb(N+1)**(arb(2)/3)
    return Ceta*y*(arb(N+1)**3*(-(beta*w)).exp() + arb(3)/2*(beta*w).gamma_upper(arb(6))/beta**6)

def phi_k(args):
    k, M, XMAX, aS = args
    ctx.prec = PREC
    h = 2*arb.pi()/M
    y = (h*k).exp()
    N = max(1, int(XMAX/float(y.mid())))     # any N >= 1 is valid: the tail bound below is for n > N
    S = arb(0)
    for n in range(1, N+1):
        S += arb(aS[n].numerator)/aS[n].denominator * K_arb(n*y)
    S = S + arb(0, tail_n(y, N).upper())
    mm, me = S.mid().man_exp(); rm, re_ = S.rad().man_exp()
    return k, N, [str(mm), int(me), str(rm), int(re_)]

if __name__ == "__main__":
    M = int(sys.argv[1]) if len(sys.argv) > 1 else 640
    XMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 160
    nproc = int(sys.argv[3]) if len(sys.argv) > 3 else 4
    import math
    kmax = int(math.log(XMAX)*M/(2*math.pi)) + 1
    t0 = time.time()
    aS = aS_table(XMAX+1)
    from multiprocessing import Pool
    jobs = [(k, M, XMAX, aS) for k in range(kmax+1)]
    with Pool(nproc) as pool:
        res = pool.map(phi_k, jobs, chunksize=4)
    table = {str(k): v for k, N, v in res}
    out = dict(M=M, XMAX=XMAX, kmax=kmax, prec=PREC, eta=str(ETA), table=table,
               aS_first=[str(aS[n]) for n in range(1, 13)])
    js = json.dumps(out, sort_keys=True) + "\n"
    fn = DATA + "/sym2_phi_M%d.json" % M
    open(fn, "w").write(js)
    print("wrote", fn, " kmax =", kmax, " sha256 =", hashlib.sha256(js.encode()).hexdigest(), " time %.1f s" % (time.time()-t0))
