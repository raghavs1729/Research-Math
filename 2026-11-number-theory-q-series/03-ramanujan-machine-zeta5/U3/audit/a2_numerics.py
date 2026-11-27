"""Audit of U3: direct numerics (independent of verify/ and of the pickled constants).
 1. PCF from the recurrence (normalised sequences w3, v3), N up to 2^17, with Richardson in 1/N.
 2. S0, S1, Z0 by own partial sums + Richardson, cross-checked against mpmath hyper (S0, Z0) and
    a second route for S1 (derivative of the 5F4 in the parameter).
 3. Closed form, the L-free third form (as stated and with the sign fixed), Remark 1 numbers, sign claims,
    e2 telescoping, Remark 2 elliptic integral.
"""
from mpmath import agm, mp, mpf, pi, zeta, log, hyper, diff, quad, ellipk, ellipe, nstr
mp.dps = 60
half = mpf(1)/2
def richardson(Ns, vals):
    # polynomial extrapolation in h = 1/N to h = 0 (Neville)
    hs = [mpf(1)/N for N in Ns]; P = list(vals)
    for lev in range(1, len(P)):
        for i in range(len(P)-1, lev-1, -1):
            P[i] = (hs[i-lev]*P[i] - hs[i]*P[i-1])/(hs[i-lev] - hs[i])
    return P[-1], P[-2]
# 1. PCF
def at3(m):
    T = m + half; return 2*T*(T**4 + 5*T**2 + mpf(25)/16)
Ns = [1000*2**j for j in range(8)]
w0, w1 = mpf(1), at3(0); v0, v1 = mpf(0), mpf(1)/8
ratios = []; mm = 1; target = iter(Ns); nxt = next(target)
while True:
    w2 = (at3(mm)*w1 - mm**5*w0)/(mm+1)**5; v2 = (at3(mm)*v1 - mm**5*v0)/(mm+1)**5
    w0, w1, v0, v1 = w1, w2, v1, v2; mm += 1
    if mm == nxt:
        ratios.append(w1/v1)
        try: nxt = next(target)
        except StopIteration: break
print("1. PCF partial (N=%d): %s" % (Ns[-1], nstr(ratios[-1], 45)))
print("   successive differences:", [nstr(ratios[i+1]-ratios[i], 3) for i in range(len(ratios)-1)])
print("   w3_N/N at N=%d: %s" % (Ns[-1], nstr(w1/Ns[-1], 15)))
pcf_r, pcf_r2 = richardson(Ns, ratios)
print("   Richardson: %s   (previous level diff %s)" % (nstr(pcf_r, 45), nstr(pcf_r - pcf_r2, 3)))
PCF = pcf_r
# 2. constants via own sums
def sums(Nmax_list):
    A = mpf(1); k = 0; s0 = s1 = z0 = mpf(0); H2k = Hk = mpf(0); L2 = log(2); out = []
    it = iter(Nmax_list); nxt = next(it)
    while True:
        s0 += A; s1 += A*(H2k - Hk - L2); z0 += A/(k+1)
        # advance
        A *= ((k+half)/(k+1))**4
        H2k += mpf(1)/(2*k+1) + mpf(1)/(2*k+2); Hk += mpf(1)/(k+1); k += 1
        if k == nxt:
            out.append((s0, s1, z0))
            try: nxt = next(it)
            except StopIteration: return out
KN = [500*2**j for j in range(9)]
out = sums(KN)
S0, _a = richardson(KN, [o[0] for o in out]); S1, _b = richardson(KN, [o[1] for o in out]); Z0, _c = richardson(KN, [o[2] for o in out])
print("2. own sums+Richardson: S0 = %s (lev diff %s)" % (nstr(S0, 40), nstr(S0-_a, 2)))
print("                        S1 = %s (lev diff %s)" % (nstr(S1, 40), nstr(S1-_b, 2)))
print("                        Z0 = %s (lev diff %s)" % (nstr(Z0, 40), nstr(Z0-_c, 2)))
mp.dps = 40
S0h = hyper([half]*4, [1, 1, 1], 1); Z0h = hyper([half]*4 + [1], [1, 1, 1, 2], 1)
Yf = lambda e: hyper([1] + [half+e]*4, [1+e]*4, 1)
S1h = (diff(Yf, 0) - 8*log(2)*Yf(0))/8
print("   mpmath hyper: S0 - own = %s, Z0 - own = %s, S1(derivative route) - own = %s"
      % (nstr(S0h-S0, 3), nstr(Z0h-Z0, 3), nstr(S1h-S1, 3)))
mp.dps = 60
print("   stated: S0 = 1.1186363871641870683, S1 = -0.7090013674619294375, Z0 = 1.0433551893343198886")
# 3. formula
K = 4 + pi**2*(Z0 - 2*S0)
inv = mpf(93)/4*zeta(5) - 4 + pi**6/64*S0*K/(1 + S1*K)
inv2 = mpf(93)/4*zeta(5) + pi**6*S0/(64*S1)
third_stated = inv2 - 4 + pi**6*S0/(64*S1*(1 + S1*K))
third_fixed = inv2 - 4 - pi**6*S0/(64*S1*(1 + S1*K))
print("3. K = %s, 1+S1K = %s, 2S0-Z0 = %s" % (nstr(K, 20), nstr(1+S1*K, 20), nstr(2*S0-Z0, 20)))
print("   1/PCF(direct) = %s" % nstr(1/PCF, 40))
print("   formula       = %s ; diff = %s" % (nstr(inv, 40), nstr(1/PCF - inv, 3)))
print("   1/formula     = %s" % nstr(1/inv, 40))
print("   stated value    22.841008565767985615909667462437")
print("   third form as stated: %s (diff %s);  sign-fixed: diff %s" % (nstr(third_stated, 20), nstr(third_stated-inv, 3), nstr(third_fixed-inv, 3)))
print("   1/U2 = %s, U2 = %s" % (nstr(inv2, 25), nstr(1/inv2, 25)))
# e2 telescoping: sum A_k/(k-1/2) = -8/pi^2 - Z0
X = hyper([half, half, half, -half], [1, 1, 1], 1)   # sum A_k (1/2)/(1/2 - k)... check: (-1/2)_k/(1/2)_k = -1/(2k-1)
tel = -2*X                                          # sum A_k/(k-1/2) = -2 sum A_k/(1-2k)
print("   sum A_k/(k-1/2) + 8/pi^2 + Z0 = %s" % nstr(tel + 8/pi**2 + Z0, 3))
mp.dps = 30
def _f(kk):
    if kk == 0 or kk == 1: return mpf(0)
    with mp.workdps(2*mp.dps + 40):
        Kp = pi/(2*agm(1, kk))                      # K'(k), stable for small k
        m1 = 1 - kk**2
        val = Kp*(ellipe(kk**2) - m1*ellipk(kk**2))/kk**2
    return +val
J = quad(_f, [0, half, 1])
print("   Remark 2: (16/pi^3) int K'(E-(1-k^2)K)/k^2 dk - Z0 = %s" % nstr(16/pi**3*J - Z0, 3))
