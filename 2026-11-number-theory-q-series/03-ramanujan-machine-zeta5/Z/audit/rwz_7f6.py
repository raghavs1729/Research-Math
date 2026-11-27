# Consistency of Theorem Z with Rogers-Wan-Zucker (31) / Wan (25):
# L(f,1) = (pi^2/32) * 7F6(5/4,1/2^6;1/4,1^5;1) = (pi^2/32) sum (4k+1) a_k^3, hence S_1 = -(pi^2/16) sum (4k+1) a_k^3.
from mpmath import mp, mpf, pi, nsum, inf, gamma, sqrt, harmonic, log, hyper, nstr
mp.dps = 40
a = lambda n: (gamma(n+mpf(1)/2)/(sqrt(pi)*gamma(n+1)))**2
X = hyper([mpf(5)/4]+[mpf(1)/2]*6, [mpf(1)/4]+[1]*5, 1)
S1 = nsum(lambda n: a(n)**2*(harmonic(2*n)-harmonic(n)-log(2)), [0, inf])
print("7F6 =", nstr(X, 30), " S1 + pi^2 X/16 =", nstr(S1 + pi**2*X/16, 3))
# Borwein-Straub-Wan, "Three-step and four-step random walk integrals" (preprint 2010 / Exp. Math. 22 (2013)):
# W4(-1) = (pi/4) 7F6(...)  [Thm 4],  W4(-1) = (4/pi) sum C(2n,n)^4/4^{4n} sum_{k>2n} (-1)^{k+1}/k  [display after (29)],
# W4(-1) = (4/pi^3) int_0^1 K'(k)^2 dk  [(65)].  Since sum_{k>2n}(-1)^{k+1}/k = -h_n, the middle one is -(4/pi) S_1.
W4 = pi/4*X
print("W4(-1) =", nstr(W4, 25), "  W4(-1) + (4/pi) S1 =", nstr(W4 + 4/pi*S1, 3))
from mpmath import quad, ellipk, isfinite
def Kp(k):
    v = ellipk((1-k)*(1+k)); return v if isfinite(v) else mpf(0)
print("  W4(-1) - (4/pi^3) int K'^2 =", nstr(W4 - 4/pi**3*quad(lambda k: Kp(k)**2, [0, 1]), 3))
