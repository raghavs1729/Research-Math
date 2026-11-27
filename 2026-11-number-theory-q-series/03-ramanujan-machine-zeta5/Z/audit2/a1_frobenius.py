"""Audit 2, script 1: Zudilin's F_j(1) for (r,t)=(1/2,1/4),(1/2,1/3) straight from eq. (5) (unilateral form).
Two independent methods:
  (A) termwise Taylor coefficients in eps via polygammas, summed with mpmath.nsum (Richardson, cross-checked with Levin);
  (B) Cauchy integral over |eps| = r of F(eps) (each F(eps) by nsum), N-point trapezoid (spectrally accurate).
Also: the BILATERAL sum of eq. (5) at small eps vs the unilateral (difference O(eps^4)), and F_3/F_1 vs (csc^2 pi t + 1/3) pi^2."""
from mpmath import mp, mpf, mpc, gamma, loggamma, psi, nsum, pi, sin, exp, inf, nstr, fsum, cos, polyroots
import sys
mp.dps = 50
half = mpf(1)/2

def coeffs_termwise(t, method='richardson'):
    t = mpf(t)
    norm = 1/(pi*gamma(t)*gamma(1-t))
    pars = [half, half, t, 1-t]
    def base(n):
        return exp(sum(loggamma(p+n) for p in pars) - 4*loggamma(1+n))
    def l(k, n):  # k-th eps-derivative of log term
        return sum(psi(k-1, p+n) for p in pars) - 4*psi(k-1, 1+n)
    def c(j):
        def term(n):
            b = base(n); l1, l2, l3 = l(1,n), l(2,n), l(3,n)
            return b*[1, l1, (l1**2+l2)/2, (l1**3+3*l1*l2+l3)/6][j]
        return norm*nsum(term, [0, inf], method=method)
    return [c(j) for j in range(4)]

def F_eps(t, eps):
    t = mpf(t)
    norm = 1/(pi*gamma(t)*gamma(1-t))
    pars = [half, half, t, 1-t]
    f = lambda n: exp(sum(loggamma(p+eps+n) for p in pars) - 4*loggamma(1+eps+n))
    return norm*nsum(f, [0, inf])

def coeffs_cauchy(t, r=mpf(1)/16, N=40):
    vals = [F_eps(t, r*exp(2j*pi*k/N)) for k in range(N)]
    return [ (fsum(vals[k]*exp(-2j*pi*k*j/N) for k in range(N))/N / r**j).real for j in range(4)]

def bilateral(t, eps, M=None):
    # sum over n<0 by reflection: term(-m) = prod_rho (-sin(pi eps)/sin(pi(eps+rho))) Gamma(m-eps)/Gamma(m+1-rho-eps)
    t = mpf(t); pars = [half, half, t, 1-t]
    norm = 1/(pi*gamma(t)*gamma(1-t))
    trig = 1
    for p in pars: trig *= -sin(pi*eps)/sin(pi*(eps+p))
    neg = nsum(lambda m: exp(sum(loggamma(m-eps)-loggamma(m+1-p-eps) for p in pars)), [1, inf])
    return norm*trig*neg

if __name__ == '__main__':
    res = {}
    for t, name in [(mpf(1)/4, '1/4'), (mpf(1)/3, '1/3')]:
        A = coeffs_termwise(t)
        A2 = coeffs_termwise(t, method='levin')
        B = coeffs_cauchy(t)
        print('t =', name)
        for j in range(4):
            print('  F_%d = %s   |A-Levin| = %s   |A-Cauchy| = %s' % (j, nstr(A[j], 40), nstr(abs(A[j]-A2[j]), 3), nstr(abs(A[j]-B[j]), 3)))
        S = 1/sin(pi*t)**2
        print('  F_3/(pi^2 F_1) =', nstr(A[3]/(pi**2*A[1]), 40), '  csc^2+1/3 =', nstr(S+mpf(1)/3, 40))
        e = mpf(10)**-3
        print('  bilateral-minus-unilateral at eps=1e-3 (should be ~eps^4 ~ 1e-12):', nstr(bilateral(t, e), 5))
        res[name] = A
    # dump for other scripts
    with open(__file__.replace('a1_frobenius.py', 'F_values.txt'), 'w') as fh:
        for name, A in res.items():
            fh.write(name + ' ' + ' '.join(nstr(a, 45) for a in A) + '\n')
