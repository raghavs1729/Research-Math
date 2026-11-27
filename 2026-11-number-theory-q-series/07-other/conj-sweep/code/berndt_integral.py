# arXiv:2401.01385 (Xu, Zhao et al., "General Berndt-type integrals ..."), Conjecture:
# closed form of int_0^oo x^9/(cos x + cosh x)^6 dx in Q[Gamma(1/4)^8, 1/pi].  High-precision numerical test.
from mpmath import mp, mpf, quad, cos, cosh, gamma, pi, inf, linspace
mp.dps = 70
f = lambda x: x**9/(cos(x) + cosh(x))**6
I = quad(f, [0, 1, 2, 3, 4, 6, 8, 12, 16, 24, 40, 80, inf])
G = gamma(mpf(1)/4)
terms = [(-63, 5*2**10, 0, 0), (1071, 5**2*2**13, 8, 2), (-21, 2**12, 8, 3), (63, 2**16, 8, 4), (-21, 5**3*2**13, 16, 4),
         (3, 5*2**13, 16, 5), (-161, 5*2**19, 16, 6), (21, 2**21, 16, 8), (1, 5**2*2**19, 24, 8), (-1, 3*2**20, 24, 9),
         (69, 5*2**25, 24, 10), (-21, 5*2**24, 24, 11), (63, 5*2**27, 24, 12), (-17, 3*5**2*2**31, 32, 14),
         (13, 5*2**34, 32, 16), (3, 5**2*2**40, 40, 20)]
R = sum(mpf(a)/b * G**g / pi**p for a, b, g, p in terms)
print('integral  =', mp.nstr(I, 55)); print('closed form =', mp.nstr(R, 55)); print('difference =', mp.nstr(I - R, 5))
