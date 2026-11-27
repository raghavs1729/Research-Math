"""Second audit: lower-limit constants (Corollary 1.3), mpmath 30 digits, by direct quadrature of li
(not via the antiderivative x li(x) - li(x^2) used by the author)."""
import mpmath as mp
mp.mp.dps = 30
li = mp.li
I02 = mp.quad(li, [0, 1, 2])           # integral_0^2 li (log singularity at 1 is integrable)
print('int_0^2 li =', I02, '  c0 = li(4)-2li(2) =', li(4) - 2 * li(2))
g = lambda x: -mp.quad(li, [0, 1, x]) if x > 1 else -mp.quad(li, [0, x])
x0 = mp.findroot(g, 2.6)
print('x_0 =', x0, ' g(2.6)=', g(2.6), ' g(2.62)=', g(2.62))
# q = 1: A_0(x) = g(x) + sum_{p<=x}(x-p) on [3,5)
A0 = lambda x: g(x) + (x - 2) + (x - 3 if x >= 3 else 0)
x1 = mp.findroot(A0, 3.5)
print('x_1 =', x1, ' A0(3)=', A0(3), ' A0(4.999999)=', A0(mp.mpf('4.999999')))
print('pi(p) < li(p) at primes <= 31:', all(i + 1 < li(p) for i, p in enumerate([2,3,5,7,11,13,17,19,23,29,31])))
