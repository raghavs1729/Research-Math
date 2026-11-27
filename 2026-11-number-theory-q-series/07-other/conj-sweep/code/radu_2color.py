# arXiv:2607.16608 (Bharadwaj, Sujatha, Chandankumar), closing conjectures.
# pbar_l(n): coefficients of f_2 f_{2l} / (f_1^2 f_l^2) = 1/(phi(-q) phi(-q^l)).
from radu import *
import flint, sys
out = '../cert/radu_2color_overpartitions.txt'; open(out, 'w').write('arXiv:2607.16608 conjectures\n\n')
def r(l): return {1: -2, 2: 1, l: -2, 2*l: 1}
jobs = [('pbar_8(32n+28) = 0 mod 32', 8, 32, 28, 32), ('pbar_6(24n+17) = 0 mod 64', 6, 24, 17, 64),
        ('pbar_6(24n+23) = 0 mod 64', 6, 24, 23, 64), ('pbar_9(36n+33) = 0 mod 128', 9, 36, 33, 128),
        ('pbar_4(36n+21) = 0 mod 9 ("left to the interested reader")', 4, 36, 21, 9)]
for name, l, m, t, u in jobs:
    ok, info = certificate(name, 2*l, r(l), m, t, u, out=out, claimed=[t])
# numerical sanity for the two infinite families proved by hand in R08
D = 60000; U = 2**20
def ser(l):
    def phim(k):
        c = [0]*D; n = 0
        while k*n*n < D: c[k*n*n] = 1 if n == 0 else (2*(-1)**n) % U; n += 1
        return flint.nmod_poly(c, U)
    return (phim(1).mul_low(phim(l), D)).inverse_series_trunc(D)
bad1 = [(l, n) for l in list(range(4, 80, 4)) + list(range(1, 80, 4)) for s in [ser(l)] for n in range(3, D, 4) if int(s[n]) % 8]
print('Conj 1 check: pbar_l(4n+3) mod 8 for l = 0,1 mod 4, l < 80, 4n+3 <', D, ': failures', bad1[:5])
s9 = ser(9)
bad3 = [(a, n) for a in range(0, 12) for n in range(0, D) if (2**a)*(8*n+7) < D and int(s9[(2**a)*(8*n+7)]) % 16]
print('Conj 3 check: pbar_9(2^a(8n+7)) mod 16, a < 12, argument <', D, ': failures', bad3[:5])
for l in (2, 3, 5, 6, 7):
    s = ser(l); f = [n for n in range(3, 400, 4) if int(s[n]) % 8]
    print(f'   (control) l = {l}: pbar_l(4n+3) mod 8 nonzero at', f[:4])
