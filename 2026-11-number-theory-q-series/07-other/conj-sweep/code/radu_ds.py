# arXiv:2611.05623 (Yaqubi, Mirzavaziri), Conjectures conj:A6 and conj:deep.
# A_k(n): coefficients of (f_6/(f_1 f_2 f_3))^k.
from radu import *
import sys
out = '../cert/radu_designated_summands.txt'
open(out, 'w').write('Radu finite-check certificates, arXiv:2611.05623 Conjectures 3 and 4\n\n')
def A(k): return {1: -k, 2: -k, 3: -k, 6: k}
# sanity: A_k(3) = k(k+1)(k+14)/6
for k in (1, 2, 5):
    s = eta_series(A(k), 10, 10**9+7); assert int(s[3]) == k*(k+1)*(k+14)//6, (k, s[3])
jobs = [('Conj 3: A_6(12n+10) = 0 mod 3^5', 6, 12, 10, 3**5),
        ('Conj 4(a): A_3(12n+7) = 0 mod 3^4', 3, 12, 7, 3**4),
        ('Conj 4(a): A_3(12n+11) = 0 mod 3^4', 3, 12, 11, 3**4),
        ('Conj 4(b): A_4(6n+5) = 0 mod 3^3', 4, 6, 5, 3**3),
        ('Conj 4(b): A_7(12n+11) = 0 mod 3^3', 7, 12, 11, 3**3),
        ('Conj 4(c): A_8(18n+15) = 0 mod 3^5', 8, 18, 15, 3**5)]
for name, k, m, t, u in jobs:
    ok, info = certificate(name, 6, A(k), m, t, u, out=out)
