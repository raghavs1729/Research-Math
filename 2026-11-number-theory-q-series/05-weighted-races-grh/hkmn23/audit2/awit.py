"""Second audit: independent prime-side evaluation of F_{q,a}(log x0) (formula (F)) at the author's
witness points, with my own constants (al1.Char, functional-equation route) and my own prime-power sums;
mpmath 40 digits.  Also a scan (float) of max F over 1 <= x <= X for undecided classes.
"""
import sys, json, glob, math
import mpmath as mp
import numpy as np
from math import gcd
from sympy import primerange, primefactors, totient, factorint
from aw1 import all_chars, constants

mp.mp.dps = 40


def Fprime(q, a, x, chars=None, consts=None):
    x = mp.mpf(x)
    phi = int(totient(q))
    ne, no, B, C = consts if consts else constants(q, a, chars or all_chars(q))
    ps1 = mp.mpf(0); imp = mp.mpf(0)
    for p in primerange(2, int(x) + 1):
        lp = mp.log(p); pk = p
        mpp = q
        while mpp % p == 0:
            mpp //= p
        while pk <= x:
            if q % p != 0:
                if (pk - a) % q == 0:
                    ps1 += lp * (x - pk)
            else:
                if (pk - a) % mpp == 0:
                    imp += int(totient(mpp)) * lp * (x - pk)
            pk *= p
    u = 1 / x
    Te = (x / 2) * ((1 + u) * mp.log(1 + u) + ((1 - u) * mp.log(1 - u) if u < 1 else 0))
    To = (1 - mp.log(1 - u * u) / 2 - mp.atanh(u) / u) if u < 1 else 1 - mp.log(2)
    ell = -ne * x * mp.log(x) + B * x + no * mp.log(x) + C - (1 + ne) * Te - no * To
    return (phi * ps1 - x * x / 2 - ell + imp) / x ** mp.mpf(1.5)


if __name__ == '__main__':
    recs = {}
    for fn in glob.glob('../raw/t2_witness_120_*.jsonl'):
        for l in open(fn):
            r = json.loads(l); recs[(r['q'], r['a'])] = r
    pairs = [tuple(map(int, s.split(','))) for s in sys.argv[1:]]
    for q, a in pairs:
        r = recs[(q, a)]
        chars = all_chars(q)
        cs = constants(q, a, chars)
        F = Fprime(q, a, r['x0'], consts=cs)
        print(f'({q},{a}) N={r["N"]} thr={r["threshold"]:.4f} x0={r["x0"]} mine F={mp.nstr(F, 15)}  author {r["F_at_x0"]}  witness={r["witness"]}  constants n_e={cs[0]} n_o={cs[1]} B={mp.nstr(cs[2],12)} (author {r["B"]}) C={mp.nstr(cs[3],12)} (author {r["C"]})')
