from mpmath import mp, mpf, pi, gamma, zeta, pslq, hyper, taylor, identify
import pickle
C = {k: [mpf(x) for x in v] for k, v in pickle.load(open('u_consts.pkl', 'rb')).items()}
mp.dps = 60
h = mpf(1) / 2
Xi = lambda e: 2 * zeta(3, h + e) * (gamma(1 + e) / gamma(h + e))**4
Y = lambda e: hyper([1, h + e, h + e, h + e, h + e], [1 + e] * 4, 1)
mp.dps = 90
xi = taylor(Xi, 0, 4); eta = taylor(Y, 0, 4)
mp.dps = 50
print('xi', [mp.nstr(x, 15) for x in xi]); print('eta', [mp.nstr(x, 15) for x in eta])
L, al1, be1, al2, be2 = C['U2']
for k in range(0, 5):
    Tk = sum(xi[k - j] * eta[j] * (-1)**j for j in range(k + 1))
    Sk = sum(xi[k - j] * eta[j] for j in range(k + 1))
    print(k, 'alt', mp.nstr(Tk, 20), ' plain', mp.nstr(Sk, 20))
    for nm, T in (('alt', Tk), ('plain', Sk)):
        for tn, tgt in (('be1', be1), ('al1', al1), ('1/U2', 1 / L)):
            for a in range(-6, 7):
                r = tgt / (T * pi**a)
                # is r a simple rational?
                from fractions import Fraction
                fr = Fraction(float(r)).limit_denominator(5000)
                if abs(r - mpf(fr.numerator) / fr.denominator) < mpf(10)**-30:
                    print('   HIT', nm, k, tn, 'a=', a, fr)
