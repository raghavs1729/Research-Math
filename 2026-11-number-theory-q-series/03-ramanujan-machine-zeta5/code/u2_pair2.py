from mpmath import mp, mpf, pi, gamma, zeta, hyper, taylor, nsum, inf, psi
import pickle
C = {k: [mpf(x) for x in v] for k, v in pickle.load(open('u_consts.pkl', 'rb')).items()}
mp.dps = 60
h = mpf(1) / 2
Xi = lambda e: 2 * zeta(3, h + e) * (gamma(1 + e) / gamma(h + e))**4
Y = lambda e: hyper([1, h + e, h + e, h + e, h + e], [1 + e] * 4, 1)
mp.dps = 200
xi = taylor(Xi, 0, 4, singular=False, h=mpf(10)**-30); eta = taylor(Y, 0, 4, h=mpf(10)**-30)
mp.dps = 50
L, al1, be1, al2, be2 = C['U2']
alt3 = sum(xi[3 - j] * eta[j] * (-1)**j for j in range(4))
print('alt3', mp.nstr(alt3, 40))
print('be1/alt3', mp.nstr(be1 / alt3, 40), ' 16*be1/alt3 - 1 =', mp.nstr(16 * be1 / alt3 - 1, 5))
for k in range(0, 5):
    Tk = sum(xi[k - j] * eta[j] * (-1)**j for j in range(k + 1))
    print(k, mp.nstr(Tk, 30), 'al1/T', mp.nstr(al1 / Tk, 25), 'be1/T', mp.nstr(be1 / Tk, 25), 'al1*pi^2/T', mp.nstr(al1 * pi**2 / Tk, 25))
