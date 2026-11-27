from mpmath import mp, mpf, pi, gamma, zeta, pslq, log, taylor, diff
import pickle
mp.dps = 120
C = {k: [mpf(x) for x in v] for k, v in pickle.load(open('u_consts.pkl', 'rb')).items()}
U2 = C['U2'][0]; U3 = C['U3'][0]
Xi = lambda e: 2 * zeta(3, mpf(1) / 2 + e) * (gamma(1 + e) / gamma(mpf(1) / 2 + e))**4
mp.dps = 160
xi = taylor(Xi, 0, 5)
mp.dps = 110
xi = [+x for x in xi]
print([mp.nstr(x, 20) for x in xi])
for a in range(-4, 6):
    vec = [1 / U2] + [x * pi**a for x in xi[:4]]
    r = pslq(vec, maxcoeff=10**8, maxsteps=100000)
    print('a=', a, r)
# more general: allow log 2 normalisation shifts: Xi_c(e) = 2^(c e) Xi(e)
for cexp in (8, -8, 4, -4, 16, -16, 2, -2):
    mp.dps = 160
    xi2 = taylor(lambda e: mpf(2)**(cexp * e) * Xi(e), 0, 5)
    mp.dps = 110
    for a in range(-4, 6):
        vec = [1 / U2] + [x * pi**a for x in xi2[:4]]
        r = pslq(vec, maxcoeff=10**8, maxsteps=100000)
        if r: print('c=', cexp, 'a=', a, r)
