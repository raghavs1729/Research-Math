from mpmath import mp, mpf, pi, zeta, pslq, log, nsum, inf, binomial, hyper
import pickle
mp.dps = 60
Ls = {k: mpf(v) for k, v in pickle.load(open('f8_L.pkl', 'rb')).items()}
L1, L2 = Ls['L1'], Ls['L2']
d = pickle.load(open('u3_ls.pkl', 'rb')); b = mpf(d['b'])
C = {k: [mpf(x) for x in v] for k, v in pickle.load(open('u_consts.pkl', 'rb')).items()}
U2, U3 = C['U2'][0], C['U3'][0]
h = mpf(1) / 2
Z0 = hyper([h, h, h, h, 1], [1, 1, 1, 2], 1)
print('Z0 =', mp.nstr(Z0, 45))
bpred = (pi**2 * L1 * Z0 + 4 * L1 - h) / (32 * L2)
print('b - (pi^2 L1 Z0 + 4 L1 - 1/2)/(32 L2) =', mp.nstr(b - bpred, 5))
pred = 1 / U2 - 4 + pi**4 * L2 / (8 * L1 * (64 * L1 * L2 - 2 * pi**2 * L1 * Z0 - 8 * L1 + 1))
print('1/U3 - formula =', mp.nstr(1 / U3 - pred, 5))
for vec in ([Z0, mpf(1), L1, L2, L1 / pi, L2 / pi, L1 / pi**2, L2 / pi**2, 1 / pi, 1 / pi**2],
            [Z0, mpf(1), L1 * pi**2, L2 * pi**2, L1 / pi**2, L2 / pi**2, 1 / pi**2, 1 / pi**4, L1 / pi**4, L2 / pi**4, 1 / (pi**2 * L1), 1 / (pi**2 * L2), 1 / L1, 1 / L2],
            [Z0 * L1, Z0 * L2, mpf(1), L1, L2, 1 / pi**2, L1 / pi**2, L2 / pi**2, L1 * L2, L1**2, L2**2, L1 * L2 / pi**2, L1**2 / pi**2, L2**2 / pi**2]):
    print(pslq(vec, maxcoeff=10**6, maxsteps=200000))
