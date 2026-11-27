from mpmath import mp, mpf, pi, zeta, log, pslq, factorial
import pickle
mp.dps = 95
D = pickle.load(open('u3_coords.pkl', 'rb'))
Ls = {k: mpf(v) for k, v in pickle.load(open('f8_L.pkl', 'rb')).items()}
L1, L2 = Ls['L1'], Ls['L2']
c8 = 8 * log(2)
def hat(c):
    c = [mpf(x) for x in c]
    return [sum(c[j] * c8**(j - i) / factorial(j - i) for j in range(i, 4)) for i in range(4)]
z3 = zeta(3)
for nm in ('nu', 'r2', 'rs', 'r0a'):
    h = hat(D[nm])
    print(nm, [mp.nstr(x, 20) for x in h])
    for i in (3, 2, 1, 0):
        best = None
        for a in range(-6, 5):
            vec = [h[i] * pi**a, L1, L2, L1 * pi, L2 * pi, L1 * pi**2, L2 * pi**2, mpf(1), pi, pi**2]
            r = pslq(vec, maxcoeff=10**6, maxsteps=100000)
            if r and r[0] != 0:
                best = (a, r); break
        if best is None and i <= 0:
            for a in range(-6, 5):
                vec = [h[i] * pi**a, L1, L2, L1 * pi**2, L2 * pi**2, z3 * L1, z3 * L2, z3, mpf(1), pi**2, z3 * L1 / pi**2, z3 * L2 / pi**2, L1 * pi**4, L2 * pi**4]
                r = pslq(vec, maxcoeff=10**6, maxsteps=200000)
                if r and r[0] != 0:
                    best = (a, r, 'ext'); break
        print('   hat c_%d' % i, best, flush=True)
