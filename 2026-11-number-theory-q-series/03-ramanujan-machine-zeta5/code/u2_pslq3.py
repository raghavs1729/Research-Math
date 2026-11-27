from mpmath import mp, mpf, pi, zeta, pslq
import pickle
mp.dps = 300
C = {k: [mpf(x) for x in v] for k, v in pickle.load(open('u_consts.pkl', 'rb')).items()}
Ls = {k: mpf(v) for k, v in pickle.load(open('f8_L.pkl', 'rb')).items()}
L1, L2 = Ls['L1'], Ls['L2']
z3, z5 = zeta(3), zeta(5)
for name in ('U2', 'U3'):
    L, al1, be1, al2, be2 = C[name]
    for tname, target in (('c_q=be1', be1), ('1/U', 1 / L), ('al1', al1), ('al2', al2), ('be2', be2)):
        for W in range(-4, 9):
            X = [z3 * pi**(W - 3), z5 * pi**(W - 5), pi**W, z3**2 * pi**(W - 6)]
            vec = [target] + [x * L1 for x in X] + [x * L2 for x in X]
            r = pslq(vec, maxcoeff=10**9, maxsteps=200000)
            if r and r[0] != 0:
                print(name, tname, 'W=', W, r, flush=True)
    print(name, 'done', flush=True)
