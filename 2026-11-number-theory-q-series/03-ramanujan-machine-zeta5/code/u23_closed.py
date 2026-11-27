from mpmath import mp, mpf, pi, zeta, pslq
import pickle
mp.dps = 300
C = {k: [mpf(x) for x in v] for k, v in pickle.load(open('u_consts.pkl', 'rb')).items()}
Ls = {k: mpf(v) for k, v in pickle.load(open('f8_L.pkl', 'rb')).items()}
L1, L2 = Ls['L1'], Ls['L2']
z5, z3 = zeta(5), zeta(3)
U2 = C['U2'][0]; U3 = C['U3'][0]
cand = 1 / (mpf(93) / 4 * z5 - pi**4 / 8 * L2 / L1)
print('U2 - 8 L1/(186 z5 L1 - pi^4 L2) =', mp.nstr(U2 - cand, 5))
print('pslq 1/U2:', pslq([1 / U2, z5, pi**4 * L2 / L1], maxcoeff=10**6, maxsteps=100000))
r = L2 / L1
for name in ('U2', 'U3'):
    L, al1, be1, al2, be2 = C[name]
    print(name)
    tests = {
     '1/U lin': [1 / L, mpf(1), z5, pi**4 * r, pi**2 * r, r, z3, pi**2, pi**4, pi**2 * z3],
     '1/U mobius in r': [1 / L, mpf(1), z5, pi**4 * r, r / L, pi**4 * r / L * 0 + pi**2 * r / L, pi**2 * r, r],
     'al1': [al1, L1 / pi**2, L2 / pi**2, L1 / pi**4, L2 / pi**4, L1, L2, L1 / pi, L2 / pi, L1 / pi**3, L2 / pi**3],
     'al2': [al2, L1 / pi**2, L2 / pi**2, L1 / pi**4, L2 / pi**4, L1, L2, L1 / pi, L2 / pi, L1 / pi**3, L2 / pi**3],
     'be1': [be1, z5 * L1 / pi**2, pi**2 * L2, z5 * L2 / pi**2, pi**2 * L1, L1, L2, z5 * L1, z5 * L2, pi**4 * L1, pi**4 * L2, z5 * L1 / pi**4, z5 * L2 / pi**4, L1 / pi**2, L2 / pi**2],
     'be2': [be2, z5 * L1 / pi**2, pi**2 * L2, z5 * L2 / pi**2, pi**2 * L1, L1, L2, z5 * L1, z5 * L2, pi**4 * L1, pi**4 * L2, z5 * L1 / pi**4, z5 * L2 / pi**4, L1 / pi**2, L2 / pi**2],
    }
    for tn, vec in tests.items():
        print('  ', tn, pslq(vec, maxcoeff=10**9, maxsteps=300000), flush=True)
