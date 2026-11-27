from mpmath import mp, mpf, pi, gamma, zeta, pslq, exp, gammainc, sqrt, log
import pickle
mp.dps = 310
C = {k: [mpf(x) for x in v] for k, v in pickle.load(open('u_consts.pkl', 'rb')).items()}
N = 700
co = [0] * (N + 1); co[0] = 1
for n in range(1, N + 1):
    for (step, e) in ((2 * n, 4), (4 * n, 4)):
        if step > N: continue
        for _ in range(e):
            for i in range(N, step - 1, -1):
                co[i] -= co[i - step]
a = [0] + co[:N]
def Lval(s, k=4, lev=8, eps=1):
    A = sqrt(mpf(lev)) / (2 * pi)
    tot = mpf(0)
    for n in range(1, N + 1):
        if a[n] == 0: continue
        x = n / A
        tot += a[n] * (gammainc(s, x) / (mpf(n) / A)**s + eps * gammainc(k - s, x) / (mpf(n) / A)**(k - s))
    return tot / (A**s * gamma(s))
L1, L2 = Lval(1), Lval(2)
L4 = sum(mpf(a[n]) / mpf(n)**4 for n in range(1, N + 1))  # abs. convergent but slow; use functional-eq free formula instead
L4, L5 = Lval(4), Lval(5)
print('L4,L5', mp.nstr(L4, 25), mp.nstr(L5, 25))
U2 = C['U2'][0]; U3 = C['U3'][0]; cw2 = C['U2'][1]; cw3 = C['U3'][1]
z3, z5 = zeta(3), zeta(5)
tests = {
 'A 1/U2': [1/U2, z5, pi**2*z3, pi**5, pi**4*L1, pi**3*L2, pi*L4, L5, z3*L2, z3*pi*L1, pi**4, pi**2, mpf(1)],
 'B 1/U2': [1/U2, mpf(1), L1**2, L1*L2, L2**2, L1**2/pi, L1*L2/pi, L2**2/pi, L1**2/pi**2, L1*L2/pi**2, L2**2/pi**2, L1**2*pi, L1*L2*pi, L2**2*pi],
 'C cw/U2': [cw2/U2, mpf(1), pi, pi**2, pi**3, z3, z3/pi, z3/pi**2, z3/pi**3, L1, L2, L1/pi, L2/pi, L1*pi, L2*pi, z5/pi**2, z5/pi**3, z5/pi**4],
 'D cw^2/U2': [cw2**2/U2, mpf(1), 1/pi, 1/pi**2, 1/pi**3, 1/pi**4, z3/pi**3, z3/pi**4, z5/pi**4, z5/pi**5, cw2**2, cw2**2*pi**2, cw2**2*z3],
 'E 1/U3 vs 1/U2, cw': [1/U3, 1/U2, mpf(1), 1/cw2**2, 1/cw2, 1/(cw2*pi), 1/(cw2**2*pi), 1/(cw2**2*pi**2), 1/(pi*cw2)**2*z3],
 'F cw3': [cw3, cw2, L1/pi**2, L2/pi**2, L1/pi, L2/pi**3, mpf(1), 1/pi, 1/pi**2],
}
for name, vec in tests.items():
    print(name, pslq(vec, maxcoeff=10**10, maxsteps=300000), flush=True)
