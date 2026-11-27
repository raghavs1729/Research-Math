from mpmath import mp, mpf, pi, gamma, zeta, pslq, exp, gammainc, sqrt, log, hyper, catalan
import pickle
mp.dps = 310
C = {k: [mpf(x) for x in v] for k, v in pickle.load(open('u_consts.pkl', 'rb')).items()}
# f8 = eta(2t)^4 eta(4t)^4, weight 4 level 8
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
L1, L2, L3 = Lval(1), Lval(2), Lval(3)
print('L(f8,1..3)', mp.nstr(L1, 30), mp.nstr(L2, 30), mp.nstr(L3, 30))
y01 = hyper([mpf(1)/2]*4, [1, 1, 1], 1)
print('4F3(1) vs 16 L2/pi^2:', mp.nstr(y01, 30), mp.nstr(16 * L2 / pi**2, 30))
print('L3/L1 * 1/pi^2 ratio', mp.nstr(L3 / L1 / pi**2, 30))
cw = C['U2'][1]; U2 = C['U2'][0]; U3 = C['U3'][0]; cw3 = C['U3'][1]
w1, w2 = L1 / pi, L2 / pi**2
tests = {
 'cw lin': [cw, mpf(1), 1/pi, 1/pi**2, 1/pi**3, 1/pi**4, w1, w2, w1/pi, w2/pi, w1/pi**2, w2/pi**2],
 'cw quad': [cw, w1*w2, w1**2, w2**2, w1*w2/pi, w1**2/pi, w2**2/pi, w1*w2/pi**2, w1**2/pi**2, w2**2/pi**2, 1/pi**2, 1/pi**3],
 'cw*pi^k quad': [cw*pi, cw*pi**2, cw*pi**3, cw, w1*w2, w1**2, w2**2],
 '1/U2 lin': [1/U2, mpf(1), zeta(3), zeta(5), pi**2, pi**4, w1, w2, w1*pi, w2*pi, w1*pi**2, w2*pi**2, w1*pi**3, w2*pi**3],
 'U2 ratio': [mpf(1), U2*w1**2, U2*w2**2, U2*w1*w2, w1**2, w2**2, w1*w2, U2],
}
for name, vec in tests.items():
    print(name, pslq(vec, maxcoeff=10**10, maxsteps=300000), flush=True)
pickle.dump({'L1': str(L1), 'L2': str(L2), 'L3': str(L3)}, open('f8_L.pkl', 'wb'))
