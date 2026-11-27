from mpmath import mp, mpf, pi, zeta, log, matrix, lu_solve, qr_solve
import pickle
mp.dps = 95
D = pickle.load(open('u3_coords.pkl', 'rb'))
Ls = {k: mpf(v) for k, v in pickle.load(open('f8_L.pkl', 'rb')).items()}
L1, L2 = Ls['L1'], Ls['L2']
C = {k: [mpf(x) for x in v] for k, v in pickle.load(open('u_consts.pkl', 'rb')).items()}
U2, U3 = C['U2'][0], C['U3'][0]
cnu, cr2, cr0 = [[mpf(x) for x in D[k]] for k in ('nu', 'r2', 'r0a')]
c8 = 8 * log(2)
e = matrix([-c8, 1, 0, 0])            # hat y_1 = y_1 - log(256) y_0 in Frobenius coordinates (y_0..y_3)
A = matrix(4, 3)
for i in range(4):
    A[i, 0] = cnu[i]; A[i, 1] = cr2[i]; A[i, 2] = cr0[i]
sol, res = qr_solve(A, e)
a, b, c = sol
print('hat y_1 = a nu + b r_2 + c r_0 :  a =', mp.nstr(a, 30), ' b =', mp.nstr(b, 30), ' c =', mp.nstr(c, 30), ' residual', mp.nstr(res, 5))
print('c + 16 L1 =', mp.nstr(c + 16 * L1, 5))
lhs = 1 / U3 - 1 / U2 + 4
rhs = pi**4 / (512 * L1 * (L1 - b))
print('1/U3 - 1/U2 + 4 =', mp.nstr(lhs, 40)); print('pi^4/(512 L1 (L1 - b)) =', mp.nstr(rhs, 40))
print('U2 formula check:', mp.nstr(1 / U2 - (mpf(93) / 4 * zeta(5) - pi**4 / 8 * L2 / L1), 5))
# hat y_0 = y_0 decomposition (has a nu~ component, so not in span): skip.  sigma-type: second-order datum of other H elements
print('ell_s(hat y_1) = b =', mp.nstr(b, 60))
pickle.dump({'b': str(b), 'a': str(a)}, open('u3_ls.pkl', 'wb'))
