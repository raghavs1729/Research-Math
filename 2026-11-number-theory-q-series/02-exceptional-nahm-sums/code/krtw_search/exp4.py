import sys
from lib2 import *
N = int(sys.argv[1]) if len(sys.argv)>1 else 30
K.setN(N)
def mx_inf(e, base=2):   # (-x q^e; q^base)_inf
    return xpoch_inf(+1, 1, e, base)
def mx_fin_inv(e, n, base=2):  # 1/(-x q^e;q^base)_n
    return xpoch_inv(+1, 1, e, n, base)
def Vsum(c, e, xs=0):
    """sum_V (x q^xs)^{2V} q^{2V^2 + cV} / ((Q)_V (-x q^{e+xs}... ) ) : returns sum_V x^{2V} q^{2V^2+cV+2 xs V}/((Q)_V (-x q^e;Q)_V)"""
    tot = S(); V = 0
    while 2*V*V + c*V + 2*xs*V < K.N:
        tot = tot + (mx_fin_inv(e, V) * iqpoch(V, 2)).shift(2*V, 2*V*V + c*V + 2*xs*V)
        V += 1
    return tot
B3 = bress(3)
print("B_3 = (-xq;q)_inf * Vsum(0,2):", B3 == mx_inf(1, 1)*Vsum(0, 2))
J2 = Cchain(2)*Sinv(Cchain(1))
print("J_2 = Vsum(0,1):", J2 == Vsum(0, 1))
B3ins = bress(3, lambda Ns: Ns[2])
for e0 in (1, 2, 3):
    for c in range(0, 5):
        for e in range(1, 6):
            if B3ins == mx_inf(e0, 1)*Vsum(c, e): print("  <q^{N_3}>_{B_3} = (-xq^%d;q)_inf * sum_V x^{2V} q^{2V^2+%dV}/((Q)_V(-xq^%d;Q)_V)" % (e0, c, e))
B3ins2 = bress(3, lambda Ns: Ns[1]+Ns[2])
for e0 in (1, 2, 3):
    for c in range(0, 5):
        for e in range(1, 6):
            if B3ins2 == mx_inf(e0, 1)*Vsum(c, e): print("  <q^{N_2+N_3}>_{B_3} = (-xq^%d;q)_inf * sum_V x^{2V} q^{2V^2+%dV}/((Q)_V(-xq^%d;Q)_V)" % (e0, c, e))
