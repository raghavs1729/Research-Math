import sys
from lib2 import *
from exp4 import mx_inf, mx_fin_inv, Vsum
N = int(sys.argv[1]) if len(sys.argv)>1 else 30
K.setN(N)
j1 = Vsum(0, 1)
base = mx_inf(1, 1)*j1
binv = Sinv(base)
def show(name, Q, nd=5, nc=14):
    print(name)
    for d in sorted(Q.d)[:nd]:
        print("   x^%d:"%d, Q.d[d].coeffs()[:nc])
for i in range(0, 4):
    Qt = calK(2, i)*binv
    show("i=%d: calK_2(i) / ((-xq;q)_inf j_1(x))" % i, Qt)
    for c in range(0, 9):
        for e in range(1, 8):
            if Qt == Vsum(c, e): print("    == sum_V x^{2V} q^{2V^2+%dV}/((Q)_V(-xq^%d;Q)_V)" % (c, e))
