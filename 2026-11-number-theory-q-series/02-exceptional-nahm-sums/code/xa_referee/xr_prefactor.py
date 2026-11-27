"""Lemma 4 exactly: with Kac-normalised form, rho-bar and h^vee computed from the root data (xr_lib.Affine),
   tr D(X)/24 + m_{kL0} == -c_SW(X,A_{k-1})/24, m = (rho|rho)(1/(2(k+h^v)) - 1/(2h^v)); and h^v(g(X)) == h(X);
   and prod_j (1-q^j)^{mult j delta} == prod_i (q^{d_i};q^{d_i})_inf (to q^60)."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from fractions import Fraction as Fr
import numpy as np
from xr_lib import Affine, coxeter
types = [('A',l) for l in range(1,8)] + [('D',l) for l in range(4,9)] + [('E',6),('E',7),('E',8)] + \
        [('B',l) for l in range(3,7)] + [('C',l) for l in range(2,7)] + [('F',4),('G',2)] + [('T',l) for l in range(1,7)]
allok = True
for X, l in types:
    A = Affine(X, l)
    h = coxeter(X, l); trD = sum(A.D)
    ok_h = (A.hv == h)
    N = 60
    H = A.heis_series(N * A.DU)
    P = np.zeros(N * A.DU + 1, dtype=object); P[0] = 1
    for d in A.D:
        for j in range(1, N + 1):
            s = d * j * A.DU
            if s <= N * A.DU:
                P[s:] = P[s:] - P[:N * A.DU + 1 - s].copy()
    ok_p = all(H[i] == P[i] for i in range(N * A.DU + 1))
    ok_c = True
    for k in range(2, 13):
        m = A.rho2 * (Fr(1, 2 * (k + A.hv)) - Fr(1, 2 * A.hv))
        c = Fr(trD * (k - 1) * h, h + k)
        if Fr(trD, 24) + m != -c / 24: ok_c = False
    lem = A.rho2 / (2 * A.hv) == Fr(trD * (h + 1), 24)
    allok &= ok_h and ok_p and ok_c and lem
    print(f"{X}_{l}: kind={A.kind} h^v(g)={A.hv} h(X)={h} [{'ok' if ok_h else 'MISMATCH'}]; (rho|rho)={A.rho2}; "
          f"(rho|rho)/2h^v == trD(h+1)/24: {lem}; Heisenberg product == prod (q^d;q^d): {ok_p}; (3.3) for k=2..12: {ok_c}")
print("ALL OK" if allok else "SOME FAILURE")
