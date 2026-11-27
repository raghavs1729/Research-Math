# Brute force K_k = B_{k+1} C_k / C_1 for ALL powers of x mod q^N (X chosen so that every term of
# x-degree > X has q-exponent >= N for all series involved).
import sys, time, direct
N = int(sys.argv[1])
for k, X in ((2, int(sys.argv[2])), (3, int(sys.argv[3]))):
    t = time.time()
    Kk = direct.K(k, X, N); Bk = direct.B(k + 1, X, N); Ck = direct.C(k, X, N); C1 = direct.C(1, X, N)
    lhs = Kk * C1; rhs = Bk * Ck
    # top x-degrees must vanish mod q^N (confirms that X covers all x-powers)
    top = [i for i in range(X + 1) if not Kk.c[i].is_zero()]
    print("k=%d N=%d X=%d: K*C1 == B*C : %s ; max x-degree of K with nonzero coeff mod q^N: %d (%.1fs)"
          % (k, N, X, lhs == rhs, max(top), time.time() - t))
