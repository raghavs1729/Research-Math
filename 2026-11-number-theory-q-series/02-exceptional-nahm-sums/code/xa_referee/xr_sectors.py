"""Sector-by-sector check:  sum_n mult(k L0 + mu - n delta) q^n * prod_j (1-q^j)^{mult j delta}
      ==  q^{(mu|mu)/2k} * F_{X,k,[mu]}(q)        (Nahm sector sum with Sun-Wang's C(X), D(X))
Multiplicities by Freudenthal (xr_lib).  Usage: python3 xr_sectors.py X l k Nmax_units
"""
import sys, itertools, time
from fractions import Fraction as Fr
import numpy as np
sys.path.insert(0, __import__('os').path.dirname(__file__))
from xr_lib import Affine, Mults, nahm_sectors, inv_poch_series, ip


def run(X, l, k, Nmax, out):
    t0 = time.time()
    A = Affine(X, l)
    DU = A.DU
    M = Mults(A, k, Nmax)
    H = A.heis_series(Nmax)
    secs = nahm_sectors(A.G, A.D, k, Nmax, DU)
    # Nahm sector series (without the q^{|mu|^2/2k} factor): exponents E0 rational; store as dict
    def sector_series(c, mu):
        s = np.zeros(Nmax + 1, dtype=object)
        sh = Fr(ip(A.G, mu, mu), 2 * k)
        for p, E0, mup in secs.get(c, []):
            e = (sh + E0) * DU
            assert e.denominator == 1, ("exponent not in (1/DU)Z", X, l, k, mu, p, e)
            e = int(e)
            if e > Nmax:
                continue
            ser = np.zeros(Nmax + 1, dtype=object); ser[e] = 1
            for i in range(l):
                for s_ in range(1, k):
                    n = p[i * (k - 1) + s_ - 1]
                    if n:
                        ip_ = inv_poch_series(A.D[i], n, Nmax, DU)
                        ser = np.convolve(ser, ip_)[:Nmax + 1]
            s += ser
        return s
    nfail = 0; nchk = 0
    classes = list(itertools.product(range(k), repeat=l))
    for c in classes:
        reps = [tuple(c), tuple(x - k if i == 0 else x for i, x in enumerate(c))]
        # a third representative: shift by k b_l
        reps.append(tuple(x + k if i == l - 1 else x for i, x in enumerate(c)))
        for mu in reps:
            st = np.array(M.string(mu), dtype=object)
            lhs = np.convolve(st, H)[:Nmax + 1]
            rhs = sector_series(c, mu)
            nchk += 1
            if any(lhs[i] != rhs[i] for i in range(Nmax + 1)):
                nfail += 1
                out.write(f"  FAIL class {c} rep {mu}\n    lhs {list(lhs[:12])}\n    rhs {list(rhs[:12])}\n")
    out.write(f"{X}_{l} ({A.kind}, r={A.r}) k={k}: h^vee={A.hv} (rho|rho)={A.rho2} |alpha0|^2={A.a0norm}; "
              f"depth <= {Nmax}/{DU}; classes {len(classes)}, reps checked {nchk}, failures {nfail}; "
              f"#dominant mult entries {len(M.table)}; {time.time()-t0:.1f}s\n")
    out.flush()
    return nfail


if __name__ == '__main__':
    X, l, k, N = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
    run(X, l, k, N, sys.stdout)
