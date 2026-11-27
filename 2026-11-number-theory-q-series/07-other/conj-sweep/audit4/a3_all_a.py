"""AUDIT-4 part 3 (own code, no monotonicity used): for every n in [N0,N1] and EVERY 1 <= a < n/2 decide exactly
      R_n(a) <,=,> 1   via the integer inequality (from the source's definition G/D, cleared of denominators):
      (n-2a+1) * a * 2(2n-1) * C(n,a)^2   vs   n (n-1) C(2n,n).
C(n,a)^2 is carried by the exact recurrence q_{a+1} = q_a (n-a)^2/(a+1)^2 and cross-checked at a = amax against
flint's independent binomial.  Output per n: sign of rho(n)-1, argmax a, number of a with R>=1.
usage: a3_all_a.py N0 N1 NPROC
"""
import sys, math
from multiprocessing import Pool
from flint import fmpz

def one(n):
    amax = (n-1)//2
    rhs = fmpz(n*(n-1))*fmpz.bin_uiui(2*n, n)
    k = 2*(2*n-1)
    q = fmpz(n*n)                      # C(n,1)^2
    ngt = 0; neq = 0
    best_a = 0; best_num = fmpz(0)     # argmax of (n-2a+1) a q
    for a in range(1, amax+1):
        if a > 1:
            q = q*((n-a+1)*(n-a+1))
            d = a*a
            assert q % d == 0
            q = q//d
        num = q*((n-2*a+1)*a)
        if num > best_num: best_num = num; best_a = a
        lhs = num*k
        if lhs > rhs: ngt += 1
        elif lhs == rhs: neq += 1
    assert q == fmpz.bin_uiui(n, amax)**2, n
    return (n, ngt, neq, best_a, amax)

if __name__ == "__main__":
    N0, N1, NP = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    # spot check of the integer criterion against the source definition with Fractions
    from fractions import Fraction as F
    for n in list(range(3, 80))+[495, 496, 497]:
        for a in range(1, (n-1)//2+1):
            Rdef = F(n-2*a+1, n-a)*math.comb(n, a)*math.comb(n-2, a-1)/F(math.comb(2*n-2, n-1), n)
            lhs = (n-2*a+1)*a*2*(2*n-1)*math.comb(n, a)**2; rhs = n*(n-1)*math.comb(2*n, n)
            assert (Rdef > 1) == (lhs > rhs) and (Rdef == 1) == (lhs == rhs)
    ns = list(range(N1, N0-1, -1))
    with Pool(NP) as p:
        res = sorted(p.imap_unordered(one, ns, chunksize=8))
    with open("a3_all_a_%d_%d.tsv" % (N0, N1), "w") as f:
        f.write("n\t#a_with_R>1\t#a_with_R=1\targmax_a\tamax\n")
        for r in res: f.write("\t".join(map(str, r))+"\n")
    gt = [r[0] for r in res if r[1] > 0]; eq = [r[0] for r in res if r[2] > 0]
    tot = sum(r[4] for r in res)
    print("range n=%d..%d ; pairs (n,a) tested exactly: %d" % (N0, N1, tot))
    print("n with rho(n) > 1: count=%d ; largest=%s ; smallest=%s" % (len(gt), max(gt) if gt else None, min(gt) if gt else None))
    print("n with some R_n(a) = 1:", eq)
    print("all n in [5,375] have rho>1:", all(n in set(gt) for n in range(max(5, N0), min(375, N1)+1)))
    print("n in 376..495 with rho(n) <= 1:", [n for n in range(376, 496) if N0 <= n <= N1 and n not in set(gt)])
    print("n in 3..4 with rho>1:", [n for n in gt if n < 5])
    print("rho(n) < 1 for ALL n in [496,%d]: %s" % (N1, all(r[1] == 0 and r[2] == 0 for r in res if r[0] >= 496)))
    print("max number of a with R_n(a)>1 over n:", max(res, key=lambda r: r[1])[:2][::-1])
