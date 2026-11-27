# minimal (preperiod, period, modulus) of the difference sequence of S_{1,d-1,2d-1}
import sys
from sumfree3 import seq
def analyse(d, nterms=None):
    n = nterms or 60*d + 400
    S = seq((1, d-1, 2*d-1), n)
    D = [S[i+1]-S[i] for i in range(n-1)]
    best = None
    for p in range(1, 8*d+10):
        # largest suffix with period p
        k = len(D) - p
        while k > 0 and D[k-1] == D[k-1+p]: k -= 1
        if k < len(D)//2 and k < len(D) - 4*p - 50:   # periodic over at least the second half of the computed terms
            best = (k, p, S[k+p] - S[k]); break
    return best   # k = number of difference terms before periodicity = preperiod length in the paper's sense
if __name__ == '__main__':
    res = {}
    for d in range(int(sys.argv[1]), int(sys.argv[2])+1):
        if d == 2: continue
        res[d] = analyse(d)
    exc = {d: v for d, v in res.items() if v != (d+3, d, 8*d-3)}
    print('(preperiod, period, modulus) for d=3..12:', {d: res[d] for d in res if d <= 12})
    print('d with (preperiod, period, modulus) != (d+3, d, 8d-3):', exc)
    print('d with modulus == 8d-1:', [d for d, v in res.items() if v and v[2] == 8*d-1])
