# arXiv:2609.16843, Conjecture conj:except: the 3-sumfree sequence S_{1,d-1,2d-1} has preperiod length d+3 and
# period length d modulo 8d-1 (d >= 3, d != 3, 5).  Simulation with big-integer bitsets.
import sys
def seq(start, nterms):
    S = list(start); el = 0; two = 0; three = 0
    def add(x):
        nonlocal el, two, three
        three |= two << x; two |= el << x; el |= 1 << x
    for x in S: add(x)
    x = S[-1]
    while len(S) < nterms:
        x += 1
        if not (three >> x) & 1:
            S.append(x); add(x)
    return S
def analyse(d, periods=40):
    n = (periods + 3)*d + 20
    S = seq((1, d-1, 2*d-1), n)
    # minimal preperiod p (number of initial terms) with S[i+d] - S[i] = 8d-1 for all i >= p (0-based), in range
    ok = [S[i+d] - S[i] == 8*d-1 for i in range(n-d)]
    p = n-d
    while p > 0 and ok[p-1]: p -= 1
    return p, n
if __name__ == '__main__':
    lo, hi = int(sys.argv[1]), int(sys.argv[2])
    exc = []
    for d in range(lo, hi+1):
        if d == 2: continue
        p, n = analyse(d)
        if p != d+3: exc.append((d, p, n))
    print(f'd in [{lo},{hi}]: exceptions to "preperiod d+3, then s(i+d) = s(i)+8d-1" (d, observed preperiod, terms computed):', exc)
