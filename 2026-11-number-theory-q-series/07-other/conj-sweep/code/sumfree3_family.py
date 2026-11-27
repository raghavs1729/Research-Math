# Boundary family f + g = d of 3-sumfree sequences S_{f,g,g+d} (excluded from the source's Theorem thm:gend>1, which needs f+g >= d+1).
import sys
from sumfree3 import seq
def runs(S):
    out = []; a = b = S[0]
    for x in S[1:]:
        if x == b+1: b = x
        else: out.append((a, b)); a = b = x
    out.append((a, b)); return out
def analyse(f, g, d, nper=30):
    n = nper*d + 200
    S = seq((f, g, g+d), n); D = [S[i+1]-S[i] for i in range(n-1)]
    for p in range(1, 12*d+50):
        k = len(D) - p
        while k > 0 and D[k-1] == D[k-1+p]: k -= 1
        if k < len(D)//2 and k < len(D) - 4*p - 50: return (k, p, S[k+p]-S[k])
    return None
if __name__ == '__main__':
    FM, DM = int(sys.argv[1]), int(sys.argv[2])
    for f in range(1, FM+1):
        good = []; bad = []
        for d in range(2*f+1, DM+1):
            g = d - f
            if g <= f: continue
            r = analyse(f, g, d)
            (good if r == (d+3, d, 8*d-1-2*f) else bad).append((d, r))
        print(f'f={f}: pattern (d+3, d, 8d-1-2f) holds for d in {[d for d, _ in good][:3]}..{good[-1][0] if good else None} ({len(good)} values); exceptions: {bad[:8]}')
    print(runs(seq((2, 8, 18), 80))[:9]); print(runs(seq((3, 9, 21), 90))[:9])
