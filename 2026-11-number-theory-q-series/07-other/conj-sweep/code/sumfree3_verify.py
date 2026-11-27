# Check of the Theorem in R06 against the greedy simulation: for 6 <= d <= DMAX,
# S_{1,d-1,2d-1} = {1,d-1} U [2d-1,3d-2] U {4d-1} U U_j [10d-5+j(8d-3), 11d-6+j(8d-3)]   (first NT terms).
import sys
from sumfree3 import seq
DMAX = int(sys.argv[1]); bad = []
for d in range(6, DMAX+1):
    nt = 25*d + 50
    S = seq((1, d-1, 2*d-1), nt)
    T = [1, d-1] + list(range(2*d-1, 3*d-1)) + [4*d-1]
    j = 0
    while len(T) < nt:
        T += list(range(10*d-5 + j*(8*d-3), 11*d-5 + j*(8*d-3))); j += 1
    if S != T[:nt]: bad.append(d)
print(f'6 <= d <= {DMAX}: greedy sequence equals the claimed set on the first 25d+50 terms; mismatches: {bad}')
for d in (3, 4, 5):
    S = seq((1, d-1, 2*d-1), 40); print(f'd={d}: first terms {S[:26]}; differences {[S[i+1]-S[i] for i in range(25)]}')
