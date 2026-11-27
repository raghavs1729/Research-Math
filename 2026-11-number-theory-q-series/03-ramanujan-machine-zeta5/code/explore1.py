from fractions import Fraction as Fr
from math import comb
# n^2 v_n = (3n^2+3n+1) v_{n-1} - 2 n^2 v_{n-2},  two solutions
def sols(N, a, b, den, init):
    v=[Fr(x) for x in init]   # v_{-1}, v_0
    for n in range(1,N+1):
        v.append((a(n)*v[-1]+b(n)*v[-2])/den(n))
    return v[1:]
N=14
for init in ([0,1],[1,0]):
    v=sols(N, lambda n:3*n*n+3*n+1, lambda n:-2*n*n, lambda n:n*n, init)
    print([str(x) for x in v[:8]])
    c=[sum((-1)**(k-j)*comb(k,j)*v[j] for j in range(k+1)) for k in range(N+1)]
    print([str(x) for x in c[:10]])
    print([str(c[k+1]/c[k]) if c[k]!=0 else None for k in range(N)])
