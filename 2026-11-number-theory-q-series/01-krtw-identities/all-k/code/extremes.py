import sys
from tower import tower
for M in range(1, 5):
    for L in range(0, 3):
        g = tower(L, M)
        print('L=%d M=%d  e=L   : %s' % (L, M, g[L].factor_str()))
        print('          e=L+2M: %s' % (g[L+2*M].factor_str()))
        sys.stdout.flush()
