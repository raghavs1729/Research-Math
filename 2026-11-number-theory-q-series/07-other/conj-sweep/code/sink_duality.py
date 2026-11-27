# arXiv:2601.18715 (Additive sink subtraction): conjectured duality  v_wall([x]_p) = sigma(v_sink([-s_2 - x]_p)),
# S = {m, m+delta, 2m+delta}, sigma: 0->2, 1->1, 2->0, 3->3, p = p(m,delta) from the Main Theorem.
import sys
from math import gcd
def seqs(S, N):
    w = [0]*N; s = [0]*(N+1)              # w[x], x >= 0 (wall);  s[x], x >= 1 (sink), s[y] = 0 for y <= 0
    for x in range(N):
        o = {w[x-t] for t in S if x-t >= 0}; m = 0
        while m in o: m += 1
        w[x] = m
    for x in range(1, N+1):
        o = {(s[x-t] if x-t >= 1 else 0) for t in S}; m = 0
        while m in o: m += 1
        s[x] = m
    return w, s
def p_formula(m, dl):
    d = dl % (2*m)
    return 3*m + 2*dl - d if d <= m else m*(m + 2*dl + d)//gcd(m, d)
sig = {0: 2, 1: 1, 2: 0, 3: 3}
MM, DD = int(sys.argv[1]), int(sys.argv[2])
stats = {'sinkper': 0, 'wallper': 0, 'tot': 0}; fails = {0: [], 1: []}
for m in range(1, MM+1):
    for dl in range(1, DD+1):
        S = (m, m+dl, 2*m+dl); p = p_formula(m, dl); N = 4*p + 4*S[2] + 10
        w, s = seqs(S, N); stats['tot'] += 1
        sink_ok = all(s[x] == s[x+p] for x in range(1, N-p)); wall_ok = all(w[x] == w[x+p] for x in range(0, N-p))
        stats['sinkper'] += sink_ok; stats['wallper'] += wall_ok
        if not (sink_ok and wall_ok): continue
        for conv in (0, 1):      # sink residue representative in [1,p] (conv 0) or sink index shifted by one (conv 1)
            ok = True
            for x in range(p):
                y = (-S[1] - x + conv) % p
                if y == 0: y = p
                if w[x] != sig.get(s[y], -1): ok = False; break
            if not ok: fails[conv].append((m, dl))
print('pairs (m,delta):', stats['tot'], '; sink purely periodic with p:', stats['sinkper'], '; wall nim-values purely periodic with p:', stats['wallper'])
print('duality fails (convention [.]_p in 1..p):', len(fails[0]), fails[0][:10])
print('duality fails (sink index shifted by +1, i.e. v_sink(2p - s_2 - x + 1) as in the commented conjecture):', len(fails[1]), fails[1][:10])
