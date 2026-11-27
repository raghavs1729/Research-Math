"""Referee: (k, p, m) for S_{1,d-1,2d-1}, d = 3..12, from the C simulation (sumfree3_ref), and a PROOF of ultimate
periodicity for d = 3, 4, 5 via Lemma 'lemma' of arXiv:2609.16843 (t = 3): if the characteristic sequence has Q
consecutive identical non-null blocks of length M from position K+1 on and (Q - 6) M >= 2K, the blocks repeat forever.
Also: is there any triple (f, g, g+d) at all with (k, p, m) = (d+3, d, 8d-1)?"""
import subprocess, sys
BIN = sys.argv[1]
def seq(f, g, h, X): return list(map(int, subprocess.run([BIN, str(f), str(g), str(h), str(X)], capture_output=True, text=True).stdout.split()))
def kpm(S, X, maxm):
    """minimal (k, p, m) such that the characteristic sequence on [1, X] is m-periodic from position K+1, with
    lemma hypothesis; returns (k, p, m, K, Q) or None"""
    C = [0] * (X + 1)
    for s in S: C[s] = 1
    for m in range(1, maxm + 1):
        # smallest K with C[x] == C[x+m] for all K < x <= X - m
        K = X - m
        while K > 0 and C[K] == C[K + m]: K -= 1
        Q = (X - K) // m
        if (Q - 6) * m >= 2 * K and any(C[K + 1:K + 1 + m]):
            return (sum(C[1:K + 1]), sum(C[K + 1:K + 1 + m]), m, K, Q)
    return None
for d in range(3, 13):
    X = 400 * d
    S = seq(1, d - 1, 2 * d - 1, X)
    r = kpm(S, X, 30 * d)
    k, p, m, K, Q = r
    D = [S[i + 1] - S[i] for i in range(len(S) - 1)]
    print(f'd={d:2d}: (k,p,m) = ({k},{p},{m})  [(d+3, d, 8d-3) = ({d+3},{d},{8*d-3}); 8d-1 = {8*d-1}]  lemma: K={K}, Q={Q} blocks, (Q-6)M >= 2K: {(Q-6)*m >= 2*K}')
    if d <= 5: print('      differences:', D[:k + 2 * p])
# finite refutation for d = 4 that needs no periodicity proof: preperiod 7 / period 4 would force D[i+4] = D[i] for i >= 8
S = seq(1, 3, 7, 200); D = [None] + [S[i + 1] - S[i] for i in range(len(S) - 1)]
print('d=4: D_8..D_12 =', D[8:13], ' D_12 == D_8 ?', D[12] == D[8], '(period length 4 after 7 entries is impossible)')
# any triple with the printed parameters?
hits = []; alld = {}
for d in range(4, 13):
    for g in range(2, 2 * d + 3):
        for f in range(1, g):
            X = 150 * (g + d + f)
            S = seq(f, g, g + d, X)
            r = kpm(S, X, 12 * (g + d + f))
            if r and r[:3] == (d + 3, d, 8 * d - 1): hits.append((f, g, g + d))
            if r and f + g == d: alld[(f, g, g + d)] = r[:3]
print('triples (f,g,g+d), 4 <= d <= 12, f < g <= 2d+2, with (k,p,m) = (d+3, d, 8d-1):', hits)
print('family f+g = d: (k,p,m) vs (d+3, d, 8d-1-2f):', all(v == (h - g + 3, h - g, 8 * (h - g) - 1 - 2 * f) for (f, g, h), v in alld.items() if h - g >= 6 and f == 1),
      {t: v for t, v in alld.items() if t[2] - t[1] in (9, 10)})
