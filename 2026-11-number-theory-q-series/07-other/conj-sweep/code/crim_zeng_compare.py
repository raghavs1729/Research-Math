# Compare the 17-state automaton of R07 with the criterion of Z. Zeng's preprint (as summarised in NOVELTY-2.md):
#   w = E^{a_1} N^{b_1} ... E^{a_k} N^{b_k};  L = number of leading pairs with a_i = 1 and b_i odd;
#   R = number of trailing pairs with a_i odd and b_i = 1;  M = min(k, L+R);  losing iff |w| even and M even.
import json, itertools
D = json.load(open('../cert/crim_dfa.json')); delta = D['delta']; acc = D['acc']
def dfa(w):
    q = 0
    for x in w: q = delta[f'{q}{x}']
    return acc[q]
def zeng(w):
    if w == '': return True
    pairs = []; i = 0
    while i < len(w):
        a = 0
        while i < len(w) and w[i] == 'E': a += 1; i += 1
        b = 0
        while i < len(w) and w[i] == 'N': b += 1; i += 1
        pairs.append((a, b))
    k = len(pairs); L = 0
    while L < k and pairs[L][0] == 1 and pairs[L][1] % 2 == 1: L += 1
    R = 0
    while R < k and pairs[k-1-R][0] % 2 == 1 and pairs[k-1-R][1] == 1: R += 1
    return len(w) % 2 == 0 and min(k, L+R) % 2 == 0
bad = 0; cnt = 0
for Ln in range(2, 23):
    for mid in itertools.product('EN', repeat=Ln-2):
        w = 'E' + ''.join(mid) + 'N'; cnt += 1
        if dfa(w) != zeng(w): bad += 1
print(f'all {cnt} canonical boundary words of length <= 22: automaton of R07 vs Zeng\'s criterion: {bad} disagreements')
