"""CRIM: further conjectures of arXiv:2606.16828 decided with the proved automaton (cert/crim_dfa.json).
 (1) "A meld of losing thick hooks is losing."     meld = concatenation of boundary words.
 (2) conj:bloby: if lambda and the hook Gamma_{r,c} are losing and Gamma_{r,c} <| lambda is defined, it is losing.
Exact: orbits of the DFA under E^k, N^k are eventually periodic, so finitely many exponents represent all.
"""
import json, itertools
from collections import deque
D = json.load(open('../cert/crim_dfa.json')); n = D['n']; acc = D['acc']
d = lambda q, x: D['delta'][f'{q}{x}']
def run(q, w):
    for x in w: q = d(q, x)
    return q
# exponents: for every state the sequence q, qx, qxx, ... enters a cycle; find K, PER with q x^(k+PER) = q x^k for all k >= K, all q
def orbit_bound(x):
    K = 0; PER = 1
    for q in range(n):
        seq = [q]
        while True:
            nx = d(seq[-1], x)
            if nx in seq:
                i = seq.index(nx); K = max(K, i); per = len(seq) - i
                from math import gcd
                PER = PER*per//gcd(PER, per); break
            seq.append(nx)
    return K, PER
KE, PE = orbit_bound('E'); KN, PN = orbit_bound('N')
print('E-orbits: preperiod <=', KE, 'period', PE, '; N-orbits: preperiod <=', KN, 'period', PN)
assert PE % 2 == 0 or PE == 1; assert PN % 2 == 0 or PN == 1
ME = KE + 2*max(PE, 2) + 2; MN = KN + 2*max(PN, 2) + 2      # exponents up to these bounds represent all (value class and parity)
def Emap(q, k): return run(q, 'E'*k)
# ---------- (1) melds of losing thick hooks
# thick hook Gamma^{a,b}_{r,c} = [c^a, b^(r-a)], r > a >= 2, c > b >= 2; losing iff r = c mod 2 (their Theorem thm:thick-hook).
# word: E^b N^(r-a) E^(c-b) N^a ;  exponents (b, s, t, a) with b >= 2, s >= 1, t >= 1, a >= 2, s + a = b + t mod 2.
H = [(b, s, t, a) for b in range(2, ME+2) for s in range(1, MN+2) for t in range(1, ME+2) for a in range(2, MN+2) if (s + a - b - t) % 2 == 0]
# sanity: single losing thick hooks are accepted, thick hooks of different parity are not
assert all(acc[run(0, 'E'*b + 'N'*s + 'E'*t + 'N'*a)] for b, s, t, a in H)
assert not any(acc[run(0, 'E'*b + 'N'*s + 'E'*t + 'N'*a)] for b in range(2, 8) for s in range(1, 7) for t in range(1, 7) for a in range(2, 8) if (s+a-b-t) % 2)
R = set(); dq = deque([0]); first = {}
while dq:
    q = dq.popleft()
    for h in H:
        t_ = run(q, 'E'*h[0] + 'N'*h[1] + 'E'*h[2] + 'N'*h[3])
        if t_ not in R:
            R.add(t_); first[t_] = (first.get(q, ()) + (h,)); dq.append(t_)
badm = [q for q in R if not acc[q]]
print('(1) states reachable from the start by melds of losing thick hooks:', sorted(R), '; non-accepting among them:', badm)
if badm:
    hs = first[badm[0]]
    print('    counterexample: meld of thick hooks with (b, r-a, c-b, a) =', hs)
# ---------- (2) conj:bloby
viol = None
for a in range(0, MN+2):
    s = run(0, 'E' + 'N'*a)
    # BFS over canonical words w (start with E): pair (state from 0, state from s)
    init = (d(0, 'E'), d(s, 'E'), 'E')
    seen = {init: 'E'}; dq = deque([init])
    while dq and viol is None:
        c = dq.popleft(); p, t, last = c
        if last == 'N' and acc[p]:             # w canonical and losing
            for b in range(a % 2, ME+2, 2):
                if not acc[run(t, 'E'*b + 'N')]: viol = (a, b, seen[c]); break
        for x in 'EN':
            c2 = (d(p, x), d(t, x), x)
            if c2 not in seen: seen[c2] = seen[c] + x; dq.append(c2)
print('(2) conj:bloby violation (a = r-1-k, b = c-1-lambda_1, boundary word of lambda):', viol)
