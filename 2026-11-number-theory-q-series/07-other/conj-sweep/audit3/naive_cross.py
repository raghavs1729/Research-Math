# AUDIT-3: naive reference implementation (definition verbatim: triple loop over distinct previous entries)
# cross-checked against the C program sf3 on small cases, including all pairs with 2f+1 <= d < 6.
import subprocess, itertools, os
HERE = os.path.dirname(os.path.abspath(__file__))
def naive(f, g, h, X):
    S = [f, g, h]
    for x in range(h+1, X+1):
        Sx = set(S); found = False
        for i in range(len(S)):
            a = S[i]
            if 3*a >= x: break
            for j in range(i+1, len(S)):
                b = S[j]; c = x - a - b
                if c <= b: break
                if c in Sx: found = True; break
            if found: break
        if not found: S.append(x)
    return S
def csim(f, g, h, X):
    out = subprocess.run([os.path.join(HERE, 'sf3'), 'seq', str(f), str(g), str(h), str(X)], capture_output=True, text=True, check=True).stdout
    return [int(t) for t in out.split()]
cases = [(1,2,3,400), (1,2,4,400), (3,7,13,400), (4,7,8,400), (1,2,14,600)]
for f in range(1, 7):
    for d in range(2*f+1, 16):
        cases.append((f, d-f, 2*d-f, 10*d + 4*(8*d-2*f-1)))
bad = 0
for (f, g, h, X) in cases:
    if naive(f, g, h, X) != csim(f, g, h, X): bad += 1; print('MISMATCH', f, g, h)
print(f'naive Python vs sf3: {len(cases)} sequences compared, mismatches: {bad}')
