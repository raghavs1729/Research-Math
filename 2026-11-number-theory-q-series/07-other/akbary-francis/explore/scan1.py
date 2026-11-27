import json, numpy as np, sys
q=13
d = json.load(open('quick_zeros_q13_T100.json'))
g = np.array(sorted(sum(d.values(), [])))
w = 2/(0.25+g*g)
print(len(g), 'partial mass', w.sum())
order = np.argsort(-w)
cum = np.cumsum(w[order])
for i in [0,1,2,3,4,5,9,19,29,39,59,99,199,399,len(g)-1]:
    print(i+1, g[order][i], w[order][i], cum[i])
F=6.45484
print('tail from table F:', F - w.sum())
# scan u in [0, U]
def B(u):
    return (w[None,:]*np.cos(np.outer(u, g))).sum(1)
best=(0,0)
for k in range(0, 2000):
    u = np.linspace(k*50, (k+1)*50, 5001)
    b = B(u)
    i = b.argmin()
    if b[i] < best[0]:
        best = (b[i], u[i]); print('new best', best, flush=True)
