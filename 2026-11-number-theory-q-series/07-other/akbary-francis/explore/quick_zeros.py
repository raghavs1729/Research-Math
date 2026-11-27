# quick UNCERTIFIED exploration: approximate zeros of all L(s,chi), chi mod q, 0<gamma<T (double precision)
import sys, json, math
from flint import *
q = int(sys.argv[1]); T = float(sys.argv[2])
ctx.prec = 64
G = dirichlet_group(q)
out = {}
def zeros_of(Z, T, h=1/64):
    zs = []
    n = int(T/h)
    prev = float(Z(0.0)); 
    for i in range(1, n+1):
        t = i*h
        cur = float(Z(t))
        if prev*cur < 0:
            a, b, fa = t-h, t, prev
            for _ in range(50):
                m = (a+b)/2; fm = float(Z(m))
                if fa*fm <= 0: b = m
                else: a, fa = m, fm
            zs.append((a+b)/2)
        prev = cur
    return zs
for n in range(1, q):
    if math.gcd(n, q) != 1: continue
    chi = dirichlet_char(q, n)
    if chi.is_principal():
        Z = lambda t: acb(0.5, t).zeta().real if False else None
        # zeta zeros
        zs = [float(z.imag) for z in acb.zeta_zeros(1, 2000) if float(z.imag) < T]
        out['zeta'] = zs
        continue
    assert chi.conductor() == q
    Z = lambda t, chi=chi: chi.hardy_z(acb(t)).real.mid()
    out[str(n)] = zeros_of(Z, T)
    print(n, chi.order(), chi.parity(), len(out[str(n)]), out[str(n)][:3], flush=True)
json.dump(out, open(f'quick_zeros_q{q}_T{int(T)}.json', 'w'))
