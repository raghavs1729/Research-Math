"""Second audit: the low zero mod 19 of the Remark in section 4.2 (own characters, Hurwitz route)."""
from azero import PC
from achars import primitive_characters, order
import mpmath as mp
for v in primitive_characters(19):
    ch = PC(v, 19)
    z0 = ch.Zt(0)
    if abs(float(z0.real.mid())) < 0.1 * 19 ** 0:   # |Lambda| small: look closer
        f = lambda t: float(ch.Zt(t).real.mid())
        ts = [x / 1000 for x in range(-60, 61, 2)]
        vals = [f(t) for t in ts]
        sc = [(ts[i], ts[i + 1]) for i in range(len(ts) - 1) if vals[i] * vals[i + 1] < 0]
        if sc:
            g = mp.findroot(f, sum(sc[0]) / 2)
            print('chi(2) =', v[2], 'order', order(v), 'parity', ch.a, ' zero at gamma =', mp.nstr(g, 8), ' sign changes', sc)
