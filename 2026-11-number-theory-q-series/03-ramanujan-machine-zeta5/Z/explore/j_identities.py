from mpmath import *
mp.dps = 30
Kc = lambda th: pi / (2 * agm(1, cos(th)))        # K(sin th)
Kpc = lambda th: pi / (2 * agm(1, sin(th)))      # K'(sin th)
h = pi / 2
pts = [0, h / 10**4, h / 100, h / 4, h / 2, 3 * h / 4, h - h / 100, h - h / 10**4, h]
def Q(f):
    return quad(lambda u: f(u), pts)
W = lambda th: (sin(th) + cos(th)) / sqrt(2 * sin(th) * cos(th))       # w(k) dk = W dth
J1l = Q(lambda th: Kc(th) * Kpc(th) * W(th)); J1r = 4 * Q(lambda th: Kc(th) * Kpc(th) * cos(th))
J2l = Q(lambda th: Kc(th)**2 * W(th)); J2r = 4 * Q(lambda th: Kpc(th)**2 * cos(th))
print("J1:", J1l, J1r, J1l / J1r)
print("J2:", J2l, J2r, J2l / J2r)
# q-expansion of G2(tau) = th2 th4 th3^4 (th2^2 + th4^2) in q4 = e^{i pi tau/4}
N = 240     # powers of q4 up to N  (th3 = sum q^(n^2) = sum q4^(4 n^2); th2 = sum q^((n+1/2)^2) = sum q4^((2n+1)^2))
def ser(f):
    c = [0] * (N + 1)
    for e, v in f: 
        if e <= N: c[e] += v
    return c
th3 = ser([(4 * n * n, 1) for n in range(-10, 11)]); th4 = ser([(4 * n * n, (-1)**n) for n in range(-10, 11)]); th2 = ser([((2 * n + 1)**2, 1) for n in range(-10, 10)])
def mul(a, b):
    r = [0] * (N + 1)
    for i, x in enumerate(a):
        if x:
            for j in range(N + 1 - i):
                if b[j]: r[i + j] += x * b[j]
    return r
add = lambda a, b: [x + y for x, y in zip(a, b)]
G2 = mul(mul(th2, th4), mul(mul(mul(th3, th3), mul(th3, th3)), add(mul(th2, th2), mul(th4, th4))))
print("G2 nonzero exponents/coeffs:", [(e, c) for e, c in enumerate(G2) if c][:30])
G1 = mul(mul(th2, th2), mul(mul(th3, th3), mul(mul(th4, th4), mul(th4, th4))))
print("G1:", [(e, c) for e, c in enumerate(G1) if c][:20])
