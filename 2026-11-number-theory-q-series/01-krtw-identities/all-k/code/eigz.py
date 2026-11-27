"""Eigenfunctions Phi_{N,M}(R;x)/phi_0(R;x) as rational functions of z = q^R: numerator polynomial P_{N,M}(z),
P(z) = sum_e c_e prod_{j=0}^{e-1}(1 - z q^{-j}) prod_{i=e+1}^{N+M} (1 - x z q^i),  c_e = gamma^{N-M,M}_e(x)."""
import sys
import flint
from tower import tower

C3 = flint.fmpz_mpoly_ctx.get(('x', 'q', 'z'), 'lex')
x3, q3, z3 = C3.gens()
ONE3 = C3.from_dict({(0, 0, 0): 1})


def lift(p):
    """fmpz_mpoly in (x,q) -> (x,q,z)"""
    d = {}
    for (i, j), c in p.to_dict().items():
        d[(i, j, 0)] = c
    return C3.from_dict(d)


def Pz(N, M):
    g = tower(N - M, M)
    # common denominator
    den = ONE3
    for e, v in g.items():
        dv = lift(v.d)
        den = den * dv / den.gcd(dv)
    tot = C3.from_dict({})
    for e, v in g.items():
        num = lift(v.n) * (den / lift(v.d))
        f = ONE3
        for j in range(e):
            # (1 - z q^{-j}) = q^{-j}(q^j - z): multiply everything by q^{binom}: track separately
            f = f * (q3 ** j - z3)
        # we multiplied by q^{sum j} = q^{e(e-1)/2}; compensate by multiplying the others with q^{(NM)} later
        shiftexp = e * (e - 1) // 2
        g2 = ONE3
        for i in range(e + 1, N + M + 1):
            g2 = g2 * (ONE3 - x3 * z3 * q3 ** i)
        Emax = (N + M) * (N + M - 1) // 2
        tot = tot + num * f * g2 * q3 ** (Emax - shiftexp)
    return tot, den   # P(z) * q^{Emax} = tot/den


if __name__ == '__main__':
    for (N, M) in [(1, 1), (2, 1), (3, 1), (2, 2), (3, 2), (4, 2), (3, 3)]:
        tot, den = Pz(N, M)
        g = tot.gcd(den)
        tot, den = tot / g, den / g
        c, fl = tot.factor()
        print('N=%d M=%d:' % (N, M))
        print('   num factors:', c, [(str(f), e) for f, e in fl])
        c2, fl2 = den.factor()
        print('   den factors:', c2, [(str(f), e) for f, e in fl2])
        sys.stdout.flush()
