# arXiv:2503.17571 (Berkovich, Dhar), Conjectures conj41 and conj26.
#  conj41:  (1/(q;q^2)_L) ( (1-q^(2L-2))/(1-q^2) + q^3 [L-2,2]_{q^2} )  >=  q [L-2,2]_q (-q^3;q)_{L-4}   coefficientwise
#  conj26:  generating function of Theorem thm24 minus that of Theorem thm25 has nonnegative coefficients.
import sys
from flint import fmpz_poly as Z
def trunc(p, D): return Z(p.coeffs()[:D])
def x(k): 
    c = [0]*(k+1); c[k] = 1; return Z(c)
def qbin(n, k, step=1):           # Gaussian binomial [n,k]_{q^step} as a polynomial
    if k < 0 or k > n: return Z([0])
    num = Z([1]); den = Z([1])
    for i in range(1, k+1):
        num *= (Z([1]) - x(step*(n-k+i))); den *= (Z([1]) - x(step*i))
    q_, r_ = divmod(num, den); assert r_ == Z([0]); return q_
def geo(top, step):                # (1 - q^top)/(1 - q^step), top multiple of step (top = 0 gives 0)
    return Z([1 if i % step == 0 else 0 for i in range(top)]) if top > 0 else Z([0])
def inv_odd(L, D):                 # 1/(q;q^2)_L mod q^D
    c = [0]*D; c[0] = 1
    for i in range(L):
        k = 2*i+1
        for n in range(k, D): c[n] += c[n-k]
    return Z(c)
def test(L):
    rhs41 = x(1) * qbin(L-2, 2) * (Z([1]) if L < 4 else Z([1]))
    if L >= 4:
        prod = Z([1])
        for i in range(L-4): prod *= (Z([1]) + x(3+i))
        rhs41 = x(1) * qbin(L-2, 2) * prod
    else:
        rhs41 = Z([0])                       # [L-2,2]_q = 0 for L < 4
    D = max(rhs41.degree(), 0) + 8*L + 60
    io = inv_odd(L, D)
    lhs41 = trunc(io * (geo(2*L-2, 2) + x(3)*qbin(L-2, 2, 2)), D)
    d41 = lhs41 - rhs41
    neg41 = [n for n, v in enumerate(d41.coeffs()) if int(v) < 0]
    # conj26
    br = x(5)*geo(2*L-2, 2) + (x(6)+x(8))*geo(4*L-4, 4) + x(8)*qbin(L-1, 2, 2) + x(16)*qbin(L-1, 2, 4)
    # (1-q^2)/(1-q^(2L)) [L,3]_{q^2} = [L-1,2]_{q^2} (1-q^2)/(1-q^6); expand 1/(1-q^6) as a series mod q^D
    t3 = qbin(L-1, 2, 2) * (Z([1]) - x(2)) * ((Z([1])+x(2))*(Z([1])+x(2*L-2)) + x(2) + x(2*L))
    q3 = trunc(t3 * Z([1 if i % 6 == 0 else 0 for i in range(D)]), D)
    br = br + x(11)*q3
    lhs26 = trunc(io * br, D); rhs26 = x(5)*rhs41            # q^6 [L-2,2]_q (-q^3;q)_{L-4}
    d26 = lhs26 - rhs26
    neg26 = [n for n, v in enumerate(d26.coeffs()) if int(v) < 0]
    return neg41, neg26, rhs41.degree()
if __name__ == '__main__':
    bad = []
    for L in range(1, int(sys.argv[1])+1):
        n41, n26, deg = test(L)
        if n41 or n26: bad.append((L, n41[:5], n26[:5]))
    print(f'L <= {sys.argv[1]}: failures (L, first n with negative difference for conj41, for conj26):', bad[:10])
