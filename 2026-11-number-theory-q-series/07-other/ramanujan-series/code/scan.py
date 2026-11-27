"""Signature scanning (Au's search procedure, arXiv:2312.14051 Sec. 3.4).

signature of a hypergeometric term u(n) (ratio u(n+1)/u(n) = z * prod (n+alpha_i)^(m_i)):
     (z, frozenset of (alpha_i mod 1 in (0,1], total multiplicity))
The summand G(n,0) = R(n,0) F(n,0) of a WZ pair has the signature of F(n,0) (R rational).
"""
from fractions import Fraction as Fr
from collections import Counter
from wz import *


def frac1(x):
    """representative in (0,1]"""
    r = x - (x.numerator // x.denominator)
    return r if r != 0 else Fr(1)


def term_signature(F, var='n', kval=0):
    """signature of n -> F(n, kval) (numeric parameters only)"""
    z = Fr(1)
    cnt = Counter()
    for b, L in F.bases:
        e = L.get(var, Fr(0))
        z *= b ** int(e)
    for L, p in F.gammas:
        cn = L.get(var, Fr(0))
        if cn == 0:
            continue
        c0 = L.get('1', Fr(0)) + L.get('k', Fr(0)) * kval
        assert all(v in (var, '1', 'k') for v in L), "symbolic parameters left"
        cn = int(cn)
        if cn > 0:
            z *= Fr(cn) ** (cn * p)
            for j in range(cn):
                cnt[frac1((c0 + j) / cn)] += p
        else:
            m = -cn
            z *= (Fr(-1) ** m / Fr(m) ** m) ** p
            for j in range(m):
                cnt[frac1((1 - c0 + j) / m)] -= p
    return z, frozenset((a, m) for a, m in cnt.items() if m != 0)


def z_of_shifts(T, params, shifts, K):
    """exponent only (independent of base point)"""
    z = Fr(1)
    for b, L in T.bases:
        e = sum(L.get(p, Fr(0)) * shifts[i] for i, p in enumerate(params)) + L.get('k', Fr(0)) * K
        if e.denominator != 1:
            return None
        z *= b ** int(e)
    for L, pw in T.gammas:
        cn = sum(L.get(p, Fr(0)) * shifts[i] for i, p in enumerate(params)) + L.get('k', Fr(0)) * K
        if cn.denominator != 1:
            return None
        cn = int(cn)
        if cn > 0:
            z *= Fr(cn) ** (cn * pw)
        elif cn < 0:
            m = -cn
            z *= (Fr(-1) ** m / Fr(m) ** m) ** pw
    return z


def sig_at(T, params, shifts, K, base, k0):
    """Pochhammer signature for base point (dict param->Fraction), k -> k0 + K n"""
    cnt = Counter()
    for L, pw in T.gammas:
        cn = sum(L.get(p, Fr(0)) * shifts[i] for i, p in enumerate(params)) + L.get('k', Fr(0)) * K
        if cn == 0:
            continue
        c0 = L.get('1', Fr(0)) + sum(L.get(p, Fr(0)) * base[i] for i, p in enumerate(params)) + L.get('k', Fr(0)) * k0
        cn = int(cn)
        if cn > 0:
            for j in range(cn):
                cnt[frac1((c0 + j) / cn)] += pw
        else:
            m = -cn
            for j in range(m):
                cnt[frac1((1 - c0 + j) / m)] -= pw
    return frozenset((a, m) for a, m in cnt.items() if m != 0)


def hyp_signature(z, upper, lower):
    """signature of z^n prod (u)_n / prod (l)_n"""
    cnt = Counter()
    for u in upper:
        cnt[frac1(Fr(u))] += 1
    for l in lower:
        cnt[frac1(Fr(l))] -= 1
    return Fr(z), frozenset((a, m) for a, m in cnt.items() if m != 0)
