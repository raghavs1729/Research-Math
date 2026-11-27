"""Extract hypergeometric signatures (z, Pochhammer classes mod 1) of the summands of the survey
entries.  The survey scripts define entries with mpmath lambdas; we re-execute their source with
exact rational arithmetic (harmonic numbers replaced by fixed rational constants, which only adds a
polynomial factor to the summand) and reconstruct the term ratio as a rational function."""
import re, sys, math
from fractions import Fraction as Fr
from collections import Counter
from scan import frac1
from wz import pari, _to_fr


class _Dummy:
    def __getattr__(self, a):
        return 0

    def __setattr__(self, a, v):
        pass


class NotRational(Exception):
    pass


def _bad(*a, **k):
    raise NotRational()


def _binom(n, k):
    n = Fr(n); k = Fr(k)
    if n.denominator != 1 or k.denominator != 1:
        raise NotRational()
    n, k = int(n), int(k)
    if k < 0:
        return Fr(0)
    if n >= 0:
        return Fr(math.comb(n, k))
    # generalized
    r = Fr(1)
    for i in range(k):
        r *= Fr(n - i, i + 1)
    return r


def load_exact(path):
    src = open(path).read()
    src = re.sub(r'^from mpmath import\s*\(.*?\)', '', src, flags=re.M | re.S)   # multi-line imports
    src = re.sub(r'^from mpmath import.*$', '', src, flags=re.M)
    src = re.sub(r'^import mpmath.*$', '', src, flags=re.M)
    # cut off the "main" part
    src = src.split("if __name__")[0]
    ns = dict(mp=_Dummy(), mpf=Fr, binomial=_binom, pi=None, zeta=_bad, log=_bad, sqrt=_bad,
              catalan=None, nsum=_bad, inf=None, polylog=_bad, gamma=_bad, rf=None, factorial=lambda n: Fr(math.factorial(int(n))),
              Fraction=Fr, __name__='survey_exact')
    def rf(a, n):
        r = Fr(1)
        for i in range(int(n)):
            r *= (Fr(a) + i)
        return r
    ns['rf'] = rf
    exec(compile(src, path, 'exec'), ns)
    return ns


def ratio_signature(vals, k0):
    """vals: exact Fractions t(k0), t(k0+1), ... ; returns (z, classes frozenset, junk) or None"""
    rs = []
    for i in range(len(vals) - 1):
        if vals[i] == 0 or vals[i + 1] == 0:
            return None
        rs.append((k0 + i, vals[i + 1] / vals[i]))
    for d in range(1, (len(rs) - 4) // 2):
        pts = rs[:2 * d + 2]
        M = pari.matrix(len(pts), 2 * d + 2)
        for i, (kk, r) in enumerate(pts):
            for j in range(d + 1):
                M[i, j] = pari(kk) ** j
                M[i, d + 1 + j] = -(pari(r.numerator) / pari(r.denominator)) * pari(kk) ** j
        ker = pari.matker(M)
        if int(pari.matsize(ker)[1]) != 1:
            continue
        k = pari('k')
        N = sum((ker[j, 0] * k ** j for j in range(d + 1)), pari(0))
        D = sum((ker[d + 1 + j, 0] * k ** j for j in range(d + 1)), pari(0))
        ok = True
        for kk, r in rs[2 * d + 2:]:
            if pari.subst(N, 'k', kk) * r.denominator != pari.subst(D, 'k', kk) * r.numerator:
                ok = False
                break
        if not ok:
            continue
        rat = N / D
        N, D = pari.numerator(rat), pari.denominator(rat)
        dn, dd = pari.poldegree(N, 'k'), pari.poldegree(D, 'k')
        if dn != dd:
            return ('nonbalanced', dn, dd)
        z = _to_fr(pari.pollead(N, 'k') / pari.pollead(D, 'k'))
        cnt = Counter()
        junk = 0
        for poly, sgn in ((N, 1), (D, -1)):
            if pari.poldegree(poly, 'k') <= 0:
                continue
            fa = pari.factor(poly)
            for i in range(int(pari.matsize(fa)[0])):
                p, e = fa[i, 0], int(fa[i, 1])
                dg = pari.poldegree(p, 'k')
                if dg == 1:
                    alpha = _to_fr(pari.polcoef(p, 0, 'k') / pari.polcoef(p, 1, 'k'))
                    cnt[frac1(alpha)] += sgn * e
                elif dg > 1:
                    junk += 1
        return z, frozenset((a, m) for a, m in cnt.items() if m != 0), junk, rat
    return None


def entry_signature(ns, e, npts=90, hconst=Fr(7, 3)):
    # constant harmonic numbers
    class Hc:
        def __call__(self, n, m=1):
            return hconst + Fr(m, 11)
        def reset(self):
            pass
    saved = {}
    for name in list(ns):
        if name in ('H', 'Hn', 'harm', 'Hk', 'har') or (name.startswith('H') and callable(ns[name]) and len(name) <= 3):
            saved[name] = ns[name]
            ns[name] = Hc()
    try:
        k0 = e.get('start', 0) + 3
        vals = [Fr(e['summand'](k)) for k in range(k0, k0 + npts)]
        return ratio_signature(vals, k0)
    except NotRational:
        return 'not-rational'
    except Exception as ex:
        return 'error: ' + repr(ex)[:60]
    finally:
        ns.update(saved)


def fmt_sig(sig):
    if not isinstance(sig, tuple) or len(sig) < 3 or not isinstance(sig[1], frozenset):
        return str(sig)
    z, cl, junk = sig[:3]
    up = sorted([(a, m) for a, m in cl if m > 0])
    lo = sorted([(a, -m) for a, m in cl if m < 0])
    f = lambda L: " ".join("(%s)^%d" % (a, m) if m > 1 else "(%s)" % a for a, m in L)
    return "z=%s  [%s] / [%s]" % (z, f(up), f(lo))


if __name__ == '__main__':
    ns = load_exact(sys.argv[1])
    ents = ns.get('ENTRIES') or ns.get('E')
    for e in ents:
        print(e['id'], fmt_sig(entry_signature(ns, e)), '|', e.get('src', ''))
