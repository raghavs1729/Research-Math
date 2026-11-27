"""
wz.py -- a small Wilf-Zeilberger engine (Gosper's algorithm for structured
hypergeometric terms), built on PARI (cypari2) for polynomial algebra.

A term is
    F = pref(n,k,params) * prod_i base_i^(L_i) * prod_j Gamma(M_j)^(p_j)
with L_i, M_j linear forms in k, n and parameters with rational coefficients,
integer coefficients of n and k (so shifts n->n+1, k->k+1 give rational ratios).

Main entry points
    Term(...)                    build a term
    term.subs({...})             affine substitution of variables
    wz_mate(F)                   R(n,k) with G = R F,  F(n+1,k)-F(n,k) = G(n,k+1)-G(n,k)
    check_wz(F, R)               exact rational-function verification of the WZ equation
"""
from fractions import Fraction as Fr
from collections import Counter
import cypari2

pari = cypari2.Pari()
pari.allocatemem(4 * 10**9)
pari.default("debugmem", 0)

# variable priority: k highest, then n, then parameters
_VARS = ['k', 'n', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'u', 'v', 'w', 'x', 'y', 'z']
for _v in _VARS:
    pari(_v)
VARORDER = {v: i for i, v in enumerate(_VARS)}


def P(s):
    return pari(str(s))


# ---------------------------------------------------------------- linear forms
def lf(*args, **kw):
    """lf(const, k=1, n=2, a=-1) -> linear form dict"""
    d = {}
    c = Fr(args[0]) if args else Fr(0)
    if c != 0:
        d['1'] = c
    for v, x in kw.items():
        x = Fr(x)
        if x != 0:
            d[v] = x
    return d


def _to_fr(x):
    x = pari.simplify(x)
    assert x.type() in ('t_INT', 't_FRAC'), x
    return Fr(int(pari.numerator(x)), int(pari.denominator(x)))


def lf_parse(s):
    """parse a linear form from a string like 'a+2*k-1/2' using PARI"""
    p = pari(s)
    d = {}
    rest = p
    for v in _VARS:
        cf = pari.simplify(pari.polcoef(rest, 1, v))
        if cf != 0:
            d[v] = _to_fr(cf)
            rest = pari.simplify(rest - cf * pari(v))
    c = _to_fr(rest)
    if c != 0:
        d['1'] = c
    return d


def lf_add(L, M, s=1):
    d = dict(L)
    for v, x in M.items():
        d[v] = d.get(v, Fr(0)) + s * x
        if d[v] == 0:
            del d[v]
    return d


def lf_scale(L, s):
    s = Fr(s)
    return {v: x * s for v, x in L.items()} if s != 0 else {}


def lf_subs(L, m):
    """m: var -> linear form"""
    d = {}
    for v, x in L.items():
        if v != '1' and v in m:
            d = lf_add(d, lf_scale(m[v], x))
        else:
            d = lf_add(d, {v: x})
    return d


def lf_key(L):
    return tuple(sorted(L.items(), key=lambda t: (VARORDER.get(t[0], 99), t[0])))


def lf_pari(L):
    r = pari(0)
    for v, x in L.items():
        q = pari(x.numerator) / pari(x.denominator)
        r = r + (q if v == '1' else q * pari(v))
    return r


def lf_str(L):
    return str(lf_pari(L))


def lf_monic(L):
    """return (lead coeff, monic form) w.r.t. variable order k, n, params; const forms -> (c, {})"""
    vs = [v for v in L if v != '1']
    if not vs:
        return L.get('1', Fr(0)), {}
    v0 = min(vs, key=lambda v: VARORDER.get(v, 99))
    c = L[v0]
    return c, lf_scale(L, 1 / c)


# ---------------------------------------------------------------- factored rationals
class Fac:
    """scalar * prod (monic linear form)^mult ; scalar is a python Fraction"""

    def __init__(self, scalar=Fr(1), facs=None):
        self.scalar = Fr(scalar)
        self.facs = Counter(facs or {})

    def mul_lin(self, L, mult=1):
        c, M = lf_monic(L)
        if c == 0:
            raise ZeroDivisionError("zero linear factor")
        self.scalar *= c ** mult
        if M:
            key = lf_key(M)
            self.facs[key] += mult
            if self.facs[key] == 0:
                del self.facs[key]

    def num_den(self):
        num = pari(1)
        den = pari(1)
        for key, m in self.facs.items():
            p = lf_pari(dict(key))
            if m > 0:
                num = num * p ** m
            else:
                den = den * p ** (-m)
        return num, den

    def to_pari(self):
        num, den = self.num_den()
        return pari(self.scalar.numerator) / pari(self.scalar.denominator) * num / den


def rising(fac, L, shift, power):
    """multiply fac by (Gamma(L+shift)/Gamma(L))^power, shift integer"""
    if shift > 0:
        for j in range(shift):
            fac.mul_lin(lf_add(L, {'1': Fr(j)}), power)
    elif shift < 0:
        for j in range(1, -shift + 1):
            fac.mul_lin(lf_add(L, {'1': Fr(-j)}), -power)


# ---------------------------------------------------------------- terms
class Term:
    def __init__(self, pref='1', bases=(), gammas=()):
        """bases: list of (base, linear form); gammas: list of (linear form, power);
        linear forms may be given as strings"""
        self.pref = pari(str(pref))
        self.bases = [(Fr(b), lf_parse(L) if isinstance(L, str) else dict(L)) for b, L in bases]
        self.gammas = [(lf_parse(L) if isinstance(L, str) else dict(L), int(p)) for L, p in gammas]

    def subs(self, m):
        """m: dict var -> linear form (string or dict). Simultaneous substitution."""
        m = {v: (lf_parse(L) if isinstance(L, str) else dict(L)) for v, L in m.items()}
        t = Term.__new__(Term)
        pr = self.pref
        # simultaneous substitution in pref through temporaries
        tmp = {}
        for i, v in enumerate(m):
            tv = 'tmpv%d' % i
            tmp[v] = tv
            pr = pari.subst(pr, v, pari(tv))
        for v, tv in tmp.items():
            pr = pari.subst(pr, tv, lf_pari(m[v]))
        t.pref = pr
        t.bases = [(b, lf_subs(L, m)) for b, L in self.bases]
        t.gammas = [(lf_subs(L, m), p) for L, p in self.gammas]
        return t

    def mul(self, other):
        t = Term.__new__(Term)
        t.pref = self.pref * other.pref
        t.bases = self.bases + other.bases
        t.gammas = self.gammas + other.gammas
        return t

    def shift_ratio(self, v):
        """F(v+1)/F(v) WITHOUT the prefactor: returns (Fac) ; base factors folded in scalar"""
        fac = Fac()
        for b, L in self.bases:
            e = L.get(v, Fr(0))
            assert e.denominator == 1, "non-integer base exponent shift"
            fac.scalar *= b ** int(e)
        for L, p in self.gammas:
            s = L.get(v, Fr(0))
            assert s.denominator == 1, ("non-integer gamma shift", lf_str(L), v)
            rising(fac, L, int(s), p)
        return fac

    def pref_ratio(self, v):
        return pari.subst(self.pref, v, pari(v) + 1) / self.pref

    def ratio(self, v):
        return self.shift_ratio(v).to_pari() * self.pref_ratio(v)

    # numeric evaluation with mpmath
    def evalf(self, vals, mp):
        """vals: dict var->mpf/mpc"""
        r = mp.mpf(1)
        pr = self.pref
        num = pari.numerator(pr)
        den = pari.denominator(pr)

        def ev(p):
            # evaluate a pari polynomial numerically via string -> python eval
            s = str(p).replace('^', '**')
            env = {v: vals[v] for v in vals}
            import re
            s = re.sub(r'(?<![\w.])(\d+)(?![\w.])', r'mpf(\1)', s)
            env['mpf'] = mp.mpf
            return eval(s, {'__builtins__': {}}, env)
        r = ev(num) / ev(den)
        for b, L in self.bases:
            e = sum((mp.mpf(x.numerator) / x.denominator) * (1 if v == '1' else vals[v]) for v, x in L.items())
            r *= mp.power(mp.mpf(b.numerator) / b.denominator, e)
        for L, p in self.gammas:
            x = sum((mp.mpf(c.numerator) / c.denominator) * (1 if v == '1' else vals[v]) for v, c in L.items())
            if p > 0:
                r *= mp.gamma(x) ** p
            else:
                r *= mp.rgamma(x) ** (-p)
        return r

    def __str__(self):
        s = "(%s)" % self.pref
        for b, L in self.bases:
            s += " * (%s)^(%s)" % (b, lf_str(L))
        for L, p in self.gammas:
            s += " * Gamma(%s)^%d" % (lf_str(L), p)
        return s


# ---------------------------------------------------------------- Gosper / WZ
class NoMate(Exception):
    pass


def _poly_coeffs(p, var, deg):
    return [pari.polcoef(p, i, var) for i in range(deg + 1)]


def gosper_multi(sigma0, Nlist, rho_den_fac, var='k', verbose=False, extra_deg=0):
    """
    Zeilberger-type step.  t_j(k) = T0(k) * N_j(k) / rho_den(k)  (j = 0..m-1), where
    T0(k+1)/T0(k) = sigma0 (Fac of linear factors), N_j PARI polynomials in k.
    Finds all (beta_0..beta_{m-1}; S) with  s(k+1)-s(k) = sum_j beta_j t_j(k),  s = S*T0,
    beta_j free of k.  Returns list of (betas, S) spanning the solution space.
    """
    k = pari(var)
    AB = Fac(sigma0.scalar, sigma0.facs)
    for key, m in rho_den_fac.facs.items():
        L = dict(key)
        AB.mul_lin(L, m)
        AB.mul_lin(lf_add(L, lf_scale({'1': Fr(1)}, L.get(var, 0))), -m)
    scal = pari(AB.scalar.numerator) / pari(AB.scalar.denominator)
    num = []
    den = []
    for key, m in AB.facs.items():
        L = dict(key)
        if var not in L:
            scal = scal * lf_pari(L) ** m
            continue
        assert L[var] == 1
        alpha = dict(L)
        del alpha[var]
        (num if m > 0 else den).append([alpha, abs(m)])
    cmul = pari(1)
    changed = True
    while changed:
        changed = False
        for ua in num:
            if ua[1] == 0:
                continue
            for vb in den:
                if vb[1] == 0 or ua[1] == 0:
                    continue
                diff = lf_add(ua[0], vb[0], -1)
                if all(v == '1' for v in diff) and diff.get('1', Fr(0)).denominator == 1:
                    h = int(diff.get('1', Fr(0)))
                    if h > 0:
                        mlt = min(ua[1], vb[1])
                        ua[1] -= mlt
                        vb[1] -= mlt
                        g = k + lf_pari(ua[0])
                        for jj in range(1, h + 1):
                            cmul = cmul * pari.subst(g, var, k - jj) ** mlt
                        changed = True
    a = scal
    for alpha, m in num:
        if m:
            a = a * (k + lf_pari(alpha)) ** m
    b = pari(1)
    for beta, m in den:
        if m:
            b = b * (k + lf_pari(beta)) ** m
    da_den = pari.denominator(pari.content(a))
    a = a * da_den
    b = b * da_den
    bm = pari.subst(b, var, k - 1)
    clist = [Nj * cmul for Nj in Nlist]
    da, db = pari.poldegree(a, var), pari.poldegree(bm, var)
    dc = max(pari.poldegree(c, var) for c in clist)
    la, lb = pari.pollead(a, var), pari.pollead(bm, var)
    cands = []
    if da != db or pari.simplify(la - lb) != 0:
        cands.append(dc - max(da, db))
    else:
        cands.append(dc - da + 1)
        A1 = pari.polcoef(a, da - 1, var) if da >= 1 else pari(0)
        B1 = pari.polcoef(bm, da - 1, var) if da >= 1 else pari(0)
        q = pari.simplify((B1 - A1) / la)
        if q.type() == 't_INT' and int(q) >= 0:
            cands.append(int(q))
    d = max(cands) + extra_deg
    if verbose:
        print("  gosper: deg a,b,c = %d,%d,%d  deg x = %d" % (da, db, dc, d))
    if d < 0:
        raise NoMate("negative degree bound")
    m = len(clist)
    nrows = max(max(da, db) + d + 1, dc + 1)
    M = pari.matrix(nrows, d + 1 + m)
    for j in range(d + 1):
        pj = a * (k + 1) ** j - bm * k ** j
        for i in range(nrows):
            M[i, j] = pari.polcoef(pj, i, var)
    for j in range(m):
        for i in range(nrows):
            M[i, d + 1 + j] = -pari.polcoef(clist[j], i, var)
    ker = pari.matker(M)
    ncol = int(pari.matsize(ker)[1])
    rd_num, rd_den = rho_den_fac.num_den()
    rho_den = rd_num * pari(rho_den_fac.scalar.numerator) / pari(rho_den_fac.scalar.denominator)
    sols = []
    for jc in range(ncol):
        betas = [ker[d + 1 + j, jc] for j in range(m)]
        if all(bb == 0 for bb in betas):
            continue
        x = sum((ker[i, jc] * k ** i for i in range(d + 1)), pari(0))
        # s = b(k-1) x / c * t,  c = cmul * sum beta N ;  t = T0 * sum beta N / rho_den
        S = bm * x / cmul / rho_den
        sols.append((betas, S))
    if not sols:
        raise NoMate("no solution")
    return sols


def gosper_structured(sigma0, Npoly, rho_den_fac, var='k', verbose=False, extra_deg=0):
    sols = gosper_multi(sigma0, [Npoly], rho_den_fac, var, verbose, extra_deg)
    betas, S = sols[0]
    return S / betas[0]


def wz_mate(F, verbose=False, extra_deg=0):
    """F a Term in n,k (+parameters). Returns R with G=R*F, F(n+1,k)-F(n,k)=G(n,k+1)-G(n,k)."""
    n, k = pari('n'), pari('k')
    rho0 = F.shift_ratio('n')            # gamma/base part of F(n+1,k)/F(n,k)
    pref = F.pref
    pref1 = pari.subst(pref, 'n', n + 1)
    # rho - 1 = (pref1*num - pref*den)/(pref*den)
    num, den = rho0.num_den()
    sc = pari(rho0.scalar.numerator) / pari(rho0.scalar.denominator)
    pn, pd = pari.numerator(pref), pari.denominator(pref)
    p1n, p1d = pari.numerator(pref1), pari.denominator(pref1)
    # rho = sc*num/den * (p1n/p1d)/(pn/pd);  t = F0 * pref*(rho-1) = F0 * (sc*num*p1n*pd*... )
    # t(k) = F0 * [ sc*num*pref1 - den*pref ] / den  =  F0 * N / (den * pd * p1d),
    N = sc * num * p1n * pd - den * pn * p1d
    assert pari.poldegree(pd, 'k') <= 0 and pari.poldegree(p1d, 'k') <= 0, "prefactor denominators must be k-free"
    N = N / (pd * p1d)
    # make N a polynomial in k with coefficients polynomial in n (clear k-free denominators; harmless scalar)
    cont_den = pari.denominator(pari.content(N))
    Nk = N * cont_den
    rho_den = Fac()
    for key, m in rho0.facs.items():
        if m < 0:
            rho_den.facs[key] = -m
    sigma0 = F.shift_ratio('k')
    S = gosper_structured(sigma0, Nk, rho_den, 'k', verbose=verbose, extra_deg=extra_deg)
    # s(k) = S * F0 with t = F0*Nk/rho_den ; true t has extra 1/cont_den
    R = S / cont_den / pref
    return R


def check_wz(F, R):
    """exact check of F(n+1,k)-F(n,k) = G(n,k+1)-G(n,k), G=R F, as rational functions:
       rho - 1 == R(n,k+1)*sigma - R(n,k)"""
    n, k = pari('n'), pari('k')
    rho = F.ratio('n')
    sigma = F.ratio('k')
    lhs = rho - 1
    rhs = pari.subst(R, 'k', k + 1) * sigma - R
    return lhs - rhs == 0


# ---------------------------------------------------------------- seeds (Au, arXiv:2312.14051 Sec 3.3)
def G_(s, p=1):
    return (s, p)


SEEDS = {}

SEEDS['Gauss2F1'] = (('a', 'b', 'c'), Term('1', [], [
    ('c-a', 1), ('a+k', 1), ('c-b', 1), ('b+k', 1),
    ('a', -1), ('b', -1), ('k+1', -1), ('c+k', -1), ('c-a-b', -1)]))

SEEDS['Dixon3F2'] = (('a', 'b', 'c'), Term('1', [], [
    ('a-b+1', 1), ('a-c+1', 1), ('2*a+k', 1), ('b+k', 1), ('c+k', 1), ('2*a-b-c+1', 1),
    ('a', -1), ('b', -1), ('c', -1), ('k+1', -1), ('a-b-c+1', -1), ('2*a-b+k+1', -1), ('2*a-c+k+1', -1)]))

SEEDS['Watson3F2'] = (('a', 'b', 'c'), Term('1', [(2, '-2*a-2*b+2*c')], [
    ('c-a+1/2', 1), ('2*a+k', 1), ('c-b+1/2', 1), ('2*b+k', 1), ('c+k', 1),
    ('a', -1), ('b', -1), ('k+1', -1), ('2*c+k', -1), ('c-a-b+1/2', -1), ('a+b+k+1/2', -1)]))

SEEDS['Dougall5F4'] = (('a', 'b', 'c', 'd'), Term('a+2*k', [], [
    ('a+k', 1), ('b+k', 1), ('c+k', 1), ('d+k', 1), ('a-b-c+1', 1), ('a-b-d+1', 1), ('a-c-d+1', 1),
    ('b', -1), ('c', -1), ('d', -1), ('k+1', -1), ('a-b+k+1', -1), ('a-c+k+1', -1), ('a-d+k+1', -1),
    ('a-b-c-d+1', -1)]))

SEEDS['Dougall7F6'] = (('a', 'b', 'c', 'd', 'e'), Term('a+2*k', [(-1, 'a+e')], [
    ('a+k', 1), ('b+k', 1), ('c+k', 1), ('d+k', 1), ('e+k', 1), ('a-b-c+1', 1), ('a-b-d+1', 1),
    ('a-c-d+1', 1), ('-a+b+c+e', 1), ('-a+b+d+e', 1), ('-a+c+d+e', 1), ('2*a-b-c-d-e+k+1', 1),
    ('b', -1), ('c', -1), ('d', -1), ('e', -1), ('k+1', -1), ('-a+b+e', -1), ('a-b+k+1', -1),
    ('-a+c+e', -1), ('a-c+k+1', -1), ('-a+d+e', -1), ('a-d+k+1', -1), ('a-e+k+1', -1),
    ('a-b-c-d+1', -1), ('2*a-b-c-d-e+1', -1), ('-a+b+c+d+e+k', -1)]))

SEEDS['Balanced3F2'] = (('a', 'b', 'c', 'd'), Term('1', [(-1, 'a+b+c-d')], [
    ('d-a', 1), ('a+k', 1), ('d-b', 1), ('b+k', 1), ('d-c', 1), ('c+k', 1),
    ('a', -1), ('b', -1), ('c', -1), ('k+1', -1), ('d+k', -1), ('d-a-b', -1), ('d-a-c', -1),
    ('d-b-c', -1), ('a+b+c-d+k+1', -1)]))

SEEDS['Seed1'] = (('a', 'c'), Term('1', [(2, 'a-2*c-2*k')], [
    ('-a+2*c-1', 1), ('a+2*k', 1),
    ('a', -1), ('k+1', -1), ('-a+c-1/2', -1), ('c+k', -1)]))

SEEDS['Seed2'] = (('a', 'c', 'd'), Term('1', [(2, '-2*c-2*k'), (-1, 'a+c-d')], [
    ('-a+2*d-1', 1), ('a+2*k', 1), ('d-c', 1), ('c+k', 1),
    ('a', -1), ('c', -1), ('k+1', -1), ('-a+d-1/2', -1), ('d+k', -1), ('-a-2*c+2*d-1', -1),
    ('a+c-d+k+3/2', -1)]))

SEEDS['Seed3'] = (('a', 'b', 'd'), Term('a+2*k', [(2, '2*d')], [
    ('a-b+1/2', 1), ('a+k', 1), ('b+2*k', 1), ('d+k', 1), ('2*a-b-2*d+1', 1),
    ('b', -1), ('d', -1), ('k+1', -1), ('a-b-d+1/2', -1), ('2*a-b+2*k+1', -1), ('a-d+k+1', -1)]))

SEEDS['Seed4'] = (('a', 'b'), Term('a+2*k', [(3, 'b+1')], [
    ('3*a-2*b', 1), ('a+k', 1), ('b+3*k', 1),
    ('b', -1), ('k+1', -1), ('a-b', -1), ('3*a-b+3*k+1', -1)]))

SEEDS['Seed7'] = (('a', 'b', 'd', 'e'), Term('a+2*k', [(-1, 'a+e')], [
    ('a-b+1/2', 1), ('a+k', 1), ('b+2*k', 1), ('d+k', 1), ('e+k', 1), ('2*a-b-2*d+1', 1),
    ('-a+b+e+1/2', 1), ('-2*a+b+2*d+2*e', 1), ('2*a-b-d-e+k+1/2', 1),
    ('b', -1), ('d', -1), ('e', -1), ('k+1', -1), ('a-b-d+1/2', -1), ('-2*a+b+2*e', -1),
    ('2*a-b+2*k+1', -1), ('-a+d+e', -1), ('a-d+k+1', -1), ('a-e+k+1', -1),
    ('2*a-b-d-e+1/2', -1), ('-a+b+d+e+k+1/2', -1)]))

SEEDS['Seed9'] = (('a', 'b', 'd'), Term('a+2*k', [(-1, 'd'), (2, '4*a-2*b-2*d')], [
    ('a-b+1/2', 1), ('a+k', 1), ('b+2*k', 1), ('d+2*k', 1), ('-2*a+2*b+d+1', 1),
    ('-2*a+b+2*d+1', 1), ('2*a-b-d+k', 1),
    ('b', -1), ('d', -1), ('k+1', -1), ('-a+d+1/2', -1), ('-2*a+b+d', -1), ('2*a-b+2*k+1', -1),
    ('2*a-d+2*k+1', -1), ('-a+b+d+k+1', -1)]))

SEEDS['Seed10'] = (('a', 'd'), Term('1', [(-1, 'a-d'), (3, '-a-3*k')], [
    ('-a+3*d-2', 1), ('a+3*k', 1),
    ('a', -1), ('k+1', -1), ('-2*a+3*d-3', -1), ('d+k', -1), ('a-d+k+2', -1)]))


def seed_instance(name, shifts, K, consts=None, k0=0):
    """F(n,k) = Seed(a + A n, ..., k + k0 + K n); shifts: dict param->A ; consts: dict param->value
    (Fraction/str) to specialise the parameter base point (otherwise kept symbolic)."""
    params, T = SEEDS[name]
    m = {}
    for p in params:
        base = {p: Fr(1)} if (consts is None or p not in consts) else (
            lf_parse(str(consts[p])))
        m[p] = lf_add(base, {'n': Fr(shifts.get(p, 0))})
    m['k'] = lf_add({'k': Fr(1), 'n': Fr(K)}, lf_parse(str(k0)))
    return T.subs(m)


# ---------------------------------------------------------------- Zeilberger (order 1) and automatic seeds
def zeil1(F, verbose=False):
    """F Term in n,k.  Find rho(n) (rational) and R(n,k) with
         F(n+1,k) - rho(n) F(n,k) = G(n,k+1) - G(n,k),   G = R F.
       Returns (rho, R)."""
    n, k = pari('n'), pari('k')
    rho0 = F.shift_ratio('n')
    pref = F.pref
    pref1 = pari.subst(pref, 'n', n + 1)
    num, den = rho0.num_den()
    sc = pari(rho0.scalar.numerator) / pari(rho0.scalar.denominator)
    pn, pd = pari.numerator(pref), pari.denominator(pref)
    p1n, p1d = pari.numerator(pref1), pari.denominator(pref1)
    assert pari.poldegree(pd, 'k') <= 0 and pari.poldegree(p1d, 'k') <= 0
    # F(n+1,k) = F0 * pref1*sc*num/den ; F(n,k) = F0*pref  ->  t_j = F0 * N_j/den
    N1 = sc * num * p1n / p1d
    N2 = den * pn / pd
    rho_den = Fac()
    for key, m in rho0.facs.items():
        if m < 0:
            rho_den.facs[key] = -m
    sigma0 = F.shift_ratio('k')
    sols = gosper_multi(sigma0, [N1, N2], rho_den, 'k', verbose=verbose)
    for betas, S in sols:
        if betas[0] != 0:
            rho = pari.simplify(-betas[1] / betas[0])
            R = S / betas[0] / pref
            return rho, R
    raise NoMate("no first-order recurrence")


def check_zeil1(F, rho, R):
    k = pari('k')
    return F.ratio('n') - rho - (pari.subst(R, 'k', k + 1) * F.ratio('k') - R) == 0


def gamma_factors_of(rho, var='n'):
    """rho(var) rational over Q, product of linear factors: returns (const, [(alpha, mult)]) with
       rho = const * prod (var + alpha)^mult ; raises if a non-linear factor occurs"""
    v = pari(var)
    num, den = pari.numerator(rho), pari.denominator(rho)
    const = pari.pollead(num, var) / pari.pollead(den, var)
    out = []
    for poly, sgn in ((num, 1), (den, -1)):
        if pari.poldegree(poly, var) <= 0:
            continue
        fa = pari.factor(poly)
        for i in range(int(pari.matsize(fa)[0])):
            p, e = fa[i, 0], int(fa[i, 1])
            if pari.poldegree(p, var) == 0:
                continue
            if pari.poldegree(p, var) != 1:
                raise ValueError("non-linear factor in rho: %s" % p)
            alpha = pari.polcoef(p, 0, var) / pari.polcoef(p, 1, var)
            out.append((_to_fr(alpha), sgn * e))
    return _to_fr(const), out


def auto_seed_2F1(A, B, C, x, param='a', extra_upper=(), extra_lower=()):
    """Given a one-parameter family  pFq(A(a),B(a),..;C(a),..;x) (strings, linear in a with rational
    coefficients) with a first-order recurrence in a, return (step, Term) where Term is the seed
    f(a,k) = t(s*a,k)/S(a) in the rescaled parameter (a -> s*a makes all shifts integral)."""
    forms_up = [lf_parse(e) for e in (A, B) + tuple(extra_upper)]
    forms_lo = [lf_parse(e) for e in (C,) + tuple(extra_lower)]
    from math import lcm
    s = 1
    for L in forms_up + forms_lo:
        s = lcm(s, L.get(param, Fr(0)).denominator)
    resc = {param: {param: Fr(s)}}
    forms_up = [lf_subs(L, resc) for L in forms_up]
    forms_lo = [lf_subs(L, resc) for L in forms_lo]
    xf = Fr(x)
    gam = []
    for L in forms_up:
        gam.append((lf_add(L, {'k': Fr(1)}), 1))
        gam.append((L, -1))
    for L in forms_lo:
        gam.append((lf_add(L, {'k': Fr(1)}), -1))
        gam.append((L, 1))
    gam.append(({'k': Fr(1), '1': Fr(1)}, -1))
    T = Term('1', [(xf, {'k': Fr(1)})], gam)
    F = T.subs({param: {param: Fr(1), 'n': Fr(1)}})
    rho, R = zeil1(F)
    assert check_zeil1(F, rho, R)
    # rho is a function of (param + n): substitute n -> n - param is implicit; factor in n at param=0
    rho_n = pari.subst(rho, param, 0)
    const, facs = gamma_factors_of(rho_n, 'n')
    # S(a+1) = rho(a) S(a):  S(a) = const^a prod Gamma(a+alpha)^mult
    bases = [(xf, {'k': Fr(1)})]
    neg = False
    if const < 0:
        bases.append((Fr(-1), {param: Fr(-1)}))
        const = -const
    if const != 1:
        bases.append((const, {param: Fr(-1)}))
    gam2 = list(gam)
    for alpha, mlt in facs:
        gam2.append(({param: Fr(1), '1': alpha}, -mlt))
    return s, Term('1', bases, gam2)
