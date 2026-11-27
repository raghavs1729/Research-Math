#!/usr/bin/env python3
"""
Independent check of the WZ certificates in ../certificates/*.gp  (audit2, referee script).

Nothing from ../proofs or ../code is imported.  The .gp files are parsed as plain text
(assignments  name = expr;  only -- the `print(...)` lines of the authors are ignored) and the
identities are re-stated here:

  generic        rho - 1 == R(n,k+1)*sg - R(n,k)                        (T1,T2,T3,T5,T6,T7,T11)
  T4, T8         (n+1)*rho - n == Q(n,k+1)*sg - Q(n,k)
  compatibility  rho(n,k)*sg(n+1,k) == sg(n,k)*rho(n,k+1)               (all)
  triangle row   rho(n,n-1) - 1 == -rho(n,n-1)*sg(n+1,n-1) - R(n,n-1)   (T7, T11)

Two independent evaluation paths:
  (1) exact rational arithmetic (fractions.Fraction) at 300 random rational points per identity
      (a non-zero rational function of the degrees occurring here cannot vanish at so many
       independent random points with 40-bit numerators: Schwartz-Zippel);
  (2) symbolic: python-flint fmpq_mpoly, rational functions kept as (numerator, denominator)
      pairs WITHOUT any cancellation, final numerator compared with the zero polynomial.
"""
import re, sys, os, random
from fractions import Fraction
import flint

CERT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'certificates')

# ---------------------------------------------------------------- parsing
def parse_gp(path):
    txt = open(path).read()
    txt = re.sub(r'\\\\[^\n]*', '', txt)            # GP comments
    defs, funs = {}, {}
    for stmt in txt.split(';'):
        s = stmt.strip()
        if not s or s.startswith('print'):
            continue
        m = re.match(r'^([A-Za-z_]\w*)\((\w+)\)\s*=\s*(.*)$', s, flags=re.S)
        if m:                                       # one-argument function, e.g. rf4(z)
            funs[m.group(1)] = (m.group(2), m.group(3))
            continue
        m = re.match(r'^([A-Za-z_]\w*)\s*=\s*(.*)$', s, flags=re.S)
        if not m:
            raise ValueError('cannot parse: ' + s[:60])
        if 'subst' in m.group(2):
            continue                                # authors' helper lines; we restate them ourselves
        defs[m.group(1)] = m.group(2)
    return defs, funs

def to_py(expr):
    e = expr.replace('\n', ' ')
    e = re.sub(r'\^\s*(\d+)', r'**@\1', e)          # protect integer exponents
    e = re.sub(r'(?<![\w@])(\d+)(?!\w)', r'Q(\1)', e)
    return e.replace('@', '')

class Env(dict):
    """names -> values; certificate names are evaluated lazily under the current binding"""
    def __init__(self, defs, funs, binding, Q):
        super().__init__(binding)
        self.defs, self.funs, self.Q = defs, funs, Q
    def __missing__(self, name):
        if name == 'Q':
            return self.Q
        if name in self.funs:
            arg, body = self.funs[name]
            code = to_py(body)
            return lambda z: eval(code, {'__builtins__': {}}, Env(self.defs, self.funs, {**dict(self), arg: z}, self.Q))
        if name in self.defs:
            return eval(to_py(self.defs[name]), {'__builtins__': {}}, self)
        raise KeyError(name)

def ev(cert, name, binding, Q):
    defs, funs = cert
    return Env(defs, funs, dict(binding), Q)[name]

# ---------------------------------------------------------------- rational functions over flint
class RF:
    __slots__ = ('n', 'd')
    def __init__(self, n, d):
        self.n, self.d = n, d
    @staticmethod
    def lift(x, ctx):
        if isinstance(x, RF):
            return x
        return RF(ctx.constant(flint.fmpq(Fraction(x).numerator, Fraction(x).denominator)) if not isinstance(x, flint.fmpq_mpoly) else x,
                  ctx.constant(1))
    def _c(self, o):
        return o if isinstance(o, RF) else RF.lift(o, self.n.context())
    def __add__(self, o): o = self._c(o); return RF(self.n*o.d + o.n*self.d, self.d*o.d)
    __radd__ = __add__
    def __neg__(self): return RF(-self.n, self.d)
    def __sub__(self, o): return self + (-self._c(o))
    def __rsub__(self, o): return self._c(o) - self
    def __mul__(self, o): o = self._c(o); return RF(self.n*o.n, self.d*o.d)
    __rmul__ = __mul__
    def __truediv__(self, o): o = self._c(o); return RF(self.n*o.d, self.d*o.n)
    def __rtruediv__(self, o): return self._c(o) / self
    def __pow__(self, e): return RF(self.n**int(e), self.d**int(e))
    def is_zero(self): return self.n == 0 and self.d != 0

# ---------------------------------------------------------------- the identities
#  (file, tag, second variable, extra parameters, kind)
JOBS = [
    ('T1.gp',     'T1',      'k', [],         'std'),
    ('T2_T9.gp',  'T2/T9',   'k', [],         'std'),
    ('T3.gp',     'T3',      'x', [],         'std'),
    ('T4_T10.gp', 'T4/T10',  'k', [],         'nh'),
    ('T5.gp',     'T5',      'k', [],         'std'),
    ('T6.gp',     'T6',      'k', ['a', 'b'], 'std'),
    ('T7.gp',     'T7',      'k', [],         'tri'),
    ('T8.gp',     'T8',      'k', ['s'],      'nh'),
    ('T11.gp',    'T11',     'k', [],         'tri'),
]

def identities(cert, kind, kv, B, Q):
    """returns dict name -> value that must be 0, at the binding B (values of n, kv, params)"""
    n, k = B['n'], B[kv]
    one = Q(1)
    def at(name, nn, kk):
        b = dict(B); b['n'] = nn; b[kv] = kk
        return ev(cert, name, b, Q)
    cname = 'Q' if kind == 'nh' else 'R'
    if kind == 'nh':            # the certificate is called Q in the file; avoid clash with literal wrapper
        cert = (dict(cert[0]), cert[1]); cert[0]['CERT'] = cert[0].pop('Q'); cname = 'CERT'
    sg, rho = at('sg', n, k), at('rho', n, k)
    out = {}
    if kind == 'nh':
        out['WZ'] = (n + one)*rho - n - (at(cname, n, k + one)*sg - at(cname, n, k))
    else:
        out['WZ'] = rho - one - (at(cname, n, k + one)*sg - at(cname, n, k))
    out['compat'] = rho*at('sg', n + one, k) - sg*at('rho', n, k + one)
    if kind == 'tri':           # boundary row k = n-1  (function of n only)
        r = at('rho', n, n - one)
        out['row k=n-1'] = r - one + r*at('sg', n + one, n - one) + at(cname, n, n - one)
    return out

def main():
    random.seed(20261124)
    ok_all = True
    for fn, tag, kv, params, kind in JOBS:
        cert = parse_gp(os.path.join(CERT, fn))
        names = ['n', kv] + params
        # (1) random rational points
        QF = lambda v: Fraction(v)
        bad = 0; npts = 300
        for _ in range(npts):
            B = {v: Fraction(random.randint(-2**40, 2**40), random.randint(1, 2**20)) for v in names}
            try:
                res = identities(cert, kind, kv, B, QF)
            except ZeroDivisionError:
                continue
            for nm, val in res.items():
                if val != 0:
                    bad += 1
        # (2) symbolic, no cancellation
        ctx = flint.fmpq_mpoly_ctx.get(names, 'lex')
        gens = dict(zip(names, ctx.gens()))
        QS = lambda v: RF.lift(Fraction(v), ctx)
        B = {v: RF(g, ctx.constant(1)) for v, g in gens.items()}
        res = identities(cert, kind, kv, B, QS)
        sym = {nm: val.is_zero() for nm, val in res.items()}
        degs = {nm: val.d.total_degree() for nm, val in res.items()}
        good = (bad == 0) and all(sym.values())
        ok_all &= good
        print(f'{tag:7s} random points: {npts} pts, {bad} failures | symbolic (flint, uncancelled): {sym} '
              f'(denominator degrees {degs}) -> {"OK" if good else "FAIL"}')
    print('ALL CERTIFICATES OK' if ok_all else 'SOME CERTIFICATE FAILED')
    # negative control: a certificate with one coefficient changed must be rejected by both paths
    neg_ok = True
    for fn, kind, old, new in [('T4_T10.gp', 'nh', '612', '613'), ('T11.gp', 'tri', '1310720', '1310721'),
                               ('T7.gp', 'tri', '10080', '10081'), ('T5.gp', 'std', '145', '146')]:
        defs, funs = parse_gp(os.path.join(CERT, fn))
        key = 'Q' if 'Q' in defs else ('NUM' if 'NUM' in defs else ('P' if 'P' in defs else 'R'))
        assert old in defs[key], (fn, key)
        defs[key] = defs[key].replace(old, new, 1)
        B = {v: Fraction(random.randint(1, 2**30), random.randint(1, 2**10)) for v in ('n', 'k')}
        r1 = identities((defs, funs), kind, 'k', B, lambda v: Fraction(v))['WZ'] != 0
        ctx = flint.fmpq_mpoly_ctx.get(['n', 'k'], 'lex'); g = dict(zip(['n', 'k'], ctx.gens()))
        r2 = not identities((defs, funs), kind, 'k', {v: RF(x, ctx.constant(1)) for v, x in g.items()},
                            lambda v: RF.lift(Fraction(v), ctx))['WZ'].is_zero()
        print(f'negative control {fn} ({old}->{new}): rejected by random points: {r1}, by flint: {r2}')
        neg_ok &= r1 and r2
    print('NEGATIVE CONTROLS OK' if neg_ok else 'NEGATIVE CONTROL FAILED (checker is vacuous!)')
    return 0 if (ok_all and neg_ok) else 1

if __name__ == '__main__':
    sys.exit(main())
