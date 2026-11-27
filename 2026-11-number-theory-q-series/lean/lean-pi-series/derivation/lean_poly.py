"""Helper: print a sympy polynomial in n, k (or m, k with n = m+1) as a Lean expression over the reals
(used to state 'difference = polynomial with non-negative coefficients' certificates)."""
import sympy as sp
def lean(expr, names):
    """names: dict sympy symbol -> Lean term, e.g. {n: '(n : ℝ)'}"""
    expr = sp.expand(expr)
    syms = list(names.keys())
    p = sp.Poly(expr, *syms)
    terms = []
    for mon, c in p.terms():
        assert c.is_rational
        fac = []
        for s, e in zip(syms, mon):
            if e == 1: fac.append(names[s])
            elif e > 1: fac.append('%s ^ %d' % (names[s], e))
        cs = str(c) if c.q == 1 else '(%d / %d)' % (c.p, c.q)
        if not fac: terms.append(cs)
        elif c == 1: terms.append(' * '.join(fac))
        else: terms.append(cs + ' * ' + ' * '.join(fac))
    return ' + '.join(terms)
