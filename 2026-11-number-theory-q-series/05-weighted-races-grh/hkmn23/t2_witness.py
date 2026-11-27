#!/usr/bin/env python3
"""
t2_witness.py -- prime-side witnesses for HKMN23.md, Theorems 3 and 6 (GRH-only exclusion).

Under GRH(q), for x >= 1 (y = log x >= 0), EXACTLY (HKMN23.md, Lemma E / eq. (F)):

  F_{q,a}(log x) := - sum_{chi mod q} conj(chi(a)) sum_{rho of chi*} x^{i gamma} / (rho(rho+1))
                  = x^{-3/2} * [ phi(q) psi_1(x;q,a) - x^2/2 - ell_{q,a}(x) + imp_{q,a}(x) ],

  psi_1(x;q,a)  = sum_{n <= x, n = a (q)} Lambda(n) (x - n),
  ell_{q,a}(x)  = -n_e x log x + B x + n_o log x + C - (1 + n_e) T_e(x) - n_o T_o(x),
  T_e(x) = sum_{r>=1} x^{1-2r}/(2r(2r-1)) = (x/2)[(1+u)log(1+u) + (1-u)log(1-u)],   u = 1/x   (T_e(1) = log 2),
  T_o(x) = sum_{n>=1} x^{-2n}/(2n(2n+1))  = 1 - (1/2) log(1-u^2) - atanh(u)/u              (T_o(1) = 1 - log 2),
  imp_{q,a}(x)  = sum_{p | q} phi(m_p) sum_{k>=1, p^k = a (mod m_p)} log p * (x - p^k)_+ ,  m_p = q / p^{v_p(q)},
  n_e, n_o, B, C as in w1_analytic.py (constants of the primitive characters, computed here in Arb).

No zero of any L-function is used.  If F_{q,a}(log x_0) > (2/3) N(q,a) at ONE x_0 >= 1, then under GRH(q)
alone A_1^pi(x;q,a) > 0 for arbitrarily large x (Theorem 3).

For each pair (q, a), q in the given range, the script
  1. scans x in [1, X] in floating point (x = 1, a geometric grid, and points inside the gaps between
     prime powers = a mod q), looking for max F and for the first x with F > threshold;
  2. certifies F_{q,a}(log x_0) in Arb (ball arithmetic) at the chosen rational x_0.
Output: one JSON line per pair.

Usage: python3 t2_witness.py QMIN QMAX X [shard nshards]
"""
import sys, json, math
import numpy as np
from fractions import Fraction
from flint import acb, arb, acb_series, dirichlet_char, dirichlet_group, ctx, fmpq

ctx.prec = 128
LWS18 = {(1, 1), (2, 1), (3, 1), (4, 1), (5, 1), (5, 4), (6, 1), (8, 1), (10, 1), (10, 9), (12, 1), (15, 1), (15, 4),
         (20, 1), (20, 9), (24, 1), (30, 1), (30, 19)}

_cache = {}


def ljet(chi, s0, n):
    ser = acb_series([acb(s0), acb(1)], prec=n).dirichlet_l(chi)
    return [ser[i] for i in range(n)]


def prim_consts(f, n):
    """(parity, b, second) for the primitive character f.n ; second = b0 (odd) or cm1 (even)."""
    key = (f, n)
    if key in _cache:
        return _cache[key]
    chi = dirichlet_char(f, n)
    assert chi.is_primitive()
    par = chi.parity()
    c0 = ljet(chi, 0, 3)
    cm = ljet(chi, -1, 3)
    if par == 1:
        assert cm[0].contains(0)
        val = (1, c0[1] / c0[0], cm[2] / cm[1])
    else:
        assert c0[0].contains(0)
        val = (0, c0[2] / c0[1], cm[1] / cm[0])
    _cache[key] = val
    return val


def zeta_consts():
    if "z" not in _cache:
        ser = acb_series([acb(-1), acb(1)], prec=2).zeta()
        _cache["z"] = (ser[1] / ser[0]).real
    return _cache["z"]


def phi_of(n):
    return sum(1 for k in range(1, n + 1) if math.gcd(k, n) == 1)


_sig = {}


def _signature(chi, f, G_exp, reps):
    """exact values of chi at fixed representatives of (Z/f)^*, as fractions e/G_exp mod 1"""
    return tuple(Fraction(int(chi.chi_exponent(k)), int(G_exp)) % 1 for k in reps)


def primitive_of(q, n):
    """(f, n') such that the primitive character f.n' induces the character q.n.
    Conrey labels are NOT always compatible with reduction mod f (e.g. q = 16, 27, 32, 48), so the
    primitive character is identified by its values: chi(k) = chi'(k) for all k coprime to q."""
    chi = dirichlet_char(q, n)
    f = chi.conductor()
    if f == 1:
        return 1, 1
    Gq = dirichlet_group(q).exponent()
    Gf = dirichlet_group(f).exponent()
    res = [r for r in range(1, f) if math.gcd(r, f) == 1]
    if (q, f) not in _sig:
        # representatives k = r mod f with gcd(k, q) = 1
        reps = []
        for r in res:
            k = r
            while math.gcd(k, q) != 1:
                k += f
            reps.append(k)
        table = {}
        for m in res:
            c2 = dirichlet_char(f, m)
            if c2.is_primitive():
                table[_signature(c2, f, Gf, res)] = m
        _sig[(q, f)] = (reps, table)
    reps, table = _sig[(q, f)]
    sig = _signature(chi, q, Gq, reps)
    m = table[sig]
    return f, m


def char_table(q):
    """for each chi mod q: (chi object, f, Conrey number of the inducing primitive character mod f, principal?)"""
    out = []
    for n in range(1, q + 1):
        if math.gcd(n, q) != 1:
            continue
        chi = dirichlet_char(q, n)
        f, m = primitive_of(q, n)
        out.append((chi, f, m, f == 1))
    return out


def pair_constants(q, a, table):
    ne = 0.0
    no = 0.0
    Bc = -(2 * arb.pi()).log()
    Cc = arb(zeta_consts())
    nea = arb(0); noa = arb(0)
    for chi, f, nf, principal in table:
        if principal:
            continue
        w = acb(chi(a)).conjugate()
        par, b, sec = prim_consts(f, nf)
        if par == 0:
            nea += w.real
            Bc += (w * (1 - b)).real
            Cc += (w * sec).real
        else:
            noa += w.real
            Bc -= (w * b).real
            Cc += (w * (1 + sec)).real
    ne_i, no_i = round(float(nea.mid())), round(float(noa.mid()))
    assert nea.contains(ne_i) and noa.contains(no_i) and nea.rad() < 1e-15 and noa.rad() < 1e-15
    return ne_i, no_i, Bc, Cc


def imp_terms(q, a, X):
    """list of (p^k, phi(m_p) log p) with p | q, p^k = a mod m_p, p^k <= X"""
    out = []
    qq, p = q, 2
    while qq > 1:
        if qq % p == 0:
            pe = 1
            while qq % p == 0:
                qq //= p; pe *= p
            m = q // pe
            w = phi_of(m)
            pk = p
            while pk <= X:
                if (pk - a) % m == 0:
                    out.append((pk, w, p))
                pk *= p
        p += 1
    return out


def Te_f(x):
    u = 1.0 / x
    um = np.minimum(u, 1 - 1e-16)
    return np.where(u >= 1.0, math.log(2.0), (x / 2) * ((1 + um) * np.log1p(um) + (1 - um) * np.log1p(-um)))


def To_f(x):
    u = 1.0 / x
    um = np.minimum(u, 1 - 1e-16)
    return np.where(u >= 1.0, 1 - math.log(2.0), 1 - 0.5 * np.log1p(-um * um) - np.arctanh(um) / um)


def Te_arb(x):
    if x == 1:
        return arb(2).log()
    u = 1 / x
    return (x / 2) * ((1 + u) * (1 + u).log() + (1 - u) * (1 - u).log())


def To_arb(x):
    if x == 1:
        return 1 - arb(2).log()
    u = 1 / x
    return 1 - (1 - u * u).log() / 2 - u.atanh() / u


class PP:
    """prime powers up to X with Lambda"""
    def __init__(self, X):
        X = int(X)
        s = np.ones(X + 1, dtype=bool); s[:2] = False
        for i in range(2, int(X ** 0.5) + 1):
            if s[i]:
                s[i * i::i] = False
        primes = np.nonzero(s)[0]
        ns, ps = [primes], [primes]
        for p in primes[primes <= int(X ** 0.5) + 1]:
            pk = int(p) * int(p)
            while pk <= X:
                ns.append(np.array([pk])); ps.append(np.array([p]))
                pk *= int(p)
        n = np.concatenate(ns); p = np.concatenate(ps)
        o = np.argsort(n)
        self.n = n[o].astype(np.int64); self.p = p[o].astype(np.int64)
        self.lam = np.log(self.p.astype(np.float64))
        self.X = X


def scan_pair(q, a, phi, consts, pp, X, thr):
    ne, no, B, C = consts
    Bf, Cf = float(B.mid()), float(C.mid())
    mask = (pp.n % q) == (a % q)
    n = pp.n[mask].astype(np.float64); lam = pp.lam[mask]
    c0 = np.concatenate(([0.0], np.cumsum(lam))); c1 = np.concatenate(([0.0], np.cumsum(lam * n)))
    imp = imp_terms(q, a, X)

    def F(x):
        x = np.asarray(x, dtype=np.float64)
        idx = np.searchsorted(n, x, side="right")
        psi1 = x * c0[idx] - c1[idx]
        val = phi * psi1 - x * x / 2
        ell = -ne * x * np.log(x) + Bf * x + no * np.log(x) + Cf - (1 + ne) * Te_f(x) - no * To_f(x)
        val = val - ell
        for (pk, w, p) in imp:
            val = val + w * math.log(p) * np.maximum(x - pk, 0.0)
        return val / x ** 1.5

    grid = np.exp(np.arange(0.0, math.log(X), 0.002))
    # interior points of the gaps between consecutive breakpoints
    bps = np.unique(np.concatenate((n, np.array([pk for pk, _, _ in imp], dtype=np.float64), [1.0, float(X)])))
    mids = np.concatenate([bps[:-1] + t * (bps[1:] - bps[:-1]) for t in (0.25, 0.5, 0.75)]) if len(bps) > 1 else np.array([])
    xs = np.unique(np.concatenate(([1.0], grid, bps, mids)))
    xs = xs[(xs >= 1.0) & (xs <= X)]
    fv = F(xs)
    i = int(np.argmax(fv))
    # first exceedance
    ex = np.nonzero(fv > thr + 1e-3 * max(1.0, thr))[0]
    first = float(xs[ex[0]]) if len(ex) else None
    return float(xs[i]), float(fv[i]), first, (float(fv[ex[0]]) if len(ex) else None), float(fv[0])


def certify(q, a, phi, consts, pp, x0):
    """Arb ball for F_{q,a}(log x0); x0 a Python float/int/Fraction-like -> converted exactly"""
    ne, no, B, C = consts
    x = arb(fmpq(*x0.as_integer_ratio())) if isinstance(x0, float) else arb(x0)
    xf = float(x0)
    mask = ((pp.n % q) == (a % q)) & (pp.n <= xf)
    psi1 = arb(0)
    for nn, p in zip(pp.n[mask], pp.p[mask]):
        psi1 += arb(int(p)).log() * (x - int(nn))
    val = phi * psi1 - x * x / 2
    lx = x.log()
    ell = -ne * x * lx + B * x + no * lx + C - (1 + ne) * Te_arb(x) - no * To_arb(x)
    val -= ell
    for (pk, w, p) in imp_terms(q, a, xf):
        if pk <= xf:
            val += w * arb(p).log() * (x - pk)
    return val / x ** 1.5


def main():
    QMIN, QMAX, X = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3])
    shard, nsh = (int(sys.argv[4]), int(sys.argv[5])) if len(sys.argv) > 5 else (0, 1)
    pp = PP(X)
    for q in range(QMIN, QMAX + 1):
        if q % nsh != shard:
            continue
        phi = phi_of(q)
        table = char_table(q)
        sqcount = {}
        for b in range(1, q):
            if math.gcd(b, q) == 1:
                sqcount[b * b % q] = sqcount.get(b * b % q, 0) + 1
        for a in range(1, q):
            if math.gcd(a, q) != 1 or (q, a) in LWS18:
                continue
            N = sqcount.get(a, 0)
            thr = 2.0 * N / 3
            consts = pair_constants(q, a, table)
            xbest, fbest, xfirst, ffirst, f1 = scan_pair(q, a, phi, consts, pp, X, thr)
            # choose witness: x = 1 (closed form, no primes) if it already works, else the maximum found by the scan
            if f1 > thr + 1e-3 * max(1.0, thr):
                x0 = 1
            else:
                x0 = xbest
            Fc = certify(q, a, phi, consts, pp, x0)
            F1c = certify(q, a, phi, consts, pp, 1)
            ok = bool(Fc > arb(2) * N / 3)
            rec = {"q": q, "a": a, "phi": phi, "N": N, "threshold": thr, "x0": x0 if isinstance(x0, int) else repr(x0),
                   "F_at_x0": Fc.str(15), "F_at_x0_lower": float(Fc.lower()), "witness": ok,
                   "F_at_1": F1c.str(15), "scan_max_x": xbest, "scan_max_F": fbest, "X": X,
                   "n_e": consts[0], "n_o": consts[1], "B": consts[2].str(12, radius=False), "C": consts[3].str(12, radius=False)}
            print(json.dumps(rec))
            sys.stdout.flush()


if __name__ == "__main__":
    main()
