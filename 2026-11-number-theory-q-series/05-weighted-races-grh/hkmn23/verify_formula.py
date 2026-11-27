#!/usr/bin/env python3
"""
verify_formula.py -- numerical check (NOT part of any proof) of the explicit formulas of
HKMN23.md, Lemma E, and of the constants b, b0, cm1 stored in l1const.json.

For a primitive character chi of conductor f > 1:
  odd :  psi_1(x,chi) = -Z(x) - b x + log x + 1 + b0 - sum_{n>=1} x^(-2n)/(2n(2n+1))
  even:  psi_1(x,chi) = -Z(x) - x log x + (1-b) x + cm1 - sum_{n>=1} x^(1-2n)/(2n(2n-1))
with Z(x) = sum_rho x^(rho+1)/(rho(rho+1)).
Zeros with |gamma| < T are found from sign changes of Hardy's Z (floating point, no
completeness claim); the test is that "direct - formula" is as small as the omitted tail
allows (a wrong constant term or sign shows up as an O(1) or O(x) discrepancy).
At x = 1 the discrepancy must equal  -sum_{|gamma|>T} 1/(rho(rho+1)) ~ +sum_{|gamma|>T} gamma^-2,
which is compared with the tail of V = sum 1/|rho|^2 minus the partial sum.

Also checks the unconditional closed form
  -Re sum_rho 1/(rho(rho+1)) = Re[(L'/L)(2,chi) - (L'/L)(1,chi)] + (1/2)[psi(1+a/2) - psi((1+a)/2)]
against  -(b0 - b + log 2) (odd)  resp.  -(1 - b + cm1 - log 2) (even).
"""
import json, math, sys
import mpmath as mp
from flint import acb, arb, dirichlet_char, ctx, acb_series

ctx.prec = 100
mp.mp.dps = 30
T = float(sys.argv[1]) if len(sys.argv) > 1 else 400.0
STEP = 0.02

data = json.load(open("l1const.json"))


def getc(entry, key):
    return complex(float(entry[key]["re"][0]), float(entry[key]["im"][0]))


def chi_val(f, n, m):
    c = dirichlet_char(f, n)
    if math.gcd(m, f) != 1:
        return 0j
    return complex(c(m))


def zeros(f, n, T):
    """ordinates gamma in (-T, T) of zeros 1/2+i gamma of L(s, chi_f(n,.)), by sign changes of Z."""
    c = dirichlet_char(f, n)
    zs = []
    def Z(t):
        return float(c.hardy_z(acb(t)).real)
    t = -T
    prev = Z(t)
    while t < T:
        t2 = t + STEP
        cur = Z(t2)
        if prev == 0.0 or prev * cur < 0:
            lo, hi, flo = t, t2, prev
            for _ in range(60):
                mid = 0.5 * (lo + hi)
                fm = Z(mid)
                if flo * fm <= 0:
                    hi = mid
                else:
                    lo, flo = mid, fm
            zs.append(0.5 * (lo + hi))
        t, prev = t2, cur
    return zs


def psi1_direct(f, n, x):
    """sum_{m<=x} chi(m) Lambda(m) (x-m)."""
    tot = 0j
    X = int(math.floor(x))
    sieve = list(range(X + 1))
    for p in range(2, X + 1):
        if sieve[p] == p:
            for k in range(p * p, X + 1, p):
                if sieve[k] == k:
                    sieve[k] = p
    for m in range(2, X + 1):
        p = sieve[m]
        mm = m
        while mm % p == 0:
            mm //= p
        if mm == 1:
            tot += chi_val(f, n, m) * math.log(p) * (x - m)
    return tot


def formula(f, n, x, zs, entry):
    par = entry["parity"]
    Zx = sum(x ** (1.5 + 1j * g) / ((0.5 + 1j * g) * (1.5 + 1j * g)) for g in zs)
    b = getc(entry, "b")
    if par == 1:
        b0 = getc(entry, "b0")
        tail = sum(x ** (-2 * k) / (2 * k * (2 * k + 1)) for k in range(1, 4000)) if x > 1 else 1 - math.log(2)
        return -Zx - b * x + math.log(x) + 1 + b0 - tail
    else:
        cm1 = getc(entry, "cm1")
        tail = sum(x ** (1 - 2 * k) / (2 * k * (2 * k - 1)) for k in range(1, 4000)) if x > 1 else math.log(2)
        return -Zx - x * math.log(x) + (1 - b) * x + cm1 - tail


def closed_form_c(f, n):
    chi = dirichlet_char(f, n)
    par = chi.parity()
    def ld(s):
        ser = acb_series([acb(s), acb(1)], prec=2).dirichlet_l(chi)
        return ser[1] / ser[0]
    val = (ld(2) - ld(1)).real + ((arb(1) + arb(par) / 2).digamma() - ((arb(1) + par) / 2).digamma()) / 2
    return float(val)


tests = [(3, 2), (4, 3), (5, 4), (5, 2), (8, 5), (15, 2), (20, 19), (24, 11), (7, 3)]
print(f"# T = {T}, grid step {STEP}")
for (f, n) in tests:
    e = data["chars"][f"{f}.{n}"]
    zs = zeros(f, n, T)
    V = float(e["V"][0])
    Vpart = sum(1 / (0.25 + g * g) for g in zs)
    Spart = sum(1 / abs((0.5 + 1j * g) * (1.5 + 1j * g)) for g in zs)
    print(f"chi = {f}.{n} parity {e['parity']} order {e['order']}: {len(zs)} zeros with |gamma|<{T}; "
          f"V - partial = {V - Vpart:+.3e} (tail, must be > 0 and ~ (1/pi)log(fT/2pi)/T*2 = {2/math.pi*(math.log(f*T/(2*math.pi))+1)/T:.3e}); "
          f"S(certified) - partial S = {float(e['S'][0]) - Spart:+.3e}")
    b = getc(e, "b")
    if e["parity"] == 1:
        lhs = -(getc(e, "b0") - b + math.log(2)).real
    else:
        lhs = -(1 - b + getc(e, "cm1") - math.log(2)).real
    print(f"    -Re sum 1/(rho(rho+1)): from explicit-formula constants {lhs:+.12f}; closed form at s=1,2 {closed_form_c(f, n):+.12f}; "
          f"from zeros {-sum((1/((0.5+1j*g)*(1.5+1j*g))).real for g in zs):+.6f} (+tail)")
    for x in (1.0, 2.5, 5.5, 10.0, 30.0, 100.0):
        d = psi1_direct(f, n, x)
        fo = formula(f, n, x, zs, e)
        bound = x ** 1.5 * (V - Vpart)
        print(f"    x={x:6.1f}  direct={d.real:+.6f}{d.imag:+.6f}i  formula={fo.real:+.6f}{fo.imag:+.6f}i  |diff|={abs(d-fo):.2e}  allowed (x^1.5*tail)={bound:.2e}")
    sys.stdout.flush()

# zeta
zz = [float(mp.im(mp.zetazero(k))) for k in range(1, 201)]
ez = data["chars"]["1.1"]
z0 = float(ez["zeta_logder_0"]["re"][0]); zm1 = float(ez["zeta_logder_m1"]["re"][0])
print(f"zeta: (zeta'/zeta)(0) = {z0:.12f} (log 2pi = {math.log(2*math.pi):.12f}), (zeta'/zeta)(-1) = {zm1:.12f}")
for x in (1.0, 2.5, 10.0, 100.0):
    X = int(x)
    d = 0.0
    for m in range(2, X + 1):
        # von Mangoldt
        k, mm, p = 0, m, None
        for r in range(2, m + 1):
            if mm % r == 0:
                p = r
                break
        while mm % p == 0:
            mm //= p
        if mm == 1:
            d += math.log(p) * (x - m)
    Zx = sum(2 * (x ** (1.5 + 1j * g) / ((0.5 + 1j * g) * (1.5 + 1j * g))).real for g in zz)
    tail = sum(x ** (1 - 2 * k) / (2 * k * (2 * k - 1)) for k in range(1, 4000)) if x > 1 else math.log(2)
    fo = x * x / 2 - Zx - x * z0 + zm1 - tail
    print(f"    x={x:6.1f}  psi_1 direct={d:+.6f}  formula(200 zeros)={fo:+.6f}  diff={d-fo:+.2e}")
