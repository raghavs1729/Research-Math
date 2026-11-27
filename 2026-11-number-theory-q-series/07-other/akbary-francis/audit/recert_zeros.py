#!/usr/bin/env python3
"""
recert_zeros.py q  -- independent re-certification of zeros_q<q>.txt (audit).

Independent of the audited code: does NOT use acb_dirichlet_hardy_z, dirichlet_char,
Conrey labels from flint, or anything in ../*.py.  Characters are built by hand from a
primitive root g mod q (prime q):  chi_j(g^k) = exp(2 pi i j k/(q-1)).  The file's
label n is mapped to j = ind_g(n) (Conrey convention for prime q, smallest primitive
root); if that fails, the conjugate convention j = -ind_g(n) is tried and reported.

L(s, chi) = q^{-s} sum_{m=1}^{q-1} chi(m) zeta(s, m/q)  (Arb Hurwitz zeta, as power series
in s for the derivative).  For each window (d - r, d + r), r = 1e-30, put
D = [1/2 +- r] + i [d +- r] (a complex box) and run the Krawczyk test
    K = c1 - L(c1)/c + (1 - L'(D)/c)(D - c1) subset of interior(D),
c1 = Newton-refined point, c = mid(L'(c1)).  By the integral form of the mean value
theorem (D convex, L'(D) enclosed in a convex box) the Newton map s -> s - L(s)/c maps D
into K subset D, so (Brouwer) L has a zero in D; |1 - L'(D)/c| < 1 gives uniqueness.
Under GRH the zero is 1/2 + i gamma with |gamma - d| < r.  Also checks that windows of one
label are pairwise disjoint and records |d - gamma_refined|.
"""
import sys, hashlib, time
from flint import arb, acb, acb_series, ctx

ctx.prec = 320
q = int(sys.argv[1])
path = sys.argv[2]
R = arb("1e-30")

raw = open(path, "rb").read()
print("file", path, "sha256", hashlib.sha256(raw).hexdigest())
rows = []
for line in raw.decode("ascii").splitlines():
    if line.startswith("#") or not line.strip():
        continue
    lab, d = line.split()
    rows.append((lab, d))
print("windows:", len(rows))


def primroot(q):
    for g in range(2, q):
        if len({pow(g, k, q) for k in range(q - 1)}) == q - 1:
            return g


g = primroot(q)
ind = {pow(g, k, q): k for k in range(q - 1)}


def chivals(j):
    # chi_j(m) for m = 1..q-1, exact via exp(2 pi i j ind(m)/(q-1))
    return [(2 * acb.pi() * acb(0, 1) * acb(j * ind[m] % (q - 1)) / (q - 1)).exp() for m in range(1, q)]


CH = {}
def L_series(s, j, n=2):
    """power series L(s + t) truncated to length n; j = None -> zeta."""
    S = acb_series([s, 1], prec=n)
    if j is None:
        return acb_series.zeta(S)
    if j not in CH:
        CH[j] = chivals(j)
    tot = acb_series([0], prec=n)
    for m in range(1, q):
        tot += CH[j][m - 1] * acb_series.zeta(S, acb(m) / q)
    qs = (-S * acb(q).log()).exp()
    return qs * tot


def LL(s, j):
    ser = L_series(s, j)
    c = ser.coeffs()
    return c[0], (c[1] if len(c) > 1 else acb(0))


def mid(z):
    return acb(z.real.mid(), z.imag.mid())


def krawczyk(d, j):
    c1 = acb(arb("0.5"), arb(d))
    for _ in range(3):   # Newton refinement of the centre
        v, dv = LL(c1, j)
        c1 = mid(c1 - v / dv)
    D = acb(arb("0.5", R), arb(arb(d).mid(), R))  # box [1/2 +- r] + i[d +- r]
    v, dv = LL(c1, j)
    c = mid(dv)
    _, dD = LL(D, j)
    K = c1 - v / c + (1 - dD / c) * (D - c1)
    lo_r, hi_r = arb("0.5") - R, arb("0.5") + R
    lo_i, hi_i = arb(d) - R, arb(d) + R
    inside = (K.real.lower() > lo_r.upper() and K.real.upper() < hi_r.lower()
              and K.imag.lower() > lo_i.upper() and K.imag.upper() < hi_i.lower())
    contr = abs(1 - dD / c).upper() < 1
    dev = abs(arb(c1.imag.mid()) - arb(d)).upper()
    offline = abs(arb(c1.real.mid()) - arb("0.5")).upper()
    return bool(inside and contr), float(dev), float(offline), v


t0 = time.time()
bad = []
maxdev = 0.0
maxoff = 0.0
conv = {}
bylab = {}
for i, (lab, d) in enumerate(rows):
    bylab.setdefault(lab, []).append(arb(d))
    if lab == "zeta":
        ok, dev, off, _ = krawczyk(d, None)
    else:
        n = int(lab)
        tried = conv.get(lab)
        if tried is None:
            ok, dev, off, _ = krawczyk(d, ind[n])
            if ok:
                conv[lab] = "+"
            else:
                ok2, dev, off, _ = krawczyk(d, (-ind[n]) % (q - 1))
                if ok2:
                    conv[lab] = "-"
                ok = ok2
        else:
            jj = ind[n] if tried == "+" else (-ind[n]) % (q - 1)
            ok, dev, off, _ = krawczyk(d, jj)
    maxdev = max(maxdev, dev)
    maxoff = max(maxoff, off)
    if not ok:
        bad.append((lab, d))
        print("FAIL", lab, d, flush=True)
    if i % 500 == 0:
        print(i, lab, d[:20], "dev %.1e" % dev, "%.0fs" % (time.time() - t0), flush=True)

# disjointness within a label
overlap = 0
for lab, L in bylab.items():
    L = sorted(L, key=lambda a: float(a.mid()))
    for a, b in zip(L, L[1:]):
        if not (b - a > 2 * R):
            overlap += 1
print("label conventions (+ = Conrey j = ind_g(n)):", conv)
print("per-label counts:", {k: len(v) for k, v in bylab.items()})
print("failures:", len(bad), "overlaps:", overlap)
print("max |d - gamma_refined| = %.2e ; max |Re(refined) - 1/2| = %.2e" % (maxdev, maxoff))
print("ALL WINDOWS CERTIFIED" if not bad and not overlap else "NOT ALL CERTIFIED")
print("time %.0fs" % (time.time() - t0))
