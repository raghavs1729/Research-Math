#!/usr/bin/env python3
"""
t2_assemble.py -- assemble the GRH-only exclusion results (HKMN23.md, Theorems 3-6).

Reads raw/t2_closed_400_*.jsonl (closed-form witness for a = 1, central non-vanishing),
raw/t2_witness_120_*.jsonl (prime-side witnesses), raw/t2_zeros_*.jsonl (zero-side witnesses),
l1const.json.  Writes t2_summary.out and t2_summary.json.

Also checks the explicit lower bound of Theorem 5(b):
    Phi(q) >= LB(q) := phi/2 + 1 - 2 log(4 phi (3 log q + 0.11)) - (log q)/3
                       - (phi/q^2)(zeta(2) log q + 1.81) - R(q),
  (i)   LB(q) <= Phi(q) for 3 <= q <= 400      [consistency test of the derivation],
  (ii)  LB(q) >  (2/3) N(q,1) for 400 < q <= 10^6   [floating point for the printed margins, then decided in Arb],
  (iii) the elementary inequality for q > 10^6.
"""
import json, glob, math, sys
from collections import defaultdict

out = open("t2_summary.out", "w")
def P(*a):
    s = " ".join(str(x) for x in a)
    print(s); out.write(s + "\n")

LWS_MODULI = [1, 2, 3, 4, 5, 6, 8, 10, 12, 15, 20, 24, 30]
LWS18 = {(1, 1), (2, 1), (3, 1), (4, 1), (5, 1), (5, 4), (6, 1), (8, 1), (10, 1), (10, 9), (12, 1), (15, 1), (15, 4),
         (20, 1), (20, 9), (24, 1), (30, 1), (30, 19)}
CAND = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 15, 16, 18, 20, 21, 24, 30, 40, 42, 48, 60]

# ---------- sieve phi, omega, N(q,1) ----------
QMAX = 10 ** 6
phi = list(range(QMAX + 1)); omega = [0] * (QMAX + 1); N1 = [1] * (QMAX + 1)
for p in range(2, QMAX + 1):
    if phi[p] == p:          # prime
        for k in range(p, QMAX + 1, p):
            phi[k] -= phi[k] // p
            omega[k] += 1
        pe = p
        while pe <= QMAX:
            # number of square roots of 1 mod p^e
            if p == 2:
                e = pe.bit_length() - 1
                r = 1 if e == 1 else (2 if e == 2 else 4)
            else:
                r = 2
            for k in range(pe, QMAX + 1, pe):
                if (k // pe) % p != 0:
                    N1[k] *= r
            pe *= p

# ---------- closed form Phi(q), central non-vanishing ----------
rows = {}
for f in glob.glob("raw/t2_closed_400_*.jsonl"):
    for l in open(f):
        d = json.loads(l); rows[d["f"]] = d
QC = max(rows)
assert all(d["central_nonvanishing_certified"] for d in rows.values())
nprim = sum(d["nprim"] for d in rows.values())
minL = min((d["minL_half_lower"], d["f"]) for d in rows.values() if d["nprim"] > 0)
P(f"## Central non-vanishing: L(1/2, chi) != 0 certified (Arb) for all {nprim} primitive characters of conductor 3 <= f <= {QC};")
P(f"   smallest certified lower bound for |L(1/2,chi)|: {minL[0]:.6g} at conductor {minL[1]}.  (zeta(1/2) = -1.46035... != 0.)")
from flint import arb, acb, acb_series, ctx
ctx.prec = 128
_ser = acb_series([acb(-1), acb(1)], prec=2).zeta()
C_ZETA = (2 * arb.pi()).log() + arb(2).log() - arb(1) / 2 - (_ser[1] / _ser[0]).real     # c(zeta), Arb ball
def _ball(x):
    """outward-rounded ball from [mid string, radius string]"""
    mid = arb(x[0]); rad = arb(x[1]) * arb("1.01") + abs(mid) * arb(10) ** (-28) + arb(10) ** (-40)
    return mid + (rad * arb(-1)).union(rad)
CF = {f: _ball(rows[f]["c_f"]) for f in rows}
def Phi_arb(q):
    s = arb(C_ZETA)
    for f in range(3, q + 1):
        if q % f == 0 and f % 4 != 2:
            s += CF[f]
    return s
def Phi(q):
    b = Phi_arb(q)
    return float(b.mid()), float(b.rad()) * 1.0001 + 1e-25
P("")
P(f"## Closed-form witness for a = 1:  Phi(q) = F_(q,1)(0) = sum_chi c(chi*),  q <= {QC}")
notw = []
for q in range(1, QC + 1):
    if not (Phi_arb(q) > arb(2) * N1[q] / 3):          # comparison of Arb balls: True only if certain
        notw.append(q)
    else:
        assert not (Phi_arb(q) < arb(2) * N1[q] / 3)
# for the moduli without witness the opposite inequality is certified as well
assert all(bool(Phi_arb(q) < arb(2) * N1[q] / 3) for q in notw)
P(f"   moduli q <= {QC} with Phi(q) < (2/3)N(q,1) (no witness at y = 0; both directions decided on Arb balls): {notw}")
P(f"   = the 13 moduli of Liu-Wang-Shi {LWS_MODULI} together with {sorted(set(notw) - set(LWS_MODULI))}")
assert sorted(set(notw) - set(LWS_MODULI)) == [7, 9, 14, 16, 18]
L1 = {int(q): float(v["L1"][0]) for q, v in json.load(open("l1const.json"))["moduli"].items()}
P("   q   phi  N(q,1)  (2/3)N     Phi(q)      L_1(q)     Phi/L_1   verdict under GRH alone for (q,1)")
for q in CAND + [11, 13]:
    ph, r = Phi(q)
    thr = 2 * N1[q] / 3
    if q in LWS_MODULI:
        v = "member (Theorem 1: negative for every x > 2)"
    elif ph - r > thr:
        v = "EXCLUDED (witness y_0 = 0)"
    else:
        v = "undecided (needs sup F > (2/3)N; Phi falls short by %.4f)" % (thr - ph)
    P(f"  {q:3d} {phi[q]:4d} {N1[q]:5d}   {thr:8.4f}  {ph:10.6f}  {L1[q]:10.6f}   {ph / L1[q]:.4f}   {v}")
    assert ph <= L1[q] + 1e-9
minmargin = min(((Phi(q)[0] - 2 * N1[q] / 3) / (2 * N1[q] / 3), q) for q in range(3, QC + 1) if q not in notw)
P(f"   smallest relative margin (Phi - (2/3)N)/((2/3)N) over excluded q <= {QC}: {minmargin[0]:.4f} at q = {minmargin[1]}")

# ---------- explicit lower bound ----------
Z2 = math.pi ** 2 / 6
kap = sum(math.log(j + 1) / j ** 2 for j in range(1, 10 ** 6 + 1))
kap_up = kap + (math.log(10 ** 6 + 1) + 1) / 10 ** 6
P("")
P(f"## Explicit lower bound (Theorem 5(b)).  kappa = sum log(j+1)/j^2 <= {kap_up:.6f} <= 1.81")
assert kap_up < 1.81
def LB(q):
    ph = phi[q]; lq = math.log(q)
    Lp = ph * (3 * lq + 0.11); X1 = 16 * Lp * Lp; lX = math.log(X1)
    ne, no = ph / 2 - 1, ph / 2
    Bp, Cp = ph * (2 * lq + 4.4), ph * (lq + 2.4)
    R = 2 * ((ne + omega[q] * ph) * (lX + 1) + Bp) / X1 + (no * (lX + 0.5) + Cp) / X1 ** 2 + 0.4 * (1 + ne) / X1 ** 3 + 0.1 * no / X1 ** 4
    return ph / 2 + 1 - 2 * math.log(4 * Lp) - lq / 3 - (ph / q ** 2) * (Z2 * lq + 1.81) - R, R
worst = min((Phi(q)[0] - LB(q)[0], q) for q in range(3, QC + 1))
P(f"   (i)  consistency: min over 3 <= q <= {QC} of Phi(q) - LB(q) = {worst[0]:.4f} at q = {worst[1]}  (must be >= 0)")
assert worst[0] >= 0
best = min((LB(q)[0] - 2 * N1[q] / 3, q) for q in range(QC + 1, QMAX + 1))
bestrel = min(((LB(q)[0]) / (2 * N1[q] / 3), q) for q in range(QC + 1, QMAX + 1))
maxR = max(LB(q)[1] for q in range(QC + 1, QMAX + 1))
P(f"   (ii) min over {QC} < q <= 10^6 of LB(q) - (2/3)N(q,1) = {best[0]:.4f} at q = {best[1]};  min ratio LB/((2/3)N) = {bestrel[0]:.3f} at q = {bestrel[1]};  max R(q) = {maxR:.2e}")
assert best[0] > 0
small = max((phi[q] / q ** 2) * (Z2 * math.log(q) + 1.81) + LB(q)[1] for q in range(QC + 1, QMAX + 1))
P(f"        max over {QC} < q <= 10^6 of the two small terms (phi/q^2)(zeta(2) log q + 1.81) + R(q) = {small:.4f}  (< 0.04)")
assert small < 0.04
# the same inequality decided in ball arithmetic (independent audit, item M1): exact phi, omega, N
# from the integer sieve; logarithms and the bound LB(q) as Arb balls; comparison of balls
def LB_arb(q):
    ph = arb(phi[q]); lq = arb(q).log()
    Lp = ph * (3 * lq + arb("0.11")); X1 = 16 * Lp * Lp; lX = X1.log()
    ne, no = ph / 2 - 1, ph / 2
    Bp, Cp = ph * (2 * lq + arb("4.4")), ph * (lq + arb("2.4"))
    R = 2 * ((ne + omega[q] * ph) * (lX + 1) + Bp) / X1 + (no * (lX + arb("0.5")) + Cp) / (X1 * X1) \
        + arb("0.4") * (1 + ne) / X1 ** 3 + arb("0.1") * no / X1 ** 4
    z2 = arb.pi() ** 2 / 6
    return ph / 2 + 1 - 2 * (4 * Lp).log() - lq / 3 - (ph / arb(q) ** 2) * (z2 * lq + arb("1.81")) - R
ctx.prec = 64
nfail = 0
minarb = None
for q in range(QC + 1, QMAX + 1):
    d = LB_arb(q) - arb(2) * N1[q] / 3
    if not (d > 0):
        nfail += 1
    lo = float(d.lower())
    if minarb is None or lo < minarb[0]:
        minarb = (lo, q)
ctx.prec = 128
P(f"        the same in Arb (64-bit balls, exact integer phi, omega, N): LB(q) > (2/3)N(q,1) certified for every {QC} < q <= 10^6; "
  f"failures {nfail}; smallest certified lower bound of the difference {minarb[0]:.4f} at q = {minarb[1]}")
assert nfail == 0
minphi = min((phi[q], q) for q in range(QC + 1, QMAX + 1))
P(f"        (min phi(q) for {QC} < q <= 10^6: {minphi[0]} at q = {minphi[1]};  max N(q,1)/phi(q) for q > 168 up to 10^6: {max(N1[q] / phi[q] for q in range(169, QMAX + 1)):.4f})")
g = lambda t: t / 3 + 1 - 2 * math.log(4 * t * (6 * math.log(t) + 2.2)) - (2 * math.log(t) + 0.7) / 3 - 0.04
P(f"   (iii) q > 10^6: phi >= sqrt(q/2) >= 707, N <= phi/4;  g(707) = {g(707):.2f} > 0, g increasing for t >= 20")

# ---------- prime-side witnesses ----------
W = {}
for f in glob.glob("raw/t2_witness_120_*.jsonl"):
    for l in open(f):
        d = json.loads(l); W[(d["q"], d["a"])] = d
QW = max(q for q, _ in W)
ZW = {}      # certified zero-side witnesses
ZB = {}      # best certified zero-side lower bound per pair (witness or not)
for f in glob.glob("raw/t2_zeros_*.jsonl"):
    for l in open(f):
        d = json.loads(l)
        k = (d["q"], d["a"])
        if k not in ZB or d["F_lower_bound"] > ZB[k]["F_lower_bound"]:
            ZB[k] = d
        if d["witness"]:
            if k not in ZW or d["F_lower_bound"] > ZW[k]["F_lower_bound"]:
                ZW[k] = d
P("")
P(f"## Witnesses, 3 <= q <= {QW}: all reduced classes except the 18 pairs of Liu-Wang-Shi ({len(W)} pairs)")
nons = [d for d in W.values() if d["N"] == 0]; sq = [d for d in W.values() if d["N"] > 0]
P(f"   non-square classes: {len(nons)}; prime-side witness F(log x_0) > 0 certified for {sum(d['witness'] for d in nons)};"
  f" the other {sum(not d['witness'] for d in nons)} are covered by Theorem 4 (central non-vanishing is certified for conductors <= {QC})")
P(f"   square classes: {len(sq)} (a = 1: {sum(1 for d in sq if d['a'] == 1)});  prime-side witness certified for {sum(d['witness'] for d in sq)}")
def status(q, a):
    """'P' prime-side, 'Z' zero-side, 'T' theorem 4 (non-square), 'U' undecided; with reduction 2m -> m for odd m"""
    d = W[(q, a)]
    if d["witness"]:
        return "P"
    if (q, a) in ZW:
        return "Z"
    if q % 4 == 2:
        m = q // 2
        if (m, a % m) in ZW or ((m, a % m) in W and W[(m, a % m)]["witness"]):
            return "Z" if (m, a % m) in ZW else "P"
    if d["N"] == 0:
        return "T"
    return "U"
und = defaultdict(list); cnt = defaultdict(int)
for (q, a) in sorted(W):
    s = status(q, a); cnt[s] += 1
    if s == "U":
        und[q].append(a)
P(f"   status counts: prime-side P = {cnt['P']}, zero-side Z = {cnt['Z']}, Theorem 4 only T = {cnt['T']}, undecided U = {cnt['U']}")
P(f"   undecided pairs (all are square classes): " + "; ".join(f"q={q}: a in {und[q]}" for q in sorted(und)))
u60 = [(q, a) for q in sorted(und) if q <= 60 for a in und[q]]
u120 = [(q, a) for q in sorted(und) if q > 60 for a in und[q]]
P(f"   undecided with q <= 60: {len(u60)} pairs: {u60}")
P(f"   undecided with 60 < q <= {QW}: {len(u120)} pairs (zero side not run for these moduli): {u120}")
n60 = sum(1 for (q, a) in W if q <= 60)
P(f"   pairs with 3 <= q <= 60 outside the 18 members: {n60};  excluded: {n60 - len(u60)};  undecided: {len(u60)}")
P("   best certified lower bound for sup F at the undecided pairs with q <= 60 (zero side if available, else prime-side scan maximum):")
for (q, a) in u60:
    qq, aa = (q // 2, a % (q // 2)) if q % 4 == 2 else (q, a)
    z = ZB.get((qq, aa))
    pm = W[(q, a)]["scan_max_F"]
    P(f"      ({q},{a}): (2/3)N = {W[(q,a)]['threshold']:.4f};  prime side max {pm:.4f}" + (f";  zero side F >= {z['F_lower_bound']:.4f} at y_0 = {z['y0']} (|y| <= {z['Ymax']:.0f})" if z else ""))
P("")
P("## The 23 candidate moduli of Liu-Wang-Shi (q <= 60): verdict for every reduced class under GRH(q) alone")
P("   q : members (Theorem 1) | excluded squares (witness) | excluded non-squares | undecided")
summ = {}
for q in CAND:
    if q < 3:
        P(f"  {q:3d}: members [1]"); continue
    mem = sorted(a for (qq, a) in LWS18 if qq == q)
    exs, exn, un = [], [], []
    for a in range(1, q):
        if math.gcd(a, q) != 1 or (q, a) in LWS18:
            continue
        s = status(q, a)
        if s == "U":
            un.append(a)
        elif W[(q, a)]["N"] > 0:
            exs.append(f"{a}{s}")
        else:
            exn.append(f"{a}{s}")
    summ[q] = {"members": mem, "excluded_squares": exs, "excluded_nonsquares": exn, "undecided": un}
    P(f"  {q:3d}: members {mem} | squares excluded {exs} | non-squares excluded {len(exn)} ({sum(1 for x in exn if x.endswith('P'))} with explicit prime-side witness) | undecided {un}")
tot_un = sum(len(v["undecided"]) for v in summ.values())
P(f"   total undecided pairs at the candidate moduli: {tot_un}")
P("")
P("## Zero-side witnesses (certified): (q,a): threshold, y_0, certified lower bound for F(y_0)")
for k in sorted(ZW):
    d = ZW[k]
    P(f"   ({k[0]},{k[1]}): (2/3)N = {d['threshold']:.4f}, y_0 = {d['y0']}, F_Z(y_0) = {d['F_Z_at_y0']}, remainder <= {d['remainder_upper']:.5f}, F(y_0) >= {d['F_lower_bound']:.5f}  [T = {d['T']:.0f}, {d['num_brackets']} zeros]")
json.dump({"candidates": summ, "undecided": {str(q): und[q] for q in und}, "counts": dict(cnt),
           "a1_no_witness": notw}, open("t2_summary.json", "w"), indent=1)
out.close()
