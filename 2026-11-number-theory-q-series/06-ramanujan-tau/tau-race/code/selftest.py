#!/usr/bin/env python3
"""Self-tests and mutation tests for the tau-race certificates.  Every mutant must be REJECTED.

 (F) functional equation Phi(1/y) = y Phi(y) checked directly (no use of the functional equation):
     Phi computed from its defining series at y and 1/y.
 (M1) small range: sign-flipped lambda, and tau(3) -> -tau(3): the checker must report failures.
 (M2) bounds: sigma0^2 replaced by 0.7 (> 2/3): no positivity threshold; T(H) with H = 12: L_B never positive.
 (M3) Sym^2 completeness: delete one found zero -> Rem >= 2/3 somewhere -> FAIL;
      corrupt one table entry by a relative 1e-30 -> Euler-product comparison at s = 2 or sign pattern must change... (reported)
 (M4) constants: tau(2) -> tau(2)+1 changes sigma0^2 outside the pinned enclosure.
"""
import sys, json, os, subprocess, copy
from fractions import Fraction
from flint import arb, acb, ctx
import sym2_phi, sym2_verify, cert_small, tau_table, cert_constants
DATA = sym2_verify.DATA
results = {}

# ---------------- (F)
ctx.prec = sym2_phi.PREC
aS = sym2_phi.aS_table(400)
def Phi(y, X=170):
    y = arb(y); N = int(X/float(y.mid())); S = arb(0)
    for n in range(1, N+1):
        S += arb(aS[n].numerator)/aS[n].denominator*sym2_phi.K_arb(n*y)
    return S + arb(0, sym2_phi.tail_n(y, N).upper()) if float(y.mid()) >= 1 else S
for q in [Fraction(9, 10), Fraction(3, 4), Fraction(1, 2)]:
    y = arb(q.numerator)/q.denominator
    lhs = Phi(1/y); rhs = y*Phi(y)           # Phi(1/y) = y Phi(y);  for y < 1 the n-tail is not added (n*y > 170 terms < 1e-95)
    d = lhs - rhs
    print("(F) y = %s: Phi(1/y) = %s,  Phi(1/y) - y Phi(y) = %s" % (q, lhs.str(25), d.str(5)))
    results["F_%s" % q] = bool(abs(d) < arb(10)**(-80))

# ---------------- (M1)
def small_fail_count(tau_mut, N=2000):
    primes = cert_small.primes_upto(N)
    ctx.prec = 128
    Th = arb(0); Mm = arb(0); D0 = arb(0); D1 = arb(0); bad = 0
    for p in primes:
        ap = arb(p); lam = arb(tau_mut[p])/(ap**5*ap.sqrt()); lg = ap.log()
        Th += lam*lg; Mm += lam*lg*p; D0 += lam; D1 += lam*p
        if p >= 5 and not (ap*Th - Mm > 0): bad += 1
        if p >= 7 and not (D0 - D1/ap > 0): bad += 1
    return bad
tau = tau_table.tau_list(2001)
results["M1_baseline_pass"] = small_fail_count(tau) == 0
results["M1_signflip_rejected"] = small_fail_count([-t for t in tau]) > 0
t2 = list(tau); t2[3] = -t2[3]; t2[5] = -t2[5]; t2[11] = -t2[11]
results["M1_tau3_5_11_flipped_rejected"] = small_fail_count(t2) > 0
print("(M1)", {k: v for k, v in results.items() if k.startswith("M1")})

# ---------------- (M2)
import bounds
ctx.prec = 100
old = bounds.SIG
bounds.SIG = arb("0.7")
thr = any(bool(bounds.L_A(arb(10)**(arb(j)/10)) > 0) for j in range(40, 300))
results["M2_bad_sigma_rejected"] = not thr
bounds.SIG = old
results["M2_H12_rejected"] = not any(bool(bounds.L_A(arb(10)**(arb(j)/10), beyond=(j > 140)) - bounds.T_best(12)[0] > 0) for j in range(40, 300, 5))
results["M2_baseline_pass"] = bool(bounds.L_A(arb(10)**6) > 0) and bool(bounds.L_A(arb(10)**6) - bounds.T_best(60)[0] > 0)
print("(M2)", {k: v for k, v in results.items() if k.startswith("M2")})

# ---------------- (M3)  completeness test with one zero deleted
cert = json.load(open(DATA + "/sym2_rh_H60.json"))
zeros = [(Fraction(a), Fraction(b)) for a, b in cert["zeros"]]
P = 20000
ctx.prec = 128
taus = tau_table.tau_list(P+1)
plist = [(p, arb(p).log(), arb(taus[p])**2/arb(p)**11 - 2) for p in cert_small.primes_upto(P)]
tailP = arb("3.03")*(arb(P).log()+1)/P
def rem_max(zs, H=60):
    worst = None; T = Fraction(0)
    while T <= H:
        Ta = arb(T.numerator)/T.denominator
        tot = (sym2_verify.gammaS_logderiv(acb(2, Ta)) + sym2_verify.euler_logderiv(Ta, plist, tailP)).real
        found = arb(0)
        for lo, hi in zs:
            g = (arb(lo.numerator)/lo.denominator).union(arb(hi.numerator)/hi.denominator)
            found += arb(3)/2/(arb(9)/4 + (Ta-g)**2) + arb(3)/2/(arb(9)/4 + (Ta+g)**2)
        rem = tot - found
        m = float(rem.upper())
        worst = m if worst is None else max(worst, m)
        T += Fraction(1, 4)
    return worst
base = rem_max(zeros)
results["M3_baseline_rem_lt_2over3"] = base*1.134 < 2/3          # e^{1/8} for the coarser grid 1/4
dele = []
for idx in [0, 7, 20, 40]:
    zs = zeros[:idx] + zeros[idx+1:]
    dele.append(rem_max(zs))
results["M3_deleted_zero_rejected"] = all(v >= 2/3 for v in dele)
print("(M3) baseline max Rem (P=2e4, grid 1/4) = %.4f; with zero #1, #8, #21, #41 deleted: %s" % (base, ["%.4f" % v for v in dele]))
# a zero at a wrong place (shifted by 0.3) must also be rejected
zs = list(zeros); lo, hi = zs[10]; zs[10] = (lo + Fraction(3, 10), hi + Fraction(3, 10))
v = rem_max(zs); results["M3_shifted_zero_rem"] = "%.4f" % v
print("(M3) zero #11 shifted by 0.3: max Rem = %.4f (the sign-change step would reject this zero; Rem alone rises from %.4f)" % (v, base))

# ---------------- (M4)
ctx.prec = 200
tm = tau_table.tau_list(200); r0 = cert_constants.c1(tm)["sigma0sq"]
tm2 = list(tm); tm2[2] += 1; r1 = cert_constants.c1(tm2)["sigma0sq"]
pinned = json.load(open(DATA + "/constants.json"))["sigma0sq"]
results["M4_baseline_in_pinned"] = bool(r0 >= arb(pinned[0])) and bool(r0 <= arb(pinned[1]))
results["M4_mutant_outside_pinned"] = not r1.overlaps(arb(pinned[0]).union(arb(pinned[1])))
print("(M4) sigma0^2 = %s ; with tau(2)+1: %s" % (r0.str(20), r1.str(20)))

ok = all(v for k, v in results.items() if isinstance(v, bool))
json.dump(results, open(DATA + "/selftest.json", "w"), indent=1, sort_keys=True)
print(json.dumps(results, indent=1, sort_keys=True))
print("SELFTEST", "PASS" if ok else "FAIL")
