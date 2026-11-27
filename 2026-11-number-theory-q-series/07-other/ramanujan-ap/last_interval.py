#!/usr/bin/env python3
"""Exact end of the last failure interval, in ball arithmetic (python-flint / Arb).

Input (argv): q a p n m   where p is the last failing prime of the class found by
apscan, n = pi(p;q,a), m = pi(p/e;q,a).  We
 (1) take n and m as given (they are NOT recomputed here; apscan checks its running
     counts against the pass-1 tables, and verify_counts_lucy.py recomputes them
     independently for q = 3, 4);
 (2) certify G(p) >= 0, i.e. n^2 phi log p >= e p m;
 (2') walk BACKWARDS from p to the start of the component containing p: the component
     starts at the largest prime p0 <= p of the class with G(p0-) < 0, i.e.
     (n0 - 1)^2 phi log p0 < e p0 m0 (p need not be a component start);
 (3) follow the failure interval to its end: with l_1 < l_2 < ... the primes of the
     class beyond p/e (listed by primesieve), G is continuous and strictly
     decreasing on [p, e l_1), jumps down at e l_1, etc.  The interval ends either
     at the root s of  x/log x = phi n^2/(e m')  (closed end, G(s) = 0 is a failure)
     or at a jump e l_i (open end);
 (4) certify that the next prime of the class after p lies beyond the end.
Usage: last_interval.py q a p n m      (primes near p and p/e are listed with the
primesieve command-line program)
"""
import sys, subprocess
from flint import arb, ctx
ctx.prec = 200
q, a, p, n, m = map(int, sys.argv[1:6])
def primes_in(lo, hi):
    out = subprocess.run(["primesieve", str(lo), str(hi), "-p"], capture_output=True, text=True, check=True).stdout
    return [int(t) for t in out.split()]
def in_cls(x): return q == 1 or x % q == a
nextp = [x for x in primes_in(p + 1, p + 100000) if in_cls(x)][0]
y = int((arb(p) / arb(1).exp()).floor().unique_fmpz())          # floor(p/e), exact
lags = [x for x in primes_in(y + 1, y + 200000) if in_cls(x)][:200]
phi = {1: 1, 3: 2, 4: 2, 5: 4, 8: 4, 12: 4}[q]
E = arb(1).exp()
def G(x, n, m):            # x arb
    return arb(n) ** 2 - E * x / (phi * x.log()) * m
def root(n, m):            # solve x/log x = phi n^2/(e m) for x > e, enclosure by bisection + monotonicity
    target = arb(phi) * n * n / (E * m)
    lo, hi = arb(p), arb(p) * 2
    assert (lo / lo.log() - target) <= 0 and (hi / hi.log() - target) > 0
    for _ in range(190):
        mid = (lo + hi) / 2
        mid = arb(mid.mid())
        v = mid / mid.log() - target
        if v > 0: hi = mid
        elif v < 0: lo = mid
        else: break
    return lo, hi
for l in lags:
    assert l % q == a % q or q == 1
assert all(E * l > p for l in lags), "lag primes must be beyond p/e"
g0 = G(arb(p), n, m)
print(f"class ({q},{a}): p* = {p}, n = {n}, m = {m}")
print("G(p*) =", g0.str(25), " certified >= 0:", bool(g0 >= 0))
assert g0 > 0
# (2') start of the component containing p
prev_lead = [x for x in primes_in(max(2, p - 2000000), p - 1) if in_cls(x)][::-1]
prev_lag = [x for x in primes_in(max(2, y - 1000000), y) if in_cls(x)][::-1]      # all have e*l < p
cp, cn, cm = p, n, m
li_ = 0
for pp in prev_lead + [None]:
    v = arb(cn - 1) ** 2 * phi * arb(cp).log() - E * cp * cm
    assert v > 0 or v < 0, "undecided start"
    if v < 0:
        break
    if pp is None: raise SystemExit("component start not found in the window; enlarge it")
    # previous prime pp of the class: lag primes with e*l > pp are not yet counted there
    while li_ < len(prev_lag) and E * prev_lag[li_] > pp:
        li_ += 1; cm -= 1
    assert li_ < len(prev_lag), "lag window too short"
    cp, cn = pp, cn - 1
    assert G(arb(cp), cn, cm) >= 0
print(f"last component starts at the prime {cp} (n = {cn}, m = {cm}); failing primes in it up to p*: {n - cn + 1}")
comp_start = cp
mc = m
for l in lags + [None]:
    lo, hi = root(n, mc)
    if l is None:
        raise SystemExit("ran out of lag primes; supply more")
    jump = E * l
    if hi < jump:
        print(f"end of failure set: root s in [{lo.str(30)}, {hi.str(30)}] (closed end), m' = {mc}")
        end_lo, end_hi, kind = lo, hi, "root"
        break
    assert lo > jump, "undecided"
    # the failure persists up to e*l, where m increases by one: sign of G just after the jump
    g = G(jump, n, mc + 1)
    if g < 0:
        print(f"end of failure set: jump at e*{l} = {jump.str(30)} (open end)")
        end_lo, end_hi, kind = jump, jump, "jump"
        break
    assert g > 0, "undecided"
    mc += 1
fl = int(end_lo.floor().unique_fmpz()) if kind == "root" else int(jump.floor().unique_fmpz())
assert int(end_hi.floor().unique_fmpz()) == fl
print("largest integer counterexample:", fl)
print("next prime of the class:", nextp, " beyond the end:", bool(arb(nextp) > end_hi))
assert arb(nextp) > end_hi
# integers fl and fl+1
def Gint(x):
    cnt_m = m + sum(1 for l in lags if (E * l < x))
    return G(arb(x), n, cnt_m)
print(f"G({fl}) =", Gint(fl).str(20), f"  G({fl+1}) =", Gint(fl + 1).str(20))
assert Gint(fl) >= 0 and Gint(fl + 1) < 0
