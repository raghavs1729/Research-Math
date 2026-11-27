#!/usr/bin/env python3
"""
w1_cert.py -- certified finite computations for HKMN23.md (ball arithmetic, python-flint/Arb).

For a pair (q, a), a a square mod q, phi = phi(q), R = {b mod q : b^2 = a}, N = #R:

PART A.  W_1(t) = phi * sum_{p<=t, p=a (q)} log p * (t-p)  -  (t-2)^2/2   (t >= 2).
  Let p_1 < p_2 < ... be the primes = a (mod q), theta_i = sum_{j<=i} log p_j,
  Lam_i = sum_{j<=i} p_j log p_j.  On [p_i, p_{i+1}],
      W_1(t) = phi*(theta_i t - Lam_i) - (t-2)^2/2,
  a concave quadratic whose maximum over ALL real t is
      M_i = (phi theta_i)^2/2 + 2 phi theta_i - phi Lam_i        (at t = phi theta_i + 2).
  On [2, p_1], W_1(t) = -(t-2)^2/2 <= 0.  Hence:
      M_i < 0 for every i with p_i <= P   ==>   W_1(t) < 0 for all 2 < t <= P'      (P' = the first
      prime = a (q) beyond P; in particular for 2 < t <= P).
  The script certifies M_i < 0 (Arb: the upper end of the ball is negative) for every such prime
  p_i <= PA, and records, per decade, an upper bound for sup W_1(t)/t^(3/2):
      on [p_i, p_{i+1}]:  W_1(t)/t^(3/2) <= M_i / p_{i+1}^(3/2)     (M_i < 0).

PART B.  Theta(y) = phi * sum_{p<=y, p mod q in R} log p.   Certified upper bound for
      c_1 = sup_{0 < y <= Y1} (N y - Theta(y)) / sqrt(y).
  Between consecutive primes r < r' of R the function (N y - Theta(r))/sqrt(y) is increasing
  (Theta >= 0), so the supremum is the maximum of (N r' - Theta(r'-))/sqrt(r') over the primes
  r' <= Y1 of R, and of the value at y = Y1; on (0, r_1) it is N sqrt(r_1).

Usage:  python3 w1_cert.py q a PA Y1  [prec]     -> one JSON object on stdout (plus log on stderr)
"""
import sys, json, math, time
import numpy as np
from flint import arb, ctx


def primes_upto(n):
    """odd-only numpy sieve; returns numpy array of primes <= n."""
    if n < 2:
        return np.array([], dtype=np.int64)
    size = (n - 1) // 2          # index i <-> 2i+3
    s = np.ones(size, dtype=bool)
    lim = int(math.isqrt(n))
    for i in range((lim - 1) // 2):
        if s[i]:
            p = 2 * i + 3
            s[(p * p - 3) // 2::p] = False
    return np.concatenate((np.array([2], dtype=np.int64), 2 * np.nonzero(s)[0].astype(np.int64) + 3))


def main():
    q, a, PA, Y1 = int(sys.argv[1]), int(sys.argv[2]), int(float(sys.argv[3])), int(float(sys.argv[4]))
    prec = int(sys.argv[5]) if len(sys.argv) > 5 else 128
    ctx.prec = prec
    a %= q
    assert math.gcd(a, q) == 1
    phi = sum(1 for k in range(1, q + 1) if math.gcd(k, q) == 1)
    R = [b for b in range(q) if math.gcd(b, q) == 1 and (b * b - a) % q == 0]
    N = len(R)
    assert N > 0, "a must be a square mod q"
    t0 = time.time()
    allp = primes_upto(max(PA, Y1) + 1000)
    res = {"q": q, "a": a, "phi": phi, "N": N, "R": R, "PA": PA, "Y1": Y1, "prec": prec}

    # ---------------- PART A ----------------
    pa = allp[(allp % q) == a] if q > 1 else allp
    pa_in = pa[pa <= PA]
    nxt = pa[pa > PA]
    assert len(nxt) > 0
    plist = [int(x) for x in pa_in] + [int(nxt[0])]
    theta = arb(0)
    Lam = arb(0)
    worst = None            # largest (least negative) M_i / p_{i+1}^{3/2}
    dec = {}                # decade k -> sup of M_i/p_{i+1}^{3/2} over pieces with p_i in [10^k, 10^(k+1))
    nfail = 0
    first = []
    for i in range(len(plist) - 1):
        p, pn = plist[i], plist[i + 1]
        lp = arb(p).log()
        theta += lp
        Lam += lp * p
        pt = phi * theta
        M = pt * pt / 2 + 2 * pt - phi * Lam
        if not (M < 0):
            nfail += 1
            sys.stderr.write(f"FAIL or undecided at p={p}: M={M}\n")
            continue
        g = M / (arb(pn) ** 1.5)
        gu = float(g.upper())
        k = len(str(p)) - 1
        if k not in dec or gu > dec[k][0]:
            dec[k] = (gu, p)
        if worst is None or gu > worst[0]:
            worst = (gu, p)
        if i < 6:
            first.append({"p": p, "M": M.str(12, radius=False), "t_star": (pt + 2).str(10, radius=False)})
    res["A"] = {"num_primes": len(plist) - 1, "first_prime": plist[0] if plist else None,
                "last_prime_checked": plist[-2] if len(plist) > 1 else None, "valid_up_to": plist[-1],
                "failures": nfail, "sup_W1_over_t32_upper": worst[0], "sup_at_p": worst[1],
                "per_decade_sup": {str(k): {"sup": v[0], "p": v[1]} for k, v in sorted(dec.items())},
                "first_pieces": first,
                "theta_final": theta.str(25, radius=False), "theta_final_rad": theta.rad().str(5, radius=False)}
    sys.stderr.write(f"[{q},{a}] part A done: {len(plist)-1} primes, failures {nfail}, {time.time()-t0:.1f}s\n")

    # ---------------- PART B ----------------
    Rset = np.zeros(q, dtype=bool)
    Rset[R] = True
    pr = allp[Rset[allp % q]] if q > 1 else allp
    pr = [int(x) for x in pr[pr <= Y1]]
    Theta = arb(0)
    c1 = arb(0)
    c1_at = 0
    for r in pr:
        v = (N * arb(r) - Theta) / arb(r).sqrt()      # left limit at r
        if float(v.upper()) > float(c1.upper()):
            c1, c1_at = v, r
        Theta += phi * arb(r).log()
    v = (N * arb(Y1) - Theta) / arb(Y1).sqrt()
    if float(v.upper()) > float(c1.upper()):
        c1, c1_at = v, Y1
    res["B"] = {"num_primes": len(pr), "c1_upper": float(c1.upper()) * (1 + 1e-15), "c1_ball": c1.str(20), "attained_before": c1_at,
                "final_deficit_over_sqrt": v.str(12, radius=False)}
    res["seconds"] = round(time.time() - t0, 1)
    print(json.dumps(res))


if __name__ == "__main__":
    main()
