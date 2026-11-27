"""Numerical check of the transcribed rows: value of a_0 + b_1/(a_1 + b_2/(a_2 + ...)) against the printed Moebius
expression; when they differ, PSLQ for the true relation."""
from mpmath import mp, mpf, e, sqrt, tanh, tan, besselj, pslq, nstr
from rows import ROWS, a_seq
mp.dps = 60
BASE = {"e": e, "e^2": e**2, "sqrt(e)": sqrt(e), "tanh(1/4)": tanh(mpf(1) / 4), "tan(1)": tan(mpf(1)), "J0/J1": besselj(0, 1) / besselj(1, 1)}
def value(init, lfsr, bper, N=600):
    a = a_seq(init, lfsr, N + 1); v = mpf(a[N])
    for j in range(N, 0, -1):
        v = a[j - 1] + mpf(bper[(j - 1) % len(bper)]) / v
    return v
if __name__ == '__main__':
    bad = 0
    for (rid, base, Mx, bper, init, lfsr) in ROWS:
        t = BASE[base]; x = value(init, lfsr, bper); x2 = value(init, lfsr, bper, 900)
        claim = (Mx[0][0] * t + Mx[0][1]) / (Mx[1][0] * t + Mx[1][1])
        ok = abs(x - claim) < mpf(10)**-40
        rel = None if ok else pslq([x * t, x, t, mpf(1)], maxcoeff=10**7, maxsteps=10**5)
        print("%-6s %-10s value %s  stable %s  printed-claim diff %s  %s" % (rid, base, nstr(x, 25), nstr(abs(x - x2), 2), nstr(x - claim, 3), "OK" if ok else "MISMATCH; pslq [x t, x, t, 1] = %s" % rel))
        bad += (not ok)
    print("mismatches:", bad)
