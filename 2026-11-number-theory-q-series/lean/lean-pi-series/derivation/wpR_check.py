"""Numerical sanity check for WP-R (N2-PLAN.md): log_sigma_window and prod_sigma_window.

For each N we print the smallest constant C that makes the two statements true on the window
(for the given eps, delta).  The statements are plausible iff these C stay bounded as N grows.
"""
from mpmath import mp, mpf, sqrt, log, floor

mp.dps = 40
xs = sqrt(3) / 2


def sigma(n, k):
    n, k = mpf(n), mpf(k)
    num = (2*k+3) * (n-k-1)**2 * (8*k+4*n+1) * (8*k+4*n+3) * (8*k+4*n+5) * (8*k+4*n+7)
    den = (2*k+1) * (n+k+1)**2 * (8*k-4*n+9) * (8*k-4*n+11) * (8*k-4*n+13) * (8*k-4*n+15)
    return num / den


def check(N, eps, delta):
    # statement 1: |log sigma + 24 d/N| <= eps |d|/N + C/N on |d| <= delta N, d = k - x* N
    C1 = mpf(0)
    lo, hi = int(floor((xs - delta) * N)) + 1, int(floor((xs + delta) * N))
    for k in range(lo, hi + 1):
        d = k - xs * N
        s = sigma(N, k)
        assert s > 0
        C1 = max(C1, N * abs(log(s) + 24 * d / N) - eps * abs(d))
    # statement 2: products over Ico k0 k / Ico k k0, j = k - k0, |j| <= delta N
    k0 = int(floor(xs * N))
    J = int(floor(delta * N))
    C2 = mpf(0)
    L = mpf(0)
    for j in range(1, J + 1):          # k = k0 + j, product over i in [k0, k)
        L += log(sigma(N, k0 + j - 1))
        # need  -(12+eps) j^2/N - C/N <= L <= -(12-eps) j^2/N + C/N
        C2 = max(C2, N * abs(L + 12 * mpf(j)**2 / N) - eps * j**2)
    L = mpf(0)
    for j in range(1, J + 1):          # k = k0 - j, product over i in [k, k0)
        L += log(sigma(N, k0 - j))
        # need  (12-eps) j^2/N - C/N <= L <= (12+eps) j^2/N + C/N
        C2 = max(C2, N * abs(L - 12 * mpf(j)**2 / N) - eps * j**2)
    return C1, C2


if __name__ == "__main__":
    for eps, delta in [(mpf(1), mpf("0.004")), (mpf("0.2"), mpf("0.0008")), (mpf("0.5"), mpf("0.002"))]:
        for N in (200, 2000, 20000):
            C1, C2 = check(N, eps, delta)
            print(f"eps={float(eps)} delta={float(delta)} N={N}: C(log_sigma_window)>={float(C1):.4f}"
                  f"  C(prod_sigma_window)>={float(C2):.4f}")
