/*
 * abel_scan.c -- Chebyshev's 1853 Abel sum on a geometric grid.
 *
 *   C(x)   = sum_{p not| q} chi(p) e^{-p/x}                 (chi = chi_{-4} or chi_{-3})
 *   W(x)   = - sum_p chi(p) log p e^{-p/x}                   (the log-weighted sum, sign of the bias)
 *   M_k^{+}(x), M_k^{-}(x) = sum_{chi(p)=+1 / -1} z^k e^{-z},  z = p/x,  k = 0..K
 * (so C = M_0^+ - M_0^-;  the M_k give the Taylor expansion of C(x e^h) in h, see certify_grid.py).
 *
 * Grid: x_j = 10^{j/M},  j = jlo..jhi.   Primes are sieved up to P = CUT * x_max, and a prime p
 * contributes to x_j only if p <= CUT * x_j  (CUT = 50: the neglected tail of M_0 is
 * < x e^{-50} (1 + 1/x) and is added to the printed error bound by the post-processing script).
 *
 * Arithmetic: e^{-z} is the double-precision libm value (error < 1 ulp claimed by glibc); the products
 * z^k e^{-z} are formed in long double from the double e^{-z} and accumulated in long double
 * (64-bit mantissa).  Nothing in the proofs depends on a sign decided here; see CHEBYSHEV1853.md.
 *
 * Usage: abel_scan q log10xmin log10xmax M  > out
 * Output line:  j  x  C  W  Wsq  M0+ M0- M1+ M1- ... MK+ MK-
 *   Wsq = sum_{p not| q} log p e^{-p^2/x}   (the prime-square sum, for the bias term)
 */
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef long double ld;
#define K 8
#define CUT 50.0

static int q;
static int jlo, jhi, M, NJ;
static double *xg, *invx, *cutp;          /* x_j, 1/x_j, CUT*x_j */
static ld (*Sp)[K+1], (*Sm)[K+1];        /* moments, + and - classes */
static ld *Wl, *Wsq;
static uint64_t nprimes = 0;

static inline void do_prime(uint64_t p)
{
    int s;
    if (q == 4) { if (p == 2) return; s = ((p & 3) == 1) ? 1 : -1; }
    else        { if (p == 3) return; s = ((p % 3) == 1) ? 1 : -1; }
    nprimes++;
    double dp = (double)p, lp = log(dp);
    /* first j with CUT*x_j >= p */
    int j0 = (int)floor((log10(dp / CUT)) * M) - jlo - 1;
    if (j0 < 0) j0 = 0;
    while (j0 < NJ && cutp[j0] < dp) j0++;
    for (int j = j0; j < NJ; j++) {
        double z = dp * invx[j];
        double e = exp(-z);
        ld t = e, zl = z;
        ld *S = (s > 0) ? Sp[j] : Sm[j];
        S[0] += t;
        Wl[j] -= s * (ld)lp * t;
        for (int k = 1; k <= K; k++) { t *= zl; S[k] += t; }
    }
    /* prime squares: p^2 <= CUT x */
    double p2 = dp * dp;
    if (p2 <= cutp[NJ-1]) {
        int j1 = 0;
        while (j1 < NJ && cutp[j1] < p2) j1++;
        for (int j = j1; j < NJ; j++) Wsq[j] += (ld)lp * (ld)exp(-p2 * invx[j]);
    }
}

int main(int argc, char **argv)
{
    if (argc < 5) { fprintf(stderr, "usage: %s q log10xmin log10xmax M\n", argv[0]); return 1; }
    q = atoi(argv[1]);
    double lmin = atof(argv[2]), lmax = atof(argv[3]);
    M = atoi(argv[4]);
    jlo = (int)lround(lmin * M); jhi = (int)lround(lmax * M); NJ = jhi - jlo + 1;
    xg = malloc(NJ * sizeof(double)); invx = malloc(NJ * sizeof(double)); cutp = malloc(NJ * sizeof(double));
    Sp = calloc(NJ, sizeof(*Sp)); Sm = calloc(NJ, sizeof(*Sm));
    Wl = calloc(NJ, sizeof(ld)); Wsq = calloc(NJ, sizeof(ld));
    for (int j = 0; j < NJ; j++) {
        xg[j] = (double)powl(10.0L, (ld)(j + jlo) / (ld)M);
        invx[j] = 1.0 / xg[j]; cutp[j] = CUT * xg[j];
    }
    uint64_t P = (uint64_t)(CUT * xg[NJ-1]) + 2;
    if (P < 100) P = 100;
    /* base primes */
    uint64_t sq = (uint64_t)sqrt((double)P) + 2;
    char *small = calloc(sq + 1, 1);
    uint32_t *bp = malloc((sq / 2 + 16) * sizeof(uint32_t)); size_t nbp = 0;
    for (uint64_t i = 2; i <= sq; i++) {
        if (!small[i]) { bp[nbp++] = (uint32_t)i; for (uint64_t k = i * i; k <= sq; k += i) small[k] = 1; }
    }
    /* segmented sieve over odd numbers */
    const uint64_t SEG = 1u << 19;        /* numbers per segment (odd only: SEG/2 flags) */
    char *flag = malloc(SEG / 2);
    do_prime(2);
    for (uint64_t lo = 3; lo <= P; lo += SEG) {
        uint64_t hi = lo + SEG - 2;       /* inclusive, odd */
        if (hi > P) hi = P;
        size_t n = (hi - lo) / 2 + 1;
        memset(flag, 0, n);
        for (size_t i = 1; i < nbp; i++) {
            uint64_t pp = bp[i];
            if (pp * pp > hi) break;
            uint64_t st = (lo + pp - 1) / pp * pp;
            if (st < pp * pp) st = pp * pp;
            if (!(st & 1)) st += pp;
            for (uint64_t k = (st - lo) / 2; k < n; k += pp) flag[k] = 1;
        }
        for (size_t k = 0; k < n; k++) if (!flag[k]) do_prime(lo + 2 * k);
    }
    printf("# abel_scan q=%d M=%d K=%d CUT=%g P=%llu nprimes(race)=%llu\n", q, M, K, CUT,
           (unsigned long long)P, (unsigned long long)nprimes);
    printf("# j x C W Wsq M0+ M0- ... M%d+ M%d-\n", K, K);
    for (int j = 0; j < NJ; j++) {
        printf("%d %.17g %.20Lg %.20Lg %.20Lg", j + jlo, xg[j], Sp[j][0] - Sm[j][0], Wl[j], Wsq[j]);
        for (int k = 0; k <= K; k++) printf(" %.20Lg %.20Lg", Sp[j][k], Sm[j][k]);
        printf("\n");
    }
    return 0;
}
