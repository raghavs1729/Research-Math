/*
 * fujii_scan.c -- Fujii's sums  C_alpha(x) = sum_p chi(p) exp(-(p/x)^alpha)  and
 *                 W_alpha(x) = -sum_p chi(p) log p exp(-(p/x)^alpha)   on a geometric grid (double precision).
 *
 * Grid x_j = 10^{j/M}, M = max(50, MPER*alpha) points per decade (the sum varies on the relative scale 1/alpha).
 * A prime p contributes  chi(p) exp(-(p/x_j)^alpha):
 *     (p/x_j)^alpha > 46      : neglected (< 1.1e-20 each)
 *     (p/x_j)^alpha < 1e-18   : counted as exactly chi(p)        (difference array)
 *     otherwise               : evaluated.
 * Floating point only (exploration of the truth for alpha != 1; nothing is certified here).
 *
 * Usage: fujii_scan q log10xmin log10xmax MPER alpha1 [alpha2 ...]
 * Output per alpha: summary lines (extremes of C_alpha, of N = -C_alpha log x/sqrt x and of W_alpha/sqrt x) and
 * a coarse table.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

static int q;
static uint32_t *pr; static size_t npr;

static void sieve(uint64_t P)
{
    char *s = calloc(P + 1, 1);
    npr = 0; pr = malloc(sizeof(uint32_t) * (size_t)(1.3 * P / log((double)P) + 100));
    for (uint64_t i = 2; i <= P; i++) {
        if (!s[i]) { pr[npr++] = (uint32_t)i; for (uint64_t k = i * i; k <= P; k += i) s[k] = 1; }
    }
    free(s);
}

int main(int argc, char **argv)
{
    if (argc < 6) { fprintf(stderr, "usage: %s q log10xmin log10xmax MPER alpha...\n", argv[0]); return 1; }
    q = atoi(argv[1]);
    double lmin = atof(argv[2]), lmax = atof(argv[3]);
    int MPER = atoi(argv[4]);
    double amin = 1e300;
    for (int a = 5; a < argc; a++) if (atof(argv[a]) < amin) amin = atof(argv[a]);
    uint64_t P = (uint64_t)(pow(10.0, lmax) * pow(46.0, 1.0 / amin)) + 10;
    if (P > 2000000000ULL) { fprintf(stderr, "P = %llu too large for the simple sieve (choose larger alpha or smaller x)\n", (unsigned long long)P); return 1; }
    sieve(P);
    for (int a = 5; a < argc; a++) {
        double alpha = atof(argv[a]);
        int M = (int)ceil(MPER * alpha); if (M < 50) M = 50;
        int jlo = (int)floor(lmin * M), jhi = (int)ceil(lmax * M), NJ = jhi - jlo + 1;
        double *lx = malloc(NJ * sizeof(double));
        double *C = calloc(NJ + 2, sizeof(double)), *W = calloc(NJ + 2, sizeof(double));
        double *dC = calloc(NJ + 2, sizeof(double)), *dW = calloc(NJ + 2, sizeof(double));
        for (int j = 0; j < NJ; j++) lx[j] = (double)(j + jlo) / M * log(10.0);
        double hiw = log(46.0) / alpha;        /* log(p/x) <= hiw  */
        double low = log(1e-18) / alpha;       /* log(p/x) >= low  (else weight 1) */
        double step = log(10.0) / M;
        for (size_t i = 0; i < npr; i++) {
            uint32_t p = pr[i]; int s;
            if (q == 4) { if (p == 2) continue; s = ((p & 3) == 1) ? 1 : -1; }
            else        { if (p == 3) continue; s = ((p % 3) == 1) ? 1 : -1; }
            double lp = log((double)p);
            /* j range: lp - lx_j <= hiw  <=>  lx_j >= lp - hiw ;  lp - lx_j >= low <=> lx_j <= lp - low */
            int j0 = (int)ceil((lp - hiw - lx[0]) / step - 1e-9); if (j0 < 0) j0 = 0;
            int j1 = (int)floor((lp - low - lx[0]) / step + 1e-9);
            if (j0 >= NJ) continue;
            if (j1 >= NJ - 1) j1 = NJ - 1; else { dC[j1 + 1] += s; dW[j1 + 1] -= s * lp; }
            for (int j = j0; j <= j1; j++) {
                double w = exp(-exp(alpha * (lp - lx[j])));
                C[j] += s * w; W[j] -= s * lp * w;
            }
        }
        double accC = 0, accW = 0;
        for (int j = 0; j < NJ; j++) { accC += dC[j]; accW += dW[j]; C[j] += accC; W[j] += accW; }
        /* summaries */
        double maxC = -1e300, maxC_x = 0, firstpos = 0, minN = 1e300, minN_x = 0, maxN = -1e300, maxN_x = 0;
        double minW = 1e300, minW_x = 0, firstposW = 0, maxC100 = -1e300, maxC100_x = 0;
        int npos = 0, nties = 0;   /* nties: |C| <= 1e-9 (below double resolution: intervals where pi(x;q,.) are tied and alpha is large) */
        for (int j = 0; j < NJ; j++) {
            double x = exp(lx[j]);
            if (x < 3.0) continue;                 /* below 3 the sum is e^{-(3/x)^alpha}-dominated and tiny */
            if (C[j] > maxC) { maxC = C[j]; maxC_x = x; }
            if (x >= 100 && C[j] > maxC100) { maxC100 = C[j]; maxC100_x = x; }
            if (C[j] > 1e-9) { npos++; if (firstpos == 0) firstpos = x; }
            else if (C[j] > -1e-9) nties++;
            if (W[j] <= 0 && firstposW == 0) firstposW = x;
            if (x >= 100) {
                double N = -C[j] * lx[j] / sqrt(x);
                if (N < minN) { minN = N; minN_x = x; }
                if (N > maxN) { maxN = N; maxN_x = x; }
                double w = W[j] / sqrt(x);
                if (w < minW) { minW = w; minW_x = x; }
            }
        }
        printf("alpha=%g q=%d M=%d grid=[%g,%g] primes<=%llu\n", alpha, q, M, exp(lx[0]), exp(lx[NJ-1]), (unsigned long long)P);
        printf("  max C_alpha on x>=3: %.6f at x=%.6g ;  on x>=100: %.6f at x=%.6g ;  grid points with C_alpha>1e-9: %d (first at x=%.8g) ;  with |C_alpha|<=1e-9: %d\n",
               maxC, maxC_x, maxC100, maxC100_x, npos, firstpos, nties);
        printf("  N=-C log x/sqrt x on x>=100: min %.6f at x=%.6g, max %.6f at x=%.6g ;  min W_alpha/sqrt x on x>=100: %.6f at x=%.6g (first W<=0 at x=%.8g)\n",
               minN, minN_x, maxN, maxN_x, minW, minW_x, firstposW);
        printf("  table x: C_alpha, N, W_alpha/sqrt x:");
        for (int j = 0; j < NJ; j++) {
            double l10 = (double)(j + jlo) / M;
            if (fabs(l10 - round(l10)) < 1e-9 && l10 >= 1)
                printf("  [1e%d: %.4f %.4f %.4f]", (int)round(l10), C[j], -C[j] * lx[j] / sqrt(exp(lx[j])), W[j] / sqrt(exp(lx[j])));
        }
        printf("\n");
        fflush(stdout);
        free(lx); free(C); free(W); free(dC); free(dW);
    }
    return 0;
}
