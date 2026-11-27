/* Observed range of delta(x) = (phi pi(x;q,a) - li(x)) log x / sqrt(x) at the table
 * points x = k*10^6 (exact prime counts from the pass-1 tables; li in long double
 * via the convergent series li(x) = gamma + log log x + sum_{k>=1} (log x)^k/(k k!)).
 *   apdelta tabledir K0 K1 > out     lines: D class decade_bin min max  (bin = floor(10 log10 x))
 * Descriptive statistics only (long double li, no error analysis): used for the
 * discussion of the gap, not for any theorem. */
#include <math.h>
#include "aptable.h"
static long double li_ld(long double x) {
    long double L = logl(x), term = 1, sum = 0;
    for (int k = 1; k < 400; k++) { term *= L / k; sum += term / k; if (term / k < 1e-22L * sum) break; }
    return 0.577215664901532860606512090082402431L + logl(L) + sum;
}
int main(int argc, char **argv) {
    if (argc < 4) return 1;
    const char *dir = argv[1]; long K0 = atol(argv[2]), K1 = atol(argv[3]);
    ap_init();
    static long double mn[NCLS][200], mx[NCLS][200]; static int used[200];
    for (int c = 0; c < NCLS; c++) for (int b = 0; b < 200; b++) { mn[c][b] = 1e30L; mx[c][b] = -1e30L; }
    uint64_t res[NRES] = {0};
    for (long K = 0; K < K1; K++) {
        uint16_t *buf = load_file(dir, K);
        for (uint64_t r = 0; r < ROWS; r++) {
            for (int i = 0; i < NRES; i++) res[i] += buf[r * NRES + i];
            if (K < K0) continue;
            uint64_t x = (uint64_t)K * SC + (r + 1) * G;      /* counts are of primes < x; x is not prime */
            uint64_t n[NCLS]; res_to_cls(res, x, n);
            long double L = logl((long double)x), sq = sqrtl((long double)x), li = li_ld((long double)x);
            int b = (int)floorl(10 * log10l((long double)x) + 1e-9L);
            used[b] = 1;
            for (int c = 0; c < NCLS; c++) {
                long double d = ((long double)CPHI[c] * n[c] - li) * L / sq;
                if (d < mn[c][b]) mn[c][b] = d;
                if (d > mx[c][b]) mx[c][b] = d;
            }
        }
        free(buf);
    }
    for (int c = 0; c < NCLS; c++) for (int b = 0; b < 200; b++) if (used[b]) printf("D %d %d %.4Lf %.4Lf\n", c, b, mn[c][b], mx[c][b]);
    return 0;
}
