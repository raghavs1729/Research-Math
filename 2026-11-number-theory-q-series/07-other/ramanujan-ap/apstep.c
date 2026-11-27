/* Stepping pass (Dudek-Platt Lemma 3.3 in table form) over the pass-1 tables.
 *
 *   apstep tabledir K0 K1 classmask tasks_out > log
 *
 * For every table row [a,b) (b = a + 10^6) inside superchunks K0..K1-1 and every
 * class c in classmask: for all real x in [a,b),
 *     pi_c(x)   <= N := #{p < b}            (table, exact)
 *     x/log x   >= a/log a                  (a >= e)
 *     pi_c(x/e) >= M := #{p < j*10^6}       (table, exact; j*10^6*e < a)
 * so   N^2 phi log a < e a M   implies the inequality on the whole row.
 * The comparison is made in long double with a relative guard of 1e-15.
 * Rows where it fails for some class are written to tasks_out ("a b mask",
 * consecutive rows merged) for the exact scan (apscan).
 *
 * stdout:  CUM x n_0..n_16      class counts of primes < x at superchunk ends
 *          S c bin maxrho minrho rows flagged   per class and bin (0.01 in log10 a):
 *               max/min over rows of the row bound N^2 phi log a/(e a M)
 */
#include <math.h>
#include "aptable.h"
#define E_LD 2.718281828459045235360287471352662498L
#define NBIN 1600

typedef struct { const char *dir; long K; uint16_t *buf; uint64_t row; uint64_t res[NRES]; } reader_t;
/* advance reader so that res = residue counts of primes (>=7) < x, x multiple of G, non-decreasing */
static void advance(reader_t *r, uint64_t x) {
    uint64_t target = x / G;
    while (r->row < target) {
        long K = (long)(r->row / ROWS);
        if (K != r->K) { free(r->buf); r->buf = load_file(r->dir, K); r->K = K; }
        const uint16_t *b = r->buf + (r->row % ROWS) * NRES;
        for (int i = 0; i < NRES; i++) r->res[i] += b[i];
        r->row++;
    }
}
int main(int argc, char **argv) {
    if (argc < 6) { fprintf(stderr, "usage: apstep tabledir K0 K1 classmask tasks_out\n"); return 1; }
    const char *dir = argv[1];
    long K0 = atol(argv[2]), K1 = atol(argv[3]);
    uint32_t cmask = (uint32_t)strtoul(argv[4], 0, 0);
    FILE *tf = fopen(argv[5], "w");
    ap_init();
    reader_t lead = {dir, -1, 0, 0, {0}}, lag = {dir, -1, 0, 0, {0}};
    static long double mxr[NCLS][NBIN], mnr[NCLS][NBIN];
    static uint64_t rows[NCLS][NBIN], flg[NCLS][NBIN];
    for (int c = 0; c < NCLS; c++) for (int b = 0; b < NBIN; b++) { mxr[c][b] = 0; mnr[c][b] = 1e30L; }
    uint64_t ta = 0, tb = 0; uint32_t tmask = 0;
    uint64_t nflag_rows = 0, nrows = 0;
    for (uint64_t a = (uint64_t)K0 * SC; a < (uint64_t)K1 * SC; a += G) {
        uint64_t b = a + G;
        uint64_t N[NCLS], M[NCLS];
        advance(&lead, b);
        res_to_cls(lead.res, b, N);
        if (b % SC == 0) {
            printf("CUM %llu", (unsigned long long)b);
            for (int c = 0; c < NCLS; c++) printf(" %llu", (unsigned long long)N[c]);
            printf("\n");
        }
        if (a < 1000000) {  /* first row: a = 0, no test possible; scan it exactly */
            if (tmask) fprintf(tf, "%llu %llu %u\n", (unsigned long long)ta, (unsigned long long)tb, tmask);
            ta = a; tb = b; tmask = cmask; continue;
        }
        uint64_t j = (uint64_t)(((long double)a / E_LD) / G);
        while (j > 0 && floor_e_times(j * G) >= a) j--;
        advance(&lag, j * G);
        res_to_cls(lag.res, j * G, M);
        long double la = logl((long double)a);
        int bin = (int)floorl(100 * log10l((long double)a));
        uint32_t fm = 0;
        for (int c = 0; c < NCLS; c++) if (cmask >> c & 1) {
            long double lhs = (long double)N[c] * (long double)N[c] * CPHI[c] * la;
            long double rhs = E_LD * (long double)a * (long double)M[c];
            long double rho = M[c] ? lhs / rhs : 1e30L;
            if (!(rho < 1 - 1e-15L)) { fm |= 1u << c; flg[c][bin]++; }
            if (rho > mxr[c][bin]) mxr[c][bin] = rho;
            if (rho < mnr[c][bin]) mnr[c][bin] = rho;
            rows[c][bin]++;
        }
        nrows++;
        if (fm) {
            nflag_rows++;
            if (tmask && tb == a && tb - ta < 1000 * G && (b % SC) != G) { tb = b; tmask |= fm; }
            else {
                if (tmask) fprintf(tf, "%llu %llu %u\n", (unsigned long long)ta, (unsigned long long)tb, tmask);
                ta = a; tb = b; tmask = fm;
            }
        }
    }
    if (tmask) fprintf(tf, "%llu %llu %u\n", (unsigned long long)ta, (unsigned long long)tb, tmask);
    fclose(tf);
    for (int c = 0; c < NCLS; c++) if (cmask >> c & 1)
        for (int b = 0; b < NBIN; b++) if (rows[c][b])
            printf("S %d %d %.12Lf %.12Lf %llu %llu\n", c, b, mxr[c][b], mnr[c][b], (unsigned long long)rows[c][b], (unsigned long long)flg[c][b]);
    printf("# rows tested %llu, flagged %llu\n", (unsigned long long)nrows, (unsigned long long)nflag_rows);
    return 0;
}
