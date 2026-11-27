/* Bias-free variant (experiment for RAMANUJAN-AP.md, section on Chebyshev's bias).
 *
 *   apscan_star tabledir K0 K1 Kmax > out
 *
 * Same scan as apscan, but for Riemann's weighted prime-power count in a class,
 *     Pi*(x;q,a) = sum_{p^k <= x, p^k = a mod q} 1/k,
 * i.e. the inequality  Pi*(x)^2 < (e x/(phi log x)) Pi*(x/e).  Under GRH,
 * Pi*(x;q,a) - li(x)/phi has no bias term of size sqrt(x)/log x: the squares of
 * primes, which cause Chebyshev's bias, are counted with weight 1/2.
 * Classes: q = 1, 3, 4 (indices 0..4).  Tasks: superchunks K0..K1-1.
 * G* can jump up only at primes and prime powers of the class; both are tested.
 * Arithmetic as in apscan (two-sided double screen with 1e-12 margins, long double
 * otherwise, near-ties |rho-1| < 1e-15 reported as "A"); the weighted counts are
 * sums of at most 10^6 terms 1/k in long double (error < 1e-13 absolute, far below
 * the margins since counts exceed 10^8 where it matters).  This is an experiment,
 * not an ingredient of any theorem.
 * Output: B a b c bin nevents nfail ; F a b c first last total
 */
#include <primesieve.h>
#include <omp.h>
#include <math.h>
#include "aptable.h"
#define E_LD 2.718281828459045235360287471352662498L
#define NBIN 1600
#define NC 5
static uint64_t EDGE[NBIN + 1];
typedef struct { uint64_t pos, eflo; int k; uint64_t l; } pp_t;
static pp_t *PP; static size_t NPP;
static int cmp_pp(const void *x, const void *y) { uint64_t a = ((const pp_t*)x)->pos, b = ((const pp_t*)y)->pos; return a < b ? -1 : a > b; }
static inline int in_cls_pp(uint64_t pos, uint64_t l, int c) {   /* l^k = pos in class c ?  (l must not divide q) */
    int q = CQ[c];
    if (q == 1) return 1;
    if (q % l == 0) return 0;
    return (int)(pos % q) == CA[c];
}
int main(int argc, char **argv) {
    if (argc < 5) { fprintf(stderr, "usage: apscan_star tabledir K0 K1 Kmax\n"); return 1; }
    const char *dir = argv[1];
    long K0 = atol(argv[2]), K1 = atol(argv[3]), Kmax = atol(argv[4]);
    ap_init();
    for (int j = 0; j <= NBIN; j++) {
        if (j % 100 == 0) { uint64_t v = 1; for (int k = 0; k < j / 100; k++) v *= 10; EDGE[j] = v; }
        else EDGE[j] = (uint64_t)ceill(powl(10.0L, j / 100.0L));
    }
    /* prime powers l^k <= X, k >= 2 */
    uint64_t X = (uint64_t)K1 * SC;
    {
        size_t cap = 1 << 20; PP = malloc(cap * sizeof *PP); NPP = 0;
        primesieve_iterator it; primesieve_init(&it);
        uint64_t l;
        while ((l = primesieve_next_prime(&it)), l * l <= X) {
            unsigned __int128 v = (unsigned __int128)l * l; int k = 2;
            while (v <= X) {
                if (NPP == cap) { cap *= 2; PP = realloc(PP, cap * sizeof *PP); }
                PP[NPP].pos = (uint64_t)v; PP[NPP].k = k; PP[NPP].l = l; PP[NPP].eflo = floor_e_times((uint64_t)v); NPP++;
                v *= l; k++;
            }
        }
        primesieve_free_iterator(&it);
        qsort(PP, NPP, sizeof *PP, cmp_pp);
        fprintf(stderr, "prime powers (k >= 2) up to %llu: %zu\n", (unsigned long long)X, NPP);
    }
    uint64_t (*cum)[NRES] = load_cum(dir, Kmax);
    double t0 = omp_get_wtime();
    #pragma omp parallel for schedule(dynamic, 1)
    for (long K = K0; K < K1; K++) {
        uint64_t a = (uint64_t)K * SC, b = a + SC;
        size_t ocap = 1 << 16, olen = 0; char *out = malloc(ocap);
        #define EMIT(...) do { if (olen + 512 > ocap) { ocap *= 2; out = realloc(out, ocap); } \
                               olen += snprintf(out + olen, 512, __VA_ARGS__); } while (0)
        uint64_t res[NRES], n[NCLS], m[NCLS];
        res_below(dir, cum, a, res); res_to_cls(res, a < 7 ? 7 : a, n);
        uint64_t j = (uint64_t)(((long double)a / E_LD) / G);
        while (j > 0 && floor_e_times(j * G) >= a) j--;
        res_below(dir, cum, j * G, res); res_to_cls(res, j * G, m);
        /* prime-power corrections: lead = powers < max(a,7); lag = powers < j*G */
        long double cl[NC] = {0}, cm[NC] = {0};
        size_t ip = 0, jp = 0;
        uint64_t a7 = a < 7 ? 7 : a;
        while (ip < NPP && PP[ip].pos < a7) { for (int c = 0; c < NC; c++) if (in_cls_pp(PP[ip].pos, PP[ip].l, c)) cl[c] += 1.0L / PP[ip].k; ip++; }
        while (jp < NPP && PP[jp].pos < j * G) { for (int c = 0; c < NC; c++) if (in_cls_pp(PP[jp].pos, PP[jp].l, c)) cm[c] += 1.0L / PP[jp].k; jp++; }
        primesieve_iterator lead, lag; primesieve_init(&lead); primesieve_init(&lag);
        primesieve_jump_to(&lead, a7, b + 1000);
        primesieve_jump_to(&lag, j * G, (uint64_t)(b / 2.7) + 1000);
        uint64_t lq = primesieve_next_prime(&lag), lagT = floor_e_times(lq);
        uint64_t p = primesieve_next_prime(&lead);
        int bin = 0; while (EDGE[bin + 1] <= a7) bin++;
        uint64_t next_edge = EDGE[bin + 1];
        uint64_t nev[NC] = {0}, nfail[NC] = {0}, ftot[NC] = {0}, ffirst[NC], flast[NC];
        double Khi[NC], Klo[NC]; uint64_t win_hi = 0;
        while (1) {
            /* next event: prime p or prime power PP[ip], whichever is smaller (they never coincide) */
            int is_pp = (ip < NPP && PP[ip].pos < p);
            uint64_t y = is_pp ? PP[ip].pos : p;
            if (y >= b) break;
            if (y >= next_edge) {
                for (int c = 0; c < NC; c++) { EMIT("B %llu %llu %d %d %llu %llu\n", (unsigned long long)a, (unsigned long long)b, c, bin, (unsigned long long)nev[c], (unsigned long long)nfail[c]); nev[c] = nfail[c] = 0; }
                while (EDGE[bin + 1] <= y) bin++;
                next_edge = EDGE[bin + 1];
            }
            if (y > win_hi) {
                win_hi = y + (y >> 24) + 1;
                double Lhi = (double)logl((long double)win_hi) * (1 + 1e-15), Llo = (double)logl((long double)y) * (1 - 1e-15);
                for (int c = 0; c < NC; c++) { Khi[c] = (double)(CPHI[c] * (long double)Lhi / E_LD) * (1 + 1e-12); Klo[c] = (double)(CPHI[c] * (long double)Llo / E_LD) * (1 - 1e-12); }
            }
            /* lag: primes l and prime powers with e*(.) < y */
            while (lagT < y) {
                if (lq >= 7) { const int *cc = CLS_OF_RES[RIDX[lq % 120]]; for (int k = 0; k < 6; k++) m[cc[k]]++; }
                else { for (int c = 0; c < NCLS; c++) if (CQ[c] % lq != 0 && in_class(lq, c)) m[c]++; }
                lq = primesieve_next_prime(&lag); lagT = floor_e_times(lq);
            }
            while (jp < NPP && PP[jp].eflo < y) { for (int c = 0; c < NC; c++) if (in_cls_pp(PP[jp].pos, PP[jp].l, c)) cm[c] += 1.0L / PP[jp].k; jp++; }
            int touched[NC];
            if (is_pp) { for (int c = 0; c < NC; c++) { touched[c] = in_cls_pp(y, PP[ip].l, c); if (touched[c]) cl[c] += 1.0L / PP[ip].k; } ip++; }
            else { const int *cc = CLS_OF_RES[RIDX[p % 120]]; for (int c = 0; c < NC; c++) touched[c] = 0;
                   for (int k = 0; k < 6; k++) { n[cc[k]]++; if (cc[k] < NC) touched[cc[k]] = 1; } }
            for (int c = 0; c < NC; c++) if (touched[c]) {
                nev[c]++;
                long double N = (long double)n[c] + cl[c], Mm = (long double)m[c] + cm[c];
                double xx = (double)N; xx *= xx; double yy = (double)y * (double)Mm;
                int fail = 0;
                if (xx * Khi[c] >= yy) {
                    if (xx * Klo[c] >= yy) fail = 1;
                    else {
                        long double rho = Mm > 0 ? N * N * CPHI[c] * logl((long double)y) / (E_LD * (long double)y * Mm) : 2;
                        if (fabsl(rho - 1) < 1e-15L) EMIT("A %llu %d %.6Lf %.6Lf %d\n", (unsigned long long)y, c, N, Mm, (int)(rho >= 1));
                        fail = rho >= 1;
                    }
                }
                if (fail) { nfail[c]++; if (!ftot[c]) ffirst[c] = y; flast[c] = y; ftot[c]++; }
            }
            if (!is_pp) p = primesieve_next_prime(&lead);
        }
        for (int c = 0; c < NC; c++) {
            EMIT("B %llu %llu %d %d %llu %llu\n", (unsigned long long)a, (unsigned long long)b, c, bin, (unsigned long long)nev[c], (unsigned long long)nfail[c]);
            if (ftot[c]) EMIT("F %llu %llu %d %llu %llu %llu\n", (unsigned long long)a, (unsigned long long)b, c, (unsigned long long)ffirst[c], (unsigned long long)flast[c], (unsigned long long)ftot[c]);
        }
        primesieve_free_iterator(&lead); primesieve_free_iterator(&lag);
        #pragma omp critical
        { fwrite(out, 1, olen, stdout); fflush(stdout); }
        free(out);
    }
    fprintf(stderr, "apscan_star: K=[%ld,%ld) done in %.1f s\n", K0, K1, omp_get_wtime() - t0);
    printf("# done\n");
    return 0;
}
