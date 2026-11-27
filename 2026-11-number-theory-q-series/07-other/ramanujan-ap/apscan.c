/* Pass 2: exact scan, prime by prime, of Ramanujan's inequality in residue classes.
 *
 *   apscan tabledir tasks.txt Kmax  > out
 *
 * tasks.txt: lines "a b mask" (a < b multiples of 10^6, mask = bitmask of classes).
 * Kmax: number of superchunk tables available (primes < Kmax*10^10 tabulated).
 *
 * For a class c = (q,a0) write n = pi(x;q,a0), m = pi(x/e;q,a0), phi = phi(q).
 * G_c(x) = n^2 - (e x/(phi log x)) m.  Between consecutive primes of the class G_c
 * is non-increasing in x (x/log x increases for x >= e, and m is non-decreasing),
 * so the inequality can only start to fail at a prime p of the class.  We test
 *      rho_c(p) = n^2 phi log p / (e p m)  >=  1     (failure at x = p)
 * at EVERY prime p >= 7 of the class in [a,b).
 *
 * Arithmetic.  pi(p/e) counts lag primes l with e*l <= p, i.e. floor(e*l) < p,
 * computed exactly in integers (floor_e_times).  The test itself is done in
 * two stages.  (1) Screens in double precision, each conservative by a factor
 * 1 +- 1e-12 (the rounding errors of the few operations involved are < 1e-15),
 * using Llo <= log p <= Lhi on a short window: "n^2 Khi < p m" proves rho < 1,
 * "n^2 Klo >= p m" proves rho > 1, "n^2 Kmx < p m" proves rho < current bin maximum.
 * (2) Whatever the screens do not decide is evaluated in x87 long double
 * (64-bit mantissa; the accumulated relative error of the dozen operations,
 * including logl, is below 1e-18).  Cases with |rho-1| < 1e-17 are reported as
 * "A" records (p, class, n, m, decision taken, rho-1) and must be re-decided in
 * ball arithmetic by recheck_amb.py.
 *
 * Output records (one line each):
 *   B a b c bin nprimes nfail nfail_m0 maxrho p n m nstart   per (task, class, bin of width 0.01 in log10 p);
 *       nfail = primes p of the class with G(p) >= 0; nstart = those among them with G(p-) < 0, i.e. the
 *       number of connected components of the failure set {x : G(x) >= 0} that begin in the bin
 *   F a b c first_p n m last_p n m total_fail          per (task, class) with failures
 *   E a b c n_end                                       class count of primes < b (checked against table;
 *       the lag counts are checked against the table too, at the first table point above the lag position)
 *   A p c n m decision rho-1 [start]                    near-tie, to be re-decided by recheck_amb.py
 */
#include <primesieve.h>
#include <omp.h>
#include <math.h>
#include "aptable.h"

#define E_LD 2.718281828459045235360287471352662498L
#define NBIN 1600
static uint64_t EDGE[NBIN + 1];

typedef struct { uint64_t a, b; uint32_t mask; } task_t;

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: apscan tabledir tasks.txt Kmax\n"); return 1; }
    const char *dir = argv[1];
    long Kmax = atol(argv[3]);
    ap_init();
    for (int j = 0; j <= NBIN; j++) {
        if (j % 100 == 0) { uint64_t v = 1; for (int k = 0; k < j / 100; k++) v *= 10; EDGE[j] = v; }
        else EDGE[j] = (uint64_t)ceill(powl(10.0L, j / 100.0L));
    }
    FILE *tf = fopen(argv[2], "r");
    if (!tf) { perror("tasks"); return 1; }
    size_t nt = 0, cap = 1024;
    task_t *T = malloc(cap * sizeof *T);
    unsigned long long ta, tb; unsigned tm;
    while (fscanf(tf, "%llu %llu %u", &ta, &tb, &tm) == 3) {
        if (nt == cap) { cap *= 2; T = realloc(T, cap * sizeof *T); }
        T[nt].a = ta; T[nt].b = tb; T[nt].mask = tm; nt++;
    }
    fclose(tf);
    uint64_t (*cum)[NRES] = load_cum(dir, Kmax);
    double t0 = omp_get_wtime();
    long ndone = 0;
    unsigned long long tot_prec = 0, tot_amb = 0, tot_lagchk = 0;

    #pragma omp parallel for schedule(dynamic, 1) reduction(+:tot_prec,tot_amb,tot_lagchk)
    for (size_t ti = 0; ti < nt; ti++) {
        uint64_t a = T[ti].a, b = T[ti].b; uint32_t mask = T[ti].mask;
        size_t ocap = 1 << 16, olen = 0;
        char *out = malloc(ocap);
        #define EMIT(...) do { if (olen + 512 > ocap) { ocap *= 2; out = realloc(out, ocap); } \
                               olen += snprintf(out + olen, 512, __VA_ARGS__); } while (0)
        uint64_t res[NRES], n[NCLS], m[NCLS];
        /* lead base: primes < a (with 2,3,5 if a <= 7 we start the scan at 7) */
        res_below(dir, cum, a, res);
        res_to_cls(res, a < 7 ? 7 : a, n);
        /* lag base: j*G with e*j*G < a */
        uint64_t j = (uint64_t)(((long double)a / E_LD) / G);
        while (j > 0 && floor_e_times(j * G) >= a) j--;
        res_below(dir, cum, j * G, res);
        res_to_cls(res, j * G, m);            /* primes < j*G; includes 2,3,5 iff below j*G */
        primesieve_iterator lead, lag;
        primesieve_init(&lead); primesieve_init(&lag);
        primesieve_jump_to(&lead, a < 7 ? 7 : a, b + 1000);
        primesieve_jump_to(&lag, j * G, (uint64_t)(b / 2.7) + 1000);
        uint64_t lq = primesieve_next_prime(&lag);
        uint64_t lagT = floor_e_times(lq);
        uint64_t p = primesieve_next_prime(&lead);
        /* per-bin state */
        int bin = 0;
        while (EDGE[bin + 1] <= p) bin++;
        uint64_t next_edge = EDGE[bin + 1];
        uint64_t n_bin0[NCLS], nfail[NCLS], nfail0[NCLS], nstart[NCLS], mx_p[NCLS], mx_n[NCLS], mx_m[NCLS];
        long double mx[NCLS];
        double Kmx[NCLS], Khi[NCLS], Klo[NCLS];
        uint64_t f_first[NCLS][3], f_last[NCLS][3], f_tot[NCLS];
        for (int c = 0; c < NCLS; c++) { n_bin0[c] = n[c]; nfail[c] = nfail0[c] = nstart[c] = 0; mx[c] = 0; mx_p[c] = mx_n[c] = mx_m[c] = 0;
                                         f_tot[c] = 0; Kmx[c] = INFINITY; Khi[c] = INFINITY; Klo[c] = 0; }
        uint64_t win_hi = 0; double Lhi = 0;
        unsigned long long nprec = 0, namb = 0, nlagchk = 0;
        while (p < b) {
            if (p >= next_edge) {
                for (int c = 0; c < NCLS; c++) if (mask >> c & 1) {
                    EMIT("B %llu %llu %d %d %llu %llu %llu %.18Lg %llu %llu %llu %llu\n", (unsigned long long)a, (unsigned long long)b, c, bin,
                         (unsigned long long)(n[c] - n_bin0[c]), (unsigned long long)nfail[c], (unsigned long long)nfail0[c], mx[c],
                         (unsigned long long)mx_p[c], (unsigned long long)mx_n[c], (unsigned long long)mx_m[c], (unsigned long long)nstart[c]);
                    n_bin0[c] = n[c]; nfail[c] = nfail0[c] = nstart[c] = 0; mx[c] = 0; mx_p[c] = mx_n[c] = mx_m[c] = 0; Kmx[c] = INFINITY;
                }
                while (EDGE[bin + 1] <= p) bin++;
                next_edge = EDGE[bin + 1];
            }
            if (p > win_hi) {
                /* window [p, win_hi]: Llo <= log x <= Lhi for every prime x in it */
                win_hi = p + (p >> 24) + 1;
                Lhi = (double)logl((long double)win_hi) * (1 + 1e-15);
                double Llo = (double)logl((long double)p) * (1 - 1e-15);
                for (int c = 0; c < NCLS; c++) {
                    Kmx[c] = mx[c] > 0 ? (double)(CPHI[c] * (long double)Lhi / (E_LD * mx[c])) * (1 + 1e-12) : INFINITY;
                    Khi[c] = (double)(CPHI[c] * (long double)Lhi / E_LD) * (1 + 1e-12);
                    Klo[c] = (double)(CPHI[c] * (long double)Llo / E_LD) * (1 - 1e-12);
                }
            }
            /* lag: count lag primes l with e*l < p */
            while (lagT < p) {
                if (lq >= 7) { const int *cl = CLS_OF_RES[RIDX[lq % 120]]; for (int k = 0; k < 6; k++) m[cl[k]]++; }
                else { for (int c = 0; c < NCLS; c++) if (CQ[c] % lq != 0 && in_class(lq, c)) m[c]++; }
                lq = primesieve_next_prime(&lag);
                lagT = floor_e_times(lq);
            }
            const int *cl = CLS_OF_RES[RIDX[p % 120]];
            double pd = (double)p;
            for (int k = 0; k < 6; k++) {
                int c = cl[k];
                n[c]++;
                if (!(mask >> c & 1)) continue;
                double xx = (double)n[c]; xx *= xx;
                double yy = pd * (double)m[c];
                int fail = 0;
                if (xx * Kmx[c] >= yy) {          /* possible new bin maximum of rho: evaluate */
                    nprec++;
                    if (m[c] == 0) { fail = 1; nfail0[c]++; }
                    else {
                        long double rho = (long double)n[c] * (long double)n[c] * CPHI[c] * logl((long double)p)
                                          / (E_LD * (long double)p * (long double)m[c]);
                        if (fabsl(rho - 1) < 1e-17L) { namb++; EMIT("A %llu %d %llu %llu %d %.3Le\n", (unsigned long long)p, c, (unsigned long long)n[c], (unsigned long long)m[c], (int)(rho >= 1), rho - 1); }
                        fail = rho >= 1;
                        if (rho > mx[c]) {
                            mx[c] = rho; mx_p[c] = p; mx_n[c] = n[c]; mx_m[c] = m[c];
                            Kmx[c] = (double)(CPHI[c] * (long double)Lhi / (E_LD * rho)) * (1 + 1e-12);
                        }
                    }
                } else if (xx * Khi[c] >= yy) {    /* rho >= 1 not excluded by the screen */
                    if (xx * Klo[c] >= yy) fail = 1;   /* certainly rho > 1 */
                    else {
                        nprec++;
                        long double rho = (long double)n[c] * (long double)n[c] * CPHI[c] * logl((long double)p)
                                          / (E_LD * (long double)p * (long double)m[c]);
                        if (fabsl(rho - 1) < 1e-17L) { namb++; EMIT("A %llu %d %llu %llu %d %.3Le\n", (unsigned long long)p, c, (unsigned long long)n[c], (unsigned long long)m[c], (int)(rho >= 1), rho - 1); }
                        fail = rho >= 1;
                    }
                }
                if (fail) {
                    /* does a new connected component of the failure set start at p?
                     * yes iff G(p-) < 0, i.e. (n-1)^2 phi log p < e p m  (m is the same at p-) */
                    int start = 0;
                    if (m[c] != 0) {
                        double x1 = (double)(n[c] - 1); x1 *= x1;
                        if (x1 * Khi[c] < yy) start = 1;
                        else if (x1 * Klo[c] >= yy) start = 0;
                        else {
                            nprec++;
                            long double rho1 = (long double)(n[c] - 1) * (long double)(n[c] - 1) * CPHI[c] * logl((long double)p)
                                               / (E_LD * (long double)p * (long double)m[c]);
                            if (fabsl(rho1 - 1) < 1e-17L) { namb++; EMIT("A %llu %d %llu %llu %d %.3Le start\n", (unsigned long long)p, c, (unsigned long long)(n[c] - 1), (unsigned long long)m[c], (int)(rho1 >= 1), rho1 - 1); }
                            start = rho1 < 1;
                        }
                    }
                    nstart[c] += start;
                    nfail[c]++;
                    if (!f_tot[c]) { f_first[c][0] = p; f_first[c][1] = n[c]; f_first[c][2] = m[c]; }
                    f_last[c][0] = p; f_last[c][1] = n[c]; f_last[c][2] = m[c];
                    f_tot[c]++;
                }
            }
            p = primesieve_next_prime(&lead);
        }
        for (int c = 0; c < NCLS; c++) if (mask >> c & 1) {
            EMIT("B %llu %llu %d %d %llu %llu %llu %.18Lg %llu %llu %llu %llu\n", (unsigned long long)a, (unsigned long long)b, c, bin,
                 (unsigned long long)(n[c] - n_bin0[c]), (unsigned long long)nfail[c], (unsigned long long)nfail0[c], mx[c],
                 (unsigned long long)mx_p[c], (unsigned long long)mx_n[c], (unsigned long long)mx_m[c], (unsigned long long)nstart[c]);
            if (f_tot[c]) EMIT("F %llu %llu %d %llu %llu %llu %llu %llu %llu %llu\n", (unsigned long long)a, (unsigned long long)b, c,
                 (unsigned long long)f_first[c][0], (unsigned long long)f_first[c][1], (unsigned long long)f_first[c][2],
                 (unsigned long long)f_last[c][0], (unsigned long long)f_last[c][1], (unsigned long long)f_last[c][2], (unsigned long long)f_tot[c]);
        }
        /* consistency with the pass-1 table at b */
        uint64_t chk[NCLS];
        if (b / SC <= (uint64_t)Kmax) {
            res_below(dir, cum, b, res);
            res_to_cls(res, b, chk);
            for (int c = 0; c < NCLS; c++) {
                if (chk[c] != n[c]) { fprintf(stderr, "COUNT MISMATCH task %llu class %d: %llu vs %llu\n", (unsigned long long)a, c, (unsigned long long)n[c], (unsigned long long)chk[c]); exit(4); }
                if (mask >> c & 1) EMIT("E %llu %llu %d %llu\n", (unsigned long long)a, (unsigned long long)b, c, (unsigned long long)n[c]);
            }
        }
        /* consistency of the LAG counts with the pass-1 table: every lag prime < lq has been
         * counted, so after adding the lag primes in [lq, JG) the counts must be those of the
         * table at the next table point JG */
        {
            uint64_t JG = (lq / G + 1) * G, mm[NCLS];
            if (JG / SC <= (uint64_t)Kmax) {
                memcpy(mm, m, sizeof mm);
                while (lq < JG) {
                    if (lq >= 7) { const int *cc = CLS_OF_RES[RIDX[lq % 120]]; for (int k = 0; k < 6; k++) mm[cc[k]]++; }
                    else { for (int c = 0; c < NCLS; c++) if (CQ[c] % lq != 0 && in_class(lq, c)) mm[c]++; }
                    lq = primesieve_next_prime(&lag);
                }
                res_below(dir, cum, JG, res);
                res_to_cls(res, JG, chk);
                for (int c = 0; c < NCLS; c++)
                    if (chk[c] != mm[c]) { fprintf(stderr, "LAG COUNT MISMATCH task %llu class %d: %llu vs %llu\n", (unsigned long long)a, c, (unsigned long long)mm[c], (unsigned long long)chk[c]); exit(4); }
                nlagchk++;
            }
        }
        primesieve_free_iterator(&lead); primesieve_free_iterator(&lag);
        tot_prec += nprec; tot_amb += namb; tot_lagchk += nlagchk;
        #pragma omp critical
        {
            fwrite(out, 1, olen, stdout); fflush(stdout);
            ndone++;
            if (ndone % 50 == 0 || nt <= 50) fprintf(stderr, "apscan: %ld/%zu tasks, %.0f s\n", ndone, nt, omp_get_wtime() - t0);
        }
        free(out);
    }
    fprintf(stderr, "apscan: %zu tasks done in %.1f s wall; precise evaluations %llu; ambiguous %llu\n", nt, omp_get_wtime() - t0, tot_prec, tot_amb);
    printf("# done tasks=%zu precise=%llu ambiguous=%llu lead_and_lag_count_checks_passed=%llu\n", nt, tot_prec, tot_amb, tot_lagchk);
    return 0;
}
