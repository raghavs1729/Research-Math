/* Pass 1: exact counts of primes by residue class mod 120 on a fine grid.
 *
 *   apcount K0 K1 outdir
 *
 * Superchunk K is [K*SC, (K+1)*SC), SC = 10^10, cut into rows of G = 10^6.
 * For each row and each reduced residue r mod 120 we store (uint16) the number
 * of primes p >= 7 in the row with p = r mod 120.  The primes 2, 3, 5 are never
 * stored; they are added by res_to_cls() in the readers.
 * Output: outdir/K%05d.u16 (10^4 rows x 32 uint16), written atomically.
 * Primes come from libprimesieve (primesieve_iterator).
 */
#include <primesieve.h>
#include <omp.h>
#include <string.h>
#include <time.h>
#include "apclasses.h"

#define SC 10000000000ULL
#define G  1000000ULL
#define ROWS (SC / G)

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: apcount K0 K1 outdir\n"); return 1; }
    long K0 = atol(argv[1]), K1 = atol(argv[2]);
    const char *dir = argv[3];
    ap_init();
    double t0 = omp_get_wtime();
    long done = 0;
    #pragma omp parallel for schedule(dynamic, 1)
    for (long K = K0; K < K1; K++) {
        char fn[512], tmp[512];
        snprintf(fn, sizeof fn, "%s/K%05ld.u16", dir, K);
        FILE *f = fopen(fn, "rb");
        if (f) { fclose(f); continue; }           /* already done */
        uint16_t *buf = calloc(ROWS * NRES, sizeof(uint16_t));
        uint64_t a = (uint64_t)K * SC, b = a + SC;
        primesieve_iterator it;
        primesieve_init(&it);
        primesieve_jump_to(&it, a < 7 ? 7 : a, b + 1000);
        uint64_t p = primesieve_next_prime(&it);
        for (uint64_t row = 0; row < ROWS; row++) {
            uint64_t end = a + (row + 1) * G;
            uint32_t cnt[NRES + 1];
            memset(cnt, 0, sizeof cnt);
            while (p < end) {
                cnt[RIDX[p % 120]]++;
                p = primesieve_next_prime(&it);
            }
            if (cnt[NRES]) { fprintf(stderr, "non-coprime prime?!\n"); exit(1); }
            for (int i = 0; i < NRES; i++) {
                if (cnt[i] > 65535) { fprintf(stderr, "overflow\n"); exit(1); }
                buf[row * NRES + i] = (uint16_t)cnt[i];
            }
        }
        primesieve_free_iterator(&it);
        snprintf(tmp, sizeof tmp, "%s/K%05ld.tmp%d", dir, K, omp_get_thread_num());
        f = fopen(tmp, "wb");
        if (!f || fwrite(buf, sizeof(uint16_t), ROWS * NRES, f) != ROWS * NRES) { perror("write"); exit(1); }
        fclose(f);
        rename(tmp, fn);
        free(buf);
        #pragma omp critical
        {
            done++;
            if (done % 20 == 0 || K1 - K0 <= 20)
                fprintf(stderr, "apcount: %ld/%ld superchunks, %.0f s\n", done, K1 - K0, omp_get_wtime() - t0);
        }
    }
    fprintf(stderr, "apcount: done K=[%ld,%ld) in %.1f s wall\n", K0, K1, omp_get_wtime() - t0);
    return 0;
}
