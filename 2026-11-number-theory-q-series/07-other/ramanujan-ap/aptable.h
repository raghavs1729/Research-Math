/* Reader for the pass-1 tables written by apcount (outdir/K%05d.u16). */
#ifndef APTABLE_H
#define APTABLE_H
#include <string.h>
#include "apclasses.h"
#define SC 10000000000ULL
#define G  1000000ULL
#define ROWS (SC / G)

static uint16_t *load_file(const char *dir, long K) {
    char fn[512];
    snprintf(fn, sizeof fn, "%s/K%05ld.u16", dir, K);
    FILE *f = fopen(fn, "rb");
    if (!f) { fprintf(stderr, "missing table %s\n", fn); exit(3); }
    uint16_t *buf = malloc(ROWS * NRES * sizeof(uint16_t));
    if (fread(buf, sizeof(uint16_t), ROWS * NRES, f) != ROWS * NRES) { fprintf(stderr, "short table %s\n", fn); exit(3); }
    fclose(f);
    return buf;
}
/* cumS[K][i] = number of primes 7 <= p < K*SC with p = RES[i] mod 120, K = 0..Kmax */
static uint64_t (*load_cum(const char *dir, long Kmax))[NRES] {
    uint64_t (*cum)[NRES] = calloc(Kmax + 1, sizeof *cum);
    for (long K = 0; K < Kmax; K++) {
        uint16_t *buf = load_file(dir, K);
        memcpy(cum[K + 1], cum[K], sizeof cum[K]);
        for (uint64_t r = 0; r < ROWS; r++)
            for (int i = 0; i < NRES; i++) cum[K + 1][i] += buf[r * NRES + i];
        free(buf);
    }
    return cum;
}
/* residue counts of primes 7 <= p < x, x a multiple of G */
static void res_below(const char *dir, uint64_t (*cum)[NRES], uint64_t x, uint64_t *res) {
    if (x % G) { fprintf(stderr, "res_below: %llu not a multiple of G\n", (unsigned long long)x); exit(3); }
    long K = (long)(x / SC);
    uint64_t rows = (x % SC) / G;
    memcpy(res, cum[K], NRES * sizeof(uint64_t));
    if (rows) {
        uint16_t *buf = load_file(dir, K);
        for (uint64_t r = 0; r < rows; r++)
            for (int i = 0; i < NRES; i++) res[i] += buf[r * NRES + i];
        free(buf);
    }
}
#endif
