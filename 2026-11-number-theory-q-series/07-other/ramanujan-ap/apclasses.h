/* Shared definitions: residue classes tracked, residues mod 120, exact floor(e*n).
 * Classes c = (q,a): q in {1,3,4,5,8,12}; 17 classes in all.  Every prime p >= 7
 * lies in exactly 6 of them, determined by p mod 120. */
#ifndef APCLASSES_H
#define APCLASSES_H
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#define NRES 32
#define NCLS 17
static const int RES[NRES] = {1,7,11,13,17,19,23,29,31,37,41,43,47,49,53,59,
                              61,67,71,73,77,79,83,89,91,97,101,103,107,109,113,119};
static const int CQ[NCLS] = {1, 3,3, 4,4, 5,5,5,5, 8,8,8,8, 12,12,12,12};
static const int CA[NCLS] = {0, 1,2, 1,3, 1,2,3,4, 1,3,5,7, 1,5,7,11};
static const int CPHI[NCLS]={1, 2,2, 2,2, 4,4,4,4, 4,4,4,4, 4,4,4,4};

static int RIDX[120];          /* residue -> index 0..31, or 32 (not coprime to 120) */
static int CLS_OF_RES[NRES][6];/* the 6 classes containing residue r (p >= 7)        */

static inline int in_class(uint64_t p, int c) {
    int q = CQ[c];
    if (q == 1) return 1;
    return (int)(p % q) == CA[c];
}
static void ap_init(void) {
    for (int r = 0; r < 120; r++) RIDX[r] = 32;
    for (int i = 0; i < NRES; i++) RIDX[RES[i]] = i;
    for (int i = 0; i < NRES; i++) {
        int k = 0;
        for (int c = 0; c < NCLS; c++) if (in_class(RES[i], c)) CLS_OF_RES[i][k++] = c;
        if (k != 6) { fprintf(stderr, "class table error\n"); exit(1); }
    }
}
/* class counts from residue counts (primes >= 7 only) plus the primes 2,3,5 below bound
 * (counting primes < bound). */
static void res_to_cls(const uint64_t *res, uint64_t bound, uint64_t *cls) {
    for (int c = 0; c < NCLS; c++) cls[c] = 0;
    for (int i = 0; i < NRES; i++)
        for (int k = 0; k < 6; k++) cls[CLS_OF_RES[i][k]] += res[i];
    static const uint64_t sp[3] = {2,3,5};
    for (int j = 0; j < 3; j++) if (sp[j] < bound)
        for (int c = 0; c < NCLS; c++) {
            int q = CQ[c];
            if (q % sp[j] == 0) continue;           /* p | q: in no class mod q */
            if (in_class(sp[j], c)) cls[c]++;
        }
}

/* floor(e*n) exactly, for n < 2^62.  E = floor(e * 2^126) = EHI*2^64 + ELO.
 * e*n lies in [n*E/2^126, n*(E+1)/2^126); the floor is certain unless the
 * fractional part is within n/2^126 of 1.  We abort on a superset of these cases:
 * whenever the 62 fractional bits kept are all ones. */
#define E_HI 0xADF85458A2BB4A9AULL
#define E_LO 0xAFDC5620273D3CF1ULL
static inline uint64_t floor_e_times(uint64_t n) {
    unsigned __int128 p1 = (unsigned __int128)n * E_HI;
    unsigned __int128 p0 = (unsigned __int128)n * E_LO;
    unsigned __int128 s = p1 + (p0 >> 64);
    const unsigned __int128 mask = (((unsigned __int128)1) << 62) - 1;
    if ((s & mask) == mask) { fprintf(stderr, "floor_e_times ambiguous at %llu\n", (unsigned long long)n); exit(2); }
    return (uint64_t)(s >> 62);
}
#endif
