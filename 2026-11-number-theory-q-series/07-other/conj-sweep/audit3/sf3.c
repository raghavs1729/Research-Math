/* AUDIT-3 referee program for Theorem 2 of R06-3-sumfree.md (greedy strict 3-sumfree sequences, arXiv:2609.16843).

   Algorithm ("pull" form, written independently of code/ and audit2/):
     cnt2[v] = number of unordered pairs {a,b}, a != b, of current entries with a+b = v.
     For a candidate x and an entry c < x the number of pairs {a,b} of entries, a != b, c not in {a,b}, with
     a+b = x-c is   t_c(x) = cnt2[x-c] - [x-2c is an entry, x-2c != c, x-2c > 0].
     x is a sum of three distinct entries  <=>  t_c(x) > 0 for some entry c;
     the number of such representations is r(x) = (1/3) sum_c t_c(x).
   All entries present when x is examined are < x, so "previous entries" = "entries".

   usage:
     sf3 seq f g h X            entries of S_{f,g,h} that are <= X
     sf3 rep f g h X            x, [x in S], r(x) for h < x <= X   (full representation counts)
     sf3 family FLO FHI DMAX NPER
        for FLO <= f <= FHI and 2f+1 <= d <= DMAX, with m = 8d-2f-1 and X = 10d-3f-2 + NPER*m - 1:
        (A) greedy simulation of S_{f,d-f,2d-f} on [1,X] compared with the set T of Theorem 2;
        (B) direct check of the statements (I), (II) about T itself:
            for every x in (2d-f, X]:  x in T  <=>  x is NOT a sum of three distinct elements of T.
        Prints one line per pair with a mismatch and a summary. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int *ent, nent;            /* entries in increasing order */
static unsigned char *isent;      /* membership */
static int *cnt2;                 /* pair counts */
static int XMAX;

static void reset(int X) {
    XMAX = X;
    ent = realloc(ent, sizeof(int) * (size_t)(X + 8));
    isent = realloc(isent, (size_t)(X + 8));
    cnt2 = realloc(cnt2, sizeof(int) * (size_t)(2 * X + 16));
    memset(isent, 0, (size_t)(X + 8));
    memset(cnt2, 0, sizeof(int) * (size_t)(2 * X + 16));
    nent = 0;
}
static void add_entry(int s) {
    for (int i = 0; i < nent; i++) cnt2[ent[i] + s]++;
    ent[nent++] = s; isent[s] = 1;
}
/* t_c(x) */
static inline int tc(int x, int c) {
    int v = x - c;
    if (v <= 0) return 0;
    int t = cnt2[v];
    int b = v - c;                       /* pair {c, b} */
    if (b > 0 && b != c && b <= XMAX && isent[b]) t--;
    return t;
}
static int is3sum(int x) {               /* early exit, largest c first */
    for (int i = nent - 1; i >= 0; i--) if (tc(x, ent[i]) > 0) return 1;
    return 0;
}
static long long nrep(int x) {           /* exact number of representations */
    long long s = 0;
    for (int i = 0; i < nent; i++) s += tc(x, ent[i]);
    if (s % 3) { fprintf(stderr, "internal error: sum of t_c not divisible by 3 at x=%d\n", x); exit(2); }
    return s / 3;
}
static void simulate(int f, int g, int h, int X) {
    reset(X);
    add_entry(f); add_entry(g); add_entry(h);
    for (int x = h + 1; x <= X; x++) if (!is3sum(x)) add_entry(x);
}
static int inT(int f, int d, int x) {    /* the set T of Theorem 2 */
    int m = 8 * d - 2 * f - 1;
    if (x == f || x == d - f) return 1;
    if (x >= 2 * d - f && x <= 3 * d - f - 1) return 1;
    if (x == 4 * d - f) return 1;
    if (x >= 10 * d - 3 * f - 2) { int r = (x - (10 * d - 3 * f - 2)) % m; if (r <= d - 1) return 1; }
    return 0;
}

int main(int argc, char **argv) {
    if (argc == 6 && !strcmp(argv[1], "seq")) {
        int f = atoi(argv[2]), g = atoi(argv[3]), h = atoi(argv[4]), X = atoi(argv[5]);
        simulate(f, g, h, X);
        for (int i = 0; i < nent; i++) printf("%d%c", ent[i], i + 1 < nent ? ' ' : '\n');
        return 0;
    }
    if (argc == 6 && !strcmp(argv[1], "rep")) {
        int f = atoi(argv[2]), g = atoi(argv[3]), h = atoi(argv[4]), X = atoi(argv[5]);
        reset(X); add_entry(f); add_entry(g); add_entry(h);
        for (int x = h + 1; x <= X; x++) {
            long long r = nrep(x);
            printf("%d %d %lld\n", x, r == 0, r);
            if (r == 0) add_entry(x);
        }
        return 0;
    }
    if (argc == 6 && !strcmp(argv[1], "family")) {
        int FLO = atoi(argv[2]), FHI = atoi(argv[3]), DMAX = atoi(argv[4]), NPER = atoi(argv[5]);
        long npairs = 0, nin = 0, badA_in = 0, badB_in = 0, nout = 0, badA_out = 0;
        for (int f = FLO; f <= FHI; f++) for (int d = 2 * f + 1; d <= DMAX; d++) {
            int m = 8 * d - 2 * f - 1, X = 10 * d - 3 * f - 2 + NPER * m - 1;
            int g = d - f, h = 2 * d - f, inrange = (d >= 6);
            npairs++;
            /* (A) greedy simulation */
            simulate(f, g, h, X);
            int firstA = 0, inS_A = 0;
            for (int x = 1; x <= X; x++) if (isent[x] != inT(f, d, x)) { firstA = x; inS_A = isent[x]; break; }
            /* (B) statements (I),(II) on T */
            reset(X);
            int firstB = 0, kind = 0;
            for (int x = 1; x <= X; x++) {
                int t = inT(f, d, x);
                if (x > h) {
                    int s3 = is3sum(x);          /* entries present = T ∩ [1,x) */
                    if (t && s3 && !firstB) { firstB = x; kind = 1; }   /* (I) violated  */
                    if (!t && !s3 && !firstB) { firstB = x; kind = 2; } /* (II) violated */
                }
                if (t) add_entry(x);
            }
            if (inrange) { nin++; if (firstA) badA_in++; if (firstB) badB_in++; }
            else { nout++; if (firstA) badA_out++; }
            if (firstA || firstB || !inrange)
                printf("(f,d)=(%d,%d) %s: (A) first x with [x in S] != [x in T]: %d (%s)   (B) first violation of (I)/(II) for T: %d%s\n",
                       f, d, inrange ? "in range" : "d<6",
                       firstA, firstA ? (inS_A ? "in S, not in T" : "in T, not in S") : "none",
                       firstB, kind == 1 ? " [(I): element of T is a 3-sum]" : kind == 2 ? " [(II): non-element is not a 3-sum]" : "");
        }
        printf("f = %d..%d, 2f+1 <= d <= %d, X = 10d-3f-2 + %d*m - 1 (m = 8d-2f-1):\n", FLO, FHI, DMAX, NPER);
        printf("  pairs with d >= 6: %ld; (A) S != T on [1,X]: %ld; (B) (I)/(II) violated: %ld\n", nin, badA_in, badB_in);
        printf("  pairs with 2f+1 <= d < 6: %ld; of these S != T: %ld\n", nout, badA_out);
        return 0;
    }
    fprintf(stderr, "usage: see header\n");
    return 1;
}
