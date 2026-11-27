// CRIM (arXiv:2606.16828): impartial game on Young diagrams; a move removes a row or a column.
// A partition = boundary word b_0..b_{L-1} over {E=1, N=0}, b_0 = E, b_{L-1} = N (L = rows + columns);
// removing a column/row = deleting an E/N, then trimming leading N's and trailing E's.
// rank = (#columns) - (#rows) = L mod 2 parity.  Conjecture conj:main: every P-position has even rank (L even).
// We compute all P-positions with L <= LMAX and count those with odd L.
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#define LMAXMAX 34
static uint8_t *P[LMAXMAX+1];   // bit arrays; index = middle bits (b_1..b_{L-2}), b_1 = least significant
static inline int getP(int L, uint64_t mid){ return (P[L][mid>>3]>>(mid&7))&1; }
static inline void setP(int L, uint64_t mid){ P[L][mid>>3] |= (uint8_t)(1u<<(mid&7)); }
// full word as integer: bit i = b_i.  canonicalise a word of length len (arbitrary bits) and return P-status
static inline int statusP(uint64_t w, int len){
    // strip leading N (b_0 = 0)
    while(len>0 && !(w&1)){ w>>=1; len--; }
    // strip trailing E (top bit = 1)
    while(len>0 && ((w>>(len-1))&1)){ len--; }
    if(len==0) return 1;            // empty partition: terminal, P-position
    w &= ((1ULL<<len)-1);
    uint64_t mid = (w>>1) & ((len>=2)? ((1ULL<<(len-2))-1) : 0);
    return getP(len, mid);
}
int main(int argc, char**argv){
    int LMAX = atoi(argv[1]);
    for(int L=2; L<=LMAX; L++){
        uint64_t n = 1ULL<<(L-2);
        P[L] = calloc((n>>3)+1, 1);
        uint64_t cntP = 0, firstmid = 0; int have = 0;
        for(uint64_t mid=0; mid<n; mid++){
            uint64_t w = 1ULL | (mid<<1);          // b_0 = 1 (E), b_{L-1} = 0 (N)
            int isP = 1;
            for(int i=0; i<L && isP; i++){
                if(i>0 && (((w>>i)&1) == ((w>>(i-1))&1))) continue;   // deleting within a run gives the same word
                uint64_t lo = w & ((1ULL<<i)-1), hi = w>>(i+1);
                uint64_t c = lo | (hi<<i);
                if(statusP(c, L-1)) isP = 0;
            }
            if(isP){ setP(L, mid); cntP++; if(!have){ have=1; firstmid=mid; } }
        }
        printf("L=%d: words %llu, P-positions %llu", L, (unsigned long long)n, (unsigned long long)cntP);
        if(cntP && (L&1)) printf("   <-- ODD-RANK P-POSITION, first middle bits (b_1 first) = %llu", (unsigned long long)firstmid);
        printf("\n"); fflush(stdout);
    }
    return 0;
}
