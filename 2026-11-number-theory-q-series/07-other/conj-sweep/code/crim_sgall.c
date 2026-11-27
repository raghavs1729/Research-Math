// CRIM (arXiv:2606.16828): Sprague-Grundy values of ALL partitions with rows+columns <= LMAX (boundary words),
// then tests of the Sprague-Grundy conjectures on rectairs R^k_{r,c} = [c^(r-k), c-1, ..., c-k] and staircases.
// word: b_0..b_{L-1}, E=1, N=0, b_0 = E, b_{L-1} = N; bottom row first:  E^{lam_r} N E^{lam_{r-1}-lam_r} N ... N.
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#define LMM 34
static uint8_t *SG[LMM+1];
static inline int sgw(uint64_t w, int len){
    while(len>0 && !(w&1)){ w>>=1; len--; }
    while(len>0 && ((w>>(len-1))&1)) len--;
    if(len==0) return 0;
    w &= ((1ULL<<len)-1);
    return SG[len][(w>>1) & ((len>=2)?((1ULL<<(len-2))-1):0)];
}
// boundary word of a partition given as parts lam[0] >= ... >= lam[r-1] > 0 ; returns length, word in *w
static int word(const int*lam,int r,uint64_t*w){
    uint64_t x=0; int pos=0, prev=0;
    for(int i=r-1;i>=0;i--){ for(int j=0;j<lam[i]-prev;j++){ x|=(1ULL<<pos); pos++; } pos++; prev=lam[i]; }
    *w=x; return pos;
}
static int sgpart(const int*lam,int r){ uint64_t w; int L=word(lam,r,&w); return sgw(w,L); }
int main(int argc,char**argv){
    int LMAX=atoi(argv[1]);
    for(int L=2;L<=LMAX;L++){
        uint64_t n=1ULL<<(L-2); SG[L]=malloc(n); int mx=0;
        for(uint64_t mid=0;mid<n;mid++){
            uint64_t w=1ULL|(mid<<1); uint32_t seen=0;
            for(int i=0;i<L;i++){
                if(i>0 && (((w>>i)&1)==((w>>(i-1))&1))) continue;
                uint64_t lo=w&((1ULL<<i)-1), hi=w>>(i+1);
                int v=sgw(lo|(hi<<i),L-1); seen|=(1u<<v);
            }
            int m=0; while(seen&(1u<<m)) m++;
            SG[L][mid]=m; if(m>mx) mx=m;
        }
        fprintf(stderr,"L=%d max SG %d\n",L,mx);
    }
    int lam[64]; long bad;
    // Conjecture 4: SG(R^k_{r,r}), 0<=k<r : 0 if r even or k<r-1; 1 if r in {3,5} and k=r-1; 2 otherwise
    bad=0; printf("Conj (R^k_{r,r}):");
    for(int r=1;2*r<=LMAX;r++) for(int k=0;k<r;k++){
        for(int i=0;i<r;i++) lam[i]=(i<r-k)?r:r-(i-(r-k)+1);
        int v=sgpart(lam,r), pred=(r%2==0||k<r-1)?0:((r==3||r==5)?1:2);
        if(v!=pred){ bad++; printf(" [r=%d k=%d SG=%d pred=%d]",r,k,v,pred); }
    }
    printf("  r<=%d: %ld mismatches\n",LMAX/2,bad);
    // Conjecture 5: for r>=7: SG(R^k_{r,r-1}) = 3 if k=r-2 and r odd, else 1   (0<=k<min(r,r-1)=r-1)
    bad=0; printf("Conj (R^k_{r,r-1}, r>=7):");
    for(int r=7;2*r-1<=LMAX;r++) for(int k=0;k<r-1;k++){
        int c=r-1; for(int i=0;i<r;i++) lam[i]=(i<r-k)?c:c-(i-(r-k)+1);
        int v=sgpart(lam,r), pred=(k==r-2&&r%2==1)?3:1;
        if(v!=pred){ bad++; printf(" [r=%d k=%d SG=%d pred=%d]",r,k,v,pred); }
    }
    printf("  r<=%d: %ld mismatches\n",(LMAX+1)/2,bad);
    // Conjecture 6: staircase s_n=[n,...,1]: 0 if n even; 1 if n=1,3,5; 2 if n odd >= 7
    bad=0; printf("Conj (staircases s_n):");
    for(int n=1;2*n<=LMAX;n++){ for(int i=0;i<n;i++) lam[i]=n-i; int v=sgpart(lam,n), pred=(n%2==0)?0:((n<=5)?1:2);
        if(v!=pred){ bad++; printf(" [n=%d SG=%d pred=%d]",n,v,pred);} }
    printf("  n<=%d: %ld mismatches\n",LMAX/2,bad);
    // Conjecture 7: padded staircase ps_n=[n,n,n-1,...,1]: 2 if n=1,3; 1 if n>=6 even; 5 if n=5; 3 otherwise
    bad=0; printf("Conj (padded staircases ps_n, n>=2 by definition; n=1 means [1,1]):");
    for(int n=1;2*n+1<=LMAX;n++){ lam[0]=n; for(int i=0;i<n;i++) lam[i+1]=n-i; int v=sgpart(lam,n+1);
        int pred=(n==1||n==3)?2:((n>=6&&n%2==0)?1:((n==5)?5:3));
        if(v!=pred){ bad++; printf(" [n=%d SG=%d pred=%d]",n,v,pred);} }
    printf("  n<=%d: %ld mismatches\n",(LMAX-1)/2,bad);
    return 0;
}
