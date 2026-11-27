/* Referee brute force for CRIM on boundary words (independent of code/crim.c).
   Word of length L = integer w, bit i = letter i (0 = E, 1 = N), letter 0 is the first letter.
   P[L][w] := "canon(w) is a losing position", computed from the game recursion only:
     leading N  -> P[L-1][w without it];  trailing E -> P[L-1][w without it];
     otherwise  -> no single-letter deletion v has P[L-1][v].
   Then compared with the automaton (read from stdin: n, then lines "acc E N") on ALL words of length <= LMAX.
   Optionally dumps the tables to a file (argv[2]) for the corollary tests. */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
static uint8_t *P[40];
#define GET(L,w) ((P[L][(w)>>3]>>((w)&7))&1)
#define SET(L,w) (P[L][(w)>>3]|=(uint8_t)(1u<<((w)&7)))
int main(int argc,char**argv){
  int LMAX=argc>1?atoi(argv[1]):24; int n; int acc[64],dE[64],dN[64];
  if(scanf("%d",&n)!=1) return 1;
  for(int q=0;q<n;q++) if(scanf("%d %d %d",&acc[q],&dE[q],&dN[q])!=3) return 1;
  for(int L=0;L<=LMAX;L++){ size_t sz=((size_t)1<<L)/8+1; P[L]=calloc(sz,1); }
  SET(0,0);
  long long mism=0;
  for(int L=1;L<=LMAX;L++){
    uint64_t tot=(uint64_t)1<<L; long long canonP=0,canon=0;
    for(uint64_t w=0;w<tot;w++){
      int v;
      if(w&1) v=GET(L-1,w>>1);
      else if(!((w>>(L-1))&1)) v=GET(L-1,w);
      else{
        v=1;
        for(int i=0;i<L;i++){
          uint64_t low=w&(((uint64_t)1<<i)-1), high=w>>(i+1);
          uint64_t u=low|(high<<i);
          if(GET(L-1,u)){v=0;break;}
        }
        canon++; canonP+=v;
      }
      if(v) SET(L,w);
      int q=0; for(int i=0;i<L;i++) q=((w>>i)&1)?dN[q]:dE[q];
      if(acc[q]!=v) mism++;
    }
    printf("L=%2d canonical words %10lld  losing %10lld%s\n",L,canon,canonP,(L%2&&canonP)?"  <-- ODD-LENGTH LOSING":"");
  }
  printf("automaton vs game recursion on all words (canonical or not) of length <= %d: mismatches %lld\n",LMAX,mism);
  if(argc>2){ FILE*f=fopen(argv[2],"wb"); for(int L=0;L<=LMAX;L++) fwrite(P[L],1,((size_t)1<<L)/8+1,f); fclose(f); }
  return 0;
}
