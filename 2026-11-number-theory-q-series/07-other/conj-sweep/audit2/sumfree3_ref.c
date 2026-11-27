/* Referee simulation of the greedy strict 3-sumfree sequence S_{f,g,h} (arXiv:2609.16843).
   pairs[v] = number of unordered pairs of distinct entries with sum v; is3[v] = v is a sum of 3 distinct entries.
   usage: sumfree3_ref f g h X      -> prints all entries <= X
          sumfree3_ref -T dlo dhi J -> for each d compares S_{1,d-1,2d-1} with the set T of R06 on [1, 10d-5+J(8d-3)] */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
static int *S, nS; static unsigned char *is3, *inS; static int *pairs;
static void add(int s,int X){
  for(int v=0; v+s<=X; v++) if(pairs[v]) is3[v+s]=1;       /* pairs of earlier (smaller) entries + s */
  for(int i=0;i<nS;i++) if(S[i]+s<=X) pairs[S[i]+s]++;
  S[nS++]=s; inS[s]=1;
}
static void simulate(int f,int g,int h,int X){
  S=malloc(sizeof(int)*(X+5)); nS=0; is3=calloc(X+2,1); inS=calloc(X+2,1); pairs=calloc(X+2,sizeof(int));
  add(f,X); add(g,X); add(h,X);
  for(int x=h+1;x<=X;x++) if(!is3[x]) add(x,X);
}
int main(int argc,char**argv){
  if(argc==5 && strcmp(argv[1],"-T")){
    int f=atoi(argv[1]),g=atoi(argv[2]),h=atoi(argv[3]),X=atoi(argv[4]);
    simulate(f,g,h,X); for(int i=0;i<nS;i++) printf("%d%c",S[i],i+1<nS?' ':'\n'); return 0;
  }
  int dlo=atoi(argv[2]),dhi=atoi(argv[3]),J=atoi(argv[4]); int bad=0;
  for(int d=dlo;d<=dhi;d++){
    int m=8*d-3, X=10*d-5+J*m;
    simulate(1,d-1,2*d-1,X);
    int mism=0;
    for(int x=1;x<=X;x++){
      int t = (x==1)||(x==d-1)||(x>=2*d-1&&x<=3*d-2)||(x==4*d-1);
      if(x>=10*d-5){ int r=(x-(10*d-5))%m; if(r<=d-1) t=1; }
      if(t!=inS[x]){ if(!mism) printf("d=%d first mismatch at x=%d (in S: %d, in T: %d)\n",d,x,inS[x],t); mism++; }
    }
    if(mism) bad++;
    free(S);free(is3);free(inS);free(pairs);
  }
  printf("d = %d..%d: S_{1,d-1,2d-1} compared with T on [1, 10d-5+%d(8d-3)]: %d values of d with a mismatch\n",dlo,dhi,J,bad);
  return 0;
}
