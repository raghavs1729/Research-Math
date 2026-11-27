// arXiv:2501.14640 (Impartial Chess on Integer Partitions), two closing conjectures on SG values of KING
// on generalized staircases GS(r,c,k) = [(kc)^r, ((k-1)c)^r, ..., c^r].
// Position = subpartition lambda[i,j]; KING moves (1,0),(0,1),(1,1) as long as the subpartition is nonempty,
// i.e. a king walking on the cells of the Young diagram (0-based cell (i,j) present iff j < lambda_{i+1}).
// g(i,j) = mex of g over present cells among (i+1,j),(i,j+1),(i+1,j+1);  SG(lambda) = g(0,0).
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
static int sgKing(int r,int c,int k){
    int R=r*k, C=c*k;
    unsigned char *cur=malloc(C+2), *nxt=malloc(C+2);   // rows i and i+1 ; value 9 = absent
    memset(nxt,9,C+2);
    for(int i=R-1;i>=0;i--){
        int len=c*(k-i/r);                               // lambda_{i+1}
        memset(cur,9,C+2);
        for(int j=len-1;j>=0;j--){
            int a=nxt[j], b=cur[j+1], d=nxt[j+1], m=0;
            while(m==a||m==b||m==d) m++;
            cur[j]=m;
        }
        unsigned char*t=cur;cur=nxt;nxt=t;
    }
    int v=nxt[0]; free(cur); free(nxt); return v;
}
int main(int argc,char**argv){
    int RM=atoi(argv[1]), KM=atoi(argv[2]);
    long bad1=0,n1=0,bad2=0,n2=0;
    for(int r=1;r<=RM;r+=2) for(int k=1;k<=KM;k++){
        int v=sgKing(r,r,k), m=k%(r+2), pred=(m%2==1)?0:(m!=0?1:2); n1++;
        if(v!=pred){ if(bad1<12) printf("Conj 1 FAILS: r=c=%d k=%d: SG=%d predicted %d\n",r,k,v,pred); bad1++; }
    }
    for(int r=2;r<=RM;r+=2) for(int c=2;c<=RM;c+=2) if(r!=c) for(int k=1;k<=KM;k++){
        int v=sgKing(r,c,k), pred=(k%2==1)?2:3; n2++;
        if(v!=pred){ if(bad2<12) printf("Conj 2 FAILS: r=%d c=%d k=%d: SG=%d predicted %d\n",r,c,k,v,pred); bad2++; }
    }
    printf("Conj 1: odd r<=%d, k<=%d: %ld cases, %ld failures.  Conj 2: even r!=c <=%d, k<=%d: %ld cases, %ld failures\n",RM,KM,n1,bad1,RM,KM,n2,bad2);
    return 0;
}
