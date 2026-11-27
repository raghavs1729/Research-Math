// arXiv:2609.05358 (H. Manabe), Conjecture conj:period: S={a,b,c}, 2<=a<b<c, gcd=1, c != a+b.
// If the Sprague-Grundy sequence of the subtraction game S is purely periodic, its least period is in {a+b, c+a, c+b}.
// G is purely periodic with period p  iff  G(n+p)=G(n) for 0<=n<c (then by induction for all n).
#include <stdio.h>
#include <stdlib.h>
static int gcd(int x,int y){ while(y){int t=x%y;x=y;y=t;} return x; }
int main(int argc,char**argv){
    int A0=atoi(argv[1]),A=atoi(argv[2]),B=atoi(argv[3]),C=atoi(argv[4]),MULT=atoi(argv[5]);
    static unsigned char G[4000000];
    long tested=0,pure=0,bad=0; long cnt[4]={0,0,0,0};
    for(int a=A0;a<=A;a++) for(int b=a+1;b<=B;b++) for(int c=b+1;c<=C;c++){
        if(c==a+b) continue;
        if(gcd(gcd(a,b),c)!=1) continue;
        int PM=MULT*c, N=PM+c+1;
        for(int n=0;n<N;n++){
            int m=0, x=(n>=a)?G[n-a]:9, y=(n>=b)?G[n-b]:9, z=(n>=c)?G[n-c]:9;
            while(m==x||m==y||m==z) m++;
            G[n]=m;
        }
        tested++;
        int p,found=0;
        for(p=1;p<=PM;p++){
            int ok=1;
            for(int n=0;n<c;n++) if(G[n+p]!=G[n]){ok=0;break;}
            if(ok){found=1;break;}
        }
        if(found){
            pure++;
            if(p==a+b) cnt[0]++; else if(p==c+a) cnt[1]++; else if(p==c+b) cnt[2]++;
            else { bad++; cnt[3]++; if(bad<=40) printf("COUNTEREXAMPLE a=%d b=%d c=%d least pure period %d (a+b=%d, c+a=%d, c+b=%d)\n",a,b,c,p,a+b,c+a,c+b); }
        }
    }
    printf("a in [%d,%d], b<=%d, c<=%d, periods searched up to %d*c: triples %ld, purely periodic %ld (period a+b: %ld, c+a: %ld, c+b: %ld), other: %ld\n",A0,A,B,C,MULT,tested,pure,cnt[0],cnt[1],cnt[2],cnt[3]);
    return 0;
}
