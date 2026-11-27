// sieve.c -- exact-prime data for log f(x; q, a) (Akbary-Francis (4)):
//   log f(x;q,a) = (1/phi(q)) log log(phi(q) theta(x;q,a)) + sum_{p<=x, p=a (q)} log(1-1/p) - log C(q,a)
// Segmented sieve of Eratosthenes; long double accumulation (Kahan for theta and the log sum).
// usage: sieve q a X logC  samples.out records.out small.out
//   samples.out : u, x=floor(e^u) grid (step DU in u): k=pi(x;q,a), p_k, theta, logsum, logf
//   records.out : k at which G(p_k) = phi sqrt(p_k) log(p_k) logf(p_k) sets a new maximum (after the
//                 first 10 primes; printed when it exceeds the last printed one by 1e-3), every k
//                 with logf >= 0 (failure of (3)), and a final END line with the exact maximum
//   small.out   : all k with p_k <= 2e6: k, p_k, logf
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
typedef long double LD;
static LD th=0, thc=0, ls=0, lsc=0;
static inline void kadd(LD*s,LD*c,LD x){LD y=x-*c; LD t=*s+y; *c=(t-*s)-y; *s=t;}
int main(int argc,char**argv){
  int q=atoi(argv[1]), a=atoi(argv[2]); uint64_t X=(uint64_t)atof(argv[3]); LD logC=strtold(argv[4],NULL);
  FILE*fs=fopen(argv[5],"w"),*fr=fopen(argv[6],"w"),*fm=fopen(argv[7],"w");
  int phi=q-1; // q prime
  uint64_t sq=(uint64_t)sqrt((double)X)+2;
  char*small=calloc(sq+1,1); uint32_t*pr=malloc(sizeof(uint32_t)*(sq/2+10)); int np=0;
  for(uint64_t i=2;i<=sq;i++){ if(!small[i]){pr[np++]=i; for(uint64_t j=i*i;j<=sq;j+=i)small[j]=1;} }
  const uint64_t SEG=1<<22; char*seg=malloc(SEG);
  const double DU=0.005; double unext=3.0; uint64_t xnext=(uint64_t)floor(exp(unext));
  uint64_t k=0, pk=0; LD logf=0, gmax=-1e30L, gprint=-1e30L; uint64_t nfail=0;
  for(uint64_t lo=2; lo<=X; lo+=SEG){
    uint64_t hi=lo+SEG-1; if(hi>X)hi=X;
    memset(seg,0,SEG);
    for(int i=0;i<np;i++){ uint64_t p=pr[i]; if(p*p>hi)break; uint64_t s=(lo+p-1)/p*p; if(s<p*p)s=p*p; for(uint64_t j=s;j<=hi;j+=p)seg[j-lo]=1; }
    // first n >= lo with n = a mod q
    uint64_t n=lo+((a+q-(lo%q))%q);
    for(;n<=hi;n+=q){
      // emit samples for grid points strictly below n
      while(xnext<n){ if(k>0) fprintf(fs,"%.3f %llu %llu %llu %.20Lg %.20Lg %.12Le\n",unext,(unsigned long long)xnext,(unsigned long long)k,(unsigned long long)pk,th,ls,logf);
        unext+=DU; xnext=(uint64_t)floor(exp(unext)); }
      if(seg[n-lo])continue;
      k++; pk=n; LD p=(LD)n;
      kadd(&th,&thc,logl(p)); kadd(&ls,&lsc,log1pl(-1.0L/p));
      logf=logl(logl(phi*th))/phi+ls-logC;
      if(n<=2000000) fprintf(fm,"%llu %llu %.12Le\n",(unsigned long long)k,(unsigned long long)n,logf);
      LD G=phi*sqrtl(p)*logl(p)*logf;
      if(logf>=0){nfail++; fprintf(fr,"FAIL k=%llu p=%llu logf=%.12Le G=%.8Lf\n",(unsigned long long)k,(unsigned long long)n,logf,G);}
      if(k>10 && G>gmax){gmax=G; if(G>gprint+1e-3){gprint=G; fprintf(fr,"MAX k=%llu p=%llu logf=%.12Le G=%.8Lf\n",(unsigned long long)k,(unsigned long long)n,logf,G); fflush(fr);}}
    }
  }
  while(xnext<=X){ fprintf(fs,"%.3f %llu %llu %llu %.20Lg %.20Lg %.12Le\n",unext,(unsigned long long)xnext,(unsigned long long)k,(unsigned long long)pk,th,ls,logf); unext+=DU; xnext=(uint64_t)floor(exp(unext)); }
  fprintf(fr,"END X=%llu k=%llu p_k=%llu theta=%.20Lg logsum=%.20Lg logf=%.12Le nfail=%llu Gmax=%.8Lf\n",(unsigned long long)X,(unsigned long long)k,(unsigned long long)pk,th,ls,logf,(unsigned long long)nfail,gmax);
  return 0;
}
