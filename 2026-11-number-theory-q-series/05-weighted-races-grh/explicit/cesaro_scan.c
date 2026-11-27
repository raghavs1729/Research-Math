/*
 * cesaro_scan.c -- exact scan of the Cesaro mean of Chebyshev's race and of the
 * integrated theta-race, for the explicit theorem in EXPLICIT.md.
 * (Sieve loop adapted from ../sieve/nest_scan.c; everything else is new.)
 *
 * s(p) = -chi(p):  q=4: +1 if p=3 mod 4, -1 if p=1 mod 4, 0 for p=2
 *                  q=3: +1 if p=2 mod 3, -1 if p=1 mod 3, 0 for p=3
 * D0(x) = sum s(p),  D1(x) = sum s(p) p                      (exact integers)
 * R(x)  = sum_{p<=x} s(p)(1-p/x) = D0(x) - D1(x)/x           (continuous in x)
 * N(x)  = R(x) log x / sqrt x
 * Th(t) = sum s(p) log p,  W(t) = sum s(p) p log p            (long double)
 * Th1(t)= int_0^t Th = t Th(t) - W(t)                         (continuous, piecewise linear)
 * G(t)  = Th1(t)/t^{3/2}
 * theta(y) = sum_{p<=y} log p (all primes), vo(y) = theta(y) - log(q') for y>=q'
 *            (q' = 2 for q=4, 3 for q=3): theta-function of the principal character.
 *
 * On a gap [p,p'] between consecutive primes R(x)=A-B/x is monotone and Th1 is
 * linear, and both are continuous, so
 *    min_{[p,p']} R   = min(R(p),R(p')),      min_{[p,p']} Th1 = min(Th1(p),Th1(p')).
 * Rigorous lower bounds on the gap (when the minima are >=0 and p>=11, where
 * log x/sqrt x and x^{-3/2} are decreasing):
 *    N(x) >= min(R(p),R(p')) * log p'/sqrt p',    G(x) >= min(Th1(p),Th1(p'))/p'^{3/2}.
 * Each gap is credited to the dyadic bins of both of its endpoints, so the
 * "low" columns are lower bounds over ALL real x in the bin, not only primes.
 *
 * Also:  c1  = sup_{y<=X} (y - vo(y))/sqrt y   (attained at left limits at primes)
 *        cup = sup_{y<=X} (vo(y) - y)/sqrt y   (attained at primes)
 *        c3  = sup_{y<=X} |Th(y)|/sqrt y       (attained at primes)
 *        C1B = the same sup as c1 restricted to y >= 1423 (left limits at primes > 1423)
 *        T1MIN, NGAP: true margins of the sign tests (see EXPLICIT.md, section 5)
 * Exactness: D0, D1 are exact integers; R(p) = D0 - D1/p and N are obtained from them by
 * a fixed handful of long double operations (no accumulation).  Th, W, theta are long
 * double ACCUMULATIONS (rounding analysis in EXPLICIT.md); the decisive range p <= 10^4
 * is certified separately in exact/ball arithmetic by cert_small.py.
 * Exact integer test of R(p)>0:  p*D0 - D1 > 0  (int128).
 *
 * Usage: cesaro_scan q X
 */
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef __int128 i128;
typedef long double ld;
#define NB 40

typedef struct {
    ld Rmin; uint64_t Rmin_at;      /* min of R at primes in bin */
    ld Rlow;                        /* rigorous lower bound of R on bin */
    ld Nmin; uint64_t Nmin_at; ld Nmax; uint64_t Nmax_at; ld Nlow;
    ld Gmin; uint64_t Gmin_at; ld Gmax; uint64_t Gmax_at; ld Glow;
    uint64_t cnt;
} Bin;
static Bin B[NB], ALL;

static void binit(Bin *b){ b->Rmin=b->Rlow=b->Nmin=b->Nlow=b->Gmin=b->Glow=1e300L; b->Nmax=b->Gmax=-1e300L; b->cnt=0;
    b->Rmin_at=b->Nmin_at=b->Nmax_at=b->Gmin_at=b->Gmax_at=0; }
static int bin_of(uint64_t p){ return 63-__builtin_clzll(p); }
static void pr128(i128 v, char *buf){
    char t[48]; int i=0, neg=v<0; if(neg) v=-v;
    do { t[i++]=(char)('0'+(int)(v%10)); v/=10; } while(v);
    int j=0; if(neg) buf[j++]='-'; while(i--) buf[j++]=t[i]; buf[j]=0;
}

static int q;
static int64_t D0=0; static i128 D1=0;
static ld Th=0, W=0, theta=0;
static uint64_t nprimes=0, prevp=0;
static ld Rprev=0, T1prev=0;
static ld c1=-1e300L, cup=-1e300L, c3=-1e300L, ThNmin=1e300L, ThNmax=-1e300L;
static uint64_t c1_at=0, cup_at=0, c3_at=0, ThNmin_at=0, ThNmax_at=0;
static uint64_t firstbadR=0, firstbadT=0, nbadR=0, nbadT=0;
static uint64_t Rmin5=0; static ld Rmin5v=1e300L; /* min of R over primes >= 5 (q=4) / >=5 (q=3) */
static ld c1b=-1e300L; static uint64_t c1b_at=0;      /* sup of (y - vo(y))/sqrt y over left limits at primes > 1423 */
static ld T1min=1e300L; static uint64_t T1min_at=0;   /* min of Theta1(p) (absolute) over race primes after the first */
static ld Gminall=1e300L; static uint64_t Gminall_at=0; /* min of Theta1(p)/p^{3/2} over race primes after the first */
static ld NgapA=1e300L; static uint64_t NgapA_at=0;   /* min over gaps [p,p'] with p >= 29 of min(R(p),R(p')) log p'/sqrt p' */
static ld NgapB=1e300L; static uint64_t NgapB_at=0;   /* same, gaps with p >= 11 */

static void checkpoint(uint64_t X)
{
    /* values at x = X (state holds all primes <= X) */
    char b[48]; pr128(D1,b);
    ld x=(ld)X, R=(ld)D0-(ld)D1/x, T1=x*Th-W;
    ld vo = theta - ((q==4)?logl(2.0L):(X>=3?logl(3.0L):0));
    printf("CHK q=%d x=%llu primes=%llu D0=%lld D1=%s R=%.12Lg N=%.12Lf Theta=%.12Lf Theta1=%.15Lg G=%.12Lf theta=%.9Lf vo-x=%.9Lf\n",
        q,(unsigned long long)X,(unsigned long long)nprimes,(long long)D0,b,R,R*logl(x)/sqrtl(x),Th,T1,T1/(x*sqrtl(x)),theta,vo-x);
    fflush(stdout);
}

static void event(uint64_t p, int s) /* s = -chi(p) in {-1,0,+1}; called for EVERY prime */
{
    ld x=(ld)p, lx=logl(x), sx=sqrtl(x);
    /* theta before adding p : sup of (y - vo(y))/sqrt y on the gap ending at p */
    ld qlog = (q==4)?logl(2.0L):logl(3.0L); uint64_t qp=(q==4)?2:3;
    {   ld vo = theta - (p>qp ? qlog : 0);
        ld v = (x-vo)/sx; if (v>c1){c1=v;c1_at=p;} if (p>1423 && v>c1b){c1b=v;c1b_at=p;} }
    theta += lx;
    {   ld vo = theta - (p>=qp ? qlog : 0);
        ld v = (vo-x)/sx; if (v>cup){cup=v;cup_at=p;} }
    if (s==0) return;
    /* values at x=p: continuity => may use the old or the new state; use new */
    D0 += s; D1 += (i128)s*(i128)p; Th += s*lx; W += s*x*lx; nprimes++;
    ld R = (ld)D0 - (ld)D1/x, T1 = x*Th - W;
    i128 exact = (i128)p*(i128)D0 - D1;            /* = p*R(p), exact */
    int first = (prevp==0);
    if (!first){
        if (exact<=0){ nbadR++; if(!firstbadR) firstbadR=p; }
        if (!(T1>0)){ nbadT++; if(!firstbadT) firstbadT=p; }
        if (R<Rmin5v){Rmin5v=R;Rmin5=p;}
        if (T1<T1min){T1min=T1;T1min_at=p;}
        { ld g0=T1/(x*sx); if (g0<Gminall){Gminall=g0;Gminall_at=p;} }
    }
    ld a=fabsl(Th)/sx; if (a>c3){c3=a;c3_at=p;}
    ld tn=Th/sx; if (tn<ThNmin){ThNmin=tn;ThNmin_at=p;} if (tn>ThNmax){ThNmax=tn;ThNmax_at=p;}
    ld w=lx/sx, N=R*w, x32=x*sx, G=T1/x32;
    int k=bin_of(p);
    Bin *bs[3]; int nb=0; bs[nb++]=&B[k]; bs[nb++]=&ALL;
    if (p>=11){
        for (int i=0;i<nb;i++){ Bin *bb=bs[i]; bb->cnt++;
            if (R<bb->Rmin){bb->Rmin=R;bb->Rmin_at=p;}
            if (N<bb->Nmin){bb->Nmin=N;bb->Nmin_at=p;} if (N>bb->Nmax){bb->Nmax=N;bb->Nmax_at=p;}
            if (G<bb->Gmin){bb->Gmin=G;bb->Gmin_at=p;} if (G>bb->Gmax){bb->Gmax=G;bb->Gmax_at=p;} }
        if (prevp>=11){
            ld gR = Rprev<R?Rprev:R, gT = T1prev<T1?T1prev:T1;
            ld gN = gR>=0 ? gR*w : gR*logl((ld)prevp)/sqrtl((ld)prevp);
            ld gG = gT>=0 ? gT/x32 : gT/((ld)prevp*sqrtl((ld)prevp));
            if (gN<NgapB){NgapB=gN;NgapB_at=prevp;}
            if (prevp>=29 && gN<NgapA){NgapA=gN;NgapA_at=prevp;}
            int kp=bin_of(prevp);
            Bin *gs[3]; int ng=0; gs[ng++]=&B[k]; if(kp!=k) gs[ng++]=&B[kp]; gs[ng++]=&ALL;
            for (int i=0;i<ng;i++){ Bin *bb=gs[i];
                if (gR<bb->Rlow) bb->Rlow=gR; if (gN<bb->Nlow) bb->Nlow=gN; if (gG<bb->Glow) bb->Glow=gG; }
        }
    }
    prevp=p; Rprev=R; T1prev=T1;
}

int main(int argc, char **argv)
{
    q = argc>1 ? atoi(argv[1]) : 4;
    uint64_t X = argc>2 ? (uint64_t)atof(argv[2]) : 1000000000ULL;
    if (q!=3 && q!=4) return 1;
    uint64_t r = (uint64_t)sqrtl((ld)X)+1;
    uint8_t *small = calloc(r+1,1);
    uint64_t nsp=0; uint64_t *sp = malloc((r/2+2)*sizeof(uint64_t));
    for (uint64_t c=3;c<=r;c+=2){ if(small[c]) continue; sp[nsp++]=c; for(uint64_t m2=c*c;m2<=r;m2+=2*c) small[m2]=1; }
    for (int k=0;k<NB;k++) binit(&B[k]); binit(&ALL);
    uint64_t cps[16]; int ncp=0, icp=0;
    for (uint64_t c=100;c<=X;c*=10) cps[ncp++]=c;
    if (ncp==0 || cps[ncp-1]!=X) cps[ncp++]=X;
    const uint64_t SEG = 1u<<22;
    uint8_t *seg = malloc(SEG/16+2);
    printf("# cesaro_scan q=%d X=%llu  (long double mantissa bits: %d)\n", q,(unsigned long long)X,__LDBL_MANT_DIG__);
    event(2, q==3 ? +1 : 0);
    for (uint64_t lo=3; lo<=X; lo+=SEG){
        uint64_t hi = lo+SEG-1<X ? lo+SEG-1 : X;
        uint64_t lo_odd = lo|1ULL;
        if (lo_odd>hi) break;
        uint64_t nbits = (hi-lo_odd)/2+1;
        memset(seg,0,(nbits+7)/8);
        for (uint64_t k=0;k<nsp;k++){
            uint64_t c=sp[k]; if (c*c>hi) break;
            uint64_t s=c*c;
            if (s<lo_odd){ s=((lo_odd+c-1)/c)*c; if(!(s&1)) s+=c; }
            for (uint64_t v=s; v<=hi; v+=2*c){ uint64_t i=(v-lo_odd)>>1; seg[i>>3]|=(uint8_t)(1u<<(i&7)); }
        }
        for (uint64_t b=0;b<nbits;b++){
            if (seg[b>>3]&(1u<<(b&7))) continue;
            uint64_t p=lo_odd+2*b;
            while (icp<ncp && p>cps[icp]) { checkpoint(cps[icp]); icp++; }
            int s;
            if (q==4) s = ((p&3)==3)?+1:-1;
            else s = (p==3)?0:((p%3==2)?+1:-1);
            event(p,s);
        }
    }
    while (icp<ncp) { checkpoint(cps[icp]); icp++; }

    printf("\n# exact tests over primes after the first race prime, up to X:\n");
    printf("EXACT  primes with p*D0-D1 <= 0 (i.e. R(p)<=0): %llu  first: %llu\n",(unsigned long long)nbadR,(unsigned long long)firstbadR);
    printf("FLOAT  primes with Theta1(p) <= 0: %llu  first: %llu\n",(unsigned long long)nbadT,(unsigned long long)firstbadT);
    printf("MINR   min of R over race primes after the first: %.15Lg at p=%llu\n",Rmin5v,(unsigned long long)Rmin5);
    printf("C1     sup (y - vo(y))/sqrt(y), 0<y<=X : %.12Lf  at y->%llu^-\n",c1,(unsigned long long)c1_at);
    printf("C1B    sup (y - vo(y))/sqrt(y), 1423<=y<=X: %.12Lf  at y->%llu^-\n",c1b,(unsigned long long)c1b_at);
    printf("CUP    sup (vo(y) - y)/sqrt(y), 2<=y<=X: %.12Lf  at y=%llu\n",cup,(unsigned long long)cup_at);
    printf("C3     sup |Theta(y)|/sqrt(y), y<=X    : %.12Lf  at y=%llu\n",c3,(unsigned long long)c3_at);
    printf("THN    Theta(p)/sqrt(p) over race primes: min %.12Lf at %llu, max %.12Lf at %llu\n",ThNmin,(unsigned long long)ThNmin_at,ThNmax,(unsigned long long)ThNmax_at);
    printf("T1MIN  min of Theta1(p) over race primes after the first: %.12Lf at p=%llu ; min of Theta1(p)/p^1.5: %.12Lf at p=%llu\n",
        T1min,(unsigned long long)T1min_at,Gminall,(unsigned long long)Gminall_at);
    printf("NGAP   min over gaps [p,p'] of min(R(p),R(p'))*log(p')/sqrt(p'):  p>=29: %.12Lf (gap starting at %llu) ;  p>=11: %.12Lf (gap starting at %llu)\n",
        NgapA,(unsigned long long)NgapA_at,NgapB,(unsigned long long)NgapB_at);
    printf("\n# dyadic bins [2^k,2^(k+1)), statistics over primes p>=11.  'low' = rigorous lower bound over all real x in the bin\n");
    printf("%-3s %-10s | %-13s %-12s %-13s | %-10s %-12s %-10s %-10s %-12s | %-10s %-12s %-10s %-10s %-12s\n",
        "k","#primes","R_min","at","R_low","N_min","at","N_low","N_max","at","G_min","at","G_low","G_max","at");
    for (int k=3;k<NB;k++){ Bin *b=&B[k]; if(!b->cnt) continue;
        printf("%-3d %-10llu | %-13.6Lf %-12llu %-13.6Lf | %-10.6Lf %-12llu %-10.6Lf %-10.6Lf %-12llu | %-10.6Lf %-12llu %-10.6Lf %-10.6Lf %-12llu\n",
            k,(unsigned long long)b->cnt,b->Rmin,(unsigned long long)b->Rmin_at,b->Rlow,
            b->Nmin,(unsigned long long)b->Nmin_at,b->Nlow,b->Nmax,(unsigned long long)b->Nmax_at,
            b->Gmin,(unsigned long long)b->Gmin_at,b->Glow,b->Gmax,(unsigned long long)b->Gmax_at); }
    { Bin *b=&ALL;
        printf("ALL %-10llu | %-13.6Lf %-12llu %-13.6Lf | %-10.6Lf %-12llu %-10.6Lf %-10.6Lf %-12llu | %-10.6Lf %-12llu %-10.6Lf %-10.6Lf %-12llu\n",
            (unsigned long long)b->cnt,b->Rmin,(unsigned long long)b->Rmin_at,b->Rlow,
            b->Nmin,(unsigned long long)b->Nmin_at,b->Nlow,b->Nmax,(unsigned long long)b->Nmax_at,
            b->Gmin,(unsigned long long)b->Gmin_at,b->Glow,b->Gmax,(unsigned long long)b->Gmax_at); }
    return 0;
}
