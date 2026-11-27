// Exploratory (uncertified, double precision) scan for deep minima of
//   B(u) = sum_{gamma>0} w_gamma cos(gamma u),  w = 2/(1/4+gamma^2)
// zeros file: lines "gamma" ; sorted internally by weight.
// usage: scan zeros.txt U0 U1 M thr1 thr2
//   windows around odd multiples of pi/gamma_1 (dominant zero), half-width HW, step DU
//   stage 1: top-M zeros; if min < thr1: full sum; print if full < thr2
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define NMAX 200000
static double g[NMAX], w[NMAX], cr[NMAX], ci[NMAX], sr[NMAX], si[NMAX];
int cmp(const void*a,const void*b){double x=*(double*)a,y=*(double*)b;return x<y?-1:x>y;}
int main(int argc,char**argv){
  FILE*f=fopen(argv[1],"r"); int n=0; double x;
  while(fscanf(f,"%lf",&x)==1) g[n++]=x;
  qsort(g,n,sizeof(double),cmp);
  for(int i=0;i<n;i++) w[i]=2.0/(0.25+g[i]*g[i]);
  double U0=atof(argv[2]),U1=atof(argv[3]); int M=atoi(argv[4]); double thr1=atof(argv[5]),thr2=atof(argv[6]);
  if(M>n)M=n;
  const double DU=0.02, HW=0.6; int NS=(int)(2*HW/DU)+1;
  for(int j=0;j<n;j++){sr[j]=cos(g[j]*DU);si[j]=sin(g[j]*DU);}
  double P=2*M_PI/g[0]; long k0=(long)ceil((U0-M_PI/g[0])/P), k1=(long)floor((U1-M_PI/g[0])/P);
  double best=0; 
  static double acc[256];
  for(long k=k0;k<=k1;k++){
    double c=M_PI/g[0]+k*P, u0=c-HW;
    for(int s=0;s<NS;s++)acc[s]=0;
    int alive=1;
    for(int j=0;j<M;j++){
      double a=cos(g[j]*u0), b=sin(g[j]*u0), wr=w[j], r=sr[j], i=si[j];
      for(int s=0;s<NS;s++){acc[s]+=wr*a; double t=a*r-b*i; b=a*i+b*r; a=t;}
    }
    double m=0; int sm=0; for(int s=0;s<NS;s++) if(acc[s]<m){m=acc[s];sm=s;}
    if(m<thr1){
      // full evaluation, fine step around sm
      double bestf=0,bu=0;
      for(int s=-10;s<=10;s++){ double u=u0+sm*DU+s*DU/10, t=0; for(int j=0;j<n;j++) t+=w[j]*cos(g[j]*u); if(t<bestf){bestf=t;bu=u;} }
      if(bestf<thr2){ printf("%.6f %.6f %.6f\n",bu,bestf,m); fflush(stdout);} 
      if(bestf<best){best=bestf; fprintf(stderr,"best %.6f at %.6f\n",best,bu);} 
    }
  }
  return 0;
}
