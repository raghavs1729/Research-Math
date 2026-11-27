// scan.c -- exploratory (UNCERTIFIED, double precision) search for deep minima of
//   S_a(u) = sum_j w_j cos(gamma_j u - phase_j),  w_j = 2/(1/4+gamma_j^2)
// input file: lines "gamma phase" (see make_scan_input.py).  Zeros are sorted by weight.
// usage: scan zeros.txt U0 U1 M thr1 thr2
//   windows of half-width HW around the minima of the heaviest term, step DU;
//   stage 1: the M heaviest zeros; if the window minimum is < thr1: all zeros on a
//   finer grid; print "u S_full S_stage1" if S_full < thr2.  Progress records on stderr.
// Candidates found here are then certified independently in Arb (certify.py).
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define NMAX 400000
static double g[NMAX], ph[NMAX], w[NMAX], sr[NMAX], si[NMAX];
static int idx[NMAX];
int cmp(const void*a,const void*b){double x=g[*(int*)a],y=g[*(int*)b];return x<y?-1:x>y;}
int main(int argc,char**argv){
  FILE*f=fopen(argv[1],"r"); int n=0; double x,y; static double g0[NMAX],p0[NMAX];
  while(fscanf(f,"%lf %lf",&x,&y)==2){ g0[n]=x;p0[n]=y;g[n]=x;idx[n]=n;n++; }
  qsort(idx,n,sizeof(int),cmp);
  for(int i=0;i<n;i++){ g[i]=g0[idx[i]]; ph[i]=p0[idx[i]]; }
  for(int i=0;i<n;i++) w[i]=2.0/(0.25+g[i]*g[i]);
  double U0=atof(argv[2]),U1=atof(argv[3]); int M=atoi(argv[4]); double thr1=atof(argv[5]),thr2=atof(argv[6]);
  if(M>n)M=n;
  const double DU=0.02, HW=0.6; int NS=(int)(2*HW/DU)+1;
  for(int j=0;j<n;j++){sr[j]=cos(g[j]*DU);si[j]=sin(g[j]*DU);}
  double P=2*M_PI/g[0], c0=(M_PI+ph[0])/g[0];
  long k0=(long)ceil((U0-c0)/P), k1=(long)floor((U1-c0)/P);
  double best=0;
  static double acc[256];
  for(long k=k0;k<=k1;k++){
    double c=c0+k*P, u0=c-HW;
    for(int s=0;s<NS;s++)acc[s]=0;
    for(int j=0;j<M;j++){
      double a=cos(g[j]*u0-ph[j]), b=sin(g[j]*u0-ph[j]), wr=w[j], r=sr[j], i=si[j];
      for(int s=0;s<NS;s++){acc[s]+=wr*a; double t=a*r-b*i; b=a*i+b*r; a=t;}
    }
    double m=0; int sm=0; for(int s=0;s<NS;s++) if(acc[s]<m){m=acc[s];sm=s;}
    if(m<thr1){
      double bestf=0,bu=0;
      for(int s=-10;s<=10;s++){ double u=u0+sm*DU+s*DU/10, t=0; for(int j=0;j<n;j++) t+=w[j]*cos(g[j]*u-ph[j]); if(t<bestf){bestf=t;bu=u;} }
      if(bestf<thr2){ printf("%.6f %.6f %.6f\n",bu,bestf,m); fflush(stdout);}
      if(bestf<best){best=bestf; fprintf(stderr,"best %.6f at %.6f\n",best,bu);}
    }
  }
  return 0;
}
