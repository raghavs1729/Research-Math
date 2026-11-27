/* psieve.c -- audit sieve (independent of ../sieve.c).
   usage: psieve q logC X u_lo u_hi u_step
   For primes p = 1 mod q up to X: A = sum log(1-1/p), th = sum log p (long double, Kahan).
   log f(x) = A + log(log(phi*th))/phi - logC  (AF (4)).  Prints
     - every k with log f > 0 (failure of (3)) for the first 20, plus the total count;
     - at checkpoints x = floor(e^u): u, x, k, G = phi*sqrt(x)*log(x)*log f(x),
       and E = (log prod_{p<=x,p=1} (1-1/p)) + log(log x)/phi  (-> log C as x -> oo).
   Build: gcc -O2 -o psieve psieve.c -lm                                            */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
typedef long double LD;
static void kadd(LD *s, LD *c, LD v) { LD y = v - *c; LD t = *s + y; *c = (t - *s) - y; *s = t; }
int main(int argc, char **argv) {
  long q = atol(argv[1]); LD logC = strtold(argv[2], 0); long long X = (long long)atof(argv[3]);
  double ulo = atof(argv[4]), uhi = atof(argv[5]), ust = atof(argv[6]);
  long phi = q - 1; /* prime q */
  long long sq = (long long)sqrtl((LD)X) + 2;
  char *small = calloc(sq + 1, 1); long long *sp = malloc(sizeof(long long) * (sq + 1)); long nsp = 0;
  for (long long i = 2; i <= sq; i++) if (!small[i]) { sp[nsp++] = i; for (long long j = i * i; j <= sq; j += i) small[j] = 1; }
  const long long SEG = 1 << 24; char *seg = malloc(SEG);
  LD A = 0, Ac = 0, th = 0, thc = 0; long long k = 0, nfail = 0;
  double u = ulo; long long nextx = (long long)floorl(expl((LD)u));
  for (long long lo = 0; lo <= X; lo += SEG) {
    long long hi = lo + SEG - 1; if (hi > X) hi = X;
    memset(seg, 0, SEG);
    for (long i = 0; i < nsp; i++) { long long p = sp[i]; if (p * p > hi) break;
      long long st = ((lo + p - 1) / p) * p; if (st < p * p) st = p * p;
      for (long long j = st; j <= hi; j += p) seg[j - lo] = 1; }
    for (long long n = (lo < 2 ? 2 : lo); n <= hi; n++) {
      while (u <= uhi && nextx < n) { /* checkpoint x in [p_k, n) */
        LD x = (LD)nextx;
        if (k > 0) { LD lf = A + logl(logl(phi * th)) / phi - logC;
          printf("CP u=%.2f x=%lld k=%lld G=%.6Lf E=%.12Lf\n", u, nextx, k, phi * sqrtl(x) * logl(x) * lf,
                 A + logl(logl(x)) / phi); }
        u += ust; nextx = (long long)floorl(expl((LD)u)); }
      if (seg[n - lo]) continue;
      if (n % q != 1) continue;
      k++; kadd(&A, &Ac, log1pl(-1.0L / n)); kadd(&th, &thc, logl((LD)n));
      LD lf = A + logl(logl(phi * th)) / phi - logC;
      if (lf > 0) { nfail++; if (nfail <= 20) printf("FAIL k=%lld p=%lld logf=%.3Le\n", k, n, lf); }
    }
  }
  printf("END X=%lld k=%lld nfail=%lld E(X)=%.12Lf\n", X, k, nfail, A + logl(logl((LD)X)) / phi);
  return 0;
}
