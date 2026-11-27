\\ Independent checks with PARI/GP (different library, different algorithms; floating point, 38 digits; no rigor claimed).
\\ Run:  gp -q < independent_check.gp
default(realprecision, 38); default(parisize, 2000000000);
print("PARI/GP version ", version());
mf = mfinit([1,12],1); f = mfbasis(mf)[1]; L = lfunmf(mf, f);
LD(s) = lfun(L, s+11/2, 1)/lfun(L, s+11/2);          \\ (L'/L)(s, Delta), analytic normalisation
LamD(s) = -log(2*Pi) + psi(s+11/2) + LD(s);          \\ Lambda'/Lambda(s)
print("== (I1) constants of L(s, Delta)");
print("sigma0^2 = 2 Lambda'/Lambda(1) = ", 2*LamD(1));
print("a0 = L'/L(0) = ", LD(0), "   a_{-1} = L'/L(-1) = ", LD(-1), "   L'/L(2) = ", LD(2), "  L(1) = ", lfun(L, 13/2));
zz = lfunzeros(L, 300);
print(#zz, " zeros with 0 < gamma < 300; first ", zz[1], ";  sum_{|gamma|<300} 1/|rho|^2 = ", sum(i=1,#zz, 2/(1/4+zz[i]^2)), ";  sum 1/|rho(rho+1)| = ", sum(i=1,#zz, 2/sqrt((1/4+zz[i]^2)*(9/4+zz[i]^2))));
print("== (I7a) residue bookkeeping for psi_1(t, Delta) at t = 1:  0 = -[Lam'/Lam(1) - Lam'/Lam(2)] - a0 + a_{-1} - 2/9 :  ", -(LamD(1)-LamD(2)) - LD(0) + LD(-1) - 2/9);
\\ direct comparison psi_1(t) with zeros to 300
tv = mfcoefs(f, 100000);     \\ tv[n+1] = tau(n)
lam(p) = tv[p+1]/p^(11/2.);
{psi1(t) = my(s = 0., q, a, b, pk, ak, akm, c);
  forprime(p = 2, floor(t),
    c = lam(p); akm = 2.; ak = c; pk = p;       \\ c_k = alpha^k + beta^k : c_0 = 2, c_1 = lam, c_{k+1} = lam c_k - c_{k-1}
    while(pk <= t, s += ak*log(p)*(t - pk); [akm, ak] = [ak, c*ak - akm]; pk *= p));
  s};
{psi1f(t) = my(s = 0.);
  for(i = 1, #zz, my(r = 1/2 + I*zz[i]); s -= 2*real(t^(r+1)/(r*(r+1))));
  s - LD(0)*t + LD(-1) - sum(n = 0, 2000, t^(-9/2-n)/((n+11/2)*(n+9/2)))};
tailb = 2*(0.11225214047 - sum(i=1,#zz, 2/(1/4+zz[i]^2)))/2;   \\ sum_{|gamma|>300} 1/|rho|^2
print("== (I7b) psi_1 direct - formula(277 zeros), allowed |diff| <= t^{3/2} * ", tailb);
forstep(i = 1, 7, 1, my(t = [1, 2.5, 5.5, 10, 30, 100, 1000][i]); print("  t = ", t, ":  direct ", psi1(t), "  diff ", psi1(t) - psi1f(t), "   bound ", t^1.5*tailb));
print("== (I4) Sym^2 Delta");
an = direuler(p = 2, 6000, my(t = ramanujantau(p)); 1/(1 - (t^2-p^11)*X + p^11*(t^2-p^11)*X^2 - p^33*X^3));
S2 = lfuncreate([an, 0, [-10,0,1], 23, 1, 1]);
print("functional equation check (bits): ", lfuncheckfeq(S2));
LS(s) = lfun(S2, s+11, 1)/lfun(S2, s+11);
gS(s) = -3/2*log(Pi) + (psi((s+1)/2) + psi((s+11)/2) + psi((s+12)/2))/2;
LamS(s) = gS(s) + LS(s);
print("L_S(1/2) = ", lfun(S2, 23/2), "  L_S(1) = ", lfun(S2, 12), "  L_S(2) = ", lfun(S2, 13));
print("L_S'/L_S(2) = ", LS(2), "   Lambda_S'/Lambda_S(2) = ", LamS(2), "   L_S'/L_S(0) = ", LS(0), "   L_S'/L_S(-2) = ", LS(-2), "  L_S'/L_S(3) = ", LS(3));
print("sum 1/|rho|^2 (Sym^2) = 2 Lambda_S'/Lambda_S(1) = ", 2*LamS(1));
\\ trivial zeros of L_S: w = -1,-3,..,-9 (simple), -12,-14,.. (simple), -11,-13,... (double); residues at s0 = w/2 <= -3/2
triv(t) = sum(j = 1, 4, t^(-(2*j+1)/2+1)/(2*((2*j+1)/2)*((2*j+1)/2-1))) + sum(w = 11, 4000, if(w%2, 2, 1)*t^(-w/2+1)/(2*(w/2)*(w/2-1)));
print("== (I7c) residue bookkeeping for J_S(t) at t = 1:  0 = -[Lam_S'/Lam_S(1) - Lam_S'/Lam_S(3)] - L_S'/L_S(0) + 2 + L_S'/L_S(-2) - triv(1) :  ", -(LamS(1) - LamS(3)) - LS(0) + 2 + LS(-2) - triv(1), "   (truncation of triv at w = 4000: ~ 1/4000)");
w = lfunzeros(S2, 60);
print(#w, " zeros of L_S with 0 < gamma < 60");
print(w);
print("== (I2) small range, tau from mfcoefs (PARI) and ramanujantau");
forprime(p = 2, 20000, if(ramanujantau(p) != tv[p+1], error("tau mismatch")));
print("ramanujantau(p) == mfcoefs for primes p <= 20000: ok");
{chk(N) = my(Th = 0., Mm = 0., D0 = 0., D1 = 0., T1, R, mG = 10., mGp, bad = 0);
  forprime(p = 2, N, my(l = lam(p), lg = log(p)); Th += l*lg; Mm += l*lg*p; D0 += l; D1 += l*p;
     T1 = p*Th - Mm; R = D0 - D1/p;
     if(p >= 5 && T1 <= 0, bad++); if(p >= 7 && R <= 0, bad++);
     if(p >= 10^4 && T1/p^1.5 < mG, mG = T1/p^1.5; mGp = p);
     if(p == 97 || p == 997 || p == 99991, print("  p = ", p, ": Theta = ", Th, "  Theta_1 = ", T1, "  D0 = ", D0, "  D1 = ", D1, "  A1 = ", D0 - Th/lg - T1/(p*lg^2), "  A2 = ", -D1 + (p*Th - T1)/lg + T1/lg^2)));
  print("  N = ", N, ": failures ", bad, ";  min Theta_1(p)/p^1.5 over primes in [1e4, N] = ", mG, " at p = ", mGp)};
chk(100000);
quit
