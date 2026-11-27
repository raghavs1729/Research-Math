\\ Test of the "symmetric chain self-fusion identity" (***):
\\  sum over distinct permutations N' of N of  W(N') = prod_i 1/(q)_{N_(i)-N_(i+1)}  (N sorted decreasingly, N_(k+1)=0)
\\  W(N) = sum_{r_1>=..>=r_k>=0, 0<=s_1<=..<=s_k, r_i+2 s_i=N_i} prod_i q^{r_i(r_i-1)/2}/((q)_{r_i-r_{i+1}} (q^2;q^2)_{s_i-s_{i-1}})
qp(n) = prod(j=1,n,1-q^j);
qp2(n) = prod(j=1,n,1-q^(2*j));
W(N) = Wrec(N, 1, 0, -1, 1);
\\ recursion over i=1..k: choose s_i >= s_{i-1}, r_i = N_i - 2 s_i >= 0 and r_i <= r_{i-1} (if i>1); factor for step i:
\\   q^{r_i(r_i-1)/2} / (q^2;q^2)_{s_i-s_{i-1}}, and 1/(q)_{r_{i-1}-r_i} for i>1; at the end 1/(q)_{r_k}
Wrec(N, i, sprev, rprev, acc) = { my(k=#N, tot=0, r);
  if(i>k, return(acc/qp(rprev)));
  for(s=sprev, N[i]\2,
    r = N[i]-2*s;
    if(i>1 && r>rprev, next);
    tot += Wrec(N, i+1, s, r, acc * q^(r*(r-1)/2) / qp2(s-sprev) / if(i>1, qp(rprev-r), 1)));
  tot }
perms(N) = { my(S=Set(vector((#N)!, j, numtoperm(#N,j-1))), L=List()); 
  for(j=1,#S, my(p=S[j], v=vector(#N,i,N[p[i]])); listput(L,v)); Set(Vec(L)) }
chk(N) = { my(Ns=vecsort(N,,4), P=perms(Ns), lhs, rhs);
  lhs = sum(j=1,#P, W(P[j]));
  rhs = 1/prod(i=1,#Ns, qp(Ns[i]-if(i<#Ns,Ns[i+1],0)));
  lhs==rhs }
bad=0; cnt=0;
forvec(N=vector(2,i,[0,9]), cnt++; if(!chk(N), bad++; print("FAIL ",N)), 1);
print("k=2: tested ",cnt," multisets, failures ",bad);
bad=0; cnt=0;
forvec(N=vector(3,i,[0,7]), cnt++; if(!chk(N), bad++; if(bad<5,print("FAIL ",N))), 1);
print("k=3: tested ",cnt," multisets, failures ",bad);
bad=0; cnt=0;
forvec(N=vector(4,i,[0,4]), cnt++; if(!chk(N), bad++; if(bad<5,print("FAIL ",N))), 1);
print("k=4: tested ",cnt," multisets, failures ",bad);
quit
