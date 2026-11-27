\\ Generic (generalized) Nahm sums as power series in t = q^(1/den), exact integer coefficients.
\\ nahm(M,B,d,N,den): sum_{n in N^r} q^{(1/2) n~ M n + B.n} / prod_i (q^{d_i};q^{d_i})_{n_i}
\\   M symmetric positive definite rational r x r, B rational vector, d positive integers;
\\   result is a series in t with q = t^den, truncated O(t^(N*den)) (i.e. O(q^N)).
\\ Pruning: exact unconstrained minimum of the quadratic form over the remaining variables
\\ (a valid lower bound for the exponent), which is convex in each coordinate.
\\ Run with: gp -q -s 4G  (state is kept in global variables NH_*)
qpoch(a,q,n) = prod(k=0,n-1,1-a*q^k);
qinf(a,q,P) = { my(v=valuation(q,t)); prod(k=0, P\v+1, 1 - a*q^k + O(t^P)); }
\\ J-type products in t: jp(a,m,P) = (t^a;t^m)_inf
jp(a,m,P) = qinf(t^a,t^m,P);

NH_lb(k) =
{ my(xv, w, val);
  if(k==0, val = 0,
    xv = vectorv(k,i,NH_x[i]);
    val = (xv~*NH_Mxx[k]*xv)/2 + sum(i=1,k,NH_B[i]*NH_x[i]));
  if(k<NH_r,
    w = if(k==0, NH_By[1], NH_Myx[k+1]*xv + NH_By[k+1]);
    val -= (w~*NH_Mi[k+1]*w)/2);
  val
}
NH_rec(k, cur) =
{ my(e, i, c2, bnd, prev);
  if(k==NH_r,
    e = NH_lb(NH_r)*NH_den;
    if(denominator(e)!=1, error("exponent not in (1/den)Z: ",e));
    if(e<NH_P, NH_S += t^e*cur; NH_cnt++);
    return);
  i = 0; c2 = cur; prev = -oo;
  while(1,
    NH_x[k+1] = i;
    bnd = NH_lb(k+1);
    if(bnd >= NH_N, if(bnd >= prev && i>0, break), NH_rec(k+1, c2));
    prev = bnd;
    i++;
    c2 = c2/(1-NH_tq^(NH_d[k+1]*i)));
  NH_x[k+1]=0;
}
nahm(M,B,d,N,den=1) =
{ my(r=#B);
  NH_r=r; NH_B=B; NH_d=d; NH_N=N; NH_den=den; NH_P=N*den; NH_tq=t^den;
  NH_Mi=vector(r); NH_Myx=vector(r); NH_Mxx=vector(r); NH_By=vector(r);
  for(k=0, r-1,
    NH_Mi[k+1] = matrix(r-k,r-k,i,j,M[k+i,k+j])^(-1);
    NH_Myx[k+1] = matrix(r-k,k,i,j,M[k+i,j]);
    NH_By[k+1] = vectorv(r-k,i,B[k+i]));
  for(k=1,r, NH_Mxx[k] = matrix(k,k,i,j,M[i,j]));
  NH_S = O(t^NH_P); NH_x = vector(r); NH_cnt=0;
  NH_rec(0, 1+O(t^NH_P));
  NH_S
}
