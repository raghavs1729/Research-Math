\\ x-refined Nahm sums: sum q^{(1/2) n M n + B.n} x^{w.n} / prod (q^{d_i};q^{d_i})_{n_i}, as polynomial in x with series coefficients in t=q^(1/den)
\r nahm.gp
NHX_rec(k, cur) =
{ my(e, i, c2, bnd, prev);
  if(k==NH_r,
    e = NH_lb(NH_r)*NH_den;
    if(e<NH_P, NH_S += t^e*cur*x^(sum(j=1,NH_r,NH_w[j]*NH_x[j])); NH_cnt++);
    return);
  i = 0; c2 = cur; prev = -oo;
  while(1,
    NH_x[k+1] = i;
    bnd = NH_lb(k+1);
    if(bnd >= NH_N, if(bnd >= prev && i>0, break), NHX_rec(k+1, c2));
    prev = bnd;
    i++;
    c2 = c2/(1-NH_tq^(NH_d[k+1]*i)));
  NH_x[k+1]=0;
}
nahmx(M,B,d,w,N,den=1) =
{ my(r=#B);
  NH_r=r; NH_B=B; NH_d=d; NH_N=N; NH_den=den; NH_P=N*den; NH_tq=t^den; NH_w=w;
  NH_Mi=vector(r); NH_Myx=vector(r); NH_Mxx=vector(r); NH_By=vector(r);
  for(k=0, r-1,
    NH_Mi[k+1] = matrix(r-k,r-k,i,j,M[k+i,k+j])^(-1);
    NH_Myx[k+1] = matrix(r-k,k,i,j,M[k+i,j]);
    NH_By[k+1] = vectorv(r-k,i,B[k+i]));
  for(k=1,r, NH_Mxx[k] = matrix(k,k,i,j,M[i,j]));
  NH_S = 0; NH_x = vector(r); NH_cnt=0;
  NHX_rec(0, 1+O(t^NH_P));
  NH_S
}
