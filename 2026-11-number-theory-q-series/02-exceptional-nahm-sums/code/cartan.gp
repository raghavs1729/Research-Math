\\ Cartan matrices (Bourbaki-like labellings as in Sun-Wang appendix where it matters little)
cartanA(r) = matrix(r,r,i,j, if(i==j,2, if(abs(i-j)==1,-1,0)));
cartanT(r) = { my(C=cartanA(r)); C[r,r]=1; C };
cartanD(r) = { my(C=cartanA(r)); C[r-1,r]=0; C[r,r-1]=0; C[r-2,r]=-1; C[r,r-2]=-1; C };
\\ E_n: chain 1-2-...-(n-1) with node n attached to node n-3  (E8: chain of 7, node 8 attached to node 5)
cartanE(n) = { my(C=cartanA(n)); C[n-1,n]=0; C[n,n-1]=0; C[n-3,n]=-1; C[n,n-3]=-1; C };
