\r nahm.gp
\r cartan.gp
N=120;
C=cartanE(8); print(matdet(C)); M=C^(-1); print(M);
f=nahm(M,vector(8),vector(8,i,1),N); print(NH_cnt);
g=1/prod(a=1,10, if(a!=5&&a!=6, jp(a,11,N), 1));
print(f-g);
C=cartanE(6); print(matdet(C)); M=C^(-1); print(M);
f6=nahm(M,vector(6),vector(6,i,1),40,3); print(NH_cnt); print(f6);
quit
