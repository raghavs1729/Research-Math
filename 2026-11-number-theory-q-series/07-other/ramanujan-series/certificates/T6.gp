\\ Theorem 6: two-parameter extension (a,b) of the pair of Theorem 1; u = k+3n+2a+b
u = k+3*n+2*a+b;
sg = -1/3*(u+1)*(k+b-a+1/2)/((u+3/2)*(k+n+b+1));
rho   = -1/3*3/4*(u+1)*(u+2)*(u+3)*(2*n+2*a+3/2)*(2*n+2*a+5/2)*(n+a+1)^2*(-(a+n+1/2))/((2*n+2*a+1)*(2*n+2*a+2)*(u+3/2)*(u+5/2)*(u+7/2)*(k+n+b+1)*(n+a+3/4)*(n+a+5/4));
R     = 3*(u+2)^2/(4*(u+2)^2-1);
print("T6 WZ: ", rho - 1 - (subst(R,k,k+1)*sg - R) == 0);
