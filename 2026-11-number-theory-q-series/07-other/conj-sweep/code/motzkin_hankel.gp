\\ independent check (PARI/GP): Hankel determinants of M(x)^r, M = Motzkin generating function
default(seriesprecision, 700);
M = (1 - x - sqrt(1 - 2*x - 3*x^2 + O(x^600)))/(2*x^2);
H(r,n) = { my(F = M^r); if(n==0, 1, matdet(matrix(n,n,i,j, polcoeff(F, i+j-2)))) };
forstep(r=3, 30, 3, print("r=", r, ": H_{rn}, n=0..4: ", vector(5, n, H(r, r*(n-1))), "  H_{rn+1}: ", vector(5, n, H(r, r*(n-1)+1))));
print("r=4: H_0..H_13 = ", vector(14, n, H(4, n-1)));
print("r=5: H_0..H_16 = ", vector(17, n, H(5, n-1)));
quit
