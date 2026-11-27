M(n)=matrix(n,n,i,j,lcm(i,j)+gcd(i,j));
for(n=1,30, d=matdet(M(n)); print(n," ",if(d==0,0,factor(d,10^5)~)));
quit
