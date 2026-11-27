\r nahmx.gp
\r cartan.gp
N=80;
C=cartanE(8); A=C^(-1);
eps=[-1,1,-1,1,-1,1,-1,1]~; w=(A*eps)~; print(w);
f=nahmx(A,vector(8),vector(8,i,1),w,N);
T=cartanT(4)^(-1);
g=nahmx(2*T,vector(4),vector(4,i,1),[1,2,3,4],N);
dd=f-g; 
for(j=0,poldegree(dd,x), c=polcoeff(dd,j,x); c=truncate(c+O(t^N)); if(c!=0, print("x^",j,": ",c)));
print("done");
quit
