\r nahm.gp
\r cartan.gp
N=70;
C=cartanE(8); A=C^(-1); T=cartanT(4)^(-1);
chk(b) = { my(B=[b[1], b[3], b[1]+b[3], b[2]+b[4], b[1]+b[2]+b[4], b[1]+b[4], b[2], b[4]], f, g);
  f=nahm(A,B,vector(8,i,1),N); g=nahm(2*T,b,vector(4,i,1),N); print(b, " ", f-g==0, " ", if(f-g!=0, f-g, "")); };
chk([0,0,0,0]); chk([1,0,0,0]); chk([0,1,0,0]); chk([0,0,1,0]); chk([0,0,0,1]); chk([2,0,1,3]); chk([1,1,1,1]); chk([0,3,1,0]);chk([-1,0,0,0]);chk([0,-1,0,1]);
quit
