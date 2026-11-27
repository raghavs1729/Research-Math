\r nahm.gp
N=200;
f=nahm(Mat(2),[0],[1],N); print(f - 1/(jp(1,5,N)*jp(4,5,N)));
f=nahm([4,2;2,2],[0,0],[1,1],N); print(f - 1/(jp(1,7,N)*jp(2,7,N)*jp(5,7,N)*jp(6,7,N)));
f=nahm([2,3;3,6],[0,0],[1,3],N); print(f - 1/(jp(1,9,N)*jp(3,9,N)*jp(6,9,N)*jp(8,9,N)));
\\ negative off-diagonal: Zagier ex (A2,T1) A=[2,-1;-1,2]: sum = (A1)_1 char = theta/eta... just print
f=nahm([2,-1;-1,2],[0,0],[1,1],60); print(f);
quit
