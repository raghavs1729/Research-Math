# Sanity tests of radu.py on known true and false congruences.
from radu import *
out='../cert/radu_sanity.txt'; open(out,'w').write('Sanity tests of the Radu engine\n\n')
P1={1:-1}
tests=[('p(5n+4) mod 5 (Ramanujan, true)',1,P1,5,4,5,True),('p(7n+5) mod 7 (true)',1,P1,7,5,7,True),
       ('p(11n+6) mod 11 (true)',1,P1,11,6,11,True),('p(25n+24) mod 25 (true)',1,P1,25,24,25,True),
       ('p(5n+3) mod 5 (false)',1,P1,5,3,5,False),('p(13n+6) mod 13 (false)',1,P1,13,6,13,False),
       ('A_6(12n+10) mod 3^6 (false: 3^5 is exact)',6,{1:-6,2:-6,3:-6,6:6},12,10,3**6,False),
       ('PD(3n+2) mod 3 (Andrews-Lewis-Lovejoy, true)',6,{1:-1,2:-1,3:-1,6:1},3,2,3,True),
       ('overpartitions pbar(40n+35) mod 40 (Chen-Xia, true)',2,{1:-2,2:1},40,35,40,True)]
for name,M,r,m,t,u,expect in tests:
    ok,info=certificate(name,M,r,m,t,u,out=out)
    assert ok==expect,(name,ok)
print('all sanity tests behaved as expected')
