from mpmath import mp, mpf, hyp2f1, quad, gamma, pi, sqrt, pslq, zeta, log, catalan
mp.dps=40
exec(open('z5_values.py').read())
v1=lambda y: hyp2f1(mpf(1)/4,mpf(1)/4,mpf(1)/2,y)
v2=lambda y: sqrt(y)*hyp2f1(mpf(3)/4,mpf(3)/4,mpf(3)/2,y)
I=quad(lambda y: v1(y)*v2(y)**3,[0,mpf(1)/2,mpf(9)/10,mpf(99)/100,1])
print(I)
print(3*gamma(mpf(1)/4)**4/(pi**4*I), U1)
