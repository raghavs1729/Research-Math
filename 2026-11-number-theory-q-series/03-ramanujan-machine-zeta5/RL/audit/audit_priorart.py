# Prior-art checks: (a) equivalence C[(n+k)/(2n(n+1))] = 1 + K((n+k)/2 / (n+1)) and
# the classical b-contiguous Kummer CF (DLMF 13.3.2 + Kummer 13.2.39) at z=-1/2, a=3/2-k, b=5/2;
# (b) classical erf CF for R2 and relation of its convergents to those of C[1/(2n)].
from mpmath import mp, mpf, hyp1f1, sqrt, pi, e, erf, nstr, fraction
from fractions import Fraction as F
mp.dps = 50
R2 = sqrt(pi*e/2)*erf(1/sqrt(2))
def convergents(b0, a, b, N):  # b0 + K(a(n)/b(n)), exact Fractions
    P0,P1,Q0,Q1 = F(1),F(b0),F(0),F(1); out=[P1/Q1]
    for n in range(1,N+1):
        P0,P1 = P1, b(n)*P1 + a(n)*P0; Q0,Q1 = Q1, b(n)*Q1 + a(n)*Q0; out.append(P1/Q1)
    return out
for k in [0,1,3,F(1,3)]:
    c1 = convergents(1, lambda n: F(n+k,2*n*(n+1)), lambda n: F(1), 40)
    c2 = convergents(1, lambda n: F(n+k,2), lambda n: F(n+1), 40)
    kk = mpf(k.numerator)/k.denominator if isinstance(k,F) else mpf(k)
    a = mpf(3)/2-kk; z = -mpf(1)/2; b=mpf(5)/2
    # classical: b M(a;b;z)/M(a;b+1;z) = b+z + (b+1-a)(-z)/(b+1+z + (b+2-a)(-z)/(b+2+z+...)) ; here X = that
    X = b*hyp1f1(a,b,z)/hyp1f1(a,b+1,z)
    print("k=",k," convergents identical:", c1==c2, "; 1+(k+1)/(2X) - C:", nstr(1+(kk+1)/(2*X)-mpf(c2[-1].numerator)/c2[-1].denominator,3),
          "(40 terms)")
# erf CF: R2 = 1/(1 - 1/(3 + 2/(5 - 3/(7 + ...))))
ce = convergents(0, lambda n: F(1) if n==1 else F((-1)**(n-1)*(n-1)), lambda n: F(2*n-1), 30)
print("erf-CF conv30 - R2:", nstr(mpf(ce[-1].numerator)/ce[-1].denominator - R2,3))
c1 = convergents(1, lambda n: F(1,2*n), lambda n: F(1), 30)
print("C[1/(2n)] convergents:", [str(x) for x in c1[:7]])
print("erf-CF   convergents:", [str(x) for x in ce[:10]])
print("common:", sorted(set(c1)&set(ce), key=float)[:8])
# (c) Perron-type CF (linear denominators) in the (a+1,b+1) direction:
#   M(a+1;b+1;z)/M(a;b;z) = b/(b - z + (a+1)z/(b+1-z + (a+2)z/(b+2-z + ...)))
from mpmath import mpc
def perron(a,b,z,N=400):
    t = b+N-z
    for n in range(N-1,0,-1):
        t = b+n-z + (a+n+1)*z/t
    return b/(b - z + (a+1)*z/t)
for (a,b,z) in [(mpf('0.7'),mpf('1.3'),mpf('0.45')),(mpf('2.2'),mpf('3.9'),mpf('1.7')),(mpf('-0.4'),mpf('2.5'),mpf('0.5'))]:
    print("Perron check a,b,z=",a,b,z, nstr(perron(a,b,z)-hyp1f1(a+1,b+1,z)/hyp1f1(a,b,z),3))
