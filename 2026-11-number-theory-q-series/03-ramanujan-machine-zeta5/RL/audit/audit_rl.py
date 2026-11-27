# Independent audit checks for THEOREM.md (RL). python3 + mpmath + sympy.
from mpmath import mp, mpf, hyp1f1, sqrt, pi, e, erf, gamma, nstr, pslq, fsum, factorial
mp.dps = 60
R2 = sqrt(pi*e/2)*erf(1/sqrt(2))
M = lambda a,b: hyp1f1(a,b,mpf(1)/2)

def C(f, N):
    # backward evaluation of 1 + f1/(1+f2/(1+...+fN))
    t = mpf(1)
    for n in range(N,0,-1):
        t = 1 + f(n)/t
    return t

def Cconv(f, N0=400):
    a, b = C(f,N0), C(f,2*N0)
    return b, abs(a-b)

print("R2 =", nstr(R2,40))
print("R2 vs sum 1/(2j+1)!!:", nstr(R2 - fsum(1/mp.fac2(2*j+1) for j in range(200)),5),
      " R2 vs M(1;3/2;1/2):", nstr(R2-M(1,mpf(3)/2),5))

# Theorem formula, second set
print("\n== second set: C[(n+k)/(2n(n+1))] vs 1+((k+1)/5) M(k+2;7/2)/M(k+1;5/2)")
for k in [0,1,2,3,4,5,6,mpf('0.3'),mpf('1.5'),mpf('2.5'),mpf('-0.5'),mpf('-0.9'),mpf(7)/3,mpf(10)]:
    k = mpf(k)
    val, err = Cconv(lambda n: (n+k)/(2*n*(n+1)))
    rhs = 1 + (k+1)/5*M(k+2,mpf(7)/2)/M(k+1,mpf(5)/2)
    print(f"k={nstr(k,4):>6}  C={nstr(val,35)}  |C-rhs|={nstr(abs(val-rhs),3)}  tailerr={nstr(err,3)}")

closed = {0:(0,1,2,-2),1:(1,0,0,1),2:(2,1,1,1),3:(15,9,8,6),4:(76,48,39,27),5:(955,615,464,312)}
print("\n== closed forms (aR2+b)/(cR2+d)")
for k,(a,b,c,d) in closed.items():
    val,_ = Cconv(lambda n: mpf(n+k)/(2*n*(n+1)))
    print(k, nstr(val,35), " diff:", nstr(val-(a*R2+b)/(c*R2+d),3))
# k=6 via PSLQ
for k in [6]:
    val,_ = Cconv(lambda n: mpf(n+k)/(2*n*(n+1)),600)
    mp.dps=45; print("k=6 PSLQ [C*R2, C, R2, 1]:", pslq([val*R2,val,R2,1],maxcoeff=10**8,maxsteps=10**6)); mp.dps=60

# Section 3 listed values
h=mpf(1)/2
print("\n== section 3 values")
chk = [("mu2",M(2,1.5),R2+h),("mu3",M(3,1.5),(5*R2+3)/4),("mu4",M(4,1.5),(19*R2+12)/12),
 ("M(1;5/2)",M(1,2.5),3*(R2-1)),("M(2;5/2)",M(2,2.5),mpf(3)/2),("M(3;5/2)",M(3,2.5),mpf(3)/4*(R2+1)),
 ("M(2;7/2)",M(2,3.5),mpf(15)/2*(3-2*R2)),("M(3;7/2)",M(3,3.5),mpf(15)/4*(R2-1)),("M(4;7/2)",M(4,3.5),5*R2/4)]
for name,a,b in chk: print(f"{name:10s} {nstr(a,25)}  diff {nstr(a-b,3)}")
# contiguous relations DLMF 13.3.1 and 13.3.4 at random parameters
from mpmath import mpc
for (a,b,z) in [(mpf('0.7'),mpf('1.3'),mpf('0.45')),(mpf('2.2'),mpf('3.9'),mpf('-1.1')),(mpc(1,2),mpf('0.6'),mpf('2.0'))]:
    r1 = (b-a)*hyp1f1(a-1,b,z)+(2*a-b+z)*hyp1f1(a,b,z)-a*hyp1f1(a+1,b,z)
    r4 = b*hyp1f1(a,b,z)-b*hyp1f1(a-1,b,z)-z*hyp1f1(a,b+1,z)
    print("13.3.1 resid", nstr(abs(r1),3), " 13.3.4 resid", nstr(abs(r4),3))

# Lemma 1 / minimal solution: J_n, m_n, recurrence, m_{-1}, m_0
print("\n== Lemma 1 / minimal solution")
for k in [mpf(0),mpf('0.3'),mpf(5),mpf('-0.7')]:
    J = lambda n: gamma(n+k+2)/gamma(n+mpf(7)/2)*M(n+k+2,n+mpf(7)/2)
    Jser = lambda n: fsum(gamma(n+k+2+j)/(2**j*factorial(j)*gamma(n+j+mpf(7)/2)) for j in range(300))
    mm = lambda n: (-h)**n*J(n)/factorial(n+1)
    res = max(abs(J(n)+2*(n+1)*J(n-1)-2*(n+k)*J(n-2))/abs(J(n-2)) for n in range(1,15))
    rec = max(abs(mm(n)-mm(n-1)-(n+k)/(2*n*(n+1))*mm(n-2))/abs(mm(n-2)) for n in range(1,15))
    print(f"k={nstr(k,3)} J series-vs-closed {nstr(abs(Jser(-1)-J(-1)),3)}; 3-term rel resid {nstr(res,3)}; m-rec resid {nstr(rec,3)}; m_-1+2J_-1={nstr(mm(-1)+2*J(-1),3)}; m_20={nstr(mm(20),3)}")
    # P,Q check of m_n = m_-1 P_n + (m_0-m_-1) Q_n
    f = lambda n: (n+k)/(2*n*(n+1))
    P=[mpf(1),mpf(1)]; Q=[mpf(0),mpf(1)]
    for n in range(1,12): P.append(P[-1]+f(n)*P[-2]); Q.append(Q[-1]+f(n)*Q[-2])
    print("   max |m_n - (m_-1 P_n + (m_0-m_-1)Q_n)| n<=11:", nstr(max(abs(mm(n)-(mm(-1)*P[n+1]+(mm(0)-mm(-1))*Q[n+1])) for n in range(-1,12)),3))

# Remark 2
print("\n== remark 2: k=3/2 ->", nstr(Cconv(lambda n:(n+mpf(1.5))/(2*n*(n+1)))[0],35),
      " k=5/2 ->", nstr(Cconv(lambda n:(n+mpf(2.5))/(2*n*(n+1)))[0]-mpf(5)/3,3), "(minus 5/3)")
# Remark 1
C0=Cconv(lambda n: mpf(1)/(2*(n+1)))[0]; C1=Cconv(lambda n: mpf(1)/(2*n))[0]
print("remark 1: C1-(1+1/(2C0))", nstr(C1-1-1/(2*C0),3))

# First set + printed values
print("\n== Table 3 printed values")
rows = [("C[(1-2n)/(4n(n+1))]", lambda n:(1-2*n)/mpf(4*n*(n+1)), (R2+1)/(2*R2)),
        ("C[(1-2n)/(4(n+1)(n+2))]", lambda n:(1-2*n)/mpf(4*(n+1)*(n+2)), (2*R2+1)/4),
        ("C[(1-2n)/(4(n+2)(n+3))]", lambda n:(1-2*n)/mpf(4*(n+2)*(n+3)), (R2+1)/(6*R2-1)),
        ("C[n/(2n(n+1))]", lambda n:mpf(n)/(2*n*(n+1)), 1/(2*R2-1)),
        ("C[(n+1)/(2n(n+1))]", lambda n:mpf(n+1)/(2*n*(n+1)), R2),
        ("C[(n+2)/(2n(n+1))]", lambda n:mpf(n+2)/(2*n*(n+1)), (2*R2+1)/(R2+1))]
for name,f,pr in rows:
    v,err = Cconv(f,1000)
    print(f"{name:26s} C={nstr(v,30)} printed={nstr(pr,15)} diff={nstr(v-pr,3)} tail={nstr(err,3)}")
v,_ = Cconv(lambda n:(1-2*n)/mpf(4*(n+2)*(n+3)),2000)
print("third row vs (R2+1)/(6R2-6):", nstr(v-(R2+1)/(6*R2-6),3))
mp.dps=40; print("PSLQ third row:", pslq([v*R2,v,R2,1],maxcoeff=10**6,maxsteps=10**6))
# alternative readings of the misprints
mp.dps=60
v0,_=Cconv(lambda n:mpf(n)/(2*n*(n+1)))
print("k=0: 1/(2R2-2) diff", nstr(v0-1/(2*R2-2),3), "; 1/(2*v0-1)?", nstr(1/(2*v0-1),15), "; v0-1=",nstr(v0-1,15))
for N in [1,2,3]:
    print(f"  shifted f_(n+{N}) reading:", nstr(Cconv(lambda n,N=N:mpf(n+N)/(2*(n+N)*(n+N+1)))[0],15))
