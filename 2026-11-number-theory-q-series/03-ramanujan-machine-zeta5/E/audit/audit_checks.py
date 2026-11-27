"""Independent referee checks for E/THEOREM.md (E1, E2, E3). Written from scratch."""
from sympy import symbols, expand, simplify, factorial, Rational, E, prod, nsimplify
from fractions import Fraction
from math import factorial as fac
import mpmath as mp

n, j = symbols('n j')
mp.mp.dps = 60

cases = {
    'E1': dict(a=n**2+6*n+7, b=-n**2*(n+3), Pi=n**3+6*n**2+11*n+8,
               N=-(j**2+3*j+1)/3, val=lambda: 6*mp.e/(2*mp.e-3), sym=6*E/(2*E-3)),
    'E2': dict(a=(n+1)*(n+3), b=n**2*(n+3), Pi=(n+3)*(n**2+7*n+11),
               N=-(j+3)/2, val=lambda: 3/(4*mp.e-10), sym=3/(4*E-10)),
    'E3': dict(a=n**2+3*n+3, b=-n**2*(n+2), Pi=1+0*n,
               N=None, val=lambda: 4*mp.e/(2*mp.e-1), sym=4*E/(2*E-1)),
}

def U(Pi, k):
    return fac(k+1)**2 * Fraction(int(Pi.subs(n, k)))

ok = True
def chk(name, cond):
    global ok
    print(('PASS ' if cond else 'FAIL ') + name)
    ok &= bool(cond)

for name, c in cases.items():
    a, b, Pi = c['a'], c['b'], c['Pi']
    A = lambda k: Fraction(int(a.subs(n, k)))
    B = lambda k: Fraction(int(b.subs(n, k)))
    # (1) polynomial identity, and the recurrence for U directly, incl. n=1 (uses U_{-1})
    chk(f'{name} poly identity', expand((n+1)**2*n**2*Pi - a*n**2*Pi.subs(n, n-1) - b*Pi.subs(n, n-2)) == 0)
    chk(f'{name} U solves recurrence n=1..60 (n=1 uses U_-1)',
        all(U(Pi, k) == A(k)*U(Pi, k-1) + B(k)*U(Pi, k-2) for k in range(1, 61)))
    chk(f'{name} U_n>0 n=-1..200', all(U(Pi, k) > 0 for k in range(-1, 201)))
    print(f'   U_-1={U(Pi,-1)}, U_0={U(Pi,0)}')
    # (2) product of (-b_i)
    for k in range(0, 15):
        P = 1
        for i in range(1, k+1):
            P *= -B(i)
        sgn = (-1)**k if name == 'E2' else 1
        shift = 2 if name == 'E3' else 3
        claimed = sgn * fac(k)**2 * fac(k+shift) // (2 if name == 'E3' else 6)
        if P != claimed:
            chk(f'{name} prod(-b) at j={k}', False)
    chk(f'{name} prod(-b_i) formula j<=14', True)
    # (3) telescoping identity
    if name == 'E1':
        Nj = c['N']; PiJ = Pi.subs(n, j); PiJm = Pi.subs(n, j-1)
        chk('E1 telescoping poly identity',
            expand((j+2)*(j+3) - ((j+1)*Nj*PiJ - Nj.subs(j, j+1)*PiJm + Rational(1, 3)*PiJ*PiJm)) == 0)
        chk('E1 Pi(j)=? summand form: T summand = (j+2)(j+3)/(6 (j+1)! Pi(j)Pi(j-1))',
            all(Fraction(fac(k)**2*fac(k+3), 6)/(U(Pi, k)*U(Pi, k-1)) ==
                Fraction((k+2)*(k+3), 6*fac(k+1)*int(Pi.subs(n, k))*int(Pi.subs(n, k-1))) for k in range(0, 20)))
    if name == 'E2':
        cm = lambda m: m**2+3*m+1
        Nj = c['N']
        chk('E2 Pi(j)=(j+3)c_{j+2}, Pi(j-1)=(j+2)c_{j+1}',
            expand(Pi.subs(n, j) - (j+3)*cm(j+2)) == 0 and expand(Pi.subs(n, j-1) - (j+2)*cm(j+1)) == 0)
        chk('E2 telescoping poly identity',
            expand(1 - ((j+1)*Nj*cm(j+2) + Nj.subs(j, j+1)*cm(j+1) + Rational(1, 2)*cm(j+1)*cm(j+2))) == 0)
        chk('E2 summand form', all(Fraction((-1)**k*fac(k)**2*fac(k+3), 6)/(U(Pi, k)*U(Pi, k-1)) ==
            Fraction((-1)**k, 6*fac(k+1)*cm(k+1)*cm(k+2)) for k in range(0, 20)))
    # (4) Lemma: exact convergents vs closed formula
    p = [Fraction(1), A(0)]; q = [Fraction(0), Fraction(1)]
    T = Fraction(0); P = Fraction(1); good = True
    Um1, U0 = U(Pi, -1), U(Pi, 0)
    for k in range(0, 41):
        if k >= 1:
            p.append(A(k)*p[-1] + B(k)*p[-2]); q.append(A(k)*q[-1] + B(k)*q[-2]); P *= -B(k)
        T += P/(U(Pi, k)*U(Pi, k-1))
        if q[-1] != U(Pi, k)*Um1*T or p[-1] != U(Pi, k)*(1/Um1 + (A(0)*Um1 - U0)*T):
            good = False
        if p[-1]/q[-1] != A(0) - U0/Um1 + 1/(Um1**2*T):
            good = False
    chk(f'{name} Lemma formulas for q_n, p_n, p_n/q_n exact n<=40', good)
    # (5) limit T and final algebra symbolically
    Tlim = {'E1': (2*E-3)/36, 'E2': (2*E-5)/(60*E), 'E3': E-Rational(1, 2)}[name]
    chk(f'{name} final algebra a0-U0/U-1+1/(U-1^2 T) == target',
        simplify(int(a.subs(n, 0)) - Rational(int(U0), int(Um1)) + 1/(int(Um1)**2*Tlim) - c['sym']) == 0)
    Tnum = mp.mpf(0); Pm = mp.mpf(1)
    for k in range(0, 200):
        if k >= 1: Pm *= -int(b.subs(n, k))
        Tnum += Pm / (mp.mpf(U(Pi, k).numerator) * mp.mpf(U(Pi, k-1).numerator))
    chk(f'{name} T series == closed form (50 digits)', abs(Tnum - mp.mpf(mp.mpmathify(str(Tlim.evalf(70))))) < mp.mpf(10)**-50)
    # (6) numerics straight from recurrence (mpmath, backward evaluation and forward ratio)
    pm, pp, qm, qp = mp.mpf(1), mp.mpf(int(a.subs(n, 0))), mp.mpf(0), mp.mpf(1)
    for k in range(1, 400):
        ak, bk = int(a.subs(n, k)), int(b.subs(n, k))
        pm, pp = pp, ak*pp + bk*pm
        qm, qp = qp, ak*qp + bk*qm
    v = pp/qp
    print(f'   {name}: PCF(400) = {mp.nstr(v, 45)}')
    print(f'   {name}: target   = {mp.nstr(c["val"](), 45)}  diff = {mp.nstr(abs(v-c["val"]()), 3)}')
    chk(f'{name} numeric agreement >= 40 digits', abs(v - c['val']()) < mp.mpf(10)**-40)

# sign misprint check for row 2: PCF with printed minus sign
a2 = lambda k: (k+1)*(k+3)
pm, pp, qm, qp = mp.mpf(1), mp.mpf(3), mp.mpf(0), mp.mpf(1)
for k in range(1, 400):
    pm, pp = pp, a2(k)*pp - (k+3)*k*k*pm
    qm, qp = qp, a2(k)*qp - (k+3)*k*k*qm
print('row 2 with printed b_n = -(n+3)n^2 :', mp.nstr(pp/qp, 30), ' vs 3/(4e-10) =', mp.nstr(3/(4*mp.e-10), 30))
print('first partial numerators E1:', [-(k+3)*k*k for k in range(1, 5)], ' E2 (+):', [(k+3)*k*k for k in range(1, 5)],
      ' E3:', [-(k+2)*k*k for k in range(1, 5)])
print('first partial denominators E1:', [k*k+6*k+7 for k in range(5)], ' E2:', [(k+1)*(k+3) for k in range(5)],
      ' E3:', [k*k+3*k+3 for k in range(5)])
print('ALL PASS' if ok else 'SOME FAIL')
