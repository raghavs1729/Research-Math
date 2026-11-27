"""Is each formula in O. David's (arXiv:2308.02567, Thm 1) Euler family?
   b(x) = -h1(x)h2(x),  f(x)a(x) = f(x-1)h1(x) + f(x+1)h2(x+1), f polynomial of degree <= 6.
   Also: equivalence (C-transform b_n/(a_n a_{n-1})) of E2 with PCF[n+2, n]."""
from sympy import symbols, expand, Poly, solve, factor, simplify, cancel
import itertools, mpmath as mp
x = symbols('x')
cases = {'E1': (x**2+6*x+7, -x**2*(x+3)), 'E2': ((x+1)*(x+3), x**2*(x+3)), 'E3': (x**2+3*x+3, -x**2*(x+2))}
factors = {'E1': [x, x, x+3], 'E2': [x, x, x+3], 'E3': [x, x, x+2]}
for name, (a, b) in cases.items():
    found = []
    fs = factors[name]
    for mask in itertools.product([0, 1], repeat=3):
        h1 = 1
        for m, fct in zip(mask, fs):
            if m: h1 *= fct
        h2 = expand(-b / h1)  # so that -h1 h2 = b
        for sgn in (1, -1):
            H1, H2 = expand(sgn*h1), expand(cancel(-b/(sgn*h1)))
            for d in range(0, 7):
                cs = symbols(f'c0:{d+1}')
                f = sum(c*x**i for i, c in enumerate(cs))
                eq = expand(f*a - f.subs(x, x-1)*H1 - f.subs(x, x+1)*H2.subs(x, x+1))
                sol = solve(Poly(eq, x).all_coeffs(), cs, dict=True)
                for s in sol:
                    fs_ = expand(f.subs(s))
                    if fs_ != 0:
                        free = [c for c in cs if c not in s]
                        fs_ = fs_.subs({c: 1 for c in free})
                        if fs_ != 0:
                            found.append((factor(H1), factor(H2), factor(fs_)))
                            break
                if found and found[-1][:2] == (factor(H1), factor(H2)):
                    break
    print(name, 'Euler-family decompositions (h1, h2, f):', found if found else 'none with deg f <= 6')

n = symbols('n')
print('E2 C-transform b_n/(a_n a_{n-1}) =', factor(cancel(n**2*(n+3)/((n+1)*(n+3)*n*(n+2)))),
      ';  PCF[n+2,n] C-transform =', factor(n/((n+2)*(n+1))))
mp.mp.dps = 50
def pcf(a, b, N=600):
    pm, p, qm, q = mp.mpf(1), mp.mpf(a(0)), mp.mpf(0), mp.mpf(1)
    for k in range(1, N):
        pm, p = p, a(k)*p + b(k)*pm
        qm, q = q, a(k)*q + b(k)*qm
    return p/q
v = pcf(lambda k: k+2, lambda k: k)
print('PCF[n+2,n] =', mp.nstr(v, 40), ' 1/(2e-5) =', mp.nstr(1/(2*mp.e-5), 40))
print('(3/2)*PCF[n+2,n] =', mp.nstr(1.5*v, 40), ' 3/(4e-10) =', mp.nstr(3/(4*mp.e-10), 40))
# E3 as Mobius transform of Cohen 5.10.8 (k=1): 2e = 3 + 6/(4 - 16/(13 - 45/(21 - ...)))
tail = pcf(lambda k: (k+3)*(k+4)+1, lambda k: -(k+2)**2*(k+4))  # tail starting at 13: a=13,21,31,...; b=-45,-96,...
c5108 = 3 + 6/(4 - 16/tail)
print('Cohen 5.10.8 (k=1) =', mp.nstr(c5108, 40), ' 2e =', mp.nstr(2*mp.e, 40))
e3 = 3 - 3/(7 - 16/tail)
print('E3 via same tail =', mp.nstr(e3, 40), ' 4e/(2e-1) =', mp.nstr(4*mp.e/(2*mp.e-1), 40))
# E2 as Mobius transform of Cohen 5.10.5: e = 5/2 + (3/4)/PCF[(n+1)(n+3), n^2(n+3)]
e2 = pcf(lambda k: (k+1)*(k+3), lambda k: k*k*(k+3))
print('Cohen 5.10.5: 5/2 + (3/4)/E2 =', mp.nstr(mp.mpf(5)/2 + mp.mpf(3)/4/e2, 40), ' e =', mp.nstr(mp.e, 40))

# Apply David's Thm 1 formula (as printed in arXiv:2308.02567v2) to each decomposition found above.
def david(a0, h1, h2, f, K=300):
    s = mp.mpf(0); P = mp.mpf(1)
    for k in range(0, K):
        if k >= 1: P *= mp.mpf(h1(k))/h2(k+1)
        s += mp.mpf(f(0)*f(1))/(f(k)*f(k+1))*P
    return a0 + mp.mpf(f(1)*h2(1))/f(0)*(1/s - 1)  # bracket spans (1/S - 1) in the printed formula
print('David Thm1, E1:', mp.nstr(david(7, lambda t: t+3, lambda t: t*t, lambda t: t**3+3*t**2+2*t+2), 40), mp.nstr(6*mp.e/(2*mp.e-3), 40))
print('David Thm1, E2:', mp.nstr(david(3, lambda t: -t-3, lambda t: t*t, lambda t: (t+2)*(t*t+5*t+5)), 40), mp.nstr(3/(4*mp.e-10), 40))
print('David Thm1, E3:', mp.nstr(david(3, lambda t: t+2, lambda t: t*t, lambda t: 1), 40), mp.nstr(4*mp.e/(2*mp.e-1), 40))
