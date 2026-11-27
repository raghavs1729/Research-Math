"""Audit of U3: numerical checks of the analytic steps (independent implementation).
 - Taylor coefficients of Y, Z at 0 (own power-series sums + Richardson), N via mpmath.
 - Connection coefficients of nu, r, r* to the Frobenius basis, SOLVED numerically at x = 1/2 from the
   local series (value + 3 theta-derivatives) and compared with U2 (5.1) / Prop. 3.1.
 - g(u), g*(u) by quadrature; Prop 3.2 (i) at non-integer u; (ii); (3.1); Prop 3.2 (iii) formulas;
   Prop 4.1 value vs direct PCF; Prop 5.1 (Casoratian = -4/(n+1)^5); Lemma 3.3 asymptotics; G(g*)_n.
 - Prop 3.1 Mellin transform at u = 0.7, 2.5.
 - Prop 6.1: T_Z at real and complex points, by the bilateral sum directly and by (6.2); zeta_i relations.
"""
import sympy as sp
from mpmath import mp, mpf, mpc, pi, gamma, zeta, log, quad, taylor, tan, matrix, lu_solve, nstr, exp
mp.dps = 36
half = mpf(1)/2
ORD = 4
def pm(a, b):
    r = [0]*ORD
    for i in range(ORD):
        if a[i] == 0: continue
        for j in range(ORD-i): r[i+j] += a[i]*b[j]
    return r
def pinv(a):
    r = [1/a[0]] + [0]*(ORD-1)
    for n_ in range(1, ORD): r[n_] = -sum(a[k]*r[n_-k] for k in range(1, n_+1))/a[0]
    return r
def richardson(Ns, vals):
    hs = [mpf(1)/N for N in Ns]; P = list(vals)
    for lev in range(1, len(P)):
        for i in range(len(P)-1, lev-1, -1):
            P[i] = (hs[i-lev]*P[i] - hs[i]*P[i-1])/(hs[i-lev] - hs[i])
    return P[-1], P[-1]-P[-2]
def taylor_sums(u0, kinds, Ns=None):
    """Taylor coeffs (order<4) at u0 of sum_k A_k(u0+e) * w_kind(k, u0+e). kinds: 'Y': w=1, 'Z': 1/(u+k+1),
       'Zt': 1/(u+k-1/2)."""
    Ns = Ns or [400*2**j for j in range(8)]
    A = [mpf(1) if not isinstance(u0, mpc) else mpc(1), 0, 0, 0]
    acc = {kd: [0]*ORD for kd in kinds}; rec = {kd: [] for kd in kinds}
    it = iter(Ns); nxt = next(it); k = 0
    while True:
        for kd in kinds:
            off = {'Y': None, 'Z': 1, 'Zt': -half}[kd]
            term = A if off is None else pm(A, pinv([u0 + k + off, 1, 0, 0]))
            acc[kd] = [x+y for x, y in zip(acc[kd], term)]
        f = pm([u0 + half + k, 1, 0, 0], pinv([u0 + 1 + k, 1, 0, 0])); f2 = pm(f, f)
        A = pm(A, pm(f2, f2)); k += 1
        if k == nxt:
            for kd in kinds: rec[kd].append(list(acc[kd]))
            try: nxt = next(it)
            except StopIteration: break
    out = {}
    for kd in kinds:
        out[kd] = [richardson(Ns, [r_[i] for r_ in rec[kd]]) for i in range(ORD)]
    return out
TS = taylor_sums(mpf(0), ['Y', 'Z'])
eta = [c[0] for c in TS['Y']]; zc = [c[0] for c in TS['Z']]
print("Taylor sums: max Richardson increments Y %s, Z %s" % (nstr(max(abs(c[1]) for c in TS['Y']), 2), nstr(max(abs(c[1]) for c in TS['Z']), 2)))
Nf = lambda e: (gamma(1+e)/gamma(half+e))**4
Nt = taylor(Nf, 0, 3)
S0 = eta[0]; S1 = (eta[1] - 8*log(2)*eta[0])/8; Z0 = zc[0]
print("S0, S1, Z0 =", nstr(S0, 25), nstr(S1, 25), nstr(Z0, 25))
etas = [-eta[i]/16 + zc[i]/32 for i in range(4)]; etas[3] += half      # Y*(u) = u^3/2 - Y/16 + Z/32
# Y/N and Z/N Taylor
invN = pinv(list(Nt))
tau = pm(eta, invN); zeta_ = pm(zc, invN)
print("U2 Prop 8: tau2 - 2pi^2 tau0 = %s, tau3 - (4pi^2/3)tau1 = %s" % (nstr(tau[2]-2*pi**2*tau[0], 3), nstr(tau[3]-4*pi**2/3*tau[1], 3)))
e0 = zeta_[0]; e1 = zeta_[1]/pi; e2 = e0 + 8
print("(6.1): zeta0 - pi^2 Z0 = %s; zeta2 - pi^2(e0+e2) = %s; zeta3 - ((4pi^3/3)e1 - 16pi^2) = %s"
      % (nstr(e0-pi**2*Z0, 3), nstr(zeta_[2]-pi**2*(e0+e2), 3), nstr(zeta_[3]-(4*pi**3/3*e1-16*pi**2), 3)))
print("   e1 = %s ; by-product 8(1+4S1+pi^2 S1 Z0)/(pi S0) = %s" % (nstr(e1, 25), nstr(8*(1+4*S1+pi**2*S1*Z0)/(pi*S0), 25)))
# ---- local series at x = 1 (s = 1-x) via the recurrence derived in a1 ----
m = sp.symbols('m'); h = sp.Rational(1, 2)
def th_apply(terms):
    out = {}
    for j, c in terms.items():
        e = m + j; out[j-1] = out.get(j-1, 0) - c*e; out[j] = out.get(j, 0) + c*e
    return out
def add(a, b, cb=1):
    o = dict(a)
    for j, c in b.items(): o[j] = o.get(j, 0) + cb*c
    return o
def times_x(a):
    o = {}
    for j, c in a.items(): o[j] = o.get(j, 0) + c; o[j+1] = o.get(j+1, 0) - c
    return o
Pw = [{0: sp.Integer(1)}]
for _ in range(4): Pw.append(th_apply(Pw[-1]))
sh = {}
for i in range(5): sh = add(sh, Pw[i], sp.binomial(4, i)*h**(4-i))
Lop = {j: sp.lambdify(m, sp.expand(c), 'mpmath') for j, c in add(Pw[4], times_x(sh), -1).items() if sp.expand(c) != 0}
MS = 190
def local(c0, c1, c2):
    c = [mpf(c0), mpf(c1), mpf(c2)]
    for pp in range(MS - 3):
        tot = 0
        for j, Pj in Lop.items():
            idx = pp - j
            if idx != pp+3 and 0 <= idx < len(c): tot += Pj(mpf(idx))*c[idx]
        c.append(-tot/Lop[-3](mpf(pp+3)))
    return c
def theta(c): return [mm*c[mm] - ((mm+1)*c[mm+1] if mm+1 < len(c) else 0) for mm in range(len(c))]
loc = {}
for name, ini in (('nu', (0, 1, mpf(3)/4)), ('r', (0, 0, 1)), ('rs', (1, 0, -mpf(1)/16))):
    c = local(*ini); ths = [c]
    for _ in range(3): ths.append(theta(ths[-1]))
    loc[name] = [t_[:MS-6] for t_ in ths]
def horner(c, s):
    r_ = 0
    for a in reversed(c): r_ = r_*s + a
    return r_
def loc_eval(name, x, i):
    return horner(loc[name][i], 1 - x)
# ---- Frobenius at 0 ----
KF = 200
Acoef = [[mpf(1), 0, 0, 0]]
for k in range(KF):
    f = pm([k + half, 1, 0, 0], pinv([mpf(k+1), 1, 0, 0])); f2 = pm(f, f); Acoef.append(pm(Acoef[-1], pm(f2, f2)))
def frob(x, imax):
    """returns F[i][j] = theta^i y_j (x), i < imax"""
    L = log(x); xe = [L**i/mp.factorial(i) for i in range(4)]
    F = [[0]*4 for _ in range(imax)]; xp = mpf(1); eps_tol = mpf(10)**(-mp.dps-8)
    for k in range(KF):
        cur = pm(Acoef[k], xe)
        for i in range(imax):
            for j in range(4): F[i][j] += cur[j]*xp
            cur = pm(cur, [mpf(k), 1, 0, 0])
        xp *= x
        if xp < eps_tol: break
    return F
x0 = half
F = frob(x0, 4)
Mx = matrix(4, 4)
for i in range(4):
    for j in range(4): Mx[i, j] = F[i][j]
conn = {}
for name in ('nu', 'r', 'rs'):
    rhs = matrix([loc_eval(name, x0, i) for i in range(4)])
    conn[name] = list(lu_solve(Mx, rhs))
pred = {'nu': [(-1)**j*Nt[3-j] for j in range(4)],
        'r': [2*(-1)**j*eta[3-j] for j in range(4)],
        'rs': [2*(-1)**j*etas[3-j] for j in range(4)]}
for name in conn:
    print("connection %-3s solved - predicted: %s" % (name, nstr(max(abs(a-b) for a, b in zip(conn[name], pred[name])), 3)))
# cross-check the solved connection at another point x = 0.3 (Frobenius) vs x = 0.3 local series (|s| = 0.7)
F3 = frob(mpf('0.3'), 2)
for name in conn:
    d = max(abs(sum(conn[name][j]*F3[i][j] for j in range(4)) - loc_eval(name, mpf('0.3'), i)) for i in range(2))
    print("   %-3s at x=0.3: Frobenius vs local series: %s" % (name, nstr(d, 3)))
# ---- mu, mu* and their Mellin integrals ----
cache = {}
def fields(x):
    if x in cache: return cache[x]
    if x < half:
        Fx = frob(x, 2)
        val = {nm: [sum(conn[nm][j]*Fx[i][j] for j in range(4)) for i in range(2)] for nm in conn}
    else:
        val = {nm: [loc_eval(nm, x, i) for i in range(2)] for nm in conn}
    mu = val['nu'][0]*val['r'][1] - val['r'][0]*val['nu'][1]
    mus = val['nu'][0]*val['rs'][1] - val['rs'][0]*val['nu'][1]
    cache[x] = (mu, mus, val['rs'][0]); return cache[x]
mp.dps = 30
def g(u, which):
    idx = {'g': 0, 'gs': 1}[which]
    return quad(lambda x: fields(x)[idx]*x**u, [0, mpf(1)/1000, mpf(1)/20, half, 1])
gs = {n_: g(n_, 'gs') for n_ in range(0, 6)}; gg = {n_: g(n_, 'g') for n_ in range(0, 6)}
atil = lambda uu: 2*(uu+half)*((uu+half)**4 + (uu+half)**2 + mpf(1)/16)
for uu in (mpf('0.37'), mpf('1.6')):
    rr = uu**5*g(uu-1, 'gs') + (uu+1)**5*g(uu+1, 'gs') - atil(uu)*g(uu, 'gs')
    print("Prop 3.2(i) residual at u=%s: %s" % (uu, nstr(rr, 3)))
c = conn
mu4s_conn = (c['nu'][2]*c['rs'][3] - c['nu'][3]*c['rs'][2])/12
mu4_conn = (c['nu'][2]*c['r'][3] - c['nu'][3]*c['r'][2])/12
mu4s_form = -(Nt[0]*(-etas[1]) + Nt[1]*etas[0])/6
print("mu4*: from solved connection %s; formula (N0 eta*1 - N1 eta*0)/6 = %s; from g*(1) - atil(0) g*(0) = -24 mu4*: %s"
      % (nstr(mu4s_conn, 20), nstr(mu4s_form, 20), nstr(-(gs[1]-atil(0)*gs[0])/24, 20)))
print("mu4 : from connection %s; 4S1/(3pi^2) = %s" % (nstr(mu4_conn, 20), nstr(4*S1/(3*pi**2), 20)))
mp.dps = 36
Xi = taylor(lambda e: 2*zeta(3, half+e)*Nf(e), 0, 3)
gs0_form = -2*sum(Xi[j]*(-1)**(3-j)*etas[3-j] for j in range(4))
mp.dps = 30
print("g*(0): quadrature %s ; -2[e^3](Xi Y*(-e)) = %s" % (nstr(gs[0], 22), nstr(gs0_form, 22)))
# (3.1)
def seq(at, a0, a1, N):
    S = [mpf(a0), mpf(a1)]
    for mm in range(1, N): S.append((at(mm)*S[mm] - mm**5*S[mm-1])/(mm+1)**5)
    return S
w = seq(atil, 1, atil(0), 40); v = seq(atil, 0, mpf(1)/8, 40)
for n_ in range(2, 6):
    rhs31 = gs[0]*w[n_] - 192*mu4s_form*v[n_]
    print("(3.1) n=%d: g*(n) quad = %s, g*(0)w_n - 192mu4* v_n = %s, diff = %s" % (n_, nstr(gs[n_], 22), nstr(rhs31, 22), mp.nstr(gs[n_] - rhs31, 3, min_fixed=0, max_fixed=0)))
# Casoratian
for n_ in range(0, 4):
    C = gg[n_]*gs[n_+1] - gg[n_+1]*gs[n_]
    print("Prop 5.1: -(n+1)^5 C[g,g*]_n/24 at n=%d: %s" % (n_, nstr(-(n_+1)**5*C/24, 22)))
# Prop 4.1 value
mp.dps = 30
ell3 = gs[0]/(192*mu4s_form) - 4
print("Prop 4.1: g*(0)/(192 mu4*) - 4 = %s  (1/PCF direct = 0.0437809038563519708389373215186)" % nstr(ell3, 25))
# Lemma 3.3 / G(g*)_n through (3.1) for large n
at3 = lambda mm: 2*(mm+half)*((mm+half)**4 + 5*(mm+half)**2 + mpf(25)/16)
w3 = seq(at3, 1, at3(0), 400); v3 = seq(at3, 0, mpf(1)/8, 400)
for n_ in (10, 50, 200, 399):
    Gn = (gs[0] - 768*mu4s_form)*w3[n_] - 192*mu4s_form*v3[n_]
    print("   G(g*)_n via (3.1), n=%d: %s   n*G = %s" % (n_, nstr(Gn, 5), nstr(n_*Gn, 6)))
for n_ in (20, 30):
    print("   g*(n) n=%d by quadrature: n*g* = %s ; g(n)*n^3 = %s" % (n_, nstr(n_*g(n_, 'gs'), 8), nstr(n_**3*g(n_, 'g'), 8)))
# Mellin transform of r*
for uu in (mpf('0.7'), mpf('2.5')):
    TSu = taylor_sums(uu, ['Y', 'Z'], Ns=[400*2**j for j in range(7)])
    Yu = TSu['Y'][0][0]; Zu = TSu['Z'][0][0]
    rhat = quad(lambda x: fields(x)[2]*x**(uu-1), [0, mpf(1)/1000, mpf(1)/20, half, 1])
    print("Prop 3.1 at u=%s: rhat* - (1/u - (2Y-Z)/(16u^4)) = %s" % (uu, nstr(rhat - (1/uu - (2*Yu - Zu)/(16*uu**4)), 3)))
# Prop 6.1
mp.dps = 30
def TZ_62(uu):
    vv = half - uu
    Zu = taylor_sums(uu, ['Z'], Ns=[400*2**j for j in range(7)])['Z'][0][0]
    Ztv = taylor_sums(vv, ['Zt'], Ns=[400*2**j for j in range(7)])['Zt'][0][0]
    return Zu/Nf(uu) - tan(pi*uu)**4*Ztv/Nf(vv)
def TZ_direct(uu, Ns=[200*2**j for j in range(7)]):
    gam = lambda z: (gamma(z+half)/gamma(z+1))**4
    # k >= 0 side and k < 0 side separately, Richardson each
    out = []
    for side in (1, -1):
        acc = 0; rec = []; it = iter(Ns); nxt = next(it); k = 0 if side == 1 else 1
        cnt = 0
        while True:
            z = uu + side*k
            acc += gam(z)/(z+1); k += 1; cnt += 1
            if cnt == nxt:
                rec.append(acc)
                try: nxt = next(it)
                except StopIteration: break
        out.append(richardson(Ns, rec)[0])
    return out[0] + out[1]
for uu in (mpf('0.2'), mpc('0.3', '0.4'), mpf('-0.35')):
    t = tan(pi*uu); R = (1+t*t)*(e0 + e1*t + e2*t*t)
    a = TZ_62(uu)
    line = "Prop 6.1 at u=%s: (6.2) - closed form = %s" % (nstr(uu, 3), nstr(abs(a - R), 3))
    if uu == mpf('0.2'):
        line += ";  direct bilateral - (6.2) = %s" % nstr(abs(TZ_direct(uu) - a), 3)
    print(line)
# e2 as limit: T_Z(u)/t^4 near u = 1/2
for d in (mpf('1e-3'), mpf('1e-5')):
    uu = half - d; t = tan(pi*uu)
    print("   T_Z(1/2-%s)/t^4 = %s  (e2 = %s)" % (nstr(d, 2), nstr(TZ_62(uu)/t**4, 15), nstr(e2, 15)))
print("a3 done")
