# invariants of Ad_{T_K} in the q-Weyl algebra  (P S = q S P,  U Y = q^2 Y U; P,S commute with Y,U)
import sympy as sp, itertools, sys
S, U = sp.symbols('S U')
qv, xv = sp.Rational(3, 7), sp.Rational(5, 11)
# element: dict {(a,b): expr in S,U}  meaning sum expr * P^a Y^b
def mul(A, B):
    C = {}
    for (a, b), f in A.items():
        for (c, d), g in B.items():
            g2 = g.subs({S: qv**a*S, U: qv**(-2*b)*U}, simultaneous=True)
            key = (a+c, b+d)
            C[key] = sp.cancel(C.get(key, 0) + f*g2)
    return {k: v for k, v in C.items() if v != 0}
def add(A, B, cB=1):
    C = dict(A)
    for k, v in B.items():
        C[k] = sp.cancel(C.get(k, 0) + cB*v)
    return {k: v for k, v in C.items() if v != 0}
one = {(0, 0): sp.Integer(1)}
gen = {'P': {(1, 0): sp.Integer(1)}, 'Y': {(0, 1): sp.Integer(1)}, 'S': {(0, 0): S}, 'U': {(0, 0): U},
       'Pi': {(-1, 0): sp.Integer(1)}, 'Yi': {(0, -1): sp.Integer(1)}}
def build(kind):
    """returns images of P, P^-1, Y, Y^-1, S, U under Ad_T"""
    if kind == 'K':
        AP = {(1, 0): 1/(1-S)}
        APi = {(-1, 0): (1-S/qv)}
        AY = {(0, 1): (1-U/qv**2)}
        AYi = {(0, -1): 1/(1-U)}
        AS = {(2, 1): xv*qv*S*(1-U/qv**2)/((1-S)*(1-qv*S))}
        # Ad(U) = Ad(Y)^-2 * x^-2 q^-2 P^-2 (1-qS)(1-S) U ; write P^-2 f(S) = f(q^-2 S) P^-2
        W = {(-2, 0): ((1-qv*S)*(1-S)*U).subs(S, S/qv**2)/(xv**2*qv**2)}
        AU = mul(mul(AYi, AYi), W)
        return AP, APi, AY, AYi, AS, AU
def image(mono, imgs):
    a, b, c, d = mono
    AP, APi, AY, AYi, AS, AU = imgs
    R = one
    for _ in range(abs(a)): R = mul(R, AP if a > 0 else APi)
    for _ in range(abs(b)): R = mul(R, AY if b > 0 else AYi)
    for _ in range(c): R = mul(R, AS)
    for _ in range(d): R = mul(R, AU)
    return R
def elem(mono):
    a, b, c, d = mono
    return {(a, b): S**c*U**d}
if __name__ == "__main__":
    imgs = build('K')
    # sanity: relations preserved: Ad(P)Ad(S) = q Ad(S)Ad(P)
    AP, APi, AY, AYi, AS, AU = imgs
    print("PS=qSP:", add(mul(AP, AS), mul(AS, AP), -qv) == {}, " UY=q^2YU:", add(mul(AU, AY), mul(AY, AU), -qv**2) == {},
          " PY=YP:", add(mul(AP, AY), mul(AY, AP), -1) == {}, " SU=US:", add(mul(AS, AU), mul(AU, AS), -1) == {},
          " PU=UP:", add(mul(AP, AU), mul(AU, AP), -1) == {}, " YS=SY:", add(mul(AY, AS), mul(AS, AY), -1) == {})
    rng_a = range(-2, 3); rng_b = range(-2, 3); rng_c = range(0, 3); rng_d = range(0, 3)
    monos = [m for m in itertools.product(rng_a, rng_b, rng_c, rng_d)]
    print("monomials:", len(monos)); sys.stdout.flush()
    # compute Ad(m) - m for each monomial, collect as vectors over basis of (a,b) and monomials S^i U^j after clearing denominators
    diffs = []
    for m in monos:
        D = add(image(m, imgs), elem(m), -1)
        diffs.append(D)
    # linear dependency: find h with sum h_m D_m = 0.  Evaluate at random numeric points (S,U) to get linear equations
    import random
    random.seed(1)
    keys = sorted(set(k for D in diffs for k in D))
    pts = [(sp.Rational(random.randint(2, 50), random.randint(51, 99)), sp.Rational(random.randint(2, 50), random.randint(51, 99))) for _ in range(6)]
    rows = []
    for key in keys:
        for (s0, u0) in pts:
            rows.append([D.get(key, sp.Integer(0)).subs({S: s0, U: u0}) if key in D else sp.Integer(0) for D in diffs])
    Mx = sp.Matrix(rows)
    ns = Mx.nullspace()
    print("nullspace dim:", len(ns))
    for v in ns:
        terms = [(monos[i], v[i]) for i in range(len(monos)) if v[i] != 0]
        print("  H =", terms)
