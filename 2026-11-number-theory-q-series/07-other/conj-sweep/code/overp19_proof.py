"""Certificate for Conjecture 2 of arXiv:2603.08510:
     pbar(4 * 17^2 * n) = 0 (mod 19)  for n = 3 (mod 8), (n/17) = -1, (n/19) = +1.

All arithmetic is in F_19 (FLINT nmod_poly).  Notation: ell = 19, p = 17,
phi = sum q^{n^2}, F = sum_{n odd} sigma(n) q^n, Lambda_{K/2} = span_Z { F^j phi^(K-4j) : 0 <= j <= K/4 }
(= forms in M_{K/2}(Gamma_0(4)) with integer coefficients; triangular basis, F^j phi^(K-4j) = q^j + ...).

 W(q) := sum_{(N/19)=+1} pbar(N) (-q)^N.
 Step 1.  A := phi^18, A_+ := sum_{(n/19)=1} r_18(n) q^n.  Power series: W * phi(q^19) = A_+ (mod 19).
 Step 2.  A_+ = (Theta^18 A + Theta^9 A)/2 (mod 19) lies in Lambda_369 (x) F_19 (Theta-lemma, E_18 = 1 mod 19).
          Compute its coordinates u_j (0 <= j <= 184) from its first 185 coefficients; check u_j = 0 for j >= 180.
          Then A_+ = phi^19 * w0 with w0 := sum_{j<180} u_j F^j phi^(719-4j), hence W = w0 (mod 19) as power series.
 Step 3.  b := T(17^2) w0 in Lambda_{719/2} (x) F_19: coordinates from b(n), n < 180, with Shimura's formula
          b(n) = w0(289 n) + (-1/17)^359 (n/17) 17^358 w0(n) + 17^717 w0(n/289).
 Step 4.  h := b + 17^358 w0;  for (n/17) = -1:  w0(289 n) = h(n).
 Step 5.  g := U_4 h in Lambda_{719/2} (x) F_19: coordinates from h(4n), n < 180.
 Step 6.  Sturm (weight 360, Gamma_1(256), applied to [g sieved to n = 3 mod 8] * phi(q^8)):
          g(n) = 0 for all n = 3 (mod 8), n <= 360*49152/12 = 1474560.
 Conclusion: for n = 3 (8), (n/17) = -1, (n/19) = 1:  pbar(4*289*n) = W(289*4n) = w0(289*4n) = h(4n) = g(n) = 0.
"""
import flint, sys, time
ell = 19; p = 17
def kron(a, q):
    a %= q
    return 0 if a == 0 else (1 if pow(a, (q-1)//2, q) == 1 else -1)
out = open('../cert/overpartition_mod19_cert.txt', 'w')
def P(*x):
    t = ' '.join(str(y) for y in x); print(t, flush=True); out.write(t + '\n')
NP = lambda c: flint.nmod_poly(c, ell)
def series(D):
    phi = [0]*D; n = 0
    while n*n < D: phi[n*n] = 1 if n == 0 else 2; n += 1
    F = [0]*D
    for d in range(1, D, 2):
        for m in range(d, D, 2*d): F[m] = (F[m] + d) % ell
    return NP(phi), NP(F)
def coords(vals, K, phin, Fn, D):
    """coordinates in the basis F^j phi^(K-4j) of the element of Lambda_{K/2}(x)F_19 with first coefficients vals."""
    dim = K//4 + 1; assert len(vals) >= dim and D >= dim
    rem = [v % ell for v in vals[:dim]]; co = []
    phi4 = phin.pow_trunc(4, dim); pw = [None]*dim; cur = phin.pow_trunc(K % 4, dim)
    for j in range(dim-1, -1, -1): pw[j] = cur; cur = cur.mul_low(phi4, dim)
    Fp = NP([1])
    for j in range(dim):
        Bj = Fp.mul_low(pw[j], dim); Fp = Fp.mul_low(Fn, dim)
        assert all(int(Bj[i]) == 0 for i in range(j)) and int(Bj[j]) == 1
        a = rem[j]; co.append(a)
        if a:
            for n in range(j, dim): rem[n] = (rem[n] - a*int(Bj[n])) % ell
    return co
def evaluate2(co, K, phin, Fn, D):
    """sum_j co_j F^j phi^(K-4j): write phi^(K-4j) = phi^(K%4) * (phi^4)^(dim-1-j);
       S = sum_j co_j F^j X^(dim-1-j) with X = phi^4, homogeneous Horner: S_0 = co_0, S_j = S_{j-1}*F ... not homogeneous,
       so use: S = X^(dim-1) * sum_j co_j (F/X)^j.  Instead compute iteratively T_j = T_{j-1} * X + co_j' ... """
    dim = K//4 + 1
    X = phin.pow_trunc(4, D)
    # T = sum_j co_j F^j X^(dim-1-j).  Recurrence over j from 0: T_0 = co_0; T_j = T_{j-1} * X + co_j F^j  (degree in X grows)
    # after dim-1 steps T = sum_j co_j F^j X^(dim-1-j).
    T = NP([co[0]]); Fp = NP([1])
    for j in range(1, dim):
        Fp = Fp.mul_low(Fn, D)
        T = T.mul_low(X, D)
        if co[j]: T = T + co[j]*Fp
    return T.mul_low(phin.pow_trunc(K % 4, D), D)
t0 = time.time()
P('Certificate: pbar(4*17^2*n) = 0 mod 19 for n = 3 mod 8, (n/17) = -1, (n/19) = +1  (arXiv:2603.08510, Conjecture 2)')
# ---------- Steps 1-2
D1 = 400
phin, Fn = series(D1)
A = phin.pow_trunc(18, D1)
Aplus = [int(A[n]) if kron(n, ell) == 1 else 0 for n in range(D1)]
u = coords(Aplus, 738, phin, Fn, D1)
P('step 2: coordinates u_j of A_+ in Lambda_369 (basis F^j phi^(738-4j)), j = 180..184:', u[180:])
assert all(x == 0 for x in u[180:]), 'A_+ not divisible by phi^19'
w0 = u[:180]
P('        A_+ = phi^19 * w0 (mod 19) with w0 = sum_{j<180} u_j F^j phi^(719-4j); u =', w0)
# sanity (not part of the proof): expansion reproduces A_+ further out
chk = evaluate2(u, 738, phin, Fn, D1)
assert all(int(chk[n]) == Aplus[n] for n in range(D1)), 'expansion mismatch'
P(f'        (sanity) the expansion reproduces A_+ for all n < {D1}')
# ---------- Step 3
lam = 359; eps = kron(-1, p)**lam; c1 = pow(p, lam-1, ell); c2 = pow(p, 2*lam-1, ell)
D3 = p*p*180 + 10
phin, Fn = series(D3)
w0s = evaluate2(w0, 719, phin, Fn, D3)
# sanity: w0 = W on this range
pb = phin.inverse_series_trunc(D3)
assert all(int(w0s[n]) == (int(pb[n]) if kron(n, ell) == 1 else 0) for n in range(D3))
P(f'        (sanity) w0(n) = W(n) = [(n/19)=1] pbar(n)(-1)^n mod 19 for all n < {D3}')
def b(n):
    v = int(w0s[p*p*n]) + eps*kron(n, p)*c1*int(w0s[n])
    if n % (p*p) == 0: v += c2*int(w0s[n//(p*p)])
    return v % ell
beta = coords([b(n) for n in range(180)], 719, phin, Fn, D3)
P(f'step 3: (-1/17)^359 = {eps:+d}, 17^358 = {c1} mod 19; coordinates of b = T(17^2) w0:', beta)
# ---------- Step 4
eta = [(beta[j] + c1*w0[j]) % ell for j in range(180)]      # h = b - (n/17) eps c1 w0 with (n/17) = -1, eps = +1
assert eps == 1
P('step 4: coordinates of h = b + 17^358 w0:', eta)
# ---------- Step 5
D5 = 4*180 + 10
phin, Fn = series(D5)
hs = evaluate2(eta, 719, phin, Fn, D5)
gco = coords([int(hs[4*n]) for n in range(180)], 719, phin, Fn, D5)
P('step 5: coordinates of g = U_4 h:', gco)
# ---------- Step 6
Bd = 360*49152//12
D6 = Bd + 9
phin, Fn = series(D6)
gs = evaluate2(gco, 719, phin, Fn, D6)
gl = [int(x) for x in gs.coeffs()] + [0]*10
bad = [n for n in range(3, Bd+1, 8) if gl[n] != 0]
P(f'step 6: Sturm bound {Bd}; g(n) mod 19 at the {len(range(3, Bd+1, 8))} values n = 3 mod 8, n <= {Bd}: ' + ('ALL ZERO' if not bad else f'NONZERO at {bad[:5]}'))
P('        (for comparison, nonzero g(n), n < 4000, by n mod 8:', {r: sum(1 for n in range(r, 4000, 8) if gl[n]) for r in range(8)}, ')')
# sanity (not part of the proof): g(n) = h(4n) and h(n) = w0(289 n) for (n/17) = -1 on a range
D7 = 200000
phin, Fn = series(D7)
hs = evaluate2(eta, 719, phin, Fn, D7); pb = phin.inverse_series_trunc(D7)
assert all(gl[n] == int(hs[4*n]) for n in range(D7//4))
assert all(int(hs[n]) == (int(pb[289*n]) if kron(289*n, ell) == 1 else 0) for n in range(D7//289) if kron(n, 17) == -1)
P(f'        (sanity) g(n) = h(4n) for n < {D7//4}; h(n) = W(289 n) for (n/17) = -1, n < {D7//289}')
P('RESULT: ' + ('PROVED' if not bad else 'NOT PROVED'), f'  [{time.time()-t0:.0f} s]')
out.close()
