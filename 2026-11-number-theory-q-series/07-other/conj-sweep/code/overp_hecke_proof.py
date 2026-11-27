"""Certificate for two overpartition congruences stated as conjectural in arXiv:2603.08510 (closing remarks):

  (A) pbar(17 * 11^2 * n) = 0 (mod 17)  for n = 3 (mod 8), (n/11) = -1;
  (B) pbar(23 * 13^2 * n) = 0 (mod 23)  for n = 5 (mod 8), (n/13) = +1.

Notation: phi = sum_{n in Z} q^{n^2},  F = sum_{n odd >0} sigma(n) q^n,  B_j = F^j phi^{K-4j}
(basis of M_{K/2}(Gamma_0(4)), Koblitz, Ch. IV, Prop. 4).  sum pbar(n) (-q)^n = 1/phi(q).

Steps verified here, all in exact integer arithmetic (reductions mod ell only at the end):
 1. alpha: T_ell(phi^(ell-1)) = sum_j alpha_j B_j  (K = ell-1), from the first K/4+1 coefficients;
    alpha_top = 0 mod ell when 4 | ell-1 (otherwise every B_j is divisible by phi), so U_ell(1/phi) = c := sum_j alpha_j F^j phi^(ell-2-4j)  (mod ell).
 2. beta:  T(p^2) c0 = sum_j beta_j B_j  (K = ell-2), c0 = integral lift of c with coefficients alpha_j mod ell,
    from the first K/4+1 coefficients of Shimura's formula
        b(n) = c0(p^2 n) + (-1/p)^lam (n/p) p^(lam-1) c0(n) + p^(2 lam-1) c0(n/p^2),  lam = (ell-3)/2.
 3. h := sum_j (beta_j - s (-1/p)^lam p^(lam-1) gamma_j) B_j; for (n/p) = s:  c0(p^2 n) = h(n).
 4. Sturm: h(n) = 0 (mod ell) for all n = r (mod 8), n <= (lam+1) * 49152 / 12
    (bound for weight lam+1 on Gamma_1(256), applied to  [sieve of h] * phi(q^8)).
"""
import flint, sys
Z = flint.fmpz_poly
def kron(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p-1)//2, p) == 1 else -1)
def series(D):
    phi = [0]*D; n = 0
    while n*n < D: phi[n*n] = 1 if n == 0 else 2; n += 1
    F = [0]*D
    for d in range(1, D, 2):            # sigma(n) for odd n
        for mlt in range(d, D, 2*d): F[mlt] += d
    return phi, F
def tr(x, D): return Z(x.coeffs()[:D])
def powers(x, emax, D):
    out = [Z([1])]
    for _ in range(emax): out.append(tr(out[-1]*x, D))
    return out
def basis_exact(K, D):
    phi, F = series(D); phiP = powers(Z(phi), K, D); FP = powers(Z(F), K//4, D)
    return [[int(v) for v in tr(FP[j]*phiP[K-4*j], D).coeffs()] + [0]*D for j in range(K//4+1)]
def expand(vals, B):
    """vals[0..dim-1] = first coefficients of a form in span(B_j); returns integer coordinates (triangular basis)."""
    dim = len(B); rem = list(vals); co = []
    for j in range(dim):
        assert all(B[j][i] == 0 for i in range(j)) and B[j][j] == 1
        a = rem[j]; co.append(a)
        for n in range(dim): rem[n] -= a*B[j][n]
    return co
def prove(ell, p, r, s, out):
    lines = []
    def P(*x):
        t = ' '.join(str(y) for y in x); lines.append(t); print(t, flush=True)
    lam = (ell-3)//2
    P(f'## pbar({ell}*{p}^2*n) = 0 mod {ell} for n = {r} mod 8, (n/{p}) = {s:+d}')
    # step 1
    K1 = ell-1; dim1 = K1//4 + 1; D1 = ell*dim1 + 5
    B1 = basis_exact(K1, D1)
    f = B1[0]                                            # phi^(ell-1)
    w1 = (ell-1)//2                                      # weight
    chi = 1 if w1 % 2 == 0 else kron(-1, ell)            # character chi_{-4}^{w1} at ell
    Tl = [f[ell*n] + (chi*ell**(w1-1)*f[n//ell] if n % ell == 0 else 0) for n in range(dim1)]
    alpha = expand(Tl, B1)
    P(f'step 1: T_{ell}(phi^{ell-1}) = sum alpha_j F^j phi^({ell-1}-4j), alpha mod {ell} = {[a % ell for a in alpha]}')
    if K1 % 4 == 0:                                      # the top basis element F^(K1/4) is not divisible by phi
        assert alpha[-1] % ell == 0, 'top coefficient must vanish mod ell (divisibility by phi)'
        gamma = [a % ell for a in alpha[:-1]]
    else:
        gamma = [a % ell for a in alpha]                 # c0 = sum gamma_j F^j phi^(ell-2-4j)
    assert len(gamma) == (ell-2)//4 + 1
    P(f'        c0 = sum gamma_j F^j phi^({ell-2}-4j), gamma = {gamma}')
    # step 2
    K2 = ell-2; dim2 = K2//4 + 1; D2 = p*p*(dim2 + 25) + 5
    B2 = basis_exact(K2, D2)
    c0 = [sum(gamma[j]*B2[j][n] for j in range(dim2)) for n in range(D2)]
    eps = kron(-1, p)**lam
    def b(n):
        v = c0[p*p*n] + eps*kron(n, p)*p**(lam-1)*c0[n]
        if n % (p*p) == 0: v += p**(2*lam-1)*c0[n//(p*p)]
        return v
    beta = expand([b(n) for n in range(dim2)], B2)
    P(f'step 2: T({p}^2) c0 = sum beta_j F^j phi^({ell-2}-4j), beta mod {ell} = {[x % ell for x in beta]}  (lam = {lam}, (-1/{p})^lam = {eps:+d})')
    # extra consistency (not needed for the proof): Shimura's formula against the expansion for 25 further n
    nchk = (D2-1)//(p*p); assert nchk >= dim2 + 25
    assert all(sum(beta[j]*B2[j][n] for j in range(dim2)) == b(n) for n in range(nchk))
    # step 3
    mu = (s*eps*pow(p, lam-1, ell)) % ell
    eta = [(beta[j] - mu*gamma[j]) % ell for j in range(dim2)]
    P(f'step 3: for (n/{p}) = {s:+d}: c0({p}^2 n) = b(n) - ({s*eps:+d})*{p}^{lam-1} c0(n) = h(n), {p}^{lam-1} = {pow(p,lam-1,ell)} mod {ell};')
    P(f'        h = sum eta_j F^j phi^({ell-2}-4j) mod {ell}, eta = {eta}')
    # step 4
    Bd = (lam+1)*49152//12
    D = Bd + 9
    phi, F = series(D)
    N = lambda c: flint.nmod_poly([x % ell for x in c], ell)
    phin, Fn = N(phi), N(F)
    h = flint.nmod_poly([0], ell)
    Fp = flint.nmod_poly([1], ell)
    for j in range(dim2):
        h = h + eta[j]*Fp.mul_low(phin.pow_trunc(K2-4*j, D), D)
        Fp = Fp.mul_low(Fn, D)
    hc = [int(h[n]) for n in range(D)]
    bad = [n for n in range(r, Bd+1, 8) if hc[n] != 0]
    cnt = len(range(r, Bd+1, 8))
    P(f'step 4: Sturm bound for weight {lam+1} on Gamma_1(256): {Bd}; h(n) mod {ell} at the {cnt} values n = {r} mod 8, n <= {Bd}: ' +
      ('ALL ZERO' if not bad else f'NONZERO at {bad[:5]}'))
    other = {rr: sum(1 for n in range(rr, 2000, 8) if hc[n] != 0) for rr in range(8)}
    P(f'        (for comparison, number of nonzero h(n), n < 2000, by n mod 8: {other})')
    P('RESULT: ' + ('PROVED' if not bad else 'NOT PROVED'))
    with open(out, 'a') as fh: fh.write('\n'.join(lines) + '\n\n')
    return not bad
if __name__ == '__main__':
    out = '../cert/overpartition_hecke_cert.txt'
    open(out, 'w').write('Certificates: overpartition congruences mod 17 and 23 (arXiv:2603.08510, closing remarks)\n\n')
    assert prove(17, 11, 3, -1, out)
    assert prove(23, 13, 5, +1, out)
