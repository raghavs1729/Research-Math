"""Sensitivity of the sector test: deliberately wrong data must FAIL."""
import sys, io, os
sys.path.insert(0, os.path.dirname(__file__))
import xr_lib, xr_sectors
from xr_lib import Affine

out = sys.stdout

def trial(name, X, l, k, N, patch):
    orig = {}
    patch(orig)
    buf = io.StringIO()
    try:
        nf = xr_sectors.run(X, l, k, N, buf)
    except AssertionError as e:
        nf = f'assertion: {str(e)[:80]}'
    finally:
        for (obj, attr), v in orig.items():
            setattr(obj, attr, v)
    out.write(f"{name}: failures = {nf}\n")

# 1. C_3 with B_3's D (i.e. transposed convention / wrong long-short assignment) in the Nahm denominators
def p1(orig):
    f = xr_lib.sw_cartan
    orig[(xr_lib, 'sw_cartan')] = f
    def g(X, l):
        C, D = f(X, l)
        if X == 'C':
            D = D[::-1]
            C = [row[::-1] for row in C[::-1]]  # relabel consistently so that G is unchanged, D moves
            C, D = f(X, l)
            D = [2] * (l - 1) + [1]
            C = [[C[i][j] for j in range(l)] for i in range(l)]
            # make CD symmetric again with the wrong D: use C^T
            C = [[f(X, l)[0][j][i] for j in range(l)] for i in range(l)]
        return C, D
    xr_lib.sw_cartan = g
    xr_sectors.Affine = xr_lib.Affine
trial("C_3 k=2 with transposed Cartan matrix (i.e. B_3 Gram) but twisted A_5^(2) data", 'C', 3, 2, 8, p1)

# 2. C_3 with untwisted imaginary multiplicities
def p2(orig):
    orig[(Affine, 'imag_mult')] = Affine.imag_mult
    Affine.imag_mult = lambda self, J: self.l
trial("C_3 k=2 with mult(j delta)=3 for all j", 'C', 3, 2, 8, p2)

# 3. C_3 with untwisted root system (long roots at all j): C_3^(1)
def p3(orig):
    orig[(Affine, 'real_beta')] = Affine.real_beta
    Affine.real_beta = lambda self, J: self.roots
trial("C_3 k=2 with C_3^(1) real roots", 'C', 3, 2, 8, p3)

# 4. G_2 k=2 with D(G_2) = diag(3,1) (orientation swapped)
def p4(orig):
    f = xr_lib.sw_cartan
    orig[(xr_lib, 'sw_cartan')] = f
    def g(X, l):
        C, D = f(X, l)
        if X == 'G':
            C, D = [[2, -3], [-1, 2]], [3, 1]
        return C, D
    xr_lib.sw_cartan = g
trial("G_2 k=2 with C=[[2,-3],[-1,2]], D=(3,1)", 'G', 2, 2, 10, p4)

# 5. T_1 k=4 with denominators (q^2;q^2) i.e. D=2
def p5(orig):
    f = xr_lib.sw_cartan
    orig[(xr_lib, 'sw_cartan')] = f
    def g(X, l):
        C, D = f(X, l)
        if X == 'T':
            D = [1] * (l - 1) + [2]
        return C, D
    xr_lib.sw_cartan = g
trial("T_1 k=4 with D=(2)", 'T', 1, 4, 20, p5)
