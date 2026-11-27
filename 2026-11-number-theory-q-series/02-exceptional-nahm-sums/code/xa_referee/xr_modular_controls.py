"""Controls for xr_modular.py: non-modular sums must show power corrections."""
import sys, os, math
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
from fractions import Fraction as Fr
from xr_modular import test
from xr_lib import sw_cartan
ts = np.linspace(0.10, 0.30, 9)
# rank one, A = 3 (not in Zagier's list): G=[[3]], k=2 means C(A_1)^{-1} = 1/2 -> use G=[[6]] so that A = 3
test("control: rank-1 Nahm sum A=3 (non-modular)", [[6]], [1], 2, Fr(1), ts, 400)
test("control: rank-1 Nahm sum A=2 (Rogers-Ramanujan, modular, c=2/5)", [[4]], [1], 2, Fr(2, 5), ts, 400)
C, D = sw_cartan('C', 3); G = [[C[i][j]*D[j] for j in range(3)] for i in range(3)]
test("control: C_3 Gram with (q;q) denominators, k=2 (D replaced by 1)", G, [1, 1, 1], 2, Fr(3), ts, 600)
C, D = sw_cartan('B', 3); G = [[C[i][j]*D[j] for j in range(3)] for i in range(3)]
test("control: B_3 Gram with D(C_3)=(1,1,2) denominators", G, [1, 1, 2], 2, Fr(15, 4), ts, 700)
