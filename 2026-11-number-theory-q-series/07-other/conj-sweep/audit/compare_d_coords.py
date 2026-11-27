"""AUDIT: compare the coordinate vectors u (=w0), beta (=b), eta (=h), g written by d.gp (d_coords.txt)
with the vectors printed in ../cert/overpartition_mod19_cert.txt."""
import re, ast
mine = {}
for line in open('d_coords.txt'):
    k, v = line.strip().split('=', 1); mine[k] = ast.literal_eval(v)
cert = open('../cert/overpartition_mod19_cert.txt').read()
def grab(pat):
    m = re.search(pat + r'\s*(\[[0-9, ]+\])', cert); return ast.literal_eval(m.group(1))
theirs = {'u': grab(r'u ='), 'beta': grab(r'coordinates of b = T\(17\^2\) w0:'), 'eta': grab(r'coordinates of h = b \+ 17\^358 w0:'), 'g': grab(r'coordinates of g = U_4 h:')}
for k in ('u', 'beta', 'eta', 'g'):
    print(k, 'length', len(mine[k]), len(theirs[k]), 'IDENTICAL' if mine[k] == theirs[k] else 'DIFFERENT')
