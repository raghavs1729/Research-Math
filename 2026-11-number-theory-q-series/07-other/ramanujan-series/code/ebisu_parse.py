"""Parse the tables of A. Ebisu, 'Special values of the hypergeometric series' (arXiv:1308.5588v3 TeX source)
and extract every evaluation  2F1(A(a,b), B(a,b); C(a,b); x)  with rational x.
Only the left-hand sides are used: the closed forms are re-derived by Zeilberger's algorithm."""
import re, sys, json
src = open(sys.argv[1]).read()
start = src.index(r'\section{Tables of special values}')
end = src.index(r'\section{Some hypergeometric identities')
body = src[start:end]
# normalise
body = body.replace(r'\mbox{$F$}', 'F').replace(r'\,', ' ').replace('\n', ' ')
pat = re.compile(r'F\s*\(([^;()]*(?:\([^()]*\)[^;()]*)*);([^;()]*(?:\([^()]*\)[^;()]*)*);([^()]*(?:\([^()]*\)[^()]*)*)\)')
out = []
seen = set()
def norm(e):
    e = e.strip().strip('{}').strip()
    e = re.sub(r'\\left|\\right', '', e)
    e = re.sub(r'(\d+/\d+)\s+([ab])', r'(\1)*\2', e)
    e = re.sub(r'(\d+)\s+([ab])', r'\1*\2', e)
    e = e.replace(' ', '')
    return e
# track current subsubsection label
pos_labels = [(m.start(), m.group(1)) for m in re.finditer(r'\\paragraph\{\(([0-9,\-]+)\)\}', body)]
def label_at(p):
    lab = None
    for q, l in pos_labels:
        if q <= p: lab = l
        else: break
    return lab
for m in pat.finditer(body):
    ab, c, x = m.group(1), m.group(2), m.group(3)
    if '\\' in x or 'sqrt' in x: 
        kind = 'alg'
    else:
        kind = 'rat'
    parts = ab.split(',')
    if len(parts) != 2: continue
    A, B, Cc, X = norm(parts[0]), norm(parts[1]), norm(c), norm(x)
    if not re.fullmatch(r'[0-9ab+\-*/()]+', A + B + Cc): continue
    key = (A, B, Cc, X)
    if key in seen: continue
    seen.add(key)
    out.append(dict(label=label_at(m.start()), A=A, B=B, C=Cc, x=X, kind=kind))
json.dump(out, open(sys.argv[2], 'w'), indent=0)
print(len(out), 'entries;', sum(1 for o in out if o['kind']=='rat'), 'rational x')
from collections import Counter
print(Counter(o['x'] for o in out if o['kind']=='rat').most_common(80))
