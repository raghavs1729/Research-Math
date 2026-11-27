import sys, pickle
from fr_search12 import svals
surv = pickle.load(open(sys.argv[1], 'rb'))
marks = [80, 160, 320, 640]
res = {}
for (A, k, B, C, E, cls, s) in surv:
    v = svals(A, k, B, C, E, marks)
    ok = all(x is not None for x in v)
    # FR-like: s does not keep decreasing: s(640) > s(160) - 0.5  (non-FR loses ~ ln 4 = 1.39)
    fr = ok and v[3] > v[1] - 0.5
    res.setdefault((tuple(cls), fr), []).append((A, k, B, C, E, v))
for key in sorted(res, key=str):
    L = res[key]
    print(key, len(L))
    for row in L[:6]:
        print("    A,k=%d,%d  b=(%d,%d,%d)  s=%s" % (row[0], row[1], row[2], row[3], row[4], ["%.3f" % x if x is not None else None for x in row[5]]))
pickle.dump(res, open(sys.argv[1].replace('stage1', 'stage2'), 'wb'))
