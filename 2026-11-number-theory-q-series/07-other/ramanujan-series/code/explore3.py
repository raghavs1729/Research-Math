"""compact classification of matched configurations"""
import sys, json, mpmath as mp, signal
from explore import *
from numerics import hsum_k, hsum_n
class TO(Exception): pass
def handler(s, f): raise TO()
signal.signal(signal.SIGALRM, handler)
def classify(r, eps='1e-9'):
    mp.mp.dps = 34
    full = [n for n in ALL if n == r['seed'] or n.startswith(r['seed'])][0]
    F = make_F(full, r['shifts'], r['K'], r['base'], r['k0'])
    R = wz_mate(F)
    W = WZNum(F, R)
    e = mp.mpf(eps)
    def safe(f, *a):
        try: return f(*a)
        except TO: raise
        except Exception as ex: return mp.nan
    SF = safe(hsum_k, W, e); SG = safe(hsum_n, W, e)
    g = safe(lambda: mp.fsum(W.Gv(n + e, 30000) for n in range(40)))
    L = safe(hsum_k, W, 40 + e)
    f = lambda x: mp.nstr(x, 9)
    return 'SF %s | SG %s | SG/eps %s | SG*eps %s | g %s | L %s | chk %s' % (f(SF), f(SG), f(SG/e), f(SG*e), f(g), f(L), f(SF + g - SG - L))
seen = set()
for l in open(sys.argv[1]):
    if not l.startswith('{'): continue
    r = json.loads(l)
    if r['id'] not in sys.argv[2:]: continue
    key = (r['seed'], tuple(r['shifts']), r['K'], tuple(r['base']), r['k0'])
    if key in seen: continue
    seen.add(key)
    signal.alarm(60)
    try: res = classify(r)
    except TO: res = 'timeout'
    except Exception as ex: res = 'failed ' + repr(ex)[:80]
    signal.alarm(0)
    print(r['id'], r['seed'][:45], r['shifts'], r['K'], r['base'], r['k0'], 'j=%d' % r['j'], '::', res, flush=True)
