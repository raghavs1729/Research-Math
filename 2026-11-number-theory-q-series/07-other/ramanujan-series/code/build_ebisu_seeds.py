import json, sys, time
from fractions import Fraction as Fr
from wz import *
S='/tmp/claude-0/-home-user-weighted-prime-races/5eb24282-ae85-5f3e-8e4f-dbc18853349d/scratchpad'
ents = json.load(open(S+'/ebisu_entries.json'))
good = []; bad = []
t0=time.time()
for e in ents:
    if e['kind'] != 'rat' or e['x']=='x': continue
    if 'b' in e['A']+e['B']+e['C']: 
        bad.append((e,'two-parameter')); continue
    try:
        s, T = auto_seed_2F1(e['A'], e['B'], e['C'], e['x'])
        # verify seed property: shift a by 1
        F = T.subs({'a': 'a+n'})
        R = wz_mate(F); assert check_wz(F, R)
        good.append(dict(label=e['label'], A=e['A'], B=e['B'], C=e['C'], x=e['x'], step=s, term=str(T)))
    except Exception as ex:
        bad.append((e, repr(ex)[:80]))
print(len(good), 'seeds built;', len(bad), 'failed', time.time()-t0)
for e, why in bad[:40]: print('  FAIL', e['label'], e['A'], e['B'], e['C'], e['x'], why)
json.dump(good, open('ebisu_seeds.json','w'), indent=0)
