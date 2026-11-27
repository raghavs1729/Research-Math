import sys, time
sys.argv = ['x', '1', '3', '/tmp/claude-0/-home-user-weighted-prime-races/5eb24282-ae85-5f3e-8e4f-dbc18853349d/scratchpad/surveyB/verify_B.py']
import run_scan as rs
names = [n for n in rs.ALL if n.startswith('Eb')][:12] + ['Seed1','Seed10','Seed4','Kummer2F1']
for nm in names:
    t=time.time(); name, hits = rs.work(nm)
    print(nm[:60], len(hits), sorted({i for h in hits for i in h[6]}), round(time.time()-t,1), flush=True)
