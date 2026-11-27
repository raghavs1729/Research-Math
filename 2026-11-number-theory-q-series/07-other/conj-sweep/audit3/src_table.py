# read-only: parse all table rows "d & ( f, g, h ) & k & p & m" from main.tex and dlistcor (python3 -I)
import sys, tarfile, re
tf = tarfile.open(sys.argv[1], 'r:*')
rx = re.compile(r'^\s*(\d+)\s*&\s*\(\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*\)\s*&\s*(-?\d+)\s*&\s*(-?\d+)\s*&\s*(-?\d+)')
for member in ('main.tex', 'dlistcor'):
    txt = tf.extractfile(member).read().decode('utf-8', 'replace').split('\n')
    rows = []
    for i, l in enumerate(txt):
        for part in l.split('\\\\'):
            mm = rx.match(part)
            if mm: rows.append((i+1,) + tuple(int(v) for v in mm.groups()))
    print(f'## {member}: {len(rows)} parsed rows; rows with inconsistent d != h-g:',
          [(r[0], r[1:5]) for r in rows if r[1] != r[4]-r[3]])
    fam = [r for r in rows if r[2] + r[3] == r[1]]
    print(f'   rows with f+g = d: {len(fam)}')
    bad = 0
    for (ln, d, f, g, h, k, p, m) in fam:
        pred = (d+3, d, 8*d-2*f-1)
        inrange = d >= max(6, 2*f+1)
        ok = (k, p, m) == pred
        if inrange and not ok: bad += 1
        print(f'   line {ln}: d={d} (f,g,h)=({f},{g},{h}) source (k,p,m)=({k},{p},{m}); Thm 2 predicts {pred}; '
              f'{"in range" if inrange else "OUT of range d<6"}; {"AGREE" if ok else "DIFFER"}')
    print('   in-range disagreements:', bad)
    # which in-range pairs (f,d) with d <= dmax are present / absent
    ds = sorted(set(r[1] for r in rows))
    print('   d values present:', ds)
    present = {(r[2], r[1]) for r in fam}
    for dmax in (10, 25):
        want = [(f, d) for d in range(6, dmax+1) for f in range(1, (d-1)//2+1)]
        print(f'   boundary pairs (f,d), d<={dmax}, in Thm-2 range: {len(want)}; absent from this file:',
              [x for x in want if x not in present])
