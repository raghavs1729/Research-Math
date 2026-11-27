# read-only: print lines of a tar member matching a regex, or a line range (python3 -I)
import sys, tarfile, re
p, member, mode = sys.argv[1], sys.argv[2], sys.argv[3]
tf = tarfile.open(p, 'r:*')
txt = tf.extractfile(member).read().decode('utf-8', 'replace').split('\n')
if mode == 'grep':
    rx = re.compile(sys.argv[4]); ctx = int(sys.argv[5]) if len(sys.argv) > 5 else 0
    hit = set()
    for i, l in enumerate(txt):
        if rx.search(l):
            for k in range(max(0, i-ctx), min(len(txt), i+ctx+1)): hit.add(k)
    prev = None
    for k in sorted(hit):
        if prev is not None and k != prev+1: print('--')
        print(f'{k+1}: {txt[k]}'); prev = k
elif mode == 'range':
    a, b = int(sys.argv[4]), int(sys.argv[5])
    for k in range(a-1, min(b, len(txt))): print(f'{k+1}: {txt[k]}')
elif mode == 'wc':
    print(len(txt), 'lines')
