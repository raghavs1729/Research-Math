# read-only inspection of the downloaded arXiv source (run with python3 -I); no extraction to disk
import sys, tarfile, gzip, io
p = sys.argv[1]
try:
    tf = tarfile.open(p, 'r:*')
    for mi in tf.getmembers():
        print(mi.name, mi.size, 'file' if mi.isfile() else 'other')
except tarfile.ReadError:
    data = gzip.open(p).read()
    print('single gzip file, bytes', len(data)); print(data[:300])
