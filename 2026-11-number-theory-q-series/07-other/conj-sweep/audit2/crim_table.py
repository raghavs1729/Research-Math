"""Parse the 17-state table printed in ../R07-crim.md (NOT the JSON).  Referee code (audit2)."""
import re, os
def load(path=os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'R07-crim.md')):
    rows = {}
    for line in open(path):
        m = re.match(r'^\|\s*(\d+)\s*\|\s*`([^`]*)`\s*\|\s*(yes|no)\s*\|\s*(\d+)\s*\|\s*(\d+)\s*\|\s*$', line)
        if m:
            q, rep, a, e, n = m.groups()
            rep = '' if rep == 'ε' else rep
            assert int(q) not in rows
            rows[int(q)] = (rep, a == 'yes', int(e), int(n))
    assert sorted(rows) == list(range(17)), sorted(rows)
    E = [rows[q][2] for q in range(17)]; N = [rows[q][3] for q in range(17)]
    ACC = [rows[q][1] for q in range(17)]; REP = [rows[q][0] for q in range(17)]
    return E, N, ACC, REP
if __name__ == '__main__':
    E, N, ACC, REP = load()
    print(17)
    for q in range(17): print(int(ACC[q]), E[q], N[q])
