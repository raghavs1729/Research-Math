# AUDIT-3, defect D1: under the LITERAL set-decomposition reading of R06's definition ("S consists of k entries a_i and
# all integers b_i + j m, lengths minimal", with no order condition a_k < b_1) the least k for T is 4, not d+3:
#   T = {f, d-f, 3d-f-1, 4d-f}  u  { b + j m : j >= 0, b in [2d-f, 3d-f-2] u {10d-3f-2} }   (p = d residue classes).
# The source's k (number of 1s in the minimal preperiod of the characteristic sequence) is d+3 (proof_check.py, Part D).
bad = []; n = 0
for f in range(1, 31):
    for d in range(max(6, 2*f+1), 151):
        m = 8*d-2*f-1; X = 10*d-3*f-2 + 12*m - 1; n += 1
        T = {f, d-f} | set(range(2*d-f, 3*d-f)) | {4*d-f}
        for j in range(12): T |= set(range(10*d-3*f-2+j*m, 11*d-3*f-2+j*m))
        A = {f, d-f, 3*d-f-1, 4*d-f}
        base = list(range(2*d-f, 3*d-f-1)) + [10*d-3*f-2]
        P = {b + j*m for b in base for j in range(14) if b + j*m <= X}
        ok = (A | P) == T and not (A & P) and len(base) == d
        # each of the four a's is forced: its class is not eventually-from-there contained in T
        forced = all((a + m) not in T for a in A)
        if not (ok and forced): bad.append((f, d))
print(f"pairs f <= 30, max(6,2f+1) <= d <= 150: {n}; decomposition with k = 4, p = d fails for: {bad}")
