#!/usr/bin/env python3
"""Which classical congruence excludes tau(p^(d-1)) = eps*l, for l in {461, 599} and every odd prime
d | l^2 - 1 (p an odd prime not in {3, 5, 7, 23, 691}).  Independent of tau_p4.gp (Python)."""
def allowed(d, m):
    R = set()
    if m in (691, 2048):
        for a in range(1, m):
            if m == 2048 and a % 2 == 0: continue
            R.add(sum(pow(a, 11 * i, m) for i in range(d)) % m)
    elif m == 729:
        for a in range(1, 729):
            if a % 3: R.add(pow(a, -610 * (d - 1), 729) * sum(pow(a, 1231 * i, 729) for i in range(d)) % 729)
    elif m == 125:
        for a in range(1, 125):
            if a % 5: R.add(pow(a, -30 * (d - 1), 125) * sum(pow(a, 71 * i, 125) for i in range(d)) % 125)
    elif m == 7:
        for a in range(1, 7): R.add(pow(a, d - 1, 7) * sum(pow(a, 3 * i, 7) for i in range(d)) % 7)
    elif m == 23:
        for a in range(1, 23):
            qr = pow(a, 11, 23) == 1
            for t in ([2, -1] if qr else [0]):
                u0, u1 = 1, t % 23
                for _ in range(2, d): u0, u1 = u1, (t * u1 - pow(a, 11, 23) * u0) % 23
                R.add(u1 % 23)
    return R
def odd_prime_divisors(n):
    out, p = [], 3
    while n % 2 == 0: n //= 2
    while p * p <= n:
        if n % p == 0:
            out.append(p)
            while n % p == 0: n //= p
        p += 2
    if n > 1: out.append(n)
    return out
MODS = [(691, "691"), (2048, "2^11"), (729, "3^6"), (125, "5^3"), (7, "7"), (23, "23")]
print("| value | d | excluded modulo | not excluded modulo |")
print("|---|---|---|---|")
for l in (461, 599):
    for eps in (1, -1):
        for d in odd_prime_divisors(l * l - 1):
            ex = [name for m, name in MODS if (eps * l) % m not in allowed(d, m)]
            ok = [name for m, name in MODS if (eps * l) % m in allowed(d, m)]
            print(f"| {eps*l} | {d} | {', '.join(ex) if ex else '**none**'} | {', '.join(ok) if ok else '-'} |")
