"""Shared helpers: pinned zero files, characters, the certified sum S_a(u)."""
import hashlib, os
from flint import arb, acb, ctx, dirichlet_char

HERE = os.path.dirname(os.path.abspath(__file__))
EPS = "1e-30"


def sha256_file(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def pinned_hash(name):
    for line in open(os.path.join(HERE, "SHA256SUMS")):
        parts = line.split()
        if len(parts) == 2 and parts[1] == name:
            return parts[0]
    raise KeyError("%s not pinned in SHA256SUMS" % name)


def load_zeros(q, check_hash=True):
    """[(label, dstr)]; refuses a file whose sha256 is not the pinned one."""
    name = "zeros_q%d.txt" % q
    path = os.path.join(HERE, name)
    if check_hash:
        h = sha256_file(path)
        if h != pinned_hash(name):
            raise ValueError("sha256 of %s is %s, not the pinned value" % (name, h))
    out = []
    for line in open(path, encoding="ascii"):
        if line.startswith("#") or not line.strip():
            continue
        lab, d = line.split()
        out.append((lab, d))
    return out


def chi_arg_table(q, a):
    """label -> exact rational e with chi(a) = exp(2 pi i e) (label 'zeta' -> 0)."""
    from fractions import Fraction
    tab = {"zeta": Fraction(0)}
    for n in range(2, q):
        chi = dirichlet_char(q, n)
        G = chi.group()
        tab[str(n)] = Fraction(int(chi.chi_exponent(a)), int(G.exponent()))
    return tab
