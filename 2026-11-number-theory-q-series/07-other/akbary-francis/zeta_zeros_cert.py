#!/usr/bin/env python3
"""
zeta_zeros_cert.py T -- zeros of the Riemann zeta function with 0 < gamma < T.

 * candidates: Arb's acb.zeta_zeros (acb_dirichlet_zeta_zeros: rigorous isolation of
   the n-th zero by Turing's method);
 * [Z] independently, for every listed decimal d the Hardy Z-function (Arb,
   character of modulus 1) has certified opposite signs at d -+ 1e-30; windows
   disjoint and increasing;
 * [N] the number K of windows below T equals Arb's certified count
   arb.zeta_nzeros(T) = N(T) (zeros in the critical strip with 0 < Im < T, with
   multiplicity).  Hence exactly one simple zero per window and no others.
Output: state/zeros_zeta.json {"T": T, "zeros": [...]}.
"""
import json, os, sys
from flint import arb, acb, ctx, dirichlet_char

HERE = os.path.dirname(os.path.abspath(__file__))
EPS = "1e-30"


def sign_of(x):
    return 1 if x > 0 else (-1 if x < 0 else 0)


def main():
    T = int(sys.argv[1])
    ctx.prec = 224
    N = arb(T).zeta_nzeros()
    K = int(N.unique_fmpz())
    chi = dirichlet_char(1, 1)
    zs = acb.zeta_zeros(1, K + 1)
    out = []
    prev = None
    e = arb(EPS)
    for z in zs[:K]:
        dstr = arb(z.imag.mid()).str(45, radius=False)
        d = arb(dstr)
        assert z.imag.overlaps(d + arb(0, 1e-40)) or abs(float((z.imag - d).mid())) < 1e-40
        sm = sign_of(chi.hardy_z(acb(d - e)).real)
        sp = sign_of(chi.hardy_z(acb(d + e)).real)
        assert sm != 0 and sm == -sp, "no certified sign change at %s" % dstr
        if prev is not None:
            assert prev + e < d - e
        prev = d
        out.append(dstr)
    assert prev + e < T and zs[K].imag > T
    json.dump({"T": T, "N(T)": K, "zeros": out}, open(os.path.join(HERE, "state", "zeros_zeta.json"), "w"))
    print("zeta: %d zeros with 0 < gamma < %d; N(T) = %d (arb.zeta_nzeros); all windows certified; first %s last %s"
          % (len(out), T, K, out[0][:20], out[-1][:20]))


if __name__ == "__main__":
    main()
