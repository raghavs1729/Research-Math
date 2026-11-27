#!/usr/bin/env python3
"""Evaluation (Arb) of the analytic lower bounds L_A, L_B of TAU-RACE.md and of the tail constant T(H).

Constants are taken as rational upper bounds of the certified enclosures in ../data/constants.json
(checked below against that file).
  L_A(t) = (2/3)(1 - eps1*[t > Y1^2]) - sigma0sq - a0 t^-1/2 - ((4/5) c_psi + C_S) t^-1/4 - 2/t - K_S t^-3/2 - 2 c_theta E_odd(t)
  E_odd(t) = (3/4) t^-1/6 + (5/6) t^-3/10 + t^-5/14 max(0, (log t/log 2 - 5)/2)
  L_B(t)  = L_A(t) - T(H)
  T(H)    = min over H' < H of 2 c [ (3/2)(log(H'/2pi)+1)/H' + SX2/(3 H'^3) + SXI/(2 H'^2) + 6 log zeta(2)/H'^2 ],
            c = sqrt2/arctan((H-H')/sqrt2), SX2 = sum x_j^2, SXI = sum 1/(24 x_j), x_j in {3/2, 13/2, 7}.
"""
import os
DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data')
import json, sys
from flint import arb, ctx
from arbutil import dec_down, dec_up
ctx.prec = 100

C = json.load(open(DATA + "/constants.json"))
SIG = arb("0.1122522"); A0 = arb("0.170658"); CS = arb("0.70598"); CPSI = arb("1.4142136"); KS = arb("0.34176")
Y1 = arb(10)**7
EPS1 = arb("0.2")/Y1.log()**2            # BKLNW Cor. 11.2, k = 2 (X0 = 3 594 641 <= Y1)
CTH = 1 + EPS1
def _chk(name, ub):
    assert arb(C[name][1]) <= ub, name
_chk("sigma0sq", SIG); _chk("a0", A0); _chk("C_S", CS); _chk("c_psi", CPSI); _chk("K_S", KS)
assert C["Y1"] == 10**7 and C["theta_lt_y"] is True

def E_odd(t):
    t = arb(t)
    k7 = (t.log()/arb(2).log() - 5)/2
    k7 = k7.max(arb(0))
    return arb(3)/4*t**(-arb(1)/6) + arb(5)/6*t**(-arb(3)/10) + t**(-arb(5)/14)*k7

def L_A(t, beyond=None):
    """beyond: None -> decide from t; True -> value on (Y1^2, oo) (with the eps1 term)"""
    t = arb(t)
    if beyond is None: beyond = bool(t > Y1*Y1)
    main = arb(2)/3*(1 - (EPS1 if beyond else 0))
    return main - SIG - A0*t**(-arb(1)/2) - (arb(4)/5*CPSI + CS)*t**(-arb(1)/4) - 2/t - KS*t**(-arb(3)/2) - 2*CTH*E_odd(t)

def T_tail(H, Hp):
    H = arb(H); Hp = arb(Hp)
    c = arb(2).sqrt()/((H-Hp)/arb(2).sqrt()).atan()
    xs = [arb(3)/2, arb(13)/2, arb(7)]
    SX2 = sum(x*x for x in xs); SXI = sum(1/(24*x) for x in xs)
    I = arb(3)/2*((Hp/(2*arb.pi())).log()+1)/Hp + SX2/(3*Hp**3) + SXI/(2*Hp**2) + 6*arb(2).zeta().log()/Hp**2
    return 2*c*I

def T_best(H):
    best = None
    for i in range(1, 400):
        Hp = arb(H)*i/400
        if not (Hp > 7): continue
        v = T_tail(H, Hp)
        if best is None or v < best[0]: best = (v, Hp)
    return best

if __name__ == "__main__":
    print("eps1 =", EPS1.str(8), " c_theta =", CTH.str(10))
    print("limit of L_A: ", dec_down(arb(2)/3*(1-EPS1) - SIG, 6), "  (2/3 - sigma0sq =", dec_down(arb(2)/3 - SIG, 6), ")")
    for e in [4, 5, 5.5, 6, 6.5, 7, 8, 9, 10, 12, 14]:
        t = arb(10)**arb(str(e))
        print("L_A(10^%s) >= %s" % (e, dec_down(L_A(t), 6)))
    print("L_A(Y1^2 +) >= ", dec_down(L_A(Y1*Y1, beyond=True), 6))
    # threshold on a grid
    import math
    lo = None
    for j in range(4000, 7001):
        t = arb(10)**(arb(j)/1000)
        if L_A(t) > 0:
            lo = j; break
    print("L_A > 0 from t = 10^%.3f = %.0f on (grid, and L_A is nondecreasing on [527, Y1^2])" % (lo/1000, 10**(lo/1000)))
    for H in [15, 20, 30, 40, 50, 60, 80, 100, 150, 200]:
        v, Hp = T_best(H)
        lim = arb(2)/3*(1-EPS1) - SIG - v
        thr = None
        for j in range(4000, 12001, 5):
            t = arb(10)**(arb(j)/1000)
            if L_A(t) - v > 0: thr = j; break
        print("H = %4d  H' = %7.3f  T(H) <= %s   limit L_B >= %s   L_B(10^7) >= %s  L_B>0 from 10^%s" % (
            H, float(Hp.mid()), dec_up(v, 6), dec_down(lim, 6), dec_down(L_A(arb(10)**7) - v, 6), (thr/1000 if thr else None)))
    # ---- numbers quoted in TAU-RACE.md
    T60, Hp60 = T_best(60)
    print("---- quoted values (H = 60, H' = %s)" % Hp60.str(8))
    print("T(60) <= %s" % dec_up(T60, 6))
    for nm, t in [("10^6", arb(10)**6), ("9*10^6", 9*arb(10)**6), ("10^7", arb(10)**7)]:
        print("L_A(%s) >= %s    L_B(%s) >= %s" % (nm, dec_down(L_A(t), 6), nm, dec_down(L_A(t) - T60, 6)))
    print("L_A(Y1^2+) >= %s   L_B(Y1^2+) >= %s" % (dec_down(L_A(Y1*Y1, beyond=True), 6), dec_down(L_A(Y1*Y1, beyond=True) - T60, 6)))
    print("lim L_A >= %s   lim L_B >= %s" % (dec_down(arb(2)/3*(1-EPS1) - SIG, 6), dec_down(arb(2)/3*(1-EPS1) - SIG - T60, 6)))
    # Rosser-Schoenfeld fallback: eps = 1/log(10^7), c_theta = 1.01624
    e_rs = 1/Y1.log()
    print("fallback (Rosser-Schoenfeld): eps = %s, lim L_A >= %s, lim L_B >= %s" % (dec_up(e_rs, 6), dec_down(arb(2)/3*(1-e_rs) - SIG, 6), dec_down(arb(2)/3*(1-e_rs) - SIG - T60, 6)))
    # fallback variant evaluated in full: eps1 -> 1/log(10^7), c_theta -> 1.01624
    def L_fb(t, beyond):
        t = arb(t); cth = arb("1.01624")
        main = arb(2)/3*(1 - (e_rs if beyond else 0))
        return main - SIG - A0*t**(-arb(1)/2) - (arb(4)/5*CPSI + CS)*t**(-arb(1)/4) - 2/t - KS*t**(-arb(3)/2) - 2*cth*E_odd(t)
    print("fallback: L_A(10^6) >= %s, L_B(10^6) >= %s, L_B(Y1^2+) >= %s" % (dec_down(L_fb(arb(10)**6, False), 6), dec_down(L_fb(arb(10)**6, False) - T60, 6), dec_down(L_fb(Y1*Y1, True) - T60, 6)))
    for nm, t in [("9*10^6", 9*arb(10)**6), ("10^7", arb(10)**7)]:
        print("fallback: L_A(%s) >= %s, L_B(%s) >= %s" % (nm, dec_down(L_fb(t, False), 6), nm, dec_down(L_fb(t, False) - T60, 6)))
    print("fallback: L_A(Y1^2+) >= %s" % dec_down(L_fb(Y1*Y1, True), 6))
