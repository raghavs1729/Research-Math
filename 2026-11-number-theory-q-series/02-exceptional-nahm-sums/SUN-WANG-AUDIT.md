# SUN-WANG-AUDIT — numerical audit of K. Sun, H. Wang, "Dynkin diagrams, generalized Nahm sums and 2d CFTs" (arXiv:2604.00847v1) and of the (E8,T1) identities of Y. Yoshida (arXiv:2608.08780v1)

Date: 2026-11-26. Scripts: `code/sw_lib.py` (exact q-series toolkit), `code/sw_audit.py`, `code/sw_yoshida.py`,
`code/sw_certify_T1.py`, `code/sw_E7T1.py`. Logs: `certificates/sw_audit.txt`, `sw_yoshida.txt`, `sw_certify_T1.txt`,
`sw_E7T1.txt`. Everything is exact integer arithmetic on truncated series; "holds mod q^N" means: all coefficients
of exponents < N agree (exponent including the prefactor q^C).

## 0. Summary

* **Misprints found (3).** (i) (3.27)/(3.28): the Ramond weights 1/32 and 5/32 are interchanged (and one identity
  needs a factor 2): correct is f_4 = χ_{R,1/32}χ^F_{R,1/16}, f_3 = 2χ_{R,5/32}χ^F_{R,1/16}. (ii) Text before (3.30): "The
  Nahm sum associated with (T1,A3) is identical to f_1" should read f_2 ((3.30) itself is correct). (iii) (3.6)
  prints the same line for r = 2n and r = 2n+1; this is harmless: the formula χ^{NS}_0 − χ^{NS}_{⌊r/2⌋+1} is correct
  for both parities (tested r = 2,…,9; the Kac labels differ, the weights do not; proved for all r in
  `proofs/T1-BC.md`, Corollary C, from Warnaar's identity).
* **No other formula fails.** All of (1.8)–(1.12), (2.1)–(2.8), (3.2), (3.5)–(3.41) and the c-column of Table 1
  were confirmed to the orders listed below (137 series comparisons in `sw_audit.txt`: 133 hold; the 4 failures
  are the two printed forms in (i), the statement (ii), and Zagier's old misprint C = 19/24, which Sun–Wang
  themselves correct).
* **Convention needed for the N = 1 formulas:** all formulas with coefficient 2 in front of a Ramond character
  ((1.8)–(1.11), (3.30)) hold when χ_R is normalised with leading coefficient 1,
  χ_R = q^{h−c/24} ∏(1+q^n)/(1−q^n) Σ_k(q^{k(kpp′+rp′−sp)/2} − q^{(kp+r)(kp′+s)/2}); the paper does not state this.
* **New in this project** (beyond the audit): `proofs/T1-BC.md` (closed formula for (T1,B_r), all r; (3.6) and
  its companions for all r as a corollary of known identities), `proofs/T1-A-even.md` ((T1,A_{2k}) = Z_k
  parafermions ⊗ U(1)_{k(2k+1)}, all k), and a conjectural formula for (E7,T1) (§4, tested only).

## 1. Conventions used in the tests

* Nahm sums: q^C Σ_n q^{½nᵀ(AD)n + B·n}/∏(q^{d_i};q^{d_i})_{n_i}, A = C(X)⊗C(Y)⁻¹, D = D(X)⊗D(Y) (`fusion2.pair`); the printed
  quadratic forms of (2.4), (3.14), (3.16), (3.37) and the matrices of (T1,C_r), (A1,C_r), (G2,·) were compared with
  A(X,Y)D(X,Y) (all agree).
* Virasoro M(p,p′): χ_{r,s} = q^{h−c/24}(q)_∞⁻¹Σ_k(q^{k(kpp′+rp′−sp)} − q^{(kp+r)(kp′+s)}). N = 1 SM(p,p′): NS as in
  `proofs/T1-BC.md` §4, R as above. "Effective" models: same characters, labelled by h − h_min (h_min over NS ∪ R).
  The weights quoted by Sun–Wang (10 NS weights of SM_eff(14,4), 9+9 weights of SM(8,6), 21+21 of SM(84,2), 2r+1
  NS and R weights of SM(8r+4,2)) are reproduced by this convention.
* U(1)_k: χ_j = η⁻¹Σ_n q^{k(n+j/2k)²}, weight j²/(4k) (this is the convention under which (3.5), (3.18), (3.24), (3.25),
  (3.35) hold; e.g. U(1)_3 has weights 0, 1/12, 1/3, 3/4).
* Z_k parafermions: η·c^ℓ_m with the Kac–Peterson formula for the string functions (checked against Ising and
  Potts characters).
* Free fermion: χ^F_{NS,0} = q^{−1/48}∏(1+q^{n−1/2}), χ^F_{R,1/16} = q^{1/24}∏(1+q^n).
* D_2A (for (3.33), (3.34)): the extension of Ising ⊗ M(5,4) by (1/2,3/2): χ_0 = [0,0]+[1/2,3/2],
  χ_{3/5} = [0,3/5]+[1/2,1/10], χ_{1/2} = [1/16,7/16] (= [0,3/2]+[1/2,0] as q-series), χ_{1/10} = [1/16,3/80]
  (= [0,1/10]+[1/2,3/5]). This is my reading of "[14, Section 5.3]"; with it both formulas hold.
* Zagier's rank-two list: order and (B,C) taken from L. Wang, arXiv:2210.10748, Tables 4, 8, 10–12 (I did not
  have Zagier's paper); consistent with Sun–Wang's own remarks (f_1 = 2f_4, B = 0 is f_2, C = 19/96).

Status codes: (a) the sum side is a theorem in the literature (reference as given by Sun–Wang unless stated),
(b) proved in this project, (c) open as far as I know. In case (a) the rewriting in terms of CFT characters is
an identity between explicit theta quotients; it was only tested numerically here unless stated.

## 2. Table

| Eq. | Content (transcribed) | Result | Order | Status |
|---|---|---|---|---|
| (1.8) | (T2,A1), A = (1/2 −1/2; −1/2 1), B = 0, C = −5/84: f_1 = NS_0+NS_{1/2}+NS_{5/2}+NS_5+2R_{1/4}+2R_{5/4} of SM_eff(84,2) | holds | q^{300} | (a) Zagier Ex. 7, L. Wang 2210.10748 |
| (1.9) | f_2 (B = (1/2,−1/2), C = 1/21) = NS_{5/14}+NS_{6/7}+NS_{13/7}−NS_{20/7}+2R_{3/28}+2R_{87/28} | holds | q^{300} | (a) same |
| (1.10) | f_3 (B = (1/2,0), C = 1/84) = NS_{1/14}−NS_{15/14}+NS_{11/7}+NS_{57/14}+2R_{23/28}−2R_{135/28} | holds | q^{300} | (a) same |
| (1.11) | (T1,D_r): f = NS_{j=1}+NS_{j=2r+1}+2R_{j=r+1} of SM_eff(8r+4,2), C = −r/(8(2r+1)) (weights 0, r/2, r/8) | holds, r = 3,…,8 | q^{100} (r ≤ 6), q^{50} (r = 7,8) | (a) Wang–Wang 2607.08606 |
| (1.12) | product form of (1.11) | holds, r = 3,…,8 | same | (a) same |
| (2.1) | Andrews–Gordon, all 1 ≤ s ≤ r+1 | holds, r ≤ 4 | q^{150} | (a) |
| (2.2), (2.3) | (T1,T_{2r}) (Melzer/BIS/Warnaar), (T1,T_{2r−1}) (Warnaar) | hold, rank ≤ 7 | q^{150} | (a) |
| (2.4) | (A1,B_r), Baker–Kanade–Russell–Sadowski Thm 3.4 | holds, r = 2,…,5 | q^{150} | (a) |
| (2.5) | Bressoud, all s | holds, r = 2,3,4 | q^{150} | (a) |
| (2.6)+(2.7) | (A1,D_r) (sum of the two parities) | holds, r = 3,…,6 | q^{150} | (a) Warnaar |
| (2.8) | (T1,C_r) = (q^{(r+1)/2},q^{(r+2)/2},q^{r+3/2};q^{r+3/2})_∞/((−q)_∞(q^{1/2};q^{1/2})_∞) | holds, r = 2,…,7 (and r ≤ 9 in `sw_certify_T1`) | q^{150} | (a) Warnaar, Discrete Math. 272 (2003), (5.14) [= (5.15) of arXiv:math/0110307v2, k ≥ 1; read in full text]; re-proved in `proofs/T1-BC.md` |
| §2, end | (T1,G2): Σ q^{i²+3ij+3j²}/((q)_i(q³;q³)_j) = 1/(q,q³,q⁶,q⁸;q⁹)_∞; (G2,T1): Wang–Wang Conj. 3.6 | hold | q^{150} | called conjectures by S–W; proved since (Mizuno 2609.09653, Xia 2609.08816, per `STATUS-MAP.md`; not re-read) |
| (3.2) | (A1,E8) = Σ q^{2m²}/(q)_{2m} (= χ_0 of M(4,3)) | holds | q^{200} | (a) Warnaar–Pearce |
| (3.5) | (A1,C_r) = (χ_0 − χ_{r+1})^{U(1)_{4(r+1)}} (labels j = 0 and j = 4(r+1)) | holds, r = 2,…,7 | q^{200} | (a) Bressoud (2.5) + triple product |
| (3.6) | (T1,C_r) = χ^{NS}_0 − χ^{NS}_{⌊r/2⌋+1} of SM_eff(4r+6,4); same line printed for r = 2n and 2n+1; "checked up to r = 7" | holds for r = 2,…,9; labels (1,r+1),(1,3r+5) (r even), (1,r+2),(1,3r+4) (r odd) | q^{200} (r ≤ 6), q^{120} (r = 7,8), q^{80} (r = 9) | (a)+(b): theorem for all r: (2.8) + dissection, `proofs/T1-BC.md` Cor. C |
| (3.7)–(3.9) | r = 2, B = (0,0),(0,1),(1,1), C = −3/56, 1/56, 9/56: NS_0−NS_2, NS_{1/14}−NS_{15/14}, NS_{3/14}−NS_{45/14} | hold | q^{200} | (a) Mizuno 2305.02267; all r: Wang–Wang 2402.06253 Thm 4.1 + Cor. C |
| (3.10)–(3.13) | r = 3, B = (0,0,0),(0,0,1),(0,1,1),(1,1,2), C = −1/18, 0, 1/9, 5/18: NS_0−NS_2, NS_{1/18}−NS_{55/18}, NS_{1/6}−NS_{7/6}, NS_{1/3}−NS_{13/3} | hold | q^{200} | (a) Wang–Wang 2402.06253 |
| (3.14), (3.15) | (A1,F4) = 1/(q²,q³,q⁴,q^{10},q^{11},q^{12};q^{14})_∞ = q^{1/28}(χ_0−χ_5)^{M(7,6)} | hold | q^{300} | (a) BKRS (1.7) |
| (3.16), (3.17) | (T1,F4) = 1/(q,q²,q⁴,q⁶,q⁸,q⁹;q^{10})_∞ = q^{1/20}χ_0^{M_eff(5,2)}(χ_0−χ_3)^{M(6,5)} | hold | q^{300} | (a) BKRS Thm 1.5 |
| (3.18) | (A1,G2) (AD = (4 6; 6 12), d = (1,3)) = (χ_0−χ_1−χ_4+χ_9)^{U(1)_{36}} (j = 0,12,24,36) | holds | q^{300} | (a) Wang–Wang 2308.14728 |
| (3.19) | (T1,E6) = χ_0^{M_eff(5,2)}(χ_0+χ_3+2χ_{2/3})^{M(6,5)} | holds | q^{200} | (b) `proofs/E8-E6.md` |
| (3.20) | (T1,E8) = χ_0^{M_eff(11,2)} | holds | q^{200} | (b) `proofs/E8-E6.md` |
| (3.21) | (E8,T1) = vacuum character of T_{10}M_eff(11,2); "1+120q+1660q²+12320q³+68210q⁴" | coefficients agree; = Hecke image T_{10} (see §3) | q^{26} (Hecke); series to q^{20} in the archived log (q^{40} in a separate run) | (c) open |
| (3.22), (3.23) | A = 1/2: χ_0+χ_{1/4}, χ_{1/20}+χ_{4/5} of M_eff(5,3) | hold | q^{500} | (a) classical (Rogers/Slater); Zagier |
| (3.24) | (T1,A2) = (χ_0+2χ_{1/3})^{U(1)_3} | holds | q^{500} | (a) Zagier Ex. 1 |
| (3.25) | (A2,T1) = (χ_0+χ_1)^{U(1)_4} = η⁻¹Σ q^{m²} | holds | q^{500} | (a) Zagier Ex. 1 |
| (3.26) | (T2,T1): f_2 = χ^{SM_eff(8,2)}_{NS,0}χ^F_{NS,0} | holds | q^{300} | (a) Zagier Ex. 3; L. Wang 2210.10748 |
| (3.27) | f_3 = χ_{R,1/32}χ^F_{R,1/16} | **fails as printed** (f_3 = 2q^{1/6}+…, right side q^{1/24}+…) | first term | corrected: f_3 = 2χ_{R,5/32}χ^F_{R,1/16}, holds mod q^{300} |
| (3.28) | f_4 = χ_{R,5/32}χ^F_{R,1/16} | **fails as printed** (f_4 = q^{1/24}+…, right side q^{1/6}+…) | first term | corrected: f_4 = χ_{R,1/32}χ^F_{R,1/16}, holds mod q^{300} |
| (3.29) | f_5 = χ_{NS,1/4}χ^F_{NS,0} (C = 19/96); "f_1 = 2f_4" | hold | q^{300} | (a) |
| p. 13 | (T1,A3)/Ex. 9: f_1 = NS_{3/14}+NS_{5/7}+2R_{5/56}, f_2 = NS_0+NS_{3/2}+2R_{3/8}, f_3 = NS_{1/14}−NS_{15/14}+2R_{53/56} of SM_eff(28,2) | hold | q^{300} | (a) Zagier Ex. 9; L. Wang |
| (3.30) | (T1,A3), C = −3/56, = NS_0+NS_{3/2}+2R_{3/8} | holds; **but** "(T1,A3) is identical to f_1" is a misprint for f_2 | q^{300} | (a) |
| (3.31), (3.32) | (A1,A2)/Ex. 10: f_1 = 2χ_{1/15}+χ^{sub}_{2/5}, f_3 = χ^{sub}_0+2χ_{2/3}, χ^{sub}_0 = χ_0+χ_3, χ^{sub}_{2/5} = χ_{2/5}+χ_{7/5} | hold | q^{300} | (a) Cao–Rosengren–Wang (see `E8-E6.md`) |
| (3.33), (3.34) | (A2,A1)/Ex. 11: f_1 = 3χ_{1/10}+χ_{3/5}, f_3 = χ_0+3χ_{1/2} of D_2A | hold (with the D_2A characters of §1) | q^{300} | (a) Vlasenko–Zwegers |
| (3.35), (3.36) | (A1,A3) = (χ_0+χ_{3/4})^{U(1)_3} = χ^{Z_4}_0+χ^{Z_4}_1+2χ^{Z_4}_{3/4} | hold | q^{300} | (a) |
| (3.37)–(3.39) | (A1,B3), B = (0,0,0), (0,−1,−1), (−1,0,1/2), C = −5/96, −1/48, 1/24: NS_0−NS_3, NS_{1/32}−NS_{33/32}, R_{3/32}−R_{67/32} of SM(8,6) | hold | q^{200} | (a) Wang–Wang 2407.21725 |
| (3.40) | (A1,A5) = χ^{Z_6}_0+χ^{Z_6}_{3/2}+2χ^{Z_6}_{5/6}+2χ^{Z_6}_{4/3} | holds | q^{60} | (a) Lepowsky–Primc |
| (3.41) | (A1,A5) = (NS_0+NS_3+2NS_{5/6})^{SM(8,6)} | holds | q^{60} | follows from (3.40) + an identity between characters (tested only) |

**Table 1.** Column c: all rows agree with (1.6) (111 instances; families for r ≤ 7, k ≤ 3). Rows without a
displayed formula:

| Row | What was tested / what is known |
|---|---|
| (A1,A_r), Z_{r+1} parafermions | Nahm sum = Σ_Q χ^{Z_{r+1}}_{0,2Q}, r ≤ 6, q^{100} (r ≤ 4), q^{66}; sector by sector for r+1 ≤ 7 (`sw_certify_T1`); (a) Lepowsky–Primc |
| (A1,B_r), CFT "?" | (2.4) |
| (A1,D_r), U(1)_r | (2.6)+(2.7) |
| (A1,E6), M_sub(7,6) | Nahm sum = (χ_0+χ_5+2χ_{4/3})^{M(7,6)}, holds mod q^{200}; conjectured by Kedem–Klassen–McCoy–Melzer; proof status not checked by me |
| (A1,E7), M(5,4) | Nahm sum = (χ_0+χ_{3/2})^{M(5,4)}, holds mod q^{200}; KKMM; Warnaar–Pearce (E7 identities) |
| (T1,A_r), CFT "?" | **even r = 2k: Z_k parafermions ⊗ U(1)_{k(2k+1)}, proved** (`proofs/T1-A-even.md`; uses Lepowsky–Primc); odd r ≥ 5: open, reduced rank-k form only |
| (T1,B_r), CFT "?" | **closed formula for all r, proved** (`proofs/T1-BC.md`, Thm B); equal to the right sides of Du–Ji–Shen–Xu 2602.14483 (1.7)/(1.9) |
| (T1,E7), CFT "?" | (1/(2η))Σ_k(q^{3k²/4}+q^{7k²/4}), proved (`proofs/E7.md`); re-tested mod q^{200} |
| (T1,G2) | Kanade–Russell mod 9 identity, see above |
| (T1,T_{2r}), (T1,T_{2r+1}) | (2.2), (2.3) |
| (A_{2r−1},T_k), osp(1|2r)_k | not tested (no explicit formula in the paper) |
| (E6,T1), (E7,T1) ("no good CFT interpretation") | (E7,T1): conjectural formula, §4; (E6,T1): first coefficients only, no formula found |

## 3. Yoshida, arXiv:2608.08780: the (E8,T1) system (matrix C(E8))

Claims ((2.68)–(2.73)): for h ∈ {0, 5/11, 8/11, 10/11, 12/11} and his charge vectors Q_h,
q^{γ_h}Σ_{m∈N⁸} q^{½mᵀC(E8)m − Q_h·m}/∏(q)_{m_i} = χ_h of T_{10}M_eff(11,2); vacuum "exact" (citing Sun–Wang (3.21), which is
numerical), the other four "through the displayed order" (q⁶).

* (2.4), (2.72): every displayed coefficient is reproduced (4, resp. 7 per series).
* The five characters of T_{10}M_eff(11,2) were computed from the Hecke formula of Harvey–Wu for conductor 33,
  b^{(10)}_i(M) = Σ_{a | gcd(10,M)} (10/a) ε_a(i) b_{π_a(i)}(10M/a²) (b_i(n) = coefficient of q^{n/33} in the M_eff(11,2)
  characters), with the signs ε_a(i) of the Galois action σ_{a⁻¹} on the S-matrix of M(11,2) (no fitted
  parameter; the formula was written from memory, not compared with Harvey–Wu's paper). **All five Nahm sums
  agree with the Hecke images mod q^{26}** (about 26 coefficients each; `certificates/sw_yoshida.txt`).
* Bilinear relation (1.9)/(1.10): Σ_h χ_h(Nahm side)·χ^∨_{σ(h)} = E_4/η⁸ with his pairing
  (0,5/11,8/11,10/11,12/11) ↔ (0,6/11,3/11,1/11,10/11) of M_eff(11,2) characters: holds mod q^{24}. This check does
  not use the Hecke formula.
* (3.75)–(3.77): Q̃^∨_h = C(E8)⁻¹Q_h reproduced exactly; (3.78) and the identification of the dual sums are
  theorems (`proofs/E8-E6.md`, Corollary 1).
* Status: the five (E8,T1) identities (including the vacuum one, = Sun–Wang (3.21)) are **not proved** by
  anybody as far as I know; they are now tested to q^{26}. Fusion does not apply (C(E8) has no defusion).

## 4. A conjectural formula for (E7,T1) (tested only)

    Σ_{n∈N⁷} q^{½nᵀC(E7)n}/∏(q)_{n_i} = ( Θ_{E6}(q) + Θ_{A6}(q) ) / ( 2 (q;q)_∞⁶ ),     Θ_L = Σ_{x∈L} q^{|x|²/2},

i.e. q^{−1/4}·(Nahm sum) = ½(χ_0^{(E6)_1} + χ_0^{(A6)_1}) (c = 6; first coefficients 1 + 63q + 609q² + 3612q³ + …,
63 = (78+48)/2). Found by analogy with (T1,E7) = (Θ_{⟨6⟩} + Θ_{⟨14⟩})/(2η) (`proofs/E7.md`): ⟨6⟩ and ⟨14⟩ are the
orthogonal complements of E6 and of A6 in the E7 lattice. **Holds mod q^{40}** (`certificates/sw_E7T1.txt`,
script `code/sw_E7T1.py`). No proof. For (E6,T1) (1 + 36q + 240q² + 1144q³ + …, 24 coefficients in the same log) I
found no formula (multiplying by (q)_∞^k, k ≤ 7, and by Rogers–Ramanujan products gives nothing sparse).

## 5. What "holds" does and does not mean

A comparison mod q^N is not a proof. For the rows marked (a) the sum = product (or sum = theta quotient)
identity is a theorem and the remaining identification with CFT characters is an identity between explicit
weight-0 theta quotients, which could be proved by a finite computation (valence formula); this was not done
here, except where a proof file is cited.
