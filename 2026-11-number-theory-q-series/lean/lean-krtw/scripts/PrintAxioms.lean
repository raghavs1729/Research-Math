/-
scripts/PrintAxioms.lean — axiom audit of the headline theorems:
    lake build && lake env lean scripts/PrintAxioms.lean
Every line must print exactly
    'KRTW.X' depends on axioms: [propext, Classical.choice, Quot.sound]
(no sorryAx, no other axiom).
-/
import KRTWSolution

#print axioms KRTW.krtw_sum_eq_melzer_sum_refined
#print axioms KRTW.krtw_sum_eq_melzer_sum
#print axioms KRTW.coeff_krtwTerm_eq_zero
#print axioms KRTW.coeff_melzerTerm_eq_zero
#print axioms KRTW.krtw_identity
#print axioms KRTW.krtw_identity_theta
#print axioms KRTW.theta_eq_triple_product
#print axioms KRTW.melzer_identity
#print axioms KRTW.coeff_qPochInf
