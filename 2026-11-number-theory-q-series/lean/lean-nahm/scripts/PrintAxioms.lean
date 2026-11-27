/-
scripts/PrintAxioms.lean — axiom audit of the headline theorems:
    lake build && lake env lean scripts/PrintAxioms.lean
Every line must print exactly
    'Nahm.X' depends on axioms: [propext, Classical.choice, Quot.sound]
(no sorryAx, no other axiom; `native_decide` would add `Lean.ofReduceBool`).
-/
import NahmSolution

#print axioms Nahm.andrews_gordon
#print axioms Nahm.andrews_gordon_field
#print axioms Nahm.coeff_agTerm_eq_zero
#print axioms Nahm.invCartanE8_eq_explicit
#print axioms Nahm.e8Form_eq_e8Exp
#print axioms Nahm.coeff_e8Term_eq_zero
#print axioms Nahm.e8_vacuum
#print axioms Nahm.e8_vacuum_field
#print axioms Nahm.e8_companion_4
#print axioms Nahm.e8_companion_3
#print axioms Nahm.e8_companion_2
#print axioms Nahm.e8_companion_1
#print axioms Nahm.invCartanE7_eq_explicit
#print axioms Nahm.e7Form_eq_e7Exp
#print axioms Nahm.coeff_e7Term_eq_zero
#print axioms Nahm.coeff_e7ThetaTerm_eq_zero
#print axioms Nahm.e7_identity
#print axioms Nahm.e7_identity_bilateral
#print axioms Nahm.e7_identity_field
#print axioms Nahm.e7_identity_bilateral_field
#print axioms Nahm.coeff_qPochInf
-- internal theorems of interest
#print axioms Nahm.e8_theorem1
#print axioms Nahm.fusion
#print axioms Nahm.ag_theta_one
#print axioms Nahm.ag_theta_two
#print axioms Nahm.jtp_gen
#print axioms Nahm.cauchy
#print axioms Nahm.lemma4
#print axioms Nahm.prop5
#print axioms Nahm.prop6
