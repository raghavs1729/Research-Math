/-
scripts/PrintAxioms.lean — axiom audit of the headline theorems:
    lake build && lake env lean scripts/PrintAxioms.lean
Every line must print exactly
    'X' depends on axioms: [propext, Classical.choice, Quot.sound]
(no sorryAx, no other axiom).
-/
import RamanujanSolution

#print axioms ramanujan_prime_counting_inequality
#print axioms ramanujan_gap_asymptotic
#print axioms primeCounting_five_term_expansion
#print axioms chebyshev_theta_sub_le
#print axioms ramanujan_of_schoenfeld_bound
#print axioms ramanujan_sharp_under_RH_reduction
