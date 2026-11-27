/-
scripts/PrintAxioms.lean — axiom audit of the headline theorems:
    lake build && lake env lean scripts/PrintAxioms.lean
Every line must print exactly
    'Tau.X' depends on axioms: [propext, Classical.choice, Quot.sound]
(no sorryAx, no Lean.ofReduceBool / Lean.trustCompiler, i.e. no native_decide).
-/
import TauSolution

#print axioms Tau.no_solution_461
#print axioms Tau.no_solution_neg599
#print axioms Tau.heckeSeq_four_ne
#print axioms Tau.tau_prime_pow_four_ne
