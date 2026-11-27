/-
scripts/PrintAxioms.lean — axiom audit of the headline theorems:
    lake build && lake env lean scripts/PrintAxioms.lean
Every line must print exactly
    'X' depends on axioms: [propext, Classical.choice, Quot.sound]
(no sorryAx, no other axiom).
-/
import PiSeriesSolution

#print axioms ramanujan_series_four_div_pi
#print axioms ramanujan_series_sixteen_div_pi
#print axioms sun_series_thirty_two_sqrt_two_div_pi
#print axioms bauer_series_two_div_pi
#print axioms sun_series_sixty_four_sqrt_two_div_pi_sq
