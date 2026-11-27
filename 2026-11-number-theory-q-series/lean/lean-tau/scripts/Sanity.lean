/-
scripts/Sanity.lean — sanity checks (not part of the proof):
    lake env lean scripts/Sanity.lean
1. `heckeSeq` reproduces known values of Ramanujan's tau: τ(2) = −24, τ(4) = −1472, τ(8) = 84480,
   τ(16) = 987136; τ(3) = 252, τ(81) = 1665188361.
2. Non-vacuity of the sieve checker: for the planted solution (U, V) = (2^11, 3^2) of
   U² − 3UV + V² = 4139089 the certificate check must FAIL (the class n ≡ 0 survives).
3. The seed search is not vacuous: `seedsOK` fails if a seed is removed.
-/
import Tau

open Tau

example : heckeSeq (-24) 2 2 = -1472 := by decide +kernel
example : heckeSeq (-24) 2 3 = 84480 := by decide +kernel
example : heckeSeq (-24) 2 4 = 987136 := by decide +kernel
example : heckeSeq 252 3 4 = 1665188361 := by decide +kernel
example : Qf ((2 : ℤ) ^ 11, 9) = 4139089 := by decide +kernel
example : certOK [(2048, 9)] [3, 7, 43, 89, 199, 967] 88 = false := by decide +kernel
example : seedsOK 461 12 50 [(-1, -23), (-1, 20), (1, -20)] = false := by decide +kernel
example : seedsOK (-599) 24 50 [(-24, -47), (-24, -25), (24, 25)] = false := by decide +kernel
