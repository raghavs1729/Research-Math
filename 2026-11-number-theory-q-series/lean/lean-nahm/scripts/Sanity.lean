/-
scripts/Sanity.lean — sanity test of the DEFINITIONS (not of the proof):
    lake build && lake env lean scripts/Sanity.lean

A computable mirror of the statement-level definitions of `Nahm/Defs.lean` (= those of
`NahmStatement.lean`), with integer coefficients truncated modulo `X^41`.  It calls the very
functions `Nahm.psum`, `Nahm.agExp` of the statement file and copies the formulas of `qPoch`, `qPoch4`,
`qPochInf`, `agSum`, `agProduct`, `e8Sum`, `e7Sum`, `e7Theta`, `e7ThetaOneSided` symbol by symbol.
Two liberties, both harmless and stated here:
  * the noncomputable exponents `e8Exp B n = ⌊½ nᵀ C(E8)⁻¹ n + B·n⌋₊` and `e7Exp n = ⌊4·½ nᵀ C(E7)⁻¹ n⌋₊`
    are computed with the explicit inverse matrices of the theorems `invCartanE8_eq_explicit`,
    `invCartanE7_eq_explicit` (proved in Lean), and `e8Form_eq_e8Exp`, `e7Form_eq_e7Exp` (proved)
    say that the floor is exact;
  * the box `{0,…,N}^r` is replaced by the tuples whose exponent is `≤ 40` (all matrix entries are
    positive, so the exponent is monotone in each entry and the others contribute nothing below `X^41`;
    for the Andrews–Gordon sum: tuples with `n_1 + ⋯ + n_{k-1} ≤ 6`, as `agExp ≥ N_1² ≥ 49` otherwise).
The `#guard`s check that the first 41 coefficients of each side of each identity are the integers
computed independently by `scripts/sanity_independent.py` (plain Python; Cartan matrices built from the
PAPER's Dynkin diagrams and inverted there with exact fractions; the paper's linear terms `B^{(i)}`, in
the paper's labelling; depth-first enumeration), plus negative controls.
-/
import Nahm.Defs

open Nahm

namespace NahmSanity

def N : Nat := 40

abbrev Ser := Array Int

def zero : Ser := Array.replicate (N + 1) 0
def one : Ser := zero.set! 0 1
/-- `X ^ e` -/
def xpow (e : Nat) : Ser := if e ≤ N then zero.set! e 1 else zero
def add (a b : Ser) : Ser := Array.zipWith (· + ·) a b
def sub (a b : Ser) : Ser := Array.zipWith (· - ·) a b
def smul (c : Int) (a : Ser) : Ser := a.map (c * ·)
def mul (a b : Ser) : Ser := Id.run do
  let mut c := zero
  for i in [0:N + 1] do
    if a.getD i 0 != 0 then
      for j in [0:N + 1 - i] do
        c := c.set! (i + j) (c.getD (i + j) 0 + a.getD i 0 * b.getD j 0)
  return c
/-- inverse of a series with constant term `1` -/
def inv (a : Ser) : Ser := Id.run do
  let mut b := zero.set! 0 1
  for n in [1:N + 1] do
    let mut s : Int := 0
    for i in [1:n + 1] do
      s := s + a.getD i 0 * b.getD (n - i) 0
    b := b.set! n (-s)
  return b
def prodL (l : List Ser) : Ser := l.foldl mul one

/-- mirror of `qPoch n = ∏ i ∈ range n, (1 - X ^ (i + 1))` -/
def qPochC (n : Nat) : Ser := prodL ((List.range n).map fun i => sub one (xpow (i + 1)))
/-- mirror of `qPoch4 n = ∏ i ∈ range n, (1 - X ^ (4 * (i + 1)))` -/
def qPoch4C (n : Nat) : Ser := prodL ((List.range n).map fun i => sub one (xpow (4 * (i + 1))))
/-- mirror of `qPochInf a d` (coefficients `≤ N`: the product over `i ≤ N`) -/
def qPochInfC (a d : Nat) : Ser := prodL ((List.range (N + 1)).map fun i => sub one (xpow (a + d * i)))

/-- tuples of length `k` with entry sum `≤ T` -/
def tuplesTot (T : Nat) : (k : Nat) → List (Fin k → Nat)
  | 0 => [Fin.elim0]
  | k + 1 => (List.range (T + 1)).flatMap fun a =>
      (tuplesTot (T - a) k).map fun t => Fin.cons a t

/-- mirror of `agSum k i` -/
def agSumC (k i : Nat) : Ser :=
  (tuplesTot 6 (k - 1)).foldl (fun acc n =>
    add acc (mul (xpow (agExp k i n)) (prodL ((List.finRange (k - 1)).map fun j => inv (qPochC (n j))))))
    zero

/-- mirror of `agProduct k i` -/
def agProductC (k i : Nat) : Ser :=
  inv (prodL (((List.range (2 * k + 1)).filter fun a => 1 ≤ a ∧ a ≠ i ∧ a ≠ 2 * k + 1 - i).map
    fun a => qPochInfC a (2 * k + 1)))

/-- `C(E8)⁻¹` (Bourbaki labelling), the matrix of `invCartanE8_eq_explicit` -/
def A8 : List (List Int) :=
  [[4, 5, 7, 10, 8, 6, 4, 2], [5, 8, 10, 15, 12, 9, 6, 3], [7, 10, 14, 20, 16, 12, 8, 4],
   [10, 15, 20, 30, 24, 18, 12, 6], [8, 12, 16, 24, 20, 15, 10, 5], [6, 9, 12, 18, 15, 12, 8, 4],
   [4, 6, 8, 12, 10, 8, 6, 3], [2, 3, 4, 6, 5, 4, 3, 2]]
/-- `2 C(E7)⁻¹` (Bourbaki labelling), twice the matrix of `invCartanE7_eq_explicit` -/
def A7x2 : List (List Int) :=
  [[4, 4, 6, 8, 6, 4, 2], [4, 7, 8, 12, 9, 6, 3], [6, 8, 12, 16, 12, 8, 4], [8, 12, 16, 24, 18, 12, 6],
   [6, 9, 12, 18, 15, 10, 5], [4, 6, 8, 12, 10, 8, 4], [2, 3, 4, 6, 5, 4, 3]]

def quad (A : List (List Int)) (n : List Nat) : Int :=
  ((A.zip n).map fun (row, x) => ((row.zip n).map fun (a, y) => a * x * y).sum).sum

/-- mirror of `e8Exp B n = ½ nᵀ C⁻¹ n + B·n` -/
def e8ExpC (B n : List Nat) : Int := quad A8 n / 2 + ((B.zip n).map fun (b, x) => (b * x : Int)).sum
/-- mirror of `e7Exp n = 4 · ½ nᵀ C⁻¹ n = nᵀ (2C⁻¹) n` -/
def e7ExpC (n : List Nat) : Int := quad A7x2 n

/-- the tuples of length `r` whose exponent (monotone in each entry) is `≤ N` -/
partial def tuplesExp (f : List Nat → Int) (r : Nat) : List (List Nat) :=
  let rec go (pre : List Nat) : List (List Nat) :=
    if pre.length = r then [pre] else
      let rec loop (x : Nat) (acc : List (List Nat)) : List (List Nat) :=
        let t := pre ++ [x] ++ List.replicate (r - pre.length - 1) 0
        if f t > N then acc else loop (x + 1) (acc ++ go (pre ++ [x]))
      loop 0 []
  go []

/-- mirror of `e8Sum B` -/
def e8SumC (B : List Nat) : Ser :=
  (tuplesExp (e8ExpC B) 8).foldl (fun acc n =>
    add acc (mul (xpow (e8ExpC B n).toNat) (prodL (n.map fun x => inv (qPochC x))))) zero

/-- mirror of `e7Sum` -/
def e7SumC : Ser :=
  (tuplesExp e7ExpC 7).foldl (fun acc n =>
    add acc (mul (xpow (e7ExpC n).toNat) (prodL (n.map fun x => inv (qPoch4C x))))) zero

/-- mirror of `e7Theta` (`Σ_{k ∈ Icc (-N) N}`) -/
def e7ThetaC : Ser :=
  ((List.range (2 * N + 1)).map fun t => (t : Int) - N).foldl (fun acc k =>
    add acc (add (xpow (3 * k.natAbs ^ 2)) (xpow (7 * k.natAbs ^ 2)))) zero

/-- mirror of `e7ThetaOneSided` -/
def e7ThetaOneSidedC : Ser :=
  (List.range (N + 1)).foldl (fun acc k =>
    add acc (add (xpow (3 * k ^ 2)) (if k = 0 then zero else xpow (7 * k ^ 2)))) zero

-- independent values (scripts/sanity_independent.py, logs/sanity_independent.txt)
def exp_ag_2_1 : List Int := [1, 0, 1, 1, 1, 1, 2, 2, 3, 3, 4, 4, 6, 6, 8, 9, 11, 12, 15, 16, 20, 22, 26, 29, 35, 38, 45, 50, 58, 64, 75, 82, 95, 105, 120, 133, 152, 167, 190, 210, 237]
def exp_ag_2_2 : List Int := [1, 1, 1, 1, 2, 2, 3, 3, 4, 5, 6, 7, 9, 10, 12, 14, 17, 19, 23, 26, 31, 35, 41, 46, 54, 61, 70, 79, 91, 102, 117, 131, 149, 167, 189, 211, 239, 266, 299, 333, 374]
def exp_ag_3_1 : List Int := [1, 0, 1, 1, 2, 2, 3, 3, 5, 6, 8, 9, 13, 14, 19, 22, 28, 32, 41, 47, 59, 68, 83, 96, 117, 134, 161, 186, 221, 254, 301, 344, 405, 464, 541, 619, 720, 820, 949, 1081, 1245]
def exp_ag_3_2 : List Int := [1, 1, 1, 2, 3, 3, 5, 6, 8, 10, 13, 16, 21, 25, 31, 38, 47, 56, 69, 82, 99, 118, 141, 166, 199, 233, 275, 322, 379, 440, 516, 598, 696, 805, 933, 1074, 1242, 1425, 1639, 1878, 2154]
def exp_ag_3_3 : List Int := [1, 1, 2, 2, 3, 4, 6, 7, 10, 12, 16, 19, 25, 30, 38, 46, 57, 68, 84, 99, 121, 143, 172, 202, 242, 283, 336, 392, 462, 537, 630, 729, 851, 982, 1140, 1312, 1518, 1741, 2006, 2295, 2635]
def exp_ag_4_2 : List Int := [1, 1, 1, 2, 3, 4, 6, 7, 10, 13, 17, 21, 28, 35, 44, 55, 69, 84, 105, 127, 156, 189, 229, 275, 333, 397, 475, 565, 673, 795, 943, 1109, 1307, 1533, 1798, 2099, 2455, 2855, 3323, 3855, 4472]
def exp_ag_5_1 : List Int := [1, 0, 1, 1, 2, 2, 4, 4, 7, 8, 11, 13, 19, 22, 30, 36, 47, 56, 73, 86, 110, 131, 163, 194, 241, 284, 348, 412, 499, 588, 709, 832, 996, 1168, 1387, 1622, 1919, 2235, 2631, 3060, 3584]
def exp_ag_5_2 : List Int := [1, 1, 1, 2, 3, 4, 6, 8, 11, 14, 19, 24, 32, 40, 51, 65, 82, 101, 127, 156, 193, 236, 289, 350, 427, 514, 620, 744, 893, 1064, 1271, 1508, 1790, 2116, 2500, 2942, 3464, 4060, 4758, 5560, 6494]
def exp_ag_5_3 : List Int := [1, 1, 2, 2, 4, 5, 8, 10, 14, 18, 25, 31, 42, 52, 68, 84, 108, 133, 168, 205, 256, 311, 384, 463, 567, 681, 826, 988, 1190, 1416, 1696, 2009, 2392, 2823, 3344, 3931, 4636, 5431, 6376, 7445, 8708]
def exp_ag_5_4 : List Int := [1, 1, 2, 3, 4, 6, 9, 11, 16, 21, 28, 36, 48, 60, 78, 98, 124, 154, 194, 238, 296, 362, 444, 539, 658, 793, 960, 1152, 1384, 1652, 1976, 2345, 2789, 3299, 3902, 4596, 5416, 6352, 7454, 8715, 10186]
def exp_ag_5_5 : List Int := [1, 1, 2, 3, 5, 6, 9, 12, 17, 22, 30, 38, 51, 64, 83, 104, 133, 164, 207, 254, 316, 386, 475, 576, 704, 848, 1027, 1232, 1483, 1768, 2116, 2512, 2989, 3534, 4184, 4926, 5808, 6812, 7996, 9348, 10932]
def exp_e8_1 : List Int := [1, 0, 1, 1, 2, 2, 4, 4, 7, 8, 11, 13, 19, 22, 30, 36, 47, 56, 73, 86, 110, 131, 163, 194, 241, 284, 348, 412, 499, 588, 709, 832, 996, 1168, 1387, 1622, 1919, 2235, 2631, 3060, 3584]
def exp_e8_2 : List Int := [1, 1, 1, 2, 3, 4, 6, 8, 11, 14, 19, 24, 32, 40, 51, 65, 82, 101, 127, 156, 193, 236, 289, 350, 427, 514, 620, 744, 893, 1064, 1271, 1508, 1790, 2116, 2500, 2942, 3464, 4060, 4758, 5560, 6494]
def exp_e8_3 : List Int := [1, 1, 2, 2, 4, 5, 8, 10, 14, 18, 25, 31, 42, 52, 68, 84, 108, 133, 168, 205, 256, 311, 384, 463, 567, 681, 826, 988, 1190, 1416, 1696, 2009, 2392, 2823, 3344, 3931, 4636, 5431, 6376, 7445, 8708]
def exp_e8_4 : List Int := [1, 1, 2, 3, 4, 6, 9, 11, 16, 21, 28, 36, 48, 60, 78, 98, 124, 154, 194, 238, 296, 362, 444, 539, 658, 793, 960, 1152, 1384, 1652, 1976, 2345, 2789, 3299, 3902, 4596, 5416, 6352, 7454, 8715, 10186]
def exp_e8_5 : List Int := [1, 1, 2, 3, 5, 6, 9, 12, 17, 22, 30, 38, 51, 64, 83, 104, 133, 164, 207, 254, 316, 386, 475, 576, 704, 848, 1027, 1232, 1483, 1768, 2116, 2512, 2989, 3534, 4184, 4926, 5808, 6812, 7996, 9348, 10932]
def exp_e7 : List Int := [1, 0, 0, 1, 1, 0, 0, 2, 2, 0, 0, 3, 4, 0, 0, 5, 6, 0, 0, 8, 9, 0, 0, 12, 14, 0, 0, 19, 21, 0, 0, 27, 30, 0, 0, 39, 43, 0, 0, 55, 60]
def exp_e7_theta_bilateral : List Int := [2, 0, 0, 2, 0, 0, 0, 2, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

-- Andrews–Gordon: both sides
#guard (agSumC 2 1).toList = exp_ag_2_1 ∧ (agProductC 2 1).toList = exp_ag_2_1
#guard (agSumC 2 2).toList = exp_ag_2_2 ∧ (agProductC 2 2).toList = exp_ag_2_2
#guard (agSumC 3 1).toList = exp_ag_3_1 ∧ (agProductC 3 1).toList = exp_ag_3_1
#guard (agSumC 3 2).toList = exp_ag_3_2 ∧ (agProductC 3 2).toList = exp_ag_3_2
#guard (agSumC 3 3).toList = exp_ag_3_3 ∧ (agProductC 3 3).toList = exp_ag_3_3
#guard (agSumC 4 2).toList = exp_ag_4_2 ∧ (agProductC 4 2).toList = exp_ag_4_2
#guard (agSumC 5 1).toList = exp_ag_5_1 ∧ (agProductC 5 1).toList = exp_ag_5_1
#guard (agSumC 5 2).toList = exp_ag_5_2 ∧ (agProductC 5 2).toList = exp_ag_5_2
#guard (agSumC 5 3).toList = exp_ag_5_3 ∧ (agProductC 5 3).toList = exp_ag_5_3
#guard (agSumC 5 4).toList = exp_ag_5_4 ∧ (agProductC 5 4).toList = exp_ag_5_4
#guard (agSumC 5 5).toList = exp_ag_5_5 ∧ (agProductC 5 5).toList = exp_ag_5_5
-- E8: vacuum (product of `e8_vacuum`) and the four companions (Bourbaki vectors of the statements)
#guard (e8SumC [0, 0, 0, 0, 0, 0, 0, 0]).toList = exp_e8_5
#guard (inv (prodL ([1, 2, 3, 4, 7, 8, 9, 10].map fun a => qPochInfC a 11))).toList = exp_e8_5
#guard (e8SumC [0, 1, 1, 1, 1, 0, 0, 0]).toList = exp_e8_4 ∧ (agProductC 5 4).toList = exp_e8_4
#guard (e8SumC [0, 2, 2, 2, 2, 1, 1, 0]).toList = exp_e8_3 ∧ (agProductC 5 3).toList = exp_e8_3
#guard (e8SumC [1, 3, 3, 4, 4, 2, 2, 0]).toList = exp_e8_2 ∧ (agProductC 5 2).toList = exp_e8_2
#guard (e8SumC [2, 4, 5, 7, 6, 4, 3, 1]).toList = exp_e8_1 ∧ (agProductC 5 1).toList = exp_e8_1
-- E7: left side, one-sided and bilateral right sides
#guard e7SumC.toList = exp_e7
#guard (mul (inv (qPochInfC 4 4)) e7ThetaOneSidedC).toList = exp_e7
#guard e7ThetaC.toList = exp_e7_theta_bilateral
#guard (smul 2 e7SumC).toList = (mul (inv (qPochInfC 4 4)) e7ThetaC).toList

-- negative controls (each must FAIL to be an identity)
#guard (agSumC 5 5).toList ≠ (agProductC 5 4).toList
#guard (e8SumC [0, 0, 0, 0, 0, 0, 0, 0]).toList ≠ (agProductC 5 4).toList
#guard (e8SumC [0, 0, 0, 1, 1, 1, 0, 1]).toList ≠ exp_e8_4   -- paper-labelled B^{(4)} used with Bourbaki matrix
#guard e7SumC.toList ≠ (mul (inv (qPochInfC 4 4)) e7ThetaC).toList   -- the factor 2 matters
#guard (agSumC 3 2).toList ≠ (agProductC 3 1).toList

#eval s!"sanity: all #guard checks passed; e8 vacuum coefficients {(e8SumC [0,0,0,0,0,0,0,0]).toList.take 12}, e7 {e7SumC.toList.take 12}"

end NahmSanity
