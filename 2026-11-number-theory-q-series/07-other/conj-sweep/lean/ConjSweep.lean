/-
  conj-sweep: machine-checked finite facts used in R01, R05, R06, R07.
  Lean 4 core only (no Mathlib).  `decide` proofs are checked by the kernel; `native_decide` proofs
  additionally trust the Lean compiler (axiom `Lean.ofReduceBool`).
  Check with:  lean ConjSweep.lean      (toolchain leanprover/lean4:v4.36.0)
-/
namespace ConjSweep

/-! ## R01: `(1+q)^6 (1+q^3)` is unimodal although the condition of the conjecture fails -/

def polyMul (p q : List Nat) : List Nat :=
  (List.range (p.length + q.length - 1)).map fun k =>
    (List.range (k+1)).foldl (fun acc i => acc + p.getD i 0 * q.getD (k - i) 0) 0

/-- `[a]_q = 1 + q + ... + q^(a-1)` -/
def qint (a : Nat) : List Nat := List.replicate a 1
/-- `[b]_{q^r}` -/
def qintPow (b r : Nat) : List Nat := (List.range (r*(b-1)+1)).map fun i => if i % r == 0 then 1 else 0

def unimodal (c : List Nat) : Bool :=
  (List.range c.length).any fun k =>
    ((List.range k).all fun i => c.getD i 0 ≤ c.getD (i+1) 0) &&
    ((List.range (c.length - 1 - k)).all fun j => c.getD (k+j) 0 ≥ c.getD (k+j+1) 0)

/-- `[2]_q^6 [2]_{q^3}` -/
def P01 : List Nat := polyMul ((List.replicate 6 (qint 2)).foldl polyMul [1]) (qintPow 2 3)

theorem r01_coeffs : P01 = [1, 6, 15, 21, 21, 21, 21, 15, 6, 1] := by native_decide
theorem r01_unimodal : unimodal P01 = true := by native_decide
/-- the hypothesis of the conjecture fails for r = 3, k = 6, a_i = 2, b = 2 -/
theorem r01_condition_fails : ¬ (3 ∣ 2) ∧ ¬ (2 ≤ 1 + 6 * (2 / 3)) := by decide

/-! ## R05: Hankel determinants of powers of the Motzkin generating function -/

def motzkin : Nat → Nat
  | 0 => 1
  | 1 => 1
  | n+2 => ((2*(n+2)+1) * motzkin (n+1) + (3*(n+2)-3) * motzkin n) / (n+4)

def seriesMul (D : Nat) (p q : Array Int) : Array Int :=
  (Array.range D).map fun k =>
    (List.range (k+1)).foldl (fun acc i => acc + p.getD i 0 * q.getD (k - i) 0) 0

/-- coefficients of `M(x)^r` up to `x^(D-1)` -/
def motzkinPow (r D : Nat) : Array Int :=
  let M : Array Int := (Array.range D).map fun n => (motzkin n : Int)
  (List.range r).foldl (fun acc _ => seriesMul D acc M) ((Array.range D).map fun n => if n == 0 then 1 else 0)

/-- fraction-free (Bareiss) determinant of an integer matrix -/
def detBareiss (M : Array (Array Int)) : Int := Id.run do
  let n := M.size
  if n == 0 then return 1
  let mut A := M
  let mut sign : Int := 1
  let mut prev : Int := 1
  for k in [0:n-1] do
    if A[k]![k]! == 0 then
      let mut found := false
      for i in [k+1:n] do
        if !found && A[i]![k]! != 0 then
          let tmp := A[k]!
          A := A.set! k A[i]!
          A := A.set! i tmp
          sign := -sign
          found := true
      if !found then return 0
    for i in [k+1:n] do
      for j in [k+1:n] do
        let v := (A[i]![j]! * A[k]![k]! - A[i]![k]! * A[k]![j]!) / prev
        A := A.set! i (A[i]!.set! j v)
    prev := A[k]![k]!
  return sign * A[n-1]![n-1]!

/-- `H_n(M(x)^r) = det (a_{i+j})_{0 ≤ i,j < n}` -/
def hankel (r n : Nat) : Int :=
  let a := motzkinPow r (2*n)
  detBareiss ((Array.range n).map fun i => (Array.range n).map fun j => a.getD (i+j) 0)

theorem r05_sanity_r1 : (List.range 8).map (hankel 1) = [1, 1, 1, 1, 1, 1, 1, 1] := by native_decide
theorem r05_sanity_r3 : (List.range 8).map (hankel 3) = [1, 1, 0, -4, -4, 0, 9, 9] := by native_decide
/-- the conjecture predicts `(-1)^1 * 2^8 = -256` for both -/
theorem r05_r9 : hankel 9 9 = 256 ∧ hankel 9 10 = 256 := by native_decide
/-- the conjecture predicts `H_7 + H_10 = ± 64` -/
theorem r05_r4 : hankel 4 7 = 64 ∧ hankel 4 10 = -64 := by native_decide

/-! ## R06: greedy 3-sumfree sequences `S_{1,d-1,2d-1}` for d = 6 and d = 4 -/

def isSum3 (S : List Nat) (x : Nat) : Bool :=
  S.any fun a => S.any fun b => S.any fun c => a < b && b < c && a + b + c == x

def sumfree3 (f g h bound : Nat) : List Nat := Id.run do
  let mut S := [f, g, h]
  for x in [h+1:bound+1] do
    if !isSum3 S x then S := S ++ [x]
  return S

/-- d = 6: the period is `55..60` repeated with modulus 45 = 8d - 3 (the conjecture says 8d - 1 = 47) -/
theorem r06_d6 : sumfree3 1 5 11 200 =
    [1, 5, 11, 12, 13, 14, 15, 16, 23, 55, 56, 57, 58, 59, 60, 100, 101, 102, 103, 104, 105,
     145, 146, 147, 148, 149, 150, 190, 191, 192, 193, 194, 195] := by native_decide
/-- d = 4: period of length 3 (not d = 4), modulus 27 -/
theorem r06_d4 : sumfree3 1 3 7 200 =
    [1, 3, 7, 8, 9, 10, 15, 29, 35, 36, 56, 62, 63, 83, 89, 90, 110, 116, 117, 137, 143, 144,
     164, 170, 171, 191, 197, 198] := by native_decide

/-! ## R07: the closure check for the CRIM automaton (Theorem 1 of R07) -/

def dE : Array Nat := #[1, 2, 4, 6, 2, 8, 9, 8, 5, 11, 1, 9, 15, 15, 16, 12, 14]
def dN : Array Nat := #[0, 3, 5, 7, 8, 2, 10, 3, 5, 12, 13, 14, 9, 10, 12, 12, 9]
def acc : Array Bool :=
  #[true, true, true, false, true, false, false, false, false, false, true, false, false, false, true, false, true]

/-- image of a set of states (bitmask) under a letter -/
def img (d : Array Nat) (T : Nat) : Nat :=
  (List.range 17).foldl (fun m s => if T.testBit s then m ||| (1 <<< d[s]!) else m) 0
def anyAcc (T : Nat) : Bool := (List.range 17).any fun s => T.testBit s && acc[s]!

/-- configuration of a word `w` starting with `E`: state of `w`, set of states of `w` minus one letter,
    parity of the length -/
structure Cfg where
  q : Nat
  T : Nat
  odd : Bool
deriving DecidableEq, BEq

def stepE (c : Cfg) : Cfg := ⟨dE[c.q]!, img dE c.T ||| (1 <<< c.q), !c.odd⟩
def stepN (c : Cfg) : Cfg := ⟨dN[c.q]!, img dN c.T ||| (1 <<< c.q), !c.odd⟩
/-- for a word ending in `N`: accepted iff no one-letter deletion is accepted; accepted words have even length -/
def okN (c : Cfg) : Bool := (acc[c.q]! == !anyAcc c.T) && (!acc[c.q]! || !c.odd)

def explore : Nat → List Cfg → List Cfg → Option (List Cfg)
  | 0, _, _ => none
  | _, [], vis => some vis
  | fuel+1, c :: rest, vis =>
    let cE := stepE c
    let cN := stepN c
    if !okN cN then none else
    let (rest, vis) := if vis.contains cE then (rest, vis) else (cE :: rest, cE :: vis)
    let (rest, vis) := if vis.contains cN then (rest, vis) else (cN :: rest, cN :: vis)
    explore fuel rest vis

/-- the word `E`: state `δ(0,E)`, the only deletion gives the empty word (state 0), odd length -/
def init : Cfg := ⟨dE[0]!, 1, true⟩

/-- (i), (iii), (ii) of R07 -/
theorem crim_trim : dN[0]! = 0 ∧ acc[0]! = true ∧
    ((List.range 17).all fun q => acc[dE[q]!]! == acc[q]!) = true := by native_decide
/-- (iv) and (v) of R07: all 169 reachable configurations pass the closure and parity checks -/
theorem crim_closure : (explore 2000 [init] [init]).map List.length = some 169 := by native_decide

end ConjSweep
