/-
scripts/Sanity.lean — sanity test of the DEFINITIONS (not of the proof):
    lake build && lake env lean scripts/Sanity.lean

A computable mirror of the definitions of `KRTW/Defs.lean` (= those of `KRTWStatement.lean`), with
integer coefficients truncated modulo `q^31`.  It uses the very functions `KRTW.nxt`, `KRTW.prv` of
the statement file and copies the formulas of `krtwTerm`, `krtwSum`, `melzerTerm`, `melzerSum`,
`qPochInf` symbol by symbol; the `#guard`s check that for `k = 1, 2, 3` the first 31 coefficients of

  * `krtwSum 1 k`                                    (the 2k-fold KRTW sum),
  * `melzerSum 1 k`                                  (the k-fold sum),
  * `qPochInf 1 (k+2) * qPochInf (k+1) (k+2) * qPochInf (k+2) (k+2) * (qPochInf 1 2 * qPochInf 1 1)⁻¹`

are the integers computed independently by `scripts/sanity_independent.py` (plain Python, depth-first
enumeration of the printed sums).  The only liberty taken: the box `{0,…,N}^k` of the definitions is
replaced by `{0,…,5}^k` (KRTW) resp. `{0,…,7}^k` (Melzer), because a tuple with an entry `> 5`
(resp. `> 7`) has a level of weight `q^{(r_i+s_i)²+s_i²+s_i}` with exponent `≥ 36 > 30`
(resp. `q^{M_i(M_i+1)/2}`, exponent `≥ 36`).
-/
import KRTW.Defs

open KRTW

namespace KRTWSanity

/-- Truncation order: series are stored as their coefficients of `q^0, …, q^N`. -/
def N : Nat := 30

abbrev Ser := Array Int

def zero : Ser := Array.replicate (N + 1) 0

def one : Ser := zero.set! 0 1

/-- `X ^ e` -/
def xpow (e : Nat) : Ser := if e ≤ N then zero.set! e 1 else zero

def add (a b : Ser) : Ser := Array.zipWith (· + ·) a b

def sub (a b : Ser) : Ser := Array.zipWith (· - ·) a b

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

/-- mirror of `qPoch n = ∏ i ∈ range n, (1 - X ^ (i + 1))` -/
def qPochC (n : Nat) : Ser := (List.range n).foldl (fun acc i => mul acc (sub one (xpow (i + 1)))) one

/-- mirror of `q2Poch n = ∏ i ∈ range n, (1 - X ^ (2 * (i + 1)))` -/
def q2PochC (n : Nat) : Ser :=
  (List.range n).foldl (fun acc i => mul acc (sub one (xpow (2 * (i + 1))))) one

/-- mirror of `krtwTerm 1 r s` -/
def krtwTermC {k : Nat} (r s : Fin k → Nat) : Ser :=
  (List.finRange k).foldl (fun acc i =>
    mul acc (mul (mul (xpow ((r i + s i) ^ 2 + s i ^ 2 + s i)) (inv (qPochC (r i - nxt r i))))
      (inv (q2PochC (s i - prv s i))))) one

/-- all tuples of length `k` with entries in `{0, …, B}` -/
def tuples (B : Nat) : (k : Nat) → List (Fin k → Nat)
  | 0 => [Fin.elim0]
  | k + 1 => (List.range (B + 1)).flatMap fun a => (tuples B k).map fun t => Fin.cons a t

/-- mirror of `krtwSum 1 k` (coefficients of `q^0, …, q^30`) -/
def krtwSumC (k : Nat) : Ser :=
  (tuples 5 k).foldl (fun acc r => (tuples 5 k).foldl (fun acc s =>
    if (∀ i, nxt r i ≤ r i) ∧ (∀ i, prv s i ≤ s i) then add acc (krtwTermC r s) else acc) acc) zero

/-- mirror of `melzerTerm 1 M` -/
def melzerTermC {k : Nat} (M : Fin k → Nat) : Ser :=
  (List.finRange k).foldl (fun acc i =>
    mul acc (mul (xpow (M i * (M i + 1) / 2)) (inv (qPochC (M i - nxt M i))))) one

/-- mirror of `melzerSum 1 k` -/
def melzerSumC (k : Nat) : Ser :=
  (tuples 7 k).foldl (fun acc M =>
    if ∀ i, nxt M i ≤ M i then add acc (melzerTermC M) else acc) zero

/-- mirror of `qPochInf a d` (for `a, d ≥ 1`, the factors with `i > N` are `≡ 1 mod q^{N+1}`) -/
def qPochInfC (a d : Nat) : Ser :=
  (List.range (N + 1)).foldl (fun acc i => mul acc (sub one (xpow (a + d * i)))) one

/-- mirror of the product side of identity (i) -/
def prodC (k : Nat) : Ser :=
  mul (mul (mul (qPochInfC 1 (k + 2)) (qPochInfC (k + 1) (k + 2))) (qPochInfC (k + 2) (k + 2)))
    (inv (mul (qPochInfC 1 2) (qPochInfC 1 1)))

/-- mirror of `thetaSum k` -/
def thetaC (k : Nat) : Ser :=
  (List.range (N + 1)).foldl (fun acc n =>
    let t := mul (xpow (((k + 2) * n ^ 2 + k * n) / 2)) (sub one (xpow (2 * n + 1)))
    if n % 2 == 0 then add acc t else sub acc t) zero

/-- Output of `scripts/sanity_independent.py` (independent Python computation). -/
def expected : Nat → Ser
  | 1 => #[1, 1, 1, 2, 2, 3, 4, 5, 6, 8, 10, 12, 15, 18, 22, 27, 32, 38, 46, 54, 64, 76, 89, 104,
      122, 142, 165, 192, 222, 256, 296]
  | 2 => #[1, 1, 2, 3, 4, 6, 9, 12, 16, 22, 29, 38, 50, 64, 82, 105, 132, 166, 208, 258, 320, 395,
      484, 592, 722, 876, 1060, 1280, 1539, 1846, 2210]
  | 3 => #[1, 1, 2, 4, 5, 8, 12, 17, 24, 34, 46, 62, 84, 111, 146, 192, 248, 320, 411, 522, 662, 836,
      1048, 1310, 1632, 2022, 2498, 3077, 3774, 4616, 5632]
  | _ => #[]

#guard krtwSumC 1 == expected 1
#guard krtwSumC 2 == expected 2
#guard krtwSumC 3 == expected 3
#guard melzerSumC 1 == expected 1
#guard melzerSumC 2 == expected 2
#guard melzerSumC 3 == expected 3
#guard prodC 1 == expected 1
#guard prodC 2 == expected 2
#guard prodC 3 == expected 3
-- theta form: (sum) * (q;q²)_∞ (q;q)_∞ = thetaSum = triple product
#guard mul (krtwSumC 1) (mul (qPochInfC 1 2) (qPochInfC 1 1)) == thetaC 1
#guard mul (krtwSumC 2) (mul (qPochInfC 1 2) (qPochInfC 1 1)) == thetaC 2
#guard mul (krtwSumC 3) (mul (qPochInfC 1 2) (qPochInfC 1 1)) == thetaC 3
#guard thetaC 3 == mul (mul (qPochInfC 1 5) (qPochInfC 4 5)) (qPochInfC 5 5)
-- negative controls: a wrong product (modulus k+3 instead of k+2) and a wrong sum are detected
#guard prodC 2 != expected 3
#guard krtwSumC 2 != melzerSumC 3

#eval IO.println s!"sanity OK: krtwSum, melzerSum and the product agree with the independent values for k = 1, 2, 3 mod q^31; k = 3: {krtwSumC 3}"

end KRTWSanity
