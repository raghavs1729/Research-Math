/-
# The kernel-evaluated programs of the finite check: a bit-mask sieve and the integer scan

Everything here is plain `Nat`/`Bool` code that the kernel evaluates on literals
(`decide +kernel`).  The hot paths are written with `Nat.rec` and the raw GMP-accelerated
kernel operations (`Nat.add`, `Nat.pow`, `Nat.shiftRight`, …; no type-class notation), and every
intermediate value that is used more than once is bound by a `fun` (a β-redex), so that the kernel
evaluates it once and every loop iteration forces its accumulators.  (Measured: about
30–40 µs per iteration, against 300–700 µs for the same code written with the equation compiler,
`Nat.testBit` and notation.)

The correctness proofs are in `AllxLib.Finite.Sieve` and `AllxLib.Finite.Scan`.
-/

namespace Allx

/-- `geo q J = 1 + 2^q + 2^(2q) + ⋯ + 2^((J-1)q)`: the bits `0, q, …, (J-1)q`. -/
def geo (q J : Nat) : Nat :=
  Nat.div (Nat.sub (Nat.pow 2 (Nat.mul q J)) 1) (Nat.sub (Nat.pow 2 q) 1)

/-- The first proper multiple of `q` that is `≥ n0` (for `q ≥ 1`). -/
def firstMul (n0 q : Nat) : Nat :=
  (fun m => Nat.mul q (bif Nat.ble 2 m then m else 2)) (Nat.div (Nat.sub (Nat.add n0 q) 1) q)

/-- Bit `i` is set iff `n0 + i` is one of the first `K / q + 1` proper multiples of `q`
that are `≥ n0`. -/
def mask (n0 K q : Nat) : Nat :=
  Nat.shiftLeft (geo q (Nat.add (Nat.div K q) 1)) (Nat.sub (firstMul n0 q) n0)

/-- `acc ||| mask q ||| mask (q+1) ||| ⋯` (`fuel` terms). The `bif` forces the accumulator. -/
noncomputable def sieveLoop (n0 K fuel : Nat) : Nat → Nat → Nat :=
  Nat.rec (motive := fun _ => Nat → Nat → Nat) (fun _ acc => acc)
    (fun _ ih q acc =>
      (fun acc' => bif Nat.beq acc' 0 then ih (Nat.add q 1) 0 else ih (Nat.add q 1) acc')
        (Nat.lor acc (mask n0 K q))) fuel

/-- Bit `i < K` of `sieve n0 K B` is set iff `n0 + i` has a proper divisor `2 ≤ q < B`. -/
noncomputable def sieve (n0 K B : Nat) : Nat := sieveLoop n0 K (Nat.sub B 2) 2 0

/-- One step of the binary lifting used to find `a` with `2^a ≤ P < 2^(a+1)`.
(The result is always a `Nat.add` application, never the bare argument `a`: the kernel caches
the weak head normal form of a term only at the entry of `whnf`, so a function returning an
unevaluated argument would make the kernel re-evaluate the whole chain of earlier values.) -/
def lift (P a d : Nat) : Nat :=
  Nat.add a (bif Nat.ble (Nat.pow 2 (Nat.add a d)) P then d else 0)

/-- Candidate for `⌊log₂ P⌋`, starting from a lower bound `a` (correct if the true value is
`< a + 512`; the scan checks the result, so nothing about this function needs proof). -/
def liftAll (P a : Nat) : Nat :=
  bif Nat.ble (Nat.pow 2 (Nat.add a 1)) P then
    lift P (lift P (lift P (lift P (lift P (lift P (lift P (lift P (lift P a 256) 128) 64) 32) 16)
      8) 4) 2) 1
  else Nat.add a 0

/-- The result type of the scan: success flag and final state `(a, up, um, vp, vm)`. -/
abbrev Res : Type := Bool × Nat × Nat × Nat × Nat × Nat

/-- The failure value of the scan. -/
def scanFail : Res := (false, 0, 0, 0, 0, 0)

/-- The scan over `n, n+1, …, n+fuel-1`. State: `a` (the last exponent `⌊256 log₂ p⌋`),
`up, um` (positive and negative part of `U = ∑_{p<n} lo(p)`), `vp, vm` (the same for
`V = ∑_{p<n} lo(p) (n-1-p)`). `C` is the composite mask of the chunk starting at `n0`. -/
noncomputable def scan (C n0 fuel : Nat) : Nat → Nat → Nat → Nat → Nat → Nat → Res :=
  Nat.rec (motive := fun _ => Nat → Nat → Nat → Nat → Nat → Nat → Res)
    (fun _ a up um vp vm => (true, a, up, um, vp, vm))
    (fun _ ih n a up um vp vm =>
      (fun vp' vm' =>
        bif (bif Nat.ble 4 n then Nat.ble (Nat.add vm' 1) vp' else true) then
          bif (bif Nat.ble 2 n then Nat.beq (Nat.land (Nat.shiftRight C (Nat.sub n n0)) 1) 0
              else false) then
            (fun P => (fun a' =>
              bif (bif Nat.ble (Nat.pow 2 a') P then
                    Nat.ble (Nat.add P 1) (Nat.pow 2 (Nat.add a' 1)) else false) then
                (fun r =>
                  bif Nat.beq r 3 then ih (Nat.add n 1) a' (Nat.add up a') um vp' vm'
                  else bif Nat.beq r 1 then ih (Nat.add n 1) a' up (Nat.add (Nat.add um a') 1) vp' vm'
                  else ih (Nat.add n 1) a' up um vp' vm') (Nat.mod n 4)
              else scanFail) (liftAll P a)) (Nat.pow n 256)
          else ih (Nat.add n 1) a up um vp' vm'
        else scanFail) (Nat.add vp up) (Nat.add vm um)) fuel

/-- The check of one chunk `[n0, n0 + K)`, with sieving bound `B` (`B² ≥ n0 + K`). -/
noncomputable def runChunk (n0 K B a up um vp vm : Nat) : Res :=
  bif Nat.ble (Nat.add n0 K) (Nat.mul B B) then scan (sieve n0 K B) n0 K n0 a up um vp vm
  else scanFail

end Allx
