# lean-tau: a Lean 4 + Mathlib proof of Theorem 1 of `open-misc/tau/PROOF.md`

The equations `x^22 − 3 x^11 y^2 + y^4 = 461` and `x^22 − 3 x^11 y^2 + y^4 = −599` have no solution in
integers. This is the Diophantine statement behind "τ(p^4) is never 461 or −599", where τ is
Ramanujan's tau function.

The proof has no `sorry` and no `native_decide`. It depends only on the standard axioms
`propext`, `Classical.choice` and `Quot.sound` (`logs/print_axioms.txt`). All finite computations
are checked by the kernel with `decide +kernel`.

## 1. What is proved (`TauStatement.lean`, which imports only Mathlib)

```lean
namespace Tau

def heckeSeq (t p : ℤ) : ℕ → ℤ
  | 0 => 1
  | 1 => t
  | k + 2 => t * heckeSeq t p (k + 1) - p ^ 11 * heckeSeq t p k

theorem no_solution_461 (x y : ℤ) : x ^ 22 - 3 * x ^ 11 * y ^ 2 + y ^ 4 ≠ 461
theorem no_solution_neg599 (x y : ℤ) : x ^ 22 - 3 * x ^ 11 * y ^ 2 + y ^ 4 ≠ -599
theorem heckeSeq_four_ne (t p : ℤ) : heckeSeq t p 4 ≠ 461 ∧ heckeSeq t p 4 ≠ -599
theorem tau_prime_pow_four_ne (τ : ℕ → ℤ) (h1 : τ 1 = 1) (p : ℕ) (_hp : p.Prime)
    (hrec : ∀ k : ℕ, τ (p ^ (k + 2)) = τ p * τ (p ^ (k + 1)) - (p : ℤ) ^ 11 * τ (p ^ k)) :
    τ (p ^ 4) ≠ 461 ∧ τ (p ^ 4) ≠ -599
```

`no_solution_461` and `no_solution_neg599` together are Theorem 1. `heckeSeq_four_ne` is the same
result with `x = p` and `y = t`, because `heckeSeq t p 4 = t^4 − 3 p^11 t^2 + p^22`. In
`tau_prime_pow_four_ne`, primality of `p` is never used (hence the binder name `_hp`).

## 2. What is NOT machine-checked

- **That Ramanujan's τ satisfies the hypotheses of `tau_prime_pow_four_ne`.** No modular-forms
  library is used: τ is an arbitrary function `ℕ → ℤ`. The facts τ(1) = 1 and
  τ(p^{k+2}) = τ(p) τ(p^{k+1}) − p^11 τ(p^k) for primes p (Hecke) are assumed, not proved. As a sanity
  check, `scripts/Sanity.lean` evaluates `heckeSeq` against known values: τ(4), τ(8), τ(16) and τ(81).
- **The reduction from τ(n) = ±l to n = p^{d−1}** (Balakrishnan–Craig–Ono–Tsai, and
  Balakrishnan–Craig–Ono; PROOF.md §4) is **not formalized**. Neither are the congruences of
  PROOF.md §4 that leave d = 5 as the only exponent for 461 and for −599. Corollary 2 of PROOF.md
  (τ(n) ∉ {±461, ±599}) is therefore **not** proved here. Only the Diophantine part, Theorem 1, is.

## 3. The proof (library `Tau/`, about 420 lines)

| file | content |
|---|---|
| `Tau/Basic.lean` | the form `Qf (U, V) = U² − 3UV + V²`; the map `T (U, V) = (V, 3V − U)` as an `Equiv.Perm (R × R)` over any commutative ring `R`; invariance of `Qf` under `T ^ n` for `n : ℤ`; `T ^ n` commutes with every ring homomorphism (in particular reduction `ℤ → ZMod q`) |
| `Tau/Descent.lean` | Lemma A (seed bound) and the finite seed search |
| `Tau/Sieve.lean` | periodicity modulo `q` and the congruence sieve, with a computable checker `certOK` and its soundness theorem `sieve` |
| `Tau/Main.lean` | the two certificates (`decide +kernel`), Theorem 1, the Hecke-sequence lemmas |
| `Tau/Defs.lean` | `heckeSeq` (inlined verbatim in `TauStatement.lean`) |

**Descent (`Tau/Descent.lean`).** Take a point p of an orbit of `T` on `ℤ × ℤ`. Choose m with
`|(T^m p).1|` minimal along the orbit (`Nat.find`). Put `z = Z_m`, `a = Z_{m−1}`, `b = Z_{m+1}`. Then
`a + b = 3z`, `ab = z² − α`, `|a| ≥ |z|` and `|b| ≥ |z|`. From
`(|a| − |z|)(|b| − |z|) ≥ 0` and `|a| + |b| ≥ |a + b| = 3|z|` we get
`|z² − α| = |a||b| ≥ 2z²` (`key_ineq`). This is a little sharper than Lemma A of PROOF.md:

- for α = 461 it gives `3z² ≤ 461`, so `|z| ≤ 12`;
- for α = −599 it gives `z² ≤ 599`, so `|z| ≤ 24`.

Then `(2b − 3z)² = 5z² + 4α` bounds `|b| ≤ 50`. A kernel-checked search over the box
(`seedsOK`, about 2,500 to 5,000 points) shows that the only solutions in the box are the four seeds
of PROOF.md:

- α = 461: `(±1, …)`, namely (−1, −23), (−1, 20), (1, −20), (1, 23);
- α = −599: (−24, −47), (−24, −25), (24, 25), (24, 47).

So every solution of `Qf (U, V) = α` is `T ^ n s` for some seed `s` and some `n : ℤ` (`orbit_461`,
`orbit_M599`).

**Sieve (`Tau/Sieve.lean`).** The arithmetic is done on natural numbers:

- `stepN q (a, b) = (b, (3b + (q − a)) % q)`, iterated by `iterN`;
- `isPow e q a` is a brute-force test of "a is an e-th power mod q".

This way the kernel evaluates everything with its GMP-backed `Nat` operations. `castP_iterN` links
`iterN` to `T (ZMod q) ^ r`.

`certOK S Qs M` checks the following. For every seed `s` and every `r < M`, some `q ∈ Qs` satisfies:
`q > 1`, the orbit of `s mod q` returns after `M` steps, and at step `r` either the first coordinate
is not an 11th power mod `q` or the second is not a square mod `q`.

`sieve` proves that `certOK S Qs M = true` implies `T ^ n s ≠ (x^11, y^2)` for every `s ∈ S`,
`n : ℤ` and `x, y : ℤ`. Reducing modulo `q` turns `T ^ n` into `T(ZMod q) ^ n`, which equals
`T(ZMod q) ^ (n mod M)` on the seed's orbit (`zpow_apply_emod`), and both directions n < 0 and n ≥ 0
are covered. **The moduli q need not be prime.** The certificates are those of PROOF.md:

- α = 461: M = 88 and q ∈ {3, 7, 43, 89, 199, 967};
- α = −599: M = 22 and q = 89.

The checker does not use the witness table of `certificate_tau_p4.txt`. It simply tries the moduli in
order. It also re-verifies the periods instead of taking them from a table.

**Assembly.** `x^22 − 3x^11y^2 + y^4 = Qf (x^11, y^2)`, followed by `orbit_*` and `sieve`.

## 4. Pins and how to check

Lean `leanprover/lean4:v4.33.0-rc2` and Mathlib `51e6992efd06126df61a496bebf8f49482a4e129`, the pins
of `../lean-krtw`.

    curl -sSfL https://raw.githubusercontent.com/leanprover/elan/master/elan-init.sh | sh -s -- -y --default-toolchain none
    export PATH="$HOME/.elan/bin:$PATH"
    cd monotonicity/lean-tau
    lake exe cache get          # Mathlib binaries; never build Mathlib from source
    lake build                  # "declaration uses sorry" comes from TauStatement.lean only (deliberate)
    lake env lean scripts/PrintAxioms.lean
    ./scripts/compare_statements.sh
    lake env lean scripts/Sanity.lean

`lake build` of the library takes about 30 s, and the two kernel certificate checks take a few
seconds. The outputs are saved in `logs/` (`print_axioms.txt`, `compare_statements.txt`,
`sanity.txt`).

`scripts/compare_statements.sh` is the script of `../lean-krtw`, adapted. It checks four things:

- `TauStatement.lean` imports only Mathlib;
- its definition block is byte for byte that of `Tau/Defs.lean`;
- each statement, docstring included, occurs byte for byte in `TauSolution.lean`;
- the solution side has no `sorry`, `axiom` or `native_decide`.

It also checks that `#check` and `#print` give identical output from both modules.
`TauStatement.lean` is generated by `scripts/gen_statement.py`.

`scripts/Sanity.lean` checks three things:

- `heckeSeq` against known values of τ;
- non-vacuity of the checker: for the planted solution (U, V) = (2^11, 3^2) of
  `Qf = 4139089`, `certOK` returns `false`;
- non-vacuity of the seed search: dropping one seed makes `seedsOK` return `false`.
