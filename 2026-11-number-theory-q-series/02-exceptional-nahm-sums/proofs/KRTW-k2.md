# The first Kanade–Russell–Tsuchioka–Warnaar Andrews–Gordon-type identity: an elementary proof for k = 2 (a KNOWN case), and a reduction conjecture for all k

**Correction (referee, `REFEREE-E7-KRTW.md` K1).** The k = 2 case is NOT new: it already follows from KRTW's own
results — Theorem 4.10 (multisum for D⁽²⁾_{2Λ₀}(z,q), all z), eq. (4.7) (C⁽¹⁾_{2Λ₀}(1,q) = D⁽²⁾_{2Λ₀}(1/q,q), the left side of (KRTW-1))
and the Meurman–Primc theorem (their Conjecture 4.1 for n = 1), which gives the product. What follows is a
different, elementary proof of that known case, plus the x-refined identity (E) and the conjecture of §4. An
earlier version of this file and of RESULTS/STATUS-MAP wrongly called the case open (confusing it with
Conjecture 2.2, which is proved only for k ≤ 1).

Source. S. Kanade, M. C. Russell, S. Tsuchioka, S. O. Warnaar, *Remarks on the conjectures of Capparelli,
Meurman, Primc and Primc*, arXiv:2404.03851 (v2, 20 Apr 2026); Selecta Math. 32 (2026), Paper No. 100. After Conjecture 4.9 the authors display four
"conjectural Andrews–Gordon-type identities"; the first one is (r_{k+1} = s_0 := 0, 1/(q)_m := 0 for m<0)

    Σ_{r_1,…,r_k ≥ 0; s_1,…,s_k ≥ 0}  ∏_{i=1}^k  q^{(r_i+s_i)²+s_i²+s_i} / ( (q;q)_{r_i−r_{i+1}} (q²;q²)_{s_i−s_{i−1}} )
        =  (q, q^{k+1}, q^{k+2}; q^{k+2})_∞ / ( (q;q²)_∞ (q;q)_∞ ).                                     (KRTW-1)

It is the z = 1 case of the multisum side of their Conjecture 2.2 (a Hall–Littlewood sum of base q² equals
this multisum), which they prove only for k = 1 (the sum = product identity itself is known for k = 1, 2; see the correction
above). What is re-proved here:

**Theorem 3.** (KRTW-1) holds for k = 2:

    Σ_{r_1≥r_2≥0, 0≤s_1≤s_2} q^{(r_1+s_1)²+s_1²+s_1+(r_2+s_2)²+s_2²+s_2} / ( (q)_{r_1−r_2}(q)_{r_2}(q²;q²)_{s_1}(q²;q²)_{s_2−s_1} )
        = (−q;q)_∞ (−q²;q²)_∞ .

(For k = 2 the right side of (KRTW-1) is (q,q³,q⁴;q⁴)_∞/((q;q²)_∞(q)_∞) = (q⁴;q⁴)_∞/((q;q²)_∞(q²;q²)_∞) = (−q;q)_∞(−q²;q²)_∞.)

Status: complete proof; the chain of 14 elementary moves was found by computer search and is certified by
`code/certify_KRTW_k2.py` (log: `certificates/certify_KRTW_k2.txt`, every move checked in exact arithmetic).
The general case k ≥ 3 is **not** proved (see §4).

## 1. Generalized Nahm data and three lemmas

Data (M,B,d): M symmetric r×r rational, B ∈ Q^r, d ∈ {1,2,…}^r;
F_{M,B,d}(q) = Σ_{n∈N^r} q^{½nᵀMn+B·n} / ∏_i (q^{d_i};q^{d_i})_{n_i}, defined when (M,B) is *admissible* (finitely many
n with exponent ≤ E, for every E), as an element of Z((q^{1/D})).

**Lemma F (fusion, base Q = q^d).** If d_a = d_b = d (a ≠ b), adjoin a species c with d_c = d, row_c = row_a + row_b,
M_cc = M_aa+M_bb+2M_ab, B_c = B_a+B_b, and replace M_ab by M_ab + d. Then F is unchanged.
(Lemma 2 of `E8-E6.md` with the Durfee identity in base Q: the extra exponent is d·m_a m_b.)

**Lemma S (self-fusion).** For a ≥ 0: 1/(Q;Q)_a = Σ_{m+2c=a} Q^{m(m−1)/2} / ( (Q;Q)_m (Q²;Q²)_c ).
*Proof.* Coefficient of z^a in 1/(z;Q)_∞ = (−z;Q)_∞ · 1/(z²;Q²)_∞, using Euler's two identities
Σ z^m Q^{m(m−1)/2}/(Q;Q)_m = (−z;Q)_∞ and Σ z^{2c}/(Q²;Q²)_c = 1/(z²;Q²)_∞. ∎
Consequence: a species a with (d, M_aa = α, B_a = β) may be replaced by two species m, c with
d_m = d, M_mm = α + d, B_m = β − d/2; d_c = 2d, M_cc = 4α, B_c = 2β; M_mc = 2α; M_{m,k} = M_{a,k}, M_{c,k} = 2M_{a,k} for
the other species k. (Substitute n_a = m + 2c: ½α(m+2c)² + β(m+2c) + d·m(m−1)/2.)

**Lemma P (split).** For n ≥ 0: 1/(Q;Q)_n = Σ_{s+t=n} Q^t / ( (Q²;Q²)_s (Q²;Q²)_t ).
*Proof.* Coefficient of z^n in 1/(z;Q)_∞ = 1/(z;Q²)_∞ · 1/(zQ;Q²)_∞. ∎
Consequence: a species n with (d, α, β) may be replaced by two species s, t with d_s = d_t = 2d, identical rows
(M_ss = M_tt = M_st = α, M_{s,k} = M_{t,k} = M_{n,k}), B_s = β, B_t = β + d.

**Admissibility, in both directions.** In each lemma the new exponent is the old exponent evaluated at the
image point (n = m + c·v, n_a = m + 2c, n = s + t respectively) plus a non-negative term (d·m_a m_b,
d·m(m−1)/2, d·t), and the map to the image point has finite fibres; so admissibility of the old data implies
that of the new data, and the regrouping of the summable family proves the equality (as in Lemma 2 of
`E8-E6.md`). Conversely, each old point n has a preimage on which the extra term vanishes
(c = min(n_a,n_b); m ∈ {0,1}; t = 0), chosen injectively, so admissibility of the new data implies that of
the old data. Hence each lemma may be used in either direction: if two data are related by one of the three
replacements and one of them is admissible, both are and their sums agree. (Some intermediate matrices
below are not positive definite — after a split two rows coincide, and the state after step 5 is even
indefinite — but they are admissible by this argument.)

## 2. Proof of Theorem 3

In the variables a_1 = r_1−r_2, a_2 = r_2 (d = 1), b_1 = s_1, b_2 = s_2−s_1 (d = 2) the left side is F_{M,B,d} with

    M = ( 2 2 2 0 ; 2 4 4 2 ; 2 4 8 4 ; 0 2 4 4 ),   B = (0,0,2,1),   d = (1,1,2,2)

(positive definite). Apply, in this order (indices refer to the current list of species, starting from 0;
new species are appended; removed species are deleted from the list):

    1 self-fuse 0        2 fuse (3,4)          3 self-defuse (m=1,c=5)   4 split 0          5 fuse (0,3)
    6 defuse (0,5 → 4)   7 fuse (3,4)          8 merge (5,6)             9 fuse (0,4)      10 merge (0,4)
    11 self-defuse (0,5) 12 defuse (0,1 → 4)   13 self-defuse (0,3)      14 self-defuse (1,2)

(defuse (a,b → c): species c is recognised as the fusion of a and b and removed, Lemma F backwards;
self-defuse (m,c): Lemma S backwards; merge (s,t): Lemma P backwards). All fifteen data sets are printed in
`certificates/certify_KRTW_k2.txt`; for every backward step the script checks in exact rational arithmetic
that the forward replacement applied to the new data reproduces the old data. The final data are

    M = ( 1 1 ; 1 2 ),  B = (1/2, 1),  d = (1,1),

i.e. the sum Σ_{n_1,n_2≥0} q^{(N_1²+N_2²+N_1+N_2)/2}/((q)_{n_1}(q)_{n_2}) with N_1 = n_1+n_2, N_2 = n_2. Finally

    Σ_{n_1} q^{n_1(n_1+1)/2 + n_1 n_2}/(q)_{n_1} = (−q^{1+n_2};q)_∞ = (−q;q)_∞/(−q;q)_{n_2}      (Euler),
    Σ_{n_2} q^{n_2²+n_2}/((q)_{n_2}(−q;q)_{n_2}) = Σ_{n_2} q^{n_2²+n_2}/(q²;q²)_{n_2} = (−q²;q²)_∞   (Euler, base q²),

which gives (−q;q)_∞(−q²;q²)_∞. ∎

Numerical confirmation: the 4-fold sum equals (−q)_∞(−q²;q²)_∞ mod q^60 (`code/test_melzer.py`), and mod q^30
at every positive definite stage of the chain (`certify_KRTW_k2.txt`).

The core of the chain (steps 1–13 after discarding the spectator b_1) is the following identity with a free
parameter x (it follows from the chain because the chain preserves the grading r_1+r_2+2s_1+2s_2, and b_1 = s_1
enters only through this grading; verified in x⁰,…,x¹⁰ mod q⁴⁰ by the referee):

    Σ_{r_1≥r_2≥0, b≥0} x^{r_1+r_2+2b} q^{r_1²+(r_2+b)²+b²+b} / ((q)_{r_1−r_2}(q)_{r_2}(q²;q²)_b)
        = Σ_{N_1≥N_2≥0} x^{N_1+N_2} q^{N_1(N_1+1)/2+N_2²} / ((q)_{N_1−N_2}(q)_{N_2}).                        (E)

(The chain preserves exactly one linear grading, x^{Σ(r_i+2s_i)}; `code/charges.py`.)

## 3. What this does and does not give

* It gives an elementary proof of the (known) k = 2 case of the first KRTW Andrews–Gordon-type identity and
  the refinement (E). It says nothing new about Conjecture 2.2 (Hall–Littlewood sum = multisum with z): the
  z-dependence z^{Σ(r_i+s_i)} is not preserved by the chain. (The Hall–Littlewood product at z = 1 is KRTW
  Proposition 4.8 with (4.12a).)
* The other three KRTW identities (linear terms Σ(r_i+2s_i), 0, and the Ω-weighted one) are not touched: a
  search of the same move class up to rank 7 found no reduction.

## 4. Conjecture for general k (numerical)

**Update 2026-11-26.** The conjecture (E_k) below, and with it (KRTW-1) for every k, is now proved: see
`proofs/KRTW-general.md` (Theorems 1–3; certificate `code/krtw_general_certify.py`). The text of this section is kept as it was.

**Conjecture (E_k / Melzer reduction).** For all k ≥ 1, as power series in x and q,

    Σ_{r,s} x^{Σ(r_i+2s_i)} ∏_{i=1}^k q^{(r_i+s_i)²+s_i²+s_i} / ((q)_{r_i−r_{i+1}}(q²;q²)_{s_i−s_{i−1}})
        = Σ_{N_1≥…≥N_k≥0} x^{N_1+…+N_k} q^{Σ_i N_i(N_i+1)/2} / ∏_i (q)_{N_i−N_{i+1}} ,

equivalently (peeling off b_1 = s_1 on the left and the self-fusion partner of n_k on the right) the same
statement with s_1 = 0 on the left and the last level weight q^{N_k²} on the right. Checked for k = 2,3,4
(mod q^40, q^36, q^24 in all powers of x; `code/test_Ek.py`). At x = 1 the right side is the Melzer-type sum
Σ q^{Σ N_i(N_i+1)/2}/∏(q)_{n_i}, which equals (−q)_∞(q,q^{k+1},q^{k+2};q^{k+2})_∞/(q)_∞ numerically for k ≤ 4
(`code/test_melzer.py`; for general k this evaluation should be among the Melzer/Stembridge/Bressoud–Ismail–
Stanton/Warnaar identities — I did not verify a reference). So (KRTW-1) for all k would follow from (E_k).
A "symmetric" termwise version (fixed multiset {N_i}) is false (`code/test_symchain.gp`).
