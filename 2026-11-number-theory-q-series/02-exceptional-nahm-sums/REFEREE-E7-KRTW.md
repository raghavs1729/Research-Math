# Hostile referee report on `proofs/E7.md` and `proofs/KRTW-k2.md`

Referee: independent third referee, 2026-11-26. Nothing in the authors' scripts or logs was used as evidence
(I read `code/certify_KRTW_k2.py` only to see what it does and does not check).
Versions refereed: `proofs/E7.md` md5 d20689c0c1aa51e6a88b5828672ec74d (the file was modified at 10:44 while I was
working: §6 "Theorem 5" was added; I refereed the later version); `proofs/KRTW-k2.md` md5 1f5fb54f16ac21eb1d7bc338c7c6496a.
Referee code (written from scratch, exact integer/rational arithmetic, python3 + python-flint + sympy, PARI/GP):
`/tmp/claude-0/-home-user-weighted-prime-races/3c4df68b-d770-5604-8028-1966b4f5b18d/scratchpad/referee3/code/`
Downloads (untrusted data, only parsed as text):
`/tmp/claude-0/-home-user-weighted-prime-races/3c4df68b-d770-5604-8028-1966b4f5b18d/scratchpad/referee3/dl/`

## 0. Verdicts

| | `proofs/E7.md` (Theorem 4, the C(E7)⁻¹ Nahm sum) | `proofs/KRTW-k2.md` (Theorem 3, KRTW identity (i), k = 2) |
|---|---|---|
| Correctness | **No mathematical error found.** Every matrix entry, every step of Lemmas 3–4, Propositions 5–6 and the final assembly was re-derived and re-computed; Theorem 4 was verified by brute force for all exponents ≤ 150. | **The proof is correct.** Lemmas F, S, P and the two-directional admissibility argument are right; the 14-move chain was replayed with my own implementation (every backward step verified exactly by the forward move); all 15 states have the same series to q^30, including the 8 non-positive-definite ones the authors' certificate skips. |
| Completeness | **Complete**, given Lemmas 1–2 of `E8-E6.md`, Bailey's lemma, the unit Bailey pair, Cauchy's identity and q-Chu–Vandermonde (all standard and correctly stated). Only citation locators and one unproved-but-standard modularity sentence are deficient. §6 (Theorem 5) is a sketch but a correct one. | Complete for Theorem 3. Identity (E) is asserted "with a free parameter x" without a written proof (it is true, and it does follow from the chain; see K3). |
| Novelty risk | **Low–moderate.** No prior appearance of the formula or of a proof of modularity found (searches and their limits in §C). The question itself is older than Sun–Wang: Cherednik–Feigin 2013 §3.6 (not cited in `E7.md`). | **The novelty claim is wrong.** Identity (i) for k = 2 is *already a theorem by the KRTW paper itself*: their Theorem 4.10 (first formula, proved for all z) at z = 1/q, with their (4.7) and Conjecture 4.1 for n = 1 (which, as KRTW state, is established by Meurman–Primc). See K1. What is new is at most an elementary proof and the x-refinement (E). |

I was told to assume there is an error. In `E7.md` I could not find one. In `KRTW-k2.md` the mathematics is right
and the error is in the claim of what was open.

## 1. Numbered list of defects

### Document 2 (`proofs/KRTW-k2.md`)

**K1 (major; false novelty / "k = 2 was open" is not supported by the source).**
The document says the first Andrews–Gordon-type identity is "the z = 1 case of the multisum side of their
Conjecture 2.2 …, which they prove only for k = 1. What is proved here: Theorem 3. (KRTW-1) holds for k = 2", and
`STATUS-MAP.md` row 2 / `RESULTS.md` list "(i) for k=2 PROVED" as a result of this project while noting only for
(ii)–(iv) that they "follow from KRTW Thm 4.10 + Bartholdi". But (i) for k = 2 follows from the KRTW paper even more
directly than (ii)–(iv), and unconditionally:
- KRTW Theorem 4.10 (`KRTW.tex` label `Thm_kis2`; proved by a computer-assisted functional-equation argument, for
  all z): D⁽²⁾_{2Λ₀}(z,q) = Σ_{r,s} ∏_{i=1}^{2} z^{r_i+s_i} q^{(r_i+s_i)²+s_i²+r_i+2s_i} / ((q)_{r_i−r_{i+1}}(q²;q²)_{s_i−s_{i−1}}).
  KRTW's own sentence before it: "Conjecture 4.9 for k = 2 can be completed to a full set of weights, which is
  provable using the functional equations".
- KRTW (4.7) (label `Eq_CDn2`, a consequence of the proved Proposition 4.4): C⁽¹⁾_{aΛ₀+(k−a)Λ₁}(z) = D⁽²⁾_{aΛ₀+(k−a)Λ₂}(z/q).
  With a = k = 2, z = 1: C⁽¹⁾_{2Λ₀}(1,q) = D⁽²⁾_{2Λ₀}(1/q,q) = Σ ∏ q^{(r_i+s_i)²+s_i²+s_i}/(…) = left side of (KRTW-1).
- KRTW Conjecture 4.1 (CMPP, type C_n⁽¹⁾) for n = 1: KRTW write (p. 20 of the arXiv v2 PDF, after the diagram):
  "we obtain the two-colour partition model of [Meurman–Primc, (11.2.10) & (11.2.11)]. By results of that same
  paper, this establishes Conjecture 4.1 for n = 1." (Independently, the weight kΛ₀ for all n is Primc–Trupčević
  [PT24, Thm 2.1], "Hence for this special case Conjecture 4.1 is now a theorem".) For n = 1, (k₀,k₁) = (2,0) the
  product of Conjecture 4.1 is (q⁴;q⁸)_∞(q⁸;q⁸)_∞θ(q;q⁴)/((q;q²)_∞(q)_∞) = (q,q³,q⁴;q⁴)_∞/((q;q²)_∞(q)_∞), the right side of
  (KRTW-1) for k = 2.
This is exactly how KRTW derive the four identities ("Combining this with (4.7) and Conjectures 4.1 and 4.2 gives
the following conjectural Andrews–Gordon-type identities"): for general k the only conjectural input to (i) is
Conjecture 4.9, and for k = 2 that input is Theorem 4.10. KRTW also describe the right side of (i) as "the
principally specialised character of the A₁⁽¹⁾-standard module L(kΛ₀)" (introduction). So Theorem 3 is a known
case; the claim that only k = 1 was known confuses Conjecture 2.2 (Hall–Littlewood sum = multisum, indeed proved
only for k ≤ 1) with the sum = product identity (i). The project's own file `status/A-kanade-russell.md` §9.4
already says "(i)–(iv) for k=2 are theorems [I: my assembly]; k ≥ 3: OPEN" — `KRTW-k2.md`, `STATUS-MAP.md` and
`RESULTS.md` contradict it. Required change: retitle to "a new elementary proof of the k = 2 case", cite Theorem
4.10 + (4.7) + Meurman–Primc, and remove "(i) for k=2 PROVED" as a new result. (Also: the last bullet-one sentence
of §3, "it gives Conjecture 2.2 for k = 2 at z = 1 only", is likewise already implied by KRTW: Prop. 4.8 gives the
Hall–Littlewood side at z = 1 as C⁽¹⁾_{kΛ₀}(1,q), which by (4.7) and Theorem 4.10 is the multisum.)

**K2 (minor; false side statement).** "Some intermediate matrices below are only positive *semi*definite — after a
split two rows coincide". The state after step 5 (rank 7) is **indefinite**: with the species order of the chain,
v = (2,0,0,−1,0,−2,0) gives vᵀMv = −4 (smallest eigenvalue ≈ −0.84). States 4 and 6–11 are singular positive
semidefinite; only states 0–3 and 12–14 are positive definite. This does not damage the proof (admissibility is
obtained from the lemma, not from definiteness, and all entries of every M are ≥ 0), but the sentence and the
message printed by `certify_KRTW_k2.py` ("matrix only semidefinite after SPLIT") are wrong for state 5.

**K3 (minor; unproved assertion and wrong step range).** Identity (E) is stated as a fact "which holds with a free
parameter x", supported only by a pointer to `code/charges.py`. No proof is written (a formal-variable version of
Lemmas F, S, P, as in the remark after Lemma 2 of `E8-E6.md`, plus the observation that b₁ is a spectator, is
needed). Also "(steps 1–12 after discarding the spectator b₁)" is wrong: state 12 without b₁ has rank 3
(M = (2 1 2; 1 3 2; 2 2 4), B = (0,½,1), d = (1,1,2)); the rank-2 right side of (E) is reached only after step 13
(self-defuse (0,3)), giving M = (1 1; 1 3), B = (½,½), i.e. N₁(N₁+1)/2 + N₂². b₁ (index 2) is untouched by steps
1–13 and is consumed in step 14. I confirm: the space of linear gradings preserved by all 14 moves is
one-dimensional, spanned by (a₁,a₂,b₁,b₂) ↦ a₁+2a₂+4b₁+2b₂ = r₁+r₂+2s₁+2s₂ (`k4_grading.py`), with final weights
(1,2) = N₁+N₂; and (E) holds in x⁰,…,x¹⁰ mod q^40 (`k3_E.py`). So (E) is true and provable from the chain, but the
document does not prove it.

**K4 (minor; certificate weaker than advertised).** `certify_KRTW_k2.py` performs the numerical comparison only at
positive definite states, i.e. at 7 of the 15 states; for states 4–11 the only check is the exact matrix identity.
That is sufficient for the proof, but "certified … every move checked" should not be read as a series check of
every state. (My evaluator, valid because all matrix entries are ≥ 0, confirms all 15 series to q^30.)

**K5 (minor; locators).** "the product formula for the Hall–Littlewood side at z = 1 (KRTW, discussion of
Conjecture 2.2; not re-derived here)": there is no such formula in the discussion of Conjecture 2.2; the relevant
statements are Conjecture 4.7, eq. (4.12a) with n = 1, Proposition 4.8 ("Moreover, (4.12a) holds for z = 1") and
Conjecture 4.1 for kΛ₀. KRTW is now published: Selecta Math. (N.S.) 32 (2026), no. 5, Paper No. 100 (arXiv metadata, v2 of 20 Apr 2026).

**K6 (minor; unverified reference, admitted by the authors).** §4: the Melzer-type evaluation
Σ q^{ΣN_i(N_i+1)/2}/∏(q)_{n_i} = (−q)_∞(q,q^{k+1},q^{k+2};q^{k+2})_∞/(q)_∞ is used with "I did not verify a reference".
Numerically true for k ≤ 4 to q^40 (`k6_melzer.py`); I did not locate the reference either (it is the companion,
with linear term ½ΣN_i, of the identity quoted as (1.6) in Wang–Wang 2607.08606 from Stembridge / Bressoud–Ismail–
Stanton / Warnaar). §4 is a conjecture section and is not part of Theorem 3.

No error was found in: Lemma S and its consequence (M_mm = α+d, B_m = β−d/2, M_cc = 4α, B_c = 2β, M_mc = 2α,
M_ck = 2M_ak, d_c = 2d); Lemma P and its consequence; Lemma F in base q^d; the encoding of the KRTW sum as
(M,B,d) = ((2 2 2 0; 2 4 4 2; 2 4 8 4; 0 2 4 4), (0,0,2,1), (1,1,2,2)); the final data ((1 1; 1 2), (½,1), (1,1)); the two
Euler evaluations; the product manipulations.

### Document 1 (`proofs/E7.md`)

**E1 (minor; wrong citation locators in Lemma 3).** "Andrews, Pacific J. Math. 114 (1984) … (2.12)–(2.13) for the
unit pair, Theorem 1/(3.1) for the lemma". (2.12)–(2.13) is right. But (3.1) is Andrews' *Theorem 2*, the
N → ∞, all b_i,c_i → ∞ limit, which is not what Lemma 3 uses (finite N). The finite statement is "Bailey's Lemma"
(2.3)–(2.6) with ρ₁,ρ₂ → ∞, or Theorem 1 (2.7) with k = 3 and all b_i,c_i → ∞ (which is literally (2.1) for a general
Bailey pair). "Warnaar, 50 years of Bailey's lemma, §§2–3" should be §§1–2: definition (1.3), Theorem 2.1, the
limiting form (2.2), unit pair (2.3). (Both sources downloaded and read.)

**E2 (minor; probable citation slip, not verified against the book).** Lemma 4: "Andrews, The Theory of
Partitions, (2.2.8) with z = q^c, divided by (q)_c". As I remember the book, (2.2.8) is
Σ q^{n²−n}zⁿ/((q)_n(z)_n) = 1/(z)_∞; then the specialisation that gives Σ q^{m(m+c)}/((q)_m(q)_{m+c}) = 1/(q)_∞ is
z = q^{c+1} (divided by (q)_c), not z = q^c. I could not consult the book. The identity itself is correct
(verified for c ≤ 11 mod q^80) and is also the a,b → ∞ limit of Lemma 1, so nothing depends on this.

**E3 (minor; unproved standard assertion).** "(ϑ(3τ/2)+ϑ(7τ/2))/(2η(τ)) … a modular function of weight 0 on a
congruence subgroup" is asserted without argument or reference. It is true (in the variable u = τ/24 it is
(θ(18u)+θ(42u))/(2η(24u)) with θ(u) = Σ e^{2πik²u}, a quotient of weight-½ forms on Γ₀(N) with the theta multiplier
and quadratic characters, hence invariant under some Γ₁(N)); one sentence and a reference should be added, since
this sentence is the whole content of "Conjecture 1.1 of Sun–Wang holds for (T1,E7)".

**E4 (minor; attribution).** "Source of the problem" cites only Sun–Wang (2026). The modularity of exactly this
series (the "♭-case", A♭ = ((ω_i,ω_j)) = C(R)⁻¹ for simply-laced R) was raised in Cherednik–Feigin, Adv. Math. 248
(2013), arXiv:1209.1978, §3.6, Table 1 (E7: L♭ = 6, so effective central charge 7 − 6 = 1), where it is said to be a
conditional corollary of Nahm's conjecture and that "our approach … cannot be directly extended to the ♭-case"
(I downloaded the paper and checked: no formula for E7). This is recorded in `PRIOR-ART.md`/`RESULTS.md` but not
in the proof document, whose §5 still refers to `STATUS-MAP.md`/`REFEREE-E8-E6.md` for the searches.

**E5 (minor; §6 is a sketch).** Theorem 5 is proved by "exactly as for k = 3"; the k-fold version of Lemma 3 is not
written out. It is immediate (k applications of the same step), the exponent identity
(k/4)a² + kar + (k+1)r² + r(r−1)/2 + j² + 2rj = (k/4)n² + r(r−1)/2 + u² is correct, and I verified Theorem 5 for
k = 1,…,5 for exponents ≤ 30, the k = 2 right side against Stembridge's product
(−q^{1/2};q)_∞(q^{5/2},q^{7/2},q⁶;q⁶)_∞/(q)_∞ and against the (T1,T₄) sum. No novelty search was made for Theorem 5 by
the authors or (beyond the searches of §C) by me.

**E6 (cosmetic).** The file was changed after being handed to the referee (§6 says so); §4's list of checks and
the "Status" header do not mention Theorem 5.

No error was found in: all four printed matrices of §1 (entry by entry), the identification with C(E7)⁻¹, (1.1),
Lemma 3 (statement, x = 1 convention, hypothesis x ≠ q^{−m}), Lemma 4, Proposition 5, (3.1), Proposition 6 (all
four Pochhammer conversions, the power q^r, the range min(u,n−u), the boundary cases, the use of Gasper–Rahman
(1.5.3) with b = q^{u−n}, c = q^{1−n}), the final assembly, the even/odd refinement, the mod-1 claim, the
½[(A1,D₃)+(A1,D₇)] remark, c(T1,E7) = 1 and C = −1/24.

## A. Log of checks on `proofs/E7.md`

Own code: `s1_fusion.py`, `s2_lemma3.gp`, `s2b_bailey.gp`, `s3_lemma4.py`, `s4_prop5.py`, `s4b_prop5_mid.gp`,
`s5_prop6.gp`, `s5b_cv.gp`, `s6_thm4.py`, `s6b_family.py`, `nahmgen.py`, `ps.py`.

A1. §1 (fusion). C(E7)⁻¹ for the stated labelling (chain 1–2–3–4–5–6, node 7 on node 4; det C = 2) recomputed
exactly: equals the printed 7×7 matrix entry by entry. The three fusion steps were redone with my own
implementation of Lemma 2: the printed 5×5, 6×6 and 7×7 stage matrices agree entry by entry; the reordered final
matrix is C(E7)⁻¹. ½(a,n)A⁽⁰⁾(a,n)ᵀ = ¾a² + a(N₂+N₃) + N₁²+N₂²+N₃² (symbolic). A⁽⁰⁾ is positive definite (det 4), hence
admissible; the document's argument (entries ≥ 0, diagonal > 0) is also valid. ½nᵀC⁻¹n − ¾(n₁+n₃+n₇)² is a quadratic
form with integer coefficients (symbolic), so the mod-1 claim and the even/odd refinement are correct; the
half-integer off-diagonal entries are exactly at {1,3}, {1,7}, {3,7}.

A2. Lemma 3. Definition of a Bailey pair, unit pair, limiting Bailey lemma: agree with Warnaar (1.3), (2.2), (2.3)
and Andrews (2.3)–(2.6), (2.12)–(2.13). (2.1) verified exactly as a rational function of x, q for N ≤ 6 (generic x)
and for x = q^a, 0 ≤ a ≤ 6, N ≤ 6 (a = 0 with α_r = (−1)^r q^{r(r−1)/2}(1+q^r)). The limiting Bailey lemma was verified
on a generic sequence (n ≤ 5), the full two-parameter lemma for n ≤ 1, the unit pair for n ≤ 6. x = 1 is handled
correctly because (x;q)_r/(1−x) = (xq;q)_{r−1} is a polynomial in x. The hypothesis x ≠ q^{−m} (m ≥ 1) is exactly the
non-vanishing of (xq;q)_{N+r}; adequate (the left side has no poles; the identity is one of rational functions).

A3. Lemma 4. Index shift m → m−j, the exponent identity j²+Rj+(m−j)(m+R+a) = m²+Rm+(m−j)(a−j), the range
0 ≤ j ≤ min(a,m), Lemma 1 with (a,b,c) = (a,m,j): correct; the j-sum is finite, so the interchange needs nothing.
(2.2) verified for 0 ≤ a,R ≤ 6 mod q^80; Cauchy for c ≤ 11 mod q^80; Lemma 1 for a,b ≤ 6 mod q^80.

A4. Proposition 5. Re-derived: correct, including a = 0, r = 0, the factor (q^a)_r/(1−q^a) = (q)_{a+r−1}/(q)_a and the
exponent 3r²+r²+r(r−1)/2 = (9r²−r)/2. The intermediate display T(a) = Σ_r q^{3ar+3r²}α_r Σ_m q^{(m+r)²}/((q)_m(q^{a+1})_{m+2r})
verified for a ≤ 4 mod q^50, and (2.3) against a brute-force triple sum for a ≤ 5 mod q^60. (T(0) is the modulus-9
Andrews–Gordon sum, as it should be.)

A5. §3. (3.1): exponent ¾n² + (3r²−r)/2, the a = 0 terms, c₀ = 1, summability (exponent ≥ ¾n², j ≤ a): correct.
Proposition 6: u = r+j; r ≤ min(u,n−u) is exactly {r ≥ 0, j ≥ 0, 2r+j ≤ n}; (q)_{n−r−1} = (q)_{n−1}/(q^{n−r};q)_r needs
r ≤ n−1, true; the three reflection formulas have the stated signs and powers; total sign +1, total power q^r;
(q^{1−n};q)_r ≠ 0 for r ≤ n−1; Gasper–Rahman (1.5.3) is ₂φ₁(a,q^{−n};c;q,q) = (c/a;q)_n aⁿ/(c;q)_n, used correctly (the sum
is a polynomial in b, so b = q^{u−n} is a legitimate specialisation; both numerator parameters terminating is
harmless; c = q^{1−n} is harmless because u ≤ n−1). Exact checks: q-Chu–Vandermonde for 9 rational pairs (b,c),
u ≤ 8; the specialised ₂φ₁ = 0 for n ≤ 16; c_n = 1+q^{n²} for 1 ≤ n ≤ 24; S_u = 0 (1 ≤ u ≤ n−1), S₀ = S_n = 1/(1−q^n) for
n ≤ 24; the u-regrouping for n ≤ 16.

A6. Theorem 4 by brute force over N⁷: all exponents ≤ 150 (601 coefficients in q^{1/4}, 10 399 lattice points):
equal to (q)_∞⁻¹(Σ_{k≥0}q^{3k²/4} + Σ_{k≥1}q^{7k²/4}). Exponents lie in Z ∪ (¾+Z). KKMM sums for D₃, D₄, D₅ (to q^30), D₇ (to
q^16) equal (q)_∞⁻¹Σ_k q^{rk²/4}. Theorem 5: k ≤ 5, exponents ≤ 30.

A7. Sun–Wang (fresh download of arXiv:2604.00847 and the supplied text). Conjecture 1.1: (A(X,Y),0,C(X,Y),D(X,Y))
modular, A = C(X)⊗C(Y)⁻¹, C = −c/24, c = tr D(X,Y)·h(X)/(h(X)+h(Y)); h(T_r) = 2r+1, D(T_r) = Id, h(E7) = 18,
D(E7) = Id₇. Hence A(T1,E7) = C(E7)⁻¹, c = 7·3/21 = 1, C = −1/24: as in the document. Table 1: "(T1,E7) 1 ? −";
(T1,E6) → (3.19), (T1,E8) → (3.20); no formula for (T1,E7) anywhere in the paper (the only other E7 statement is
"for Nahm sums of type (E6,T1) and (E7,T1), we have not found good CFT interpretations", which concerns the dual
pair). Wang–Wang 2607.08606 Conjecture 1.1 is the ADET folklore conjecture with the same C: confirmed.
"Conjecture 1.1 of Sun–Wang holds for (T1,E7)" does follow from Theorem 4 (with E3).

## B. Log of checks on `proofs/KRTW-k2.md`

Own code: `k1_chain.py`, `k2_numeric.py`, `k3_E.py`, `k4_grading.py`, `k5_E3.py`, `k6_melzer.py`.

B1. Source. `KRTW.tex`: the display after Conjecture 4.9 (label `Con_Shun2`; numbering reconstructed from the
source and confirmed in the PDF: Conjecture 2.2 = `Con_Shun`, 4.9 = `Con_Shun2`, Theorem 4.10 = `Thm_kis2`, (4.7) =
`Eq_CDn2`) has first identity Σ ∏ q^{(r_i+s_i)²+s_i²+s_i}/((q;q)_{r_i−r_{i+1}}(q²;q²)_{s_i−s_{i−1}}) =
(q,q^{k+1},q^{k+2};q^{k+2})_∞/((q;q²)_∞(q;q)_∞), r_{k+1} = s₀ := 0: as printed in the document. The initial (M,B,d) was
rebuilt symbolically from this exponent in the variables a₁ = r₁−r₂, a₂ = r₂, b₁ = s₁, b₂ = s₂−s₁: equals the
document's. What KRTW prove for k = 2: see K1.

B2. Lemmas. S: (z;Q)_∞(−z;Q)_∞ = (z²;Q²)_∞ and Euler's two sums; substitution n_a = m+2c gives the stated data.
P: (z;Q)_∞ = (z;Q²)_∞(zQ;Q²)_∞; n = s+t gives the stated data. F: Durfee in base Q = q^d adds d·m_a m_b. Admissibility:
forward — new exponent = old exponent at the image + a non-negative term, finite fibres; backward — the injections
n ↦ (m − c v, c = min(n_a,n_b)), n_a ↦ (n_a mod 2, ⌊n_a/2⌋), n ↦ (n,0) land on points where the extra term vanishes, so
{Q_old ≤ E} injects into {Q_new ≤ E}; then the forward statement gives equality. Correct; the backward direction
is justified, and the non-definite intermediate states are admissible.

B3. Chain. My forward moves were written from the "Consequence" paragraphs; each backward move constructs the
smaller data and requires forward(smaller) = current position by position. All 14 steps pass; final data =
((1 1; 1 2), (½,1), (1,1)). Definiteness of the 15 states: PD, PD, PD, PD, PSD, **indefinite**, PSD, PSD, PSD, PSD,
PSD, PSD, PD, PD, PD (K2). All 15 series agree with (−q;q)_∞(−q²;q²)_∞ for exponents ≤ 30.

B4. Evaluations. Final sum: Euler twice, as in the document; correct. Product: (q,q³,q⁴;q⁴)_∞/((q;q²)_∞(q)_∞) =
(q⁴;q⁴)_∞/(q)_∞ = (−q;q)_∞(−q²;q²)_∞ (three forms agree to q^30). The 4-fold sum in the original (r,s) variables equals the
product mod q^100. (E): x⁰,…,x¹⁰ mod q^40. x-refined (E₂): x⁰,…,x¹⁰ mod q^40; (E₁), (E₃): x⁰,…,x⁸ mod q^30.

## C. Source / novelty searches (2026-11-26) and their limits

Done:
- arXiv API (export.arxiv.org), 29 queries, among them: "Nahm sum" ∧ E_7 (0 hits); "Nahm" ∧ "E_7" (0); "Nahm sums" ∧
  "inverse Cartan" (1: Wang–Wang 2607.08606); "Nahm sum" ∧ "Dynkin" (2); abs "Nahm" ∧ "Dynkin" (5); all "Nahm"
  2026-04-01…2026-11-26 (27, all titles read); "Nahm sums" (49, all titles read); "Cartan matrix" ∧
  "Rogers-Ramanujan" (4); "Nahm" ∧ "Bailey" (6); "E_7" ∧ "Rogers-Ramanujan" (1: Han 2512.06480, perfect crystals);
  "(T_1,E_7)" (0); Cherednik–Feigin follow-ups; abstracts read of 2610.15587, 2610.00967, 2611.12868, 2608.08347,
  2609.04099, 2607.08606.
- Papers citing Sun–Wang: Semantic Scholar (3: 2609.09653, 2609.08816, 2607.08606) and INSPIRE (1: 2608.08780).
  All four re-downloaded as PDF and grepped for E7/E_7: no occurrence except Yoshida's "(E7)₁ × SU(2)₁ ⊂ (E8)₁".
  Wang–Wang 2607.08606 proves the (T1,D_r) case by Bailey pairs (methodologically close) and does not treat E_n.
- Papers citing Cherednik–Feigin 1209.1978: Semantic Scholar list (62 titles read; none about the ♭-case for E_n).
  Cherednik–Feigin §3.6 read in full.
- KRTW 2404.03851: source read around Conjectures 2.2, 4.1, 4.7, 4.9, Proposition 4.8, Theorem 4.10 and its proof.
Not done / not possible: full-text search (Google Scholar, MathSciNet, zbMATH, arXiv full text); Semantic Scholar
keyword search and the citation list of Nakanishi 1210.0226 (rate-limited); my INSPIRE queries returned no record
for Cherednik–Feigin or KRTW, so no INSPIRE citation lists for them; theses (Lee, Storzer); the 1990s fermionic-character literature (Kuniba–Nakanishi–
Suzuki, Bytsko–Fring, Kedem et al.) was not re-read by me — `PRIOR-ART.md` reports it contains only 2C⁻¹ sums, which
I did not verify; no search at all for prior occurrences of the rank-4 sum (1.1) or of the family of Theorem 5 in
the Bailey-lattice literature (a natural place for it: for k = 2 it is a Melzer identity); arXiv listings after
2026-11-25.
Assessment for E7: no anticipation found; the formula is not in Sun–Wang, in its four citing papers, or in
Cherednik–Feigin. Residual risk: an unindexed or full-text-only occurrence, and the possibility that Theorem 5
(hence Theorem 4 after three known fusions) is a known Bailey-chain identity in other notation.

## D. Reproduction

    cd /tmp/claude-0/-home-user-weighted-prime-races/3c4df68b-d770-5604-8028-1966b4f5b18d/scratchpad/referee3/code
    python3 s1_fusion.py; gp -q s2_lemma3.gp </dev/null; gp -q s2b_bailey.gp; python3 s3_lemma4.py; python3 s4_prop5.py
    gp -q s4b_prop5_mid.gp; gp -q s5_prop6.gp; gp -q s5b_cv.gp; python3 s6_thm4.py 150; python3 s6b_family.py
    python3 k1_chain.py; python3 k2_numeric.py 30; python3 k3_E.py; python3 k4_grading.py; python3 k5_E3.py; python3 k6_melzer.py

(`s2b_bailey.gp` and the generic-(b,c) line of `s5_prop6.gp` end with a PARI stack overflow on their last, most
expensive symbolic test; the results quoted above are those printed before that.)
