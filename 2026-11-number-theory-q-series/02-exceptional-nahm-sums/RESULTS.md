# RESULTS — rr-identities (2026-11-26)

**Headline: the Nahm sums of the inverse Cartan matrices of E6, E7, E8 are evaluated and proved modular.**
F_R(q) := Σ_{n∈N^r} q^{½nᵀC(R)⁻¹n}/∏_i(q;q)_{n_i}. These are the "♭-case" Rogers–Ramanujan series of Cherednik–Feigin
(Adv. Math. 248 (2013), §3.6: "our approach … cannot be directly extended to the ♭-case"; modularity only
conditional on Nahm's conjecture) and the pairs (T1,E_n) of Sun–Wang, arXiv:2604.00847 (1 Apr 2026), Conj. 1.1.

| | Statement | Source of the conjecture | Proof status | Depends on |
|---|---|---|---|---|
| E8 | F_{E8} = 1/(q,q²,q³,q⁴,q⁷,q⁸,q⁹,q¹⁰;q¹¹)_∞; four companions with linear terms give the other mod-11 Andrews–Gordon products; 4-parameter refinement = Andrews–Gordon 4-fold sum | Sun–Wang eq. (3.20) (numerical); companions: Yoshida 2608.08780 (3.77)–(3.80) (numerical to q⁶) | **Proved** (`proofs/E8-E6.md`, Thm 1, Cor. 1). Two adversarial reads (ours: `REFEREE-E8-E6.md`; parent session's independent one: "proof complete") | Durfee rectangle identity; Andrews–Gordon (k=5) |
| E6 | F_{E6} = G(q)(J_{7,15} − qJ_{2,15} + 2q^{2/3}J_{3,15})/J_1 (G = Rogers–Ramanujan); five companions | Sun–Wang eq. (3.19) (numerical) | **Proved** (Thm 2, Cor. 2), same two reads | Durfee; Rogers–Ramanujan; Cao–Rosengren–Wang Thm 1.4 (JCTA 2024) on Zagier's Example 10 (their proof ends with a Maple-verified theta identity) |
| E7 | F_{E7} = (1/(2(q)_∞)) Σ_{k∈Z}(q^{3k²/4} + q^{7k²/4}) (**new formula**; Sun–Wang: CFT "?") | Sun–Wang Conj. 1.1 for (T1,E7) (modularity only) | **Proved** (`proofs/E7.md`, Thm 4); adversarial read `REFEREE-E7-KRTW.md`: no mathematical error (all lemmas re-derived; brute force to exponent 150); six minor citation/wording defects, fixed | Durfee; Bailey's lemma + unit Bailey pair; Cauchy's identity; q-Chu–Vandermonde |

Further results
1. **Theorem 5 (`proofs/E7.md` §6):** for every k ≥ 1, Σ q^{(k/4)a²+a(N_2+…+N_k)+N_1²+…+N_k²}/((q)_a∏(q)_{n_i}) = (1/(2(q)_∞))Σ_m(q^{km²/4}+q^{(k+4)m²/4})
   (k=3 is E7 after fusion). k=2 is a KNOWN identity (the (T1,T_4) case of Melzer's conjecture, proved by
   Bressoud–Ismail–Stanton Thm 5.1 and Warnaar): for k=2 only the proof is new. Same proof for all k; checked
   k ≤ 6 mod q⁸⁰ (parent's referee: k ≤ 8 to q¹⁵⁰). Parent session's independent referee on E7.md: Theorems 4 and 5
   "proof complete"; its findings F1–F5 are fixed.
2. **`proofs/KRTW-k2.md` (NOT a new theorem — corrected after the referee's finding K1):** an elementary proof
   (machine-found chain of 14 exact moves) of the k=2 case of the first Andrews–Gordon-type identity of
   Kanade–Russell–Tsuchioka–Warnaar (arXiv:2404.03851): 4-fold sum = (−q;q)_∞(−q²;q²)_∞. That case already follows
   from KRTW Thm 4.10 + (4.7) + Meurman–Primc. New there: only the x-refined identity (E) and the numerical
   reduction of the open cases k ≥ 3 to a k-fold Melzer-type sum (conjecture (E_k), k ≤ 4).
3. Method: "fusion" = Durfee rectangle identity applied to two summation variables (the lemma is known:
   Cao–Wang 2606.13590 (1.5)–(1.6); quiver "linking"). New: C(E8)⁻¹ = 4 fusions of 2min(i,j); C(E6)⁻¹ = 3 fusions
   of 2C(A2)⁻¹⊕(2); C(E7)⁻¹ = 3 fusions of an Andrews–Gordon chain coupled to one extra variable.

Second round (2026-11-26), each with an adversarial read (`REFEREE-KRTW-general.md`, `REFEREE-T1-ABC.md`)
4. **Proved: KRTW identity (i) for ALL k ≥ 1** (`proofs/KRTW-general.md`, Theorem 3): the first of the four
   conjectural Andrews–Gordon-type identities of Kanade–Russell–Tsuchioka–Warnaar (arXiv:2404.03851, after
   Conj. 4.9; conjectural there for k ≥ 3). Consequence: their Conjecture 2.2 at z = 1 for all k. New step: the
   x-refined reduction of the 2k-fold sum to the k-fold Melzer-type sum (finite refinement + a q-binomial
   lemma); the evaluation of the k-fold sum is known (Andrews, Kim–Yee, …). Referee: no error. (ii)–(iv): not
   proved; hold mod q⁴⁰⁰ for k ≤ 8; reduced to a factorisation conjecture (F_k).
5. **Proved (conditional on Lepowsky–Primc and Kac–Peterson, cited): (T1,A_{2k}) for all k**
   (`proofs/T1-A-even.md`): Σ q^{½nᵀC(A_{2k})⁻¹n}/∏(q)_{n_i} = Σ_N c⁰_{2N}(τ)q^{N²/(k(2k+1))}·q^{c/24} (level-k string functions), hence
   Sun–Wang Conj. 1.1 for (T1,A_{2k}) (their Table 1: "?"). Referee: no error; novelty risk medium.
6. **Not new (referee):** (T1,B_r) — the identities of `proofs/T1-BC.md` Theorem B are Warnaar's and
   Du–Ji–Shen–Xu's; only the identification with Sun–Wang's "?" entry is ours. (T1,C_r) = Warnaar.
7. **Audit of Sun–Wang's formulas** (`SUN-WANG-AUDIT.md`, 137 comparisons to q¹⁰⁰–q⁵⁰⁰): no false identity; misprints:
   (3.27)/(3.28) have the Ramond weights 1/32 and 5/32 interchanged (as printed both fail at the leading term);
   "identical to f_1" before (3.30) should be f_2; (3.6) prints the same case twice.
8. **Tested only:** Yoshida's five (E8,T1) sectors to q²⁶; new conjecture (E7,T1) = (Θ_{E6}+Θ_{A6})/(2(q)_∞⁶) to q⁴⁰.

Status of every third-session item (what a later referee should attack)

| Item | File | Status | Depends on | Adversarial reads |
|---|---|---|---|---|
| KRTW identity (i), all k | `proofs/KRTW-general.md` | **proved** (also machine-checked in Lean by the parent: `monotonicity/lean-krtw`, theorem `KRTW.krtw_identity`) | classical q-series only | two reads + Lean |
| KRTW identities (ii), (iii) at k = 3; (F_2), (F_3) | `proofs/KRTW-F-moves.md`, Theorem A | **proved** | Durfee rectangle identity, Euler splitting; Bressoud's even-modulus identities / Jacobi triple product for the product sides | child referee (`REFEREE-KRTW-F-moves.md`) + parent referee (`REFEREE-2-KRTW-F.md`): proof complete |
| Bailey structure of B_{k+1}, J_k, and of B_{k+1}J_k, all k | same, Theorem B | **proved** (short consequences of Bailey's lemma; not claimed as new) | Bailey's lemma | child referee |
| (F_k) ⇔ Hall–Littlewood form; (ii) ⇔ KRTW Conj. 2.2 at z = q; Conj. 2.2 at z = q for k = 3 | same, Corollary C | **conditional** | Bartlett–Warnaar C_n Rogers–Selberg identity; Griffin–Ono–Warnaar Thm 1.3 (both cited, not re-derived) | child referee |
| (F_k), (ii), (iii) for k ≥ 4; (iv) for k ≥ 3 | — | **open** (numerics only) | — | — |
| Sun–Wang Conj. 1.1 for (X, A_{k−1}), X = A, D, E (incl. Wang–Wang's ½C(D_k) conjecture) | `proofs/X-A-parafermion.md`, Theorem 1 | **corollary of published theorems, not new mathematics** (it is Kuniba–Nakanishi–Suzuki's formula) | Butorac–Kožić–Primc Thm `thm38` (Math. Z. 2021); Lepowsky–Wilson; Kac–Peterson Theorem A (untwisted; cited, primary source not accessed) | child referee (`REFEREE-X-A-parafermion.md`): identities and prefactor correct |
| same for X = B_l, C_l, F_4, G_2 (incl. Sun–Wang's (C_3,A_1), (B_3,A_1)) | same, Theorem 2, Lemma 4 | **corollary of published theorems; identification and prefactor are new bookkeeping** | Okado–Takenaka (ART 2022; Correction not read); Lepowsky–Wilson; modularity: Kac §13.9 + Kac–Peterson Theorem A for the transposed algebra (cited, not accessed) **or** Proposition M | child referee: identities correct; modularity citation flagged |
| same for X = T_l (incl. (T_1, A_r) for all r) | same, Theorem 3, Corollary 5 | **conditional** | Takenaka arXiv:2205.05271 (unrefereed preprint) for the identity; modularity: Proposition M (or Kac–Peterson for A_{2l}^{(2)}, not accessed) | child referee: identity confirmed numerically (r ≤ 7), modularity = gap |
| Proposition M (string functions of every affine type are modular on a congruence subgroup; classical theta-function argument) | same, §2 | **written, NOT adversarially read** — attack step (4) (uniform multiplier of theta series with characteristics, odd lattice for A_{2l}^{(2)}) and step (7) (cusps) | Weyl–Kac formula; classical theta transformation (Schoeneberg) | none |
| 18 of the 24 rank-three signed Nahm sum identities of Gang–Kim–Park–Stubbs | `proofs/GKPS-signed.md` | **proved (elementary), NOT adversarially read** — check each change of variables and the transcription of the source data | q-binomial theorem, q-Chu–Vandermonde, Euler, Cauchy/Durfee, finite pentagonal number theorem, Rogers–Ramanujan, Andrews–Gordon, Jacobi triple product | none (certificate `code/gkps_certify.py`, 87 sums; an adversarial read was started and stopped for budget reasons, nothing came of it) |
| (X, T_k) = principal-subspace characters of g(X) | `RESULTS.md` item 14 | **reformulation only, no proof** | same papers | none |

Third session (2026-11-26), KRTW identities (ii), (iii) — `proofs/KRTW-F-moves.md`
9. **Proved: KRTW identities (ii) and (iii) for k = 3** (second and third displayed identities after Conjecture 4.9
   of arXiv:2404.03851v2; conjectural there for k ≥ 3; k = 2 followed from KRTW Thm 4.10 + Bartholdi), via the
   x-refined factorisation (F_k): K_k(x) = B_{k+1}(x)·J_k(x), proved for k = 2, 3 by explicit chains of two exact
   moves (Durfee rectangle identity; 1/(q)_n = Σ_{u+v=n} q^v/((q²;q²)_u(q²;q²)_v)) on generalized Nahm data. Exact
   certificate `code/krtw_moves_certify.py`. Consequence (with Griffin–Ono–Warnaar Thm 1.3): KRTW Conjecture 2.2
   holds at z = q for k = 3. The chains were found by computer search; no general pattern is known.
10. **Structure for all k (proved):** both factors and their product in closed form (Theorem B there):
   B_{k+1}(x) = (−xq²;q²)_∞·(base-q² Bailey chain), J_k, and
   (xq)_∞(x²q²;q²)_∞ B_{k+1}J_k = the C₂ Rogers–Selberg series L^{(0)}_k of Bartlett–Warnaar/Griffin–Ono–Warnaar at
   (√x, −√x). Hence (conditional on the cited C₂ Rogers–Selberg identity) (F_k) ⇔ K_k(x) = Σ_{λ_1≤k} q^{|λ|} P_{2λ}(±√x q^i; q)
   — a Hall–Littlewood form parallel to KRTW's Conjecture 2.2 but in the grading x^{Σ(r_i+2s_i)} and base q; and
   (ii) ⇔ Conjecture 2.2 at z = q. For even k = 2m: K_{2m} = (−xq;q)_∞ 𝔍^{(1)}_m 𝔍^{(2)}_m (product of two base-q² Andrews–Gordon
   chains with seeds 1/(−xq;q²)_σ, 1/(−xq²;q²)_σ), equivalent to (F_{2m}).
11. **Not proved:** (F_k), (ii), (iii) for k ≥ 4; (iv) for k ≥ 3. Exact facts about the obstruction (only the
   x-grading survives; a commuting rank-two q-Toda Hamiltonian; product-form spectral coefficients) are in §6 of
   the proof file.

Third session, Sun–Wang pairs (X, A_{k−1}) — `proofs/X-A-parafermion.md` (corollaries of published theorems)
12. **Sun–Wang Conjecture 1.1 for the pairs (X, A_{k−1}), k ≥ 2** — corollary of cited theorems (`proofs/X-A-parafermion.md`;
   adversarial read `REFEREE-X-A-parafermion.md`: identities and prefactor correct, attribution and the modularity
   citations corrected). The generalized Nahm sum of C(X)⊗C(A_{k−1})^{−1} is the character of the whole parafermionic
   space of the level-k vacuum module of g(X) (X^{(1)} for ADE; A_{2l−1}^{(2)}, D_{l+1}^{(2)}, E_6^{(2)}, D_4^{(3)} for C_l, B_l, F_4, G_2; A_{2l}^{(2)} for T_l):
   q^{−c/24}F = ∏_iη(d_iτ)·Σ_{μ∈Q̄/kQ̄} c^{kΛ_0}_μ(τ) (string functions).
   * X simply laced: **not new** — this is Kuniba–Nakanishi–Suzuki's 1993 formula (total parafermion character, prefactor
     q^{−c_PF/24}), proved by Georgiev / Butorac–Kožić–Primc, plus Kac–Peterson. The observation is that it settles
     Sun–Wang's conjecture for these pairs (no paper found saying so).
   * X = B_l, C_l, F_4, G_2: via Okado–Takenaka (ART 2022) for the twisted algebras; **new bookkeeping**: identification with
     Sun–Wang's pairs and Lemma 4 (Sun–Wang's c = tr(D)h(k−1)/(h+k) is exactly the Kac–Peterson anomaly plus tr(D)/24).
     Modularity via the S-transformation to the transposed untwisted algebra (Kac §13.9 + Kac–Peterson Theorem A; cited).
     Precedent: Mizuno (2305.02267 §4) for (C_2,A_1), (G_2,A_1).
   * X = T_l: identity with A_{2l}^{(2)} string functions via Takenaka (arXiv:2205.05271, preprint); **modularity for a congruence
     subgroup is an unverified citation** (Kac–Peterson for A_{2l}^{(2)}; primary source not accessible).
   Checks: `code/xa_parafermion_certify.py` (all sectors vs Weyl–Kac multiplicities, 41 cases), referee's independent
   Freudenthal code (`code/xa_referee/`), numerical modularity tests (prefactor confirmed to 10 digits).
13. **Consequences of 12** (each stated as open/conjectural in 2025–26 papers): (a) Wang–Wang's conjecture
   (arXiv:2512.07790): f_{½C(D_k),0,(1−k)/24} modular for all k ≥ 3 [simply laced case: KNS + BKP + KP]; (b) Sun–Wang's rank-three
   cases (C_3,A_1), (B_3,A_1) (unresolved per Xia 2609.08816) [twisted case]; (c) (A_{n−1},A_{k−1}) without gcd((n−1)!,k) = 1; all
   (E_n,A_{k−1}), (F_4,A_{k−1}), (G_2,A_{k−1}); (d) **(T_1, A_r) for all r**: q^{−r/(8(r+4))}Σq^{½nᵀC(A_r)^{−1}n}/(q)_n = η·Σ_{m mod r+1}(A_2^{(2)} string functions of level
   r+1) — the identity is conditional on Takenaka's preprint (verified numerically for r ≤ 7), modularity on the
   unverified citation above; for even r modularity is proved independently in `T1-A-even.md`. So the statement
   "Sun–Wang's conjecture holds for (T_1, X) for every simply laced X" is: A_even (ours, unconditional up to
   Lepowsky–Primc/Kac–Peterson), D_r (Wang–Wang 2607.08606), E_6, E_7, E_8 (ours), **A_odd conditional** as just described.
14. **Reformulation (not a proof):** by the same papers the pairs (X, T_k) are the principal-subspace characters
   ch W(kΛ_0)|_{y=1} of g(X) (matrix C(X)D(X)⊗(min(s,t))_{s,t≤k}); Sun–Wang's conjecture for (X,T_k) is the modularity of these
   (known for X = A_r). This covers the open cases (E_n,T_1), (C_3,T_1), (B_3,T_1).

15. **Gang–Kim–Park–Stubbs signed Nahm sums (arXiv:2411.06081 §4.1.3; "conjectures based on the q-series expansion"):**
   18 of the 24 rank-three examples (all versions and Wilson lines) proved in `proofs/GKPS-signed.md` by the q-binomial
   theorem / q-Chu–Vandermonde / Euler ("integrating out" the signed variable) + Rogers–Ramanujan / Andrews–Gordon /
   Jacobi triple product; certificate `code/gkps_certify.py` (87 sums). Not treated: 3-1, 3-2(3), 3-15, 3-19–3-22, 3-24.
   **Not adversarially read.**

Verification
* Exact certificates: `code/certify_E8_E6.py`, `code/certify_E7.py`, `code/certify_KRTW_k2.py` (logs in `certificates/`).
* Power series: E8 mod q⁵⁰⁰ (vacuum), q³⁰⁰ (all five; referee q²⁰¹, parent's referee q⁴⁰⁰); E6 mod q³⁰⁰/q¹⁵⁰; E7 mod q³⁰⁰,
  with each intermediate formula (T(a), the finite identity c_n = 1+q^{n²} for n ≤ 30) checked separately.
* Lean 4 (`lean/`, no sorry, axioms propext/Classical.choice/Quot.sound): Durfee identity, fusion algebra,
  E8 and E6 chains, Cartan inverses. Not in Lean: summability/regrouping, Andrews–Gordon, Bailey's lemma, E7,
  the final power-series identities.

How openness was checked: `PRIOR-ART.md` (table of sources) and `STATUS-MAP.md`. Sun–Wang state (3.19)/(3.20)
as numerical and give no formula for E7; the papers citing them, Cherednik–Feigin, the 1990s fermionic-character
papers (2C⁻¹ only), the Wang-school papers 2022–2026, the quiver literature, the OpenAI math index and
whataifound.org contain no proof. No full-text search engine was available.

Novelty risk: moderate for E8/E6 (four/three applications of a known lemma to known identities; an unindexed
independent proof is plausible); lower for E7 (the formula itself is new and the proof is less immediate), but
the same caveat on searches applies. The KRTW k=2 item is a re-proof of a known case.

Not proved (numerical only): KRTW identity (i) for k ≥ 3 and (ii)–(iv); see `STATUS-MAP.md` for the rest of
the frontier (Sun–Wang D_r companions, rank-3 Dynkin pairs, (T1,A_r) r ≥ 4, (T1,C_r), Gang–Kim–Park–Stubbs, …).

Files: `proofs/` (E8-E6.md, E7.md, KRTW-k2.md), `REFEREE-*.md`, `PRIOR-ART.md`, `STATUS-MAP.md`, `status/`,
`PROGRESS.md`, `code/`, `certificates/`, `lean/`.
