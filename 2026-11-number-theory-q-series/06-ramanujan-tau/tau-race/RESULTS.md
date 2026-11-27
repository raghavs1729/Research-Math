# tau-race: results

**Summary (15 lines).**
1. Scope: work-order item T1 only (a fixed-sign theorem at every x for a natural smoothed sum of τ(p)). Work-order item T2 (density one for Σ τ(p)p^{-6} > 0 under GRH alone) is already Sheth 2025, Cor. B — not to be confused with Theorem 2 below, which is new as far as our searches go; work-order items T3/T4 were dropped by instruction.
2. Notation: λ(p) = τ(p) p^{-11/2}; R(x) = Σ_{p≤x} λ(p)(1 − p/x) (Cesàro mean of Σ_{p≤t} λ(p)); Θ_1(t) = Σ_{p≤t} λ(p) log p (t − p).
3. **Theorem 2.** Under GRH for L(s, Δ) alone: R(x) > 0 for EVERY x > x* = 5.518054631888…, and R(x) ≤ 0 for x ≤ x*. Also R(x) > 0.15 √x/log x for every x ≥ 11.
4. **Theorem 1.** Under GRH(Δ) alone: Θ_1(t) > 0 for EVERY t > t* = 4.266791956347…, with Θ_1(t) ≥ L_B(t) t^{3/2}, L_B explicit, L_B(10^6) ≥ 0.027, lim L_B ≥ 0.369. With GRH(Sym²Δ) in addition: limit constant 0.5539.
5. **Theorem 3.** GRH(Δ) ⇔ R(x) > 0 for all x > x* ⇔ R(x) ≥ 0 for all large x (same for Θ_1). (Converse by Landau's theorem.)
6. **Theorem 4 (no hypothesis, computer-assisted).** All zeros of L(s, Sym²Δ) with |Im s| ≤ 60 are simple and on the critical line (41 with 0 < Im ≤ 60). This is what removes GRH(Sym²Δ) from Theorems 1–3.
7. Sign and sizes: the bias is POSITIVE (prime squares: λ(p)² − 2 has mean −1). Bias 2/3; zeros of L(s,Δ): σ_0² = Σ|ρ|^{-2} = 2Λ'/Λ(1) = 0.11225214047 (certified, closed form, no zero list); zeros of L(s,Sym²Δ) above height 60: ≤ 0.184605 (unconditional bound), below: lower order t^{5/4}.
8. Deligne's bound alone does NOT suffice for the prime-square term (it only gives |…| ≤ (4/3) t^{3/2} against the bias 2/3); the identity Σ_{even powers} = ψ_{Sym²}(√u) − ψ(√u) and partial RH for Sym² are what make it work.
9. Analytic range: t ≥ 10^6 (bound positive from 7.1·10^5; from 1.12·10^5 under GRH(Sym²Δ)). Exact/ball range: every prime ≤ 9 999 991 (τ(p) exact from η^24). Overlap: one to two orders of magnitude.
10. Certificates (sha256 in TAU-RACE.md §7): constants.json, small_range_10000000.json, sym2_phi_M640.json, sym2_rh_H60.json; all Arb ball arithmetic / exact integers.
11. Second implementations: PARI/GP for every constant, the residue bookkeeping of both explicit formulas, the zeros of Sym²Δ, and the small range; the functional equation of Sym²Δ checked directly to 10^{-83}.
12. Mutation tests: 4 families, all mutants rejected (selftest.py).
13. Novelty (LITERATURE.md): no all-x or explicit-threshold statement for τ, no such equivalence, and no rigorous RH check for Sym²Δ was found; Koyama–Kurokawa need DRH for an asymptotic at the critical weight; Sheth gives density one under GRH; Sarnak's letter gives the mean 1 of the limiting distribution.
14. Not proved here: anything pointwise for Σ τ(p)p^{-6} or for the un-smoothed sum; other smoothings at every x; T3/T4.
15. Review: adversarial referee (AI agent): no fatal error, 2 repairable + 10 minor items, all applied (REFEREE-REPORT.md, TAU-RACE.md §10); independent numerical audit (second AI agent, own PARI/mpmath code incl. an argument-principle zero count): 6/6 pass (AUDIT-NUMERICS.md); second referee commissioned by the parent session: no error, three citation/presentation gaps, repaired (TAU-RACE.md §10). The completeness test of Theorem 4 is a variant of the Hiary–Ireland–Kyi method (arXiv:2408.00187). No human review.

Files: `TAU-RACE.md` (statements, complete proofs, certificate table), `LITERATURE.md` (prior art with
quotations), `REFEREE-REPORT.md`, `AUDIT-NUMERICS.md`, `SOURCES-CHECK.md` (sources read for the cited inputs), `PROGRESS.md`, `code/`, `data/`.
