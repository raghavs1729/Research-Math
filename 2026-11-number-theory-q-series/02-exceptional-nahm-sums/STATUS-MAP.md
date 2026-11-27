# STATUS-MAP — conjectured Rogers–Ramanujan-type identities, status on 2026-11-24/26

Detailed, source-by-source notes (statements copied from the papers, with [F]/[A]/[I] labels = full text read /
abstract only / inferred) are in `status/A-kanade-russell.md`, `status/B-A2-cylindric.md`,
`status/C-nahm-multisums.md`, `status/D-misc.md`. This file is the ranked summary. "Open" always means: no
proof found by arXiv-API searches (to 2026-11-24), Semantic Scholar citation lists, the OpenAI math index
(374 families, no q-series item) and whataifound.org (only the Huang–Jiang–Oblomkov entries). Almost all 2026
closures are unrefereed preprints, several AI-assisted; none was verified by us.

## A. Closed (do not attack)

| Family | Closed by |
|---|---|
| Kanade–Russell mod 9 (I1–I4, I4a) | Mizuno 2609.09653 (I1–I3); Xia 2609.08816 (all five) |
| Kanade–Russell mod 12 (17 identities) | Bringmann–Jennings-Shaffer–Mahlburg 1809.06089; Rosengren 1912.03689; L. Wang 2609.07822 (last six) |
| Uncu–Zudilin / Warnaar / Hickerson / Konenkov reflections (18) | Xia 2609.24129 |
| Wang–Wang Conj. 3.6–3.8, Li–Wang Conj. 6.4 (Nahm duals) | Xia 2609.08816 |
| Nandi; A_2^(2) levels 5, 7; Capparelli; Afsharijoo companions | Takigiku–Tsuchioka; ADJM (earlier) |
| Takigiku–Tsuchioka A_13^(2), A_{2l-1}^(2) level 2 (8 identities) | L. Wang 2610.18666 |
| Huang–Jiang–Oblomkov (Conj. 1.7 and sum-to-sum) | Lau–Ono 2608.05480; Huang–Lau–Ono(–Paule) 2609.20567, 2608.15219 |
| A_2 Andrews–Schilling–Warnaar completion, cylindric KR Conj. 3.1/5.1, Warnaar 2111.07550 Conj. 2.1, 2.5, 2.8, 2.11, 7.4 | Warnaar 2303.09069; Uncu 2301.01359; Chern–Nguyen–Tran–Wakhare 2610.40079 |
| Zagier rank 2/3 examples; Mizuno's rank 2/3 tables | Wang 2210.10748, 2211.04375; Cao–Rosengren–Wang; Wang–Wang 2308.14728, 2402.06253, 2407.21725 |
| Tadpole Nahm sums (Calinescu–Milas–Penn; Milas–Wang) | Chern–Tran–Wakhare 2607.22904 |
| Cao–Wang, Shi–Wang, Li rank-3/4 conjectures | Li 2606.25866; Shi–Wang 2607.23257 |
| Sun–Wang (1.12) (C(D_r)^{-1}, B=0) | Wang–Wang 2607.08606 |
| CMPP type C_n^(1) (all levels) | Ding 2611.10654 (earlier cases Dousse–Konan, Primc–Trupčević, Russell) |
| Primc D_3^(2) level 2 | Bartholdi 2611.04020 |
| Warnaar affine Jacobi–Trudi Conj. 1.2–1.4 | Fischer–Gangl–Gutiérrez–Kumari 2610.08500 |
| Andrews–Uncu; Russell 3-colour mod 10; Garvan mock theta convolution | Chern; Chern 2607.03977; Bradley-Thrush et al. 2608.13324 |

## B. Open, ranked by (standing) x (chance); what we did

| # | Conjecture (source) | Statement type | Our status |
|---|---|---|---|
| 1 | **Sun–Wang 2604.00847 (1 Apr 2026), eq. (3.20)**: (T1,E8) Nahm sum = vacuum character of M_eff(11,2); **eq. (3.19)**: (T1,E6) Nahm sum = χ^{M_eff(5,2)}(χ_0+χ_3+2χ_{2/3})^{M(6,5)}. Both "checked to high order"; cases of their Conj. 1.1 (Dynkin-pair modularity). Yoshida 2608.08780 (3.77)–(3.80): four E8 companions, checked to q^6 | 8-fold / 6-fold sum = product | **PROVED** (`proofs/E8-E6.md`) |
| 2 | **KRTW 2404.03851** (Selecta Math. 2026), four conjectural Andrews–Gordon-type identities after Conj. 4.9; k=1,2 known; conjectural for k≥3 | sum = product | **(i) PROVED for all k** (`proofs/KRTW-general.md`; ⇒ Conj. 2.2 at z=1); **(ii), (iii) PROVED for k = 3** (`proofs/KRTW-F-moves.md`; ⇒ Conj. 2.2 at z=q, k=3), reduced for all k to a Hall–Littlewood statement; (ii), (iii) for k≥4 and (iv) for k≥3 open (tested mod q⁴⁰⁰, k≤8) |
| 3 | Sun–Wang Conj. 1.1 for (T1,E7) (no formula given; c=1) | modularity | **PROVED** with new explicit formula (1/(2(q)_∞))Σ_k(q^{3k²/4}+q^{7k²/4}) (`proofs/E7.md`) |
| 4 | Sun–Wang: r−1 companion modular Nahm sums for C(D_r)^{-1}; Wang–Wang found ⌊(r+4)/2⌋, "⌊(r−3)/2⌋ remain" | existence | not attacked |
| 5 | Sun–Wang Conj. 1.1, rank-3 pairs (C3,A1), (C3,T1), (B3,A1), (B3,T1) (per Xia, Sept 2026) | modularity, no formula | **(C3,A1), (B3,A1): modular** — total twisted parafermion characters of A_5^{(2)}, D_4^{(2)} at level 2 (Okado–Takenaka + Kac–Peterson; `proofs/X-A-parafermion.md`, Thm 2); (C3,T1), (B3,T1) open (= level-1 principal-subspace characters of the same twisted algebras) |
| 6 | Gang–Kim–Park–Stubbs 2411.06081 §4.1.3: rank-3 signed Nahm sums = Rogers–Ramanujan / mod-7 products ("conjectures based on the q-series expansion") | sum = product | **18 of 24 examples PROVED** (`proofs/GKPS-signed.md`; elementary); 3-1, 3-2(3), 3-15, 3-19–3-22, 3-24 not treated |
| 7 | KRTW Conj. 2.2, 4.9, 1.4, 4.7, 5.1, 5.2 (Hall–Littlewood / partition g.f. = multisum, bivariate) | g.f. = sum | not attacked |
| 8 | CMPP type A_{2n}^(2) (n≥2,k≥2); Primc D_{n+1}^(2) (n=2,k≥3; n≥3) | partitions = product | not attacked |
| 9 | Li–Uncu 2501.19272 Conj. 1.3, 1.6, 7.1, 7.2 (finite Andrews–Gordon/Bressoud, k≥5) | polynomial identities | not attacked |
| 10 | Huang–Jiang–Oblomkov 2608.16086 Conj. 1.8 (linear-term variants, a≥3) | sum = product | not attacked |
| 11 | Hakimoglu-Brown 2609.14940 Conj. 6.1 (mod 13), 2609.28509 Conj. 1.2 (index 31) | double sums = products | not attacked (single-author preprints; low standing) |
| 11a | Sun–Wang 2604.00847 eq. (3.6): (T1,C_r) generalized Nahm sum = Σ q^{½(N_1²+…+N_r²)}/((q)_{n_1}⋯(q)_{n_{r−1}}(q²;q²)_{n_r}) = χ₀^{NS} − χ^{NS}_{⌊r/2⌋+1} of SM_eff(4r+6,4), "checked up to r=7"; r=2,3 known (Mizuno; Wang–Wang). Not covered by Du–Ji–Shen–Xu 2602.14483 (different matrices); may be contained in Bressoud–Ismail–Stanton / Warnaar base-change results (not checked) | sum = theta/eta | not attacked |
| 12 | Wang–Wang 2512.07790 concluding conjecture (½C(D_k), k≥5) | modularity | **PROVED for all k** as a corollary of Butorac–Kožić–Primc + Kac–Peterson (`proofs/X-A-parafermion.md`, Thm 1: D_k^{(1)} level-2 parafermions, c = k−1) |
| 12b | Cao–Wang leftovers; Kanade–Russell 2203.05690 Conj. 6.5; Warnaar 2511.17034 §6 | various | not attacked |

Exact statements for every row: see the `status/` files (sections cited there).

| 13 | Sun–Wang Conj. 1.1 for (T_1, A_r), odd r ≥ 5 (Table 1: "?") | modularity | identified: = η·Σ(A_2^{(2)} string functions, level r+1), conditional on Takenaka arXiv:2205.05271 (preprint); modularity then follows from Kac–Peterson for A_2^{(2)} — **citation not verified** (`proofs/X-A-parafermion.md`, Thm 3, Cor. 5; referee) |
| 14 | Sun–Wang Conj. 1.1 for all (X, A_{k−1}) | modularity | **PROVED** (same file), conditional on the cited theorems; (X, T_k) = principal-subspace characters, open for X ≠ A |
