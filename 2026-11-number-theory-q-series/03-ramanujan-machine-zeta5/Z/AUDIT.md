# Audit of `Z/THEOREM.md` (Theorem Z: Zudilin's relations (6) for r = t = 1/2)

Hostile referee audit, 2026-11-26. I did not edit `THEOREM.md`. My scripts are in `Z/audit/`. I did not use
`Z/verify/z_checks.py`.

## Verdict

**Mathematics: complete. The novelty section needs correcting.**

I checked every step of Lemma 1, Lemma 2, §4 and the proof of U2 Proposition 8 by hand. All the constants, signs,
substitutions and interchanges are correct, and all the claimed identities hold numerically to 40 digits. I found no
gap in the proof.

The description of prior art in §0 is **not accurate**. The step the file calls "what is added here",
`∫_0^1 K'^2 dk = −π^2 S_1` (equivalently `F_1(1) = −(8/π^2) ∫_0^1 K'^2`), is already in Borwein–Straub–Wan
(2010 preprint, Exp. Math. 2013). Their three formulas for `W_4(−1)` give it directly. Rogers–Wan–Zucker (31) then
gives `L(f,1)` in terms of the same integral, and also as a 7F6 value, which the file does not quote. So
`L(f,1) = −S_1/2` follows from published results in a few lines (finding N1). The constants `−1/16, 1/32, −3/256`
for Zudilin's `F_j(1)` are not stated anywhere I could find. They follow at once from the published results above
together with Zudilin's own Remark 2 and Remark 5.

## Findings

**N1 (major; novelty and attribution).** Borwein, Straub and Wan, *Three-step and four-step random walk integrals*,
preprint 23 Sept 2010 (https://carmamaths.org/jon/walks2.pdf), Exp. Math. 22 (2013), doi:10.1080/10586458.2013.748379.
Equation numbers below are from the preprint.
* The display right after (29) ("For comparison, (27) produces") reads
  `W_4(−1) = (4/π) Σ_n C(2n,n)^4/4^{4n} Σ_{k≥2n+1} (−1)^{k+1}/k`.
  The inner sum equals `log 2 − (H_{2n} − H_n) = −h_n`, so this is `W_4(−1) = −(4/π) S_1`. It is derived from their
  analytic continuation (27), a 4F3 formula, using L'Hôpital's rule.
* Their (65) is `W_4(−1) = (4/π^3) ∫_0^1 K'(k)^2 dk`.
* Their Theorem 4 is `W_4(−1) = (π/4) 7F6(5/4, 1/2^6; 1/4, 1^5; 1)`.

Together these give `∫_0^1 K'^2 = −π^2 S_1`, which is exactly the content of §3 and §4 that the file calls new.

Rogers–Wan–Zucker (31) says
`L(f_2,3) = (π^2/4) L(f_2,1) = (1/4) ∫K^2 = (π^4/128) 7F6(5/4,1/2^6;1/4,1^5;1)`. The file quotes (31) without its
7F6 part, and that part already gives `L(f,1) = (π^2/32) Σ (4k+1) a_k^3`. RWZ cite J. Wan's paper for the 7F6, and
Wan's (25) at n = 0 is the moment formula for `∫ K'^2`. With BSW this gives `L(f,1) = −S_1/2` immediately.

Required fix: rewrite the "Status and novelty" paragraph. It should say that every analytic ingredient is published
(RWZ (30), (31); BSW 2010/2013; Zudilin, Remark 2 for `F_2 = 2π^2 F_0` and Remark 5 for `F_3 ∝ π^2 F_1`). What is
left is the observation that these combine into (6) for (1/2,1/2) with the constants −1/16, 1/32, −3/256, plus a
self-contained proof (Lemmas 1–2). The sentence "what is added here is ... the evaluation of the derivative value
F_1(1) = ... (§3)" must go or be qualified.

The same correction applies to `../U2/THEOREM.md` §0. Its claim that "neither I nor the auditor found a proof" of
`8 S_1 = −16 L(f,1)` was already false given BSW + RWZ.

**N2 (minor; numerical typo).** §5 Remark 1 gives `S_0 = 1.118636387164190...`. The correct value is
`S_0 = 1.11863638716418706834961925753...` (mpmath `hyper` and `nsum` agree to 40 digits). The value
`L(f,2) = 0.690031163123397525119...` in the same remark is correct and does equal `π^2 S_0/16` with the corrected
`S_0`, so the error is only in the printed digits of `S_0`.

**N3 (minor; quotation).** Zudilin's eq. (5) is the *bilateral* sum over `n ∈ Z`, equated to the unilateral sum
`n ≥ 0` plus `O(ε^4)`; the `n < 0` terms vanish to order 4 because of `Γ(1+ε+n)^{-4}`. The file quotes only the
unilateral form as "eq. (5)". The mathematics is unaffected, but the quotation should show both forms or say
"(5), after dropping the n < 0 terms, which are O(ε^4)".

**N4 (cosmetic).** Zudilin's Table 1 gives the LMFDB label in the old format, 8.4.1.a. The file's 8.4.a.a is the
current label for the same form. Zudilin's Remark 2 already states `F_0(1) = F_2(1)/(2π^2)` for r = t = 1/2, so the
`F_2` relation could be credited to him as well as to Prop. 8.

**Steps checked by hand. All correct; no fix needed.**
1. **(Z1).** With `Q = q^2`, `q = e^{πiσ}`:
   * `θ_2 = 2η(2σ)^2/η(σ)`; both sides have leading term `2q^{1/4}`.
   * `θ_4 = (q;q)^2/(q^2;q^2) = η(σ/2)^2/η(σ)`.
   * `θ_3 = η(σ)^5/(η(σ/2)^2 η(2σ)^2)`; the q-power exponent is 0.

   Hence `θ_2θ_3θ_4 = 2η^3` and `θ_2^4θ_3^2θ_4^2 = 16 η(2σ)^4 η(σ)^4`. At `σ = 2τ` this gives
   `f = θ_2^4θ_3^2θ_4^2/16`.
2. **f in terms of k.** `f(iy/2) = k^2 k' K^4/π^4`, using `θ_2^4 = k^2θ_3^4`, `θ_4^2 = k'θ_3^2` and
   `θ_3^8 = 16K^4/π^4`. The nome is `e^{−πK'/K}`, so `y = K'/K`. Checked numerically at k = 0.4: difference 0.
3. **dy/dk.** `dy/dk = −π/(2kk'^2K^2)` (the Legendre-relation Wronskian). Checked numerically: difference 2e-41.
4. **Fricke relation.** `f(−1/(8τ)) = (4τ/i)^2 (2τ/i)^2 f(τ) = 64τ^4 f(τ)`. Checked numerically to 1e-42.
   * The substitution `t ↦ 1/(8t)` gives `Λ(4−s) = 8^{s−2} Λ(s)`; I redid it by hand.
   * At s = 3 this gives `L(f,3) = (π^2/4) L(f,1)`.
   * The small-t decay `O(t^{−4} e^{−π/(4t)})` is correct.
5. **L-values as integrals.**
   * `L(f,1) = 2π ∫_0^∞ f(it) dt = π ∫_0^∞ f(iy/2) dy` and `L(f,2) = π^2 ∫_0^∞ f(iy/2) y dy`.
   * Orientation: k runs from 1 to 0 as y runs from 0 to ∞, and the minus sign in dy/dk absorbs the reversal.
   * The final substitution is `κ = k'`.
   * Constants: `1/(2π^2)` for L(f,1) and `1/(2π)` for L(f,2), both correct.
6. **Lemma 2.**
   * The series `K'(k) = (π/2) Σ a_m k'^{2m}` is correct.
   * The substitution `x = k^2` gives `(1/2) B(c, m+1)` with `c = (s+1)/2`.
   * The 2F1 at 1 is summed by Gauss's formula, valid for c > 0.
   * `M(2n) = (π^2/4) a_n` (using `Γ(n+1/2) = √π (1/2)_n`).
   * `M'/M = ψ(c) − ψ(c+1/2)`.
   * Differentiation under the integral is justified by domination by `C log^2(1/k) k^{s_0−δ}`.
   * Checked numerically at n = 0 and n = 3.
7. **(Z3).** DLMF 19.12.1 is `K(k) = Σ ((1/2)_m/m!)^2 k'^{2m} (ln(1/k') + d(m))` with
   `d(m) = ψ(1+m) − ψ(1/2+m)`; replacing k by k' gives the file's form. The file's statement and citation are correct.
   * `ψ(n+1/2) = −γ − 2 log 2 + 2H_{2n} − H_n` (checked at n = 0 and n = 1), so `d(n) = −2h_n`.
   * `h_n < 0`, because `H_{2n} − H_n` is an even partial sum of the alternating harmonic series.
   * Checked numerically at k = 0.3: difference 5e-41.
8. **§4.**
   * Tonelli applies because all terms are non-negative.
   * `∫KK' = (π/2) Σ a_n M(2n) = (π^3/8) S_0`.
   * `∫K'^2 = Σ a_n(−(π^2/2) a_n h_n − 2h_n (π^2/4) a_n) = −π^2 S_1`.
   * `F_1(1) = Σ 4(ψ(n+1/2) − ψ(n+1)) a_n^2 = 8 S_1`.
   * Termwise analyticity holds because the terms are uniformly O(n^{-2}).
   * The final ratios work out to −1/16, 1/32 and −3/256 exactly.
9. **U2 Proposition 8.**
   * The reflection identity `Γ(z+1/2)/Γ(z+1) = −tan(πz) Γ(−z)/Γ(1/2−z)` is correct.
   * The `k = −m` part of the bilateral sum re-indexes to `tan^4(πu) (Y/N)(1/2−u)`.
   * With `v = 1/2 − u`, `cot^4(πv) = (πv)^{−4}(1 − (4π^2/3)v^2 + ...)`, so the principal part at u = 1/2 is as
     claimed.
   * The Laurent expansions are correct:
     - `π^2 csc^2(πv) = v^{−2} + π^2/3 + ...`
     - `π^3 csc^2 cot = v^{−3} + 0·v^{−1} + ...`
     - `π^4 csc^4 = v^{−4} + (2π^2/3) v^{−2} + ...`
   * The v^{−2} coefficient of `T_0` matches.
   * `E` is entire, 1-periodic and bounded on a strip, hence constant. Its limits at ±i∞ are ∓iπc, so c = 0 and
     E = 0.
   * `T(0) = τ_0`, which gives `τ_2 = 2π^2 τ_0`.
   * The identification `π^2 F(ε) = Y(ε)/N(ε)`, with `N(0) = π^{−2}`, gives `F_j(1) = τ_j/π^2`.
   * Numerically, `F_3 − (4π^2/3)F_1` and `F_2 − 2π^2 S_0` vanish to 40 digits (table below).
10. **Zudilin's statement and Table 1.** The quotations of (6) and of Remark 2 ("None of the relations in (6) seem
    to be proved") match `p6/z.txt`. The normalisation `Γ(r)Γ(1−r)Γ(t)Γ(1−t) = π^2` is correct, and `z^{n+ε} = 1` at
    z = 1. So the file's F_j(1) are Zudilin's. Table 1, row (1/2,1/2), is `η_2^4 η_4^4`, level 8 = 2^3, LMFDB
    "8.4.1.a". Remark 5 claims a rigorous proof that F_1(1) and F_3(1)/π^2 are proportional. The file reports this
    correctly.

## Numerical table

Script: `Z/audit/audit_numerics.py`, run with mpmath at dps = 40.
* `c_n` come from the eta product to `q^400`, with Hecke checks on `c_15`, `c_9` and `c_25`.
* The L-values come from the incomplete-gamma series split at `t_0 = 1/√8`, and were cross-checked by direct
  quadrature of `f(it)`.
* `F_j` come from term-wise polygamma Taylor coefficients summed with `nsum`, cross-checked by a Cauchy integral on
  `|ε| = 1/8`. That cross-check agrees to about 1e-11, a limit set by the cross-check itself.

Second script: `Z/audit/rwz_7f6.py` checks the BSW and RWZ forms.

| quantity | value |
|---|---|
| L(f,1) | 0.354500683730964718765559891495 |
| L(f,2) | 0.690031163123397525119105420218 |
| L(f,3) | 0.874695377085079044944594728357 |
| S_0 = F_0(1) | 1.11863638716418706834961925753 |
| S_1 | −0.70900136746192943753111978299 |
| F_1(1) | −5.67201093969543550024895826392 |
| F_2(1) | 22.080997219948720803811373447 |
| F_3(1) | −74.6406721779267451686054168198 |
| 7F6(5/4,1/2^6;1/4,1^5;1) | 1.14938972408446010238938375031 |
| W_4(−1) | 0.9027285783238348241903934 |

| check | residual |
|---|---|
| L(f,1) + S_1/2 | −1.7e-41 |
| L(f,2) − π^2 S_0/16 | 1.2e-41 |
| L(f,3) − (π^2/4) L(f,1) | −1.2e-41 |
| F_1 − 8 S_1 | −9e-41 |
| F_2 − 2π^2 S_0 | 4e-40 |
| F_3 − (4π^2/3) F_1 | 0 (at 40 digits) |
| L(f,1)/F_1(1) + 1/16 | 4e-42 |
| L(f,2)/F_2(1) − 1/32 | 0 (at 40 digits) |
| L(f,3)/F_3(1) + 3/256 | 7e-43 |
| ∫K'^2 + π^2 S_1 | −9e-21 (quadrature limited) |
| ∫KK' − π^3 S_0/8 | −3e-22 |
| ∫K^2 − (1/2)∫K'^2 | 5e-21 |
| S_1 + (π^2/16)·7F6 | 0 (at 40 digits) |
| W_4(−1) + (4/π) S_1 | 1e-41 |

All three ratios are confirmed to far more than 15 digits.

## Novelty findings

* **Zudilin, SIGMA 14 (2018) 086**, https://arxiv.org/abs/1805.00544. The relations (6) are stated as
  numerical-conjectural, and Remark 2 says none are proved.
* **Rogers–Wan–Zucker**, arXiv:1303.2259, Ramanujan J. 37 (2015), https://arxiv.org/abs/1303.2259.
  * (30) gives L(f_2,2) in terms of the 4F3.
  * (31) gives `L(f_2,3) = (π^2/4) L(f_2,1) = (1/4)∫K^2 = (π^4/128) 7F6(...)`.
  * The file's account of these is accurate but omits the 7F6 part of (31).
* **J. Wan**, arXiv:1101.1132, Adv. Appl. Math. 48 (2012), https://arxiv.org/abs/1101.1132.
  * (3) is `2∫K^2 = ∫K'^2`, first proved in BSW, his ref. [9].
  * (25) gives `∫x^n K'^2` as a 7F6.
  * Wan remarks that the non-terminating 7F6 transformations write this 7F6 as two 4F3's.
* **Borwein–Straub–Wan**, *Three-step and four-step random walk integrals*, https://carmamaths.org/jon/walks2.pdf
  (preprint 2010), Exp. Math. 22 (2013), https://doi.org/10.1080/10586458.2013.748379. Their formulas for
  `W_4(−1)` give `∫K'^2 = −π^2 S_1` (finding N1). **This is the key identity of the file's §3–§4, already published.**
* **Bönisch–Klemm–Scheidegger–Zagier**, arXiv:2203.09426, https://arxiv.org/abs/2203.09426.
  * (3.11) is stated as "numerically we further find".
  * §4, (4.25)–(4.28), says "Up to a multiplicative constant this confirms two rows from (3.11)".
  * The file's account of BKSZ is accurate.
* **Bönisch–Golyshev–Klemm**, *Fibering out Calabi–Yau motives*, arXiv:2510.03939,
  https://arxiv.org/abs/2510.03939. It proves modularity of mixed periods for the quintic conifold fibre. It does not
  treat the (1/2,1/2) family or Zudilin's (6), and does not cite Zudilin's paper.
* **Papers citing arXiv:1805.00544** (Semantic Scholar list): 2607.25173 (Grove–Long–Rosen–Tu, EHMM III),
  2509.19960 (Au), 2503.05599, 2502.08760, 2512.23549, 2203.09426, 2110.02962, 1906.07384, 1806.05207
  (Osburn–Straub), 1812.07560, 1810.03659, and the Sankhyā B paper on Pearson's walk.
  * None of the abstracts, nor the full texts I read (Au; EHMM III, in a garbled HTML rendering; BGK), states
    (6) for (1/2,1/2) as a theorem.
  * EHMM III, §6 asks "Can we use the EHMM to prove other remaining conjectures in [23, 35]?" This suggests the
    authors regard the relations as open, though I could not confirm what [35] refers to.
* **Bottom line.**
  * The three equalities with Zudilin's normalisation and the constants −1/16, 1/32, −3/256: not found stated as a
    theorem anywhere.
  * `L(f,1) = −S_1/2`, equivalently `Σ_k C(2k,k)^4/256^k (H_{2k} − H_k − log 2) = −2L(f,1)`: follows in two lines
    from BSW 2010/2013 combined with RWZ 2013/2015.
  * `F_1(1) = −(8/π^2)∫K'^2`: is BSW's formula after the trivial step `F_1(1) = 8S_1`.

  So the file is right to call the novelty low, but wrong about what the new part is. Nothing analytic is new; only
  putting the pieces together is.
