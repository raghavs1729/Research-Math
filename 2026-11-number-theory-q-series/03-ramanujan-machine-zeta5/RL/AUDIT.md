# AUDIT of RL/THEOREM.md (hostile referee), 2026-11-26

## Verdict

**Complete after listed fixes.** The proof of Theorem RL is correct: every step checks by hand and to 50+ digits.
The fixes are expository. One is substantive for how the result is presented: the novelty claim. The formula is a
one-line specialization of **Perron's continued fraction for M(a+1;b+1;z)/M(a;b;z)**, a classical textbook result
(§B below). The "still unproven" second set of rows in Table 3 of arXiv:2412.12361 therefore follows directly from
known theory. Credit for the observation should go to that classical fraction. The two misprint claims check out
against the paper (v2, and also the ICLR 2025 proceedings PDF). arXiv has no v3.

## A. Proof check: numbered findings

1. **Lemma 1. Correct.** Index bookkeeping: with `A_j = Γ(n+k+j)/Γ(n+j+3/2)` one gets `J_{n−2} = ΣA_j/(2^j j!)`,
   `J_{n−1} = ΣA_{j+1}/(2^j j!)` and `J_n = ΣA_{j+2}/(2^j j!)`.
   - `(n+k+j)A_j = (n+j+3/2)A_{j+1}` is `Γ(x+1) = xΓ(x)`, and it needs `n+k+j > 0`. For `n ≥ 1` this follows from `k > −1`.
   - Coefficient of `A_J`: `4J(J−1) + 4(n+1)J − 2(n+k) = 4J(n+J+1/2) − 2(n+k+J)`. Both sides equal `4J² + 4nJ − 2n − 2k`.
   - Telescoping: substituting `J → J+1` in the second sum gives exactly the first sum.
   - Splitting is legitimate. `A_J = O(J^{k−3/2})` uniformly, so each of the three series and the two split sums
     converge absolutely.
   - Passage to `m_n`: `−(−1/2)^{n−1} = 2(−1/2)^n` and `(−1/2)^{n−2} = 4(−1/2)^n`. Both sides equal
     `2(n+k)(−1/2)^n J_{n−2}/(n+1)!`.
   - The case `n = 1` uses `J_{−1}` and `m_{−1} = (−1/2)^{−1} J_{−1}/0! = −2J_{−1}`. Both are well defined because
     `Γ(k+1)` is finite for `k > −1`.
   - Numerically, the 3-term residual is ≤ 2e−59 and the m-recurrence residual ≤ 3e−61, for k = 0, 0.3, 5, −0.7.

   Severity: none.

   *Remark (no fix needed):* the relation is DLMF 13.3.2, the b-contiguous relation, applied to
   `M(3/2−k; b; −1/2)` after Kummer's transformation 13.2.39. Indeed `J_n = e^{1/2} Γ(b−a')/Γ(b) · M(a'; b; −1/2)`
   with `a' = 3/2−k` and `b = n+7/2`, and 13.3.2 becomes `N(b+1) + 2(b−3/2)N(b) = 2(b−1−a')N(b−1)`. I re-derived this.

2. **Lemma 2. Correct.** `(n+k+2+i)/(n+7/2+i) = 1 + (k−3/2)/x` with `x ≥ 5/2`. The maximum is `2(k+1)/5` if
   `k ≥ 3/2`, and the ratio is `< 1` otherwise, so `ρ` is right. This gives `M ≤ e^{ρ/2}`, and the Gamma ratio is
   `O(n^{k−3/2})`. Severity: none.
   - *Cosmetic:* "`|m_n| ≤ c_k^n/n!`" needs a constant factor, or `c_k` large enough. All the proof uses is
     `m_n → 0`.

3. **§2. Correct.**
   - The initial values `(P_{−1}, P_0) = (1, 1)` and `(Q_{−1}, Q_0) = (0, 1)` match `C[f] = 1 + K(f_n/1)`, giving
     `P_1/Q_1 = 1+f_1`.
   - The decomposition `m_n = m_{−1}P_n + (m_0−m_{−1})Q_n` was checked at n = −1..11 (residual < 1e−57).
   - `Q_n ≥ 1` holds since `f_n > 0` for `n ≥ 1`. `m_{−1} = −2J_{−1} < 0`, so it is nonzero.
   - The limit is `1 − m_0/m_{−1} = 1 + J_0/(2J_{−1})`.
   - Gamma algebra: `Γ(k+2)Γ(5/2)/(Γ(k+1)Γ(7/2)) = (k+1)/(5/2)`, which gives `(k+1)/5`. Correct.

   Severity: none. *Note:* this direct minimal-solution argument is not optional pedantry. Perron's fraction is known
   to show *anomalous convergence* (convergence to the wrong limit) for some parameters (Gautschi 1977, §B). Here,
   with all `f_n > 0` and `m_n/Q_n → 0`, the correct limit is established.

4. **§3. Correct.**
   - `R_2 = M(1;3/2;1/2) = Σ1/(2j+1)!!` (difference 0 at 60 digits). The ODE argument is right.
   - DLMF 13.3.1 and 13.3.4 are quoted correctly. The residuals are ≤ 5e−60 at three random (a,b,z), including a
     complex `a`.
   - All nine listed values (`μ_2, μ_3, μ_4`, `M(1..3;5/2)`, `M(2..4;7/2)`) agree to 60 digits.
   - The three closed forms k = 0, 1, 2 agree.
   - My own exact rational propagation of the quoted relations (`audit/audit_table.py`) reproduces the k = 3, 4, 5
     entries exactly, already in lowest terms: `(15R+9)/(8R+6)`, `(76R+48)/(39R+27)`, `(955R+615)/(464R+312)`.
     It also gives k = 6, 7, 8: `(7122R+4626)/(3275R+2175)`, `(122633R+79989)/(53544R+35352)`,
     `(1195240R+781320)/(497441R+327453)`. PSLQ on the 60-digit value of C for k = 6 independently returns the
     k = 6 entry.

   **Fix (minor, F1):** the Theorem says "integers given by the recurrences of §3", but §3 never writes the k = 3..5
   derivation and defers to `verify/rl_checks.py`. Add one sentence: "the α..δ for k ≥ 3 are obtained by
   propagating the three displayed relations in exact arithmetic (script)". Alternatively, display
   the next values, for example `M(4;5/2;1/2) = 3(μ_4−μ_3) = (4R_2+3)/4`.

5. **Remark 1 and Remark 2. Correct.**
   - Remark 1: `C[1/(2n)] − 1 − 1/(2C[1/(2(n+1))]) = 4e−62`.
   - Remark 2: k = 3/2 gives exactly 3/2, and k = 5/2 gives 5/3 (residual 2e−61). By hand via Kummer:
     `M(9/2;7/2;1/2)/M(7/2;5/2;1/2) = (8/7)/(6/5)`, so the value is `1 + (7/10)(20/21) = 5/3`.
   - "k = 3/2, 5/2, …" correctly excludes k = 1/2.
   - The heuristic integral in §1 is the Euler integral, which is valid for k < 3/2. Correct as stated.

6. **Novelty/standing statement (F2, required).** §0 "Standing: low … cf. Lorentzen–Waadeland, Jones–Thron" is too
   vague, and "New means only not found" understates what is known.
   - Replace it with the precise specialization in §B and state: *"Theorem RL is the case a = k+1, b = 5/2,
     z = 1/2 of Perron's continued fraction for M(a+1;b+1;z)/M(a;b;z), after the equivalence transformation
     c_n = n+1. The contribution is only the identification of the Table-3 rows with this classical fraction (and
     the misprints)."*
   - Also cite Gautschi 1977 for the convergence subtlety, as in finding 3.

7. **Misprint claims (verified against the paper).**
   - raw.txt l.562–569 (arXiv v2) and the ICLR 2025 proceedings PDF (l.609–616 of its text) both print
     `C[n/(2n(n+1))] = 1/(2R2 − 1)` and `C[(1−2n)/(4(n+2)(n+3))] = (R2+1)/(6R2−1)`.
   - The true values are 1.21747475218590112715866566091… = `1/(2R_2−2)` (proved), and
     0.978316501457267418105777107274… = `(R_2+1)/(6R_2−6)`. The second is found by PSLQ as `[6,−6,−1,−1]` and
     agrees to 60 digits.
   - The printed values 0.5490… and 0.3230… are off by 0.668 and 0.655.
   - No reading under shifted indexing (`f_{n+1}`, `f_{n+2}`, `f_{n+3}`) or of `C−1` rescues them.
   - The other four printed rows are correct (difference 0 at 30 digits).
   - The THEOREM correctly labels the second misprint as numerical only.

   **Fix (minor, F3):** say "C-transform of positive terms is > 1" applies to the first misprint only. Also note
   that the same misprints are in the ICLR camera-ready version.

8. **arXiv versions.** The arXiv abstract page (fetched 2026-11-26) lists only v1 (16 Dec 2024) and v2
   (19 Jan 2025). There is no v3, and the table is not corrected anywhere online that I found.

## Numerical table (own script `audit/audit_rl.py`, mp.dps = 60, depths 400/800)

| k | C[(n+k)/(2n(n+1))] (35 digits) | abs(C − formula) |
|---|---|---|
| 0 | 1.217474752185901127158665660910746 | 0 |
| 1 | 1.410686134642447997690824711419115 | 0 |
| 2 | 1.585180341136230246073095737732012 | 0 |
| 3 | 1.7448330148508397467109997758283831 | 0 |
| 4 | 1.8924442717581361403831743199063716 | 0 |
| 5 | 2.0300949499216048422719028461431771 | 0 |
| 6 | 2.1593690850077228521136509459212261 | 0 |
| 10 | 2.6138514283642676918396319855303131 | 0 |
| 0.3 | 1.2776927060659756432979876361320771 | 0 |
| 1.5 | 1.5 | 0 |
| 2.5 | 1.6666666666666666666666666666666667 | 1.6e−61 |
| 7/3 | 1.6398907131659374423272480721892051 | 1.6e−61 |
| −0.5 | 1.112206218978170233581918770692011 | 0 |
| −0.9 | 1.0230467015853186655422690275901949 | 0 |

The closed forms for k = 0..5 agree with these values to ≤ 2e−61.

*(Methodology trap, confirmed: building `f_n` as `(n+k)/(2*n*(n+1))` with Python ints silently gives float
precision, about 1e−17. My first run showed exactly this, and the fix was `mpf(n+k)`. A grep of rl_checks.py shows its row
lambdas use `mpf(...)`. My table above uses mpf throughout.)*

## B. Prior art

**B1. The general formula is a direct specialization of a textbook formula. Novelty: nil.**

Perron's continued fraction (O. Perron, *Die Lehre von den Kettenbrüchen*, 3rd ed., Bd. II, Teubner 1957, in the
chapter on confluent hypergeometric series; analysed in W. Gautschi, "Anomalous convergence of a continued
fraction for ratios of Kummer functions", *Math. Comp.* 31 (1977), no. 140, 994–999,
https://www.ams.org/mcom/1977-31-140/S0025-5718-1977-0442204-3/; also in Cuyt–Petersen–Verdonk–Waadeland–Jones,
*Handbook of Continued Fractions for Special Functions*, Springer 2008, ch. 16, Kummer functions) reads:

    M(a+1; b+1; z) / M(a; b; z) = b / ( b − z + (a+1)z / ( b+1−z + (a+2)z / ( b+2−z + ... ))).

I verified this to 50 digits at three parameter sets (`audit/audit_priorart.py` (c)). The exact equation numbers in
Perron and in Cuyt et al. were not checked against a physical copy; the Gautschi paper is the verified pointer to
"Perron's continued fraction for ratios of Kummer functions".

Put `a = k+1`, `b = 5/2`, `z = 1/2`. Then `b−z+n = n+2` and `(a+n+1)z = (n+k+2)/2`, so

    M(k+2;7/2;1/2)/M(k+1;5/2;1/2) = (5/2) / (2 + ((k+2)/2)/(3 + ((k+3)/2)/(4 + ...))),
    1 + ((k+1)/5)·[that] = 1 + ((k+1)/2)/(2 + ((k+2)/2)/(3 + ((k+3)/2)/(4 + ...))) = 1 + K_{n≥1}( ((n+k)/2) / (n+1) ),

and the equivalence transformation with `c_n = n+1` turns this into `1 + K((n+k)/(2n(n+1)) / 1) = C[f_n]`. The
convergents of `C[(n+k)/(2n(n+1))]` and of `1 + K(((n+k)/2)/(n+1))` are identical as exact rationals (checked up
to n = 40 for k = 0, 1, 3, 1/3). So Theorem RL, *including the explicit 1F1 form*, is a substitution into Perron's
fraction. The rational/R_2 values for integer k then follow from DLMF 13.3.1 and 13.3.4, which is standard.

The only non-trivial point in the paper's framing is the limit. Gautschi shows that Perron's fraction can converge
to the wrong value for some parameters, so one needs either Perron's/Gautschi's convergence conditions or the
THEOREM's minimal-solution argument. The latter is valid here.

Equivalently, via Kummer's transformation, this is the b-contiguous (DLMF 13.3.2) fraction for
`M(a';b+1;−1/2)/M(a';b;−1/2)` with `a' = 3/2−k` (finding 1).

**B2. Is `C[1/(2n)] = Σ 1/(2j+1)!! = R_2` classical?**

- Yes, as the k = 1 case above: `C[1/(2n)] = 1 + 1/(2 + (3/2)/(3 + 2/(4 + (5/2)/(5 + ...))))`, which is Perron's
  fraction for `M(2;7/2;1/2)/M(1;5/2;1/2)`.
- It is *not* simply an equivalent form of the classical erf fraction
  `R_2 = e^{1/2}∫_0^1 e^{−t²/2}dt = 1/(1 − 1/(3 + 2/(5 − 3/(7 + ...))))`. That fraction is the z = 1/√2 case of
  `erf z = (2z/√π)e^{−z²}/(1 − 2z²/(3 + 4z²/(5 − 6z²/(7 + ...))))` (Cuyt et al. ch. 13; Jones–Thron; I verified it
  to 5e−51 at depth 30).
- The two fractions' convergents differ, sharing only 1 and 3/2: `1, 3/2, 7/5, 24/17, 213/151, …` versus
  `1, 3/2, 17/12, 55/39, 529/375, …`. One is the Gauss-type fraction (DLMF 13.5.1, Jones–Thron Thm 6.3, partial
  numerators of degree −2 after normalization, alternating sign). The other is the Perron/b-contiguous type
  (degree −1, positive).
- Ramanujan's own `√(πe/2) = Σ1/(2j+1)!! + 1/(1+1/(1+2/(1+…)))` (Berndt, *Ramanujan's Notebooks* Part II, ch. 12)
  concerns `C[n]`, the erfc part, not `C[1/(2n)]`.

**B3. Later proofs or corrections of Table 3.**

- No v3. The ICLR 2025 proceedings PDF
  (https://proceedings.iclr.cc/paper_files/paper/2025/file/6782c18960808848174cfe60742b415a-Paper-Conference.pdf)
  still says "still unproven" and has the same two misprints.
- Semantic Scholar shows 0 citations. OpenAlex shows 2 citing works:
  - "From Euler to AI: Unifying Formulas for Mathematical Constants" (arXiv:2502.17533). It harvests only the π
    formulas of 2412.12361 (its Table 15) and does not treat the √(πe) rows; I grepped it for erf, πe and the
    authors.
  - The SCML-2026 extended abstracts (RISC, doi:10.35011/risc-proceedings-scml.3). I could not retrieve the text
    (HTML landing page only), so this is not excluded.
- Also not related: arXiv:2609.33529 (Bhatnagar, on Ramanujan's elementary continued fractions, ch. 12) does not
  mention these rows.
- Conclusion: I found no explicit published proof of the "latter 4 rows". They are, however, an immediate
  specialization of Perron's classical fraction, so the honest claim is "observed to be a special case of a
  classical formula", not a new theorem.

## Required fixes (summary)

- **F2 (required):** replace the novelty/standing paragraph with the Perron specialization of §B1. Cite Perron
  Bd. II, Gautschi 1977 and Cuyt et al. ch. 16, and state that novelty is limited to the identification and the
  misprints.
- **F1 (minor):** say explicitly how the k = 3..5 integers are obtained, or display the intermediate values.
- **F3 (minor):** note that the misprints persist in the ICLR camera-ready version, and that there is no v3 as of
  2026-11-26.
- **Optional:** prove the first-set misprint `(R_2+1)/(6R_2−6)`. It should follow the same way, via the Gauss-type
  fraction (DLMF 13.5.1) or tail relations to the proved first row. As it stands it is numerical (PSLQ, 60 digits).

## Files

- `audit/audit_rl.py`: theorem formula, closed forms, §3 values, DLMF 13.3.1 and 13.3.4, Lemma 1 and §2 checks,
  remarks, the Table 3 rows, the misprints.
- `audit/audit_table.py`: exact integer table for k = 0..8 from the quoted contiguous relations.
- `audit/audit_priorart.py`: equivalence to `1+K(((n+k)/2)/(n+1))`, Perron's fraction, the erf fraction and its
  convergents.
