# Akbary–Francis: the Nicolas-type inequality in the progression 1 mod 13 fails infinitely often

Date: 2026-11-24. Working directory: `monotonicity/akbary-francis/`. Nothing here was committed by this
agent. Brief: `monotonicity/RAMANUJAN-TRIAGE.md`, Topic 3, Work order 1.

Source read in full for the parts used (PDF fetched, `source/af.pdf`, pages rendered and read as images because
`pdftotext` loses overlines and fractions): A. Akbary, F. J. Francis, "Euler's function on products of primes
in progressions", arXiv:1810.11524v1 (26 Oct 2018, the only arXiv version). Sections 1, 3, 4, 5, 6, 7, 8.
**All statement and equation numbers below refer to arXiv v1.**
Journal version: A. Akbary, F. J. Francis, "Euler's function on products of primes in a fixed arithmetic
progression", Math. Comp. 89 (2020), no. 322, 993–1026, DOI 10.1090/mcom/3463 (bibliographic data confirmed from
Crossref). **Its full text was not accessible to us** (the AMS site returned HTTP 403 to this agent and to the
auditor); only its abstract was seen, and it states the same results as the arXiv abstract. Whether the journal
text differs from v1 in numbering or content is therefore unchecked.

Status: independently audited, `AUDIT.md` (no rigor break; two wording gaps and seven minor points, all applied
here; see §10, Review history).

## 0. Verdict

**Theorem A (GRH_13).** Assume the Riemann Hypothesis and GRH for the eleven non-principal Dirichlet
L-functions mod 13 (= "GRH_13" of the paper = GRH for the Dedekind zeta function of Q(ζ_13)). Let N̄_k be the
product of the first k primes ≡ 1 (mod 13). Then

    N̄_k / ( φ(N̄_k) · (log(12 log N̄_k))^{1/12} )  <  1 / C(13, 1)

for infinitely many k, i.e. inequality (3) of Akbary–Francis with q = 13, a = 1 **fails for infinitely many k**.
Quantitatively, with f as in their (4),

    lim sup_{x→∞} 12 √x log x · log f(x; 13, 1)  ≥  0.2716 ,
    lim sup_{x→∞} Σ_χ Σ_{ρ∈Z(χ')} x^{i Im ρ} / (ρ(ρ−1))  ≥  4.2716  >  4 = 2 R_{13,1}     (their sum (7), a = 1).

So condition (8) of their Corollary 1.19 fails for (q, a) = (13, 1): there is no x₀ beyond which
log f(x; 13, 1) < 0, and their Theorem 1.22 (GRH_q ⇔ (3) for all k, proved for q ≤ 10 and q = 12, 14) does not
extend to q = 13. The same holds for every square class a mod 13 (Theorem A′, §5) and, trivially, for q = 26
(same primes, same φ, same constant).

**Corollary B (unconditional, using Akbary–Francis Theorem 1.16(a)).** Inequality (3) with q = 13, a = 1 fails
for infinitely many k. (If GRH_13 holds this is Theorem A. If GRH_13 fails, their Theorem 1.16(a) gives
log f(x; 13, 1) = Ω_±(x^{−b}) for some b < 1/2; its hypothesis "no singularity of 𝓛(s; 13, 1) on (1−b, 1)" holds
because no L(s, χ), χ mod 13, has a real zero in (0, 1): certified here, §3.1, independently by the auditor, and
also Platt via their Cor. 1.13.) So, unconditionally, (3) fails for infinitely many k, and no criterion of the
form "GRH_13 ⇔ (3) holds for all k ≥ k₀" can hold unless GRH_13 is false.

What is **proved** (given the trusted inputs of §8): Theorem A, A′, Corollary B.
What is **computed** in support: certified zeros, the Arb certificate at one point u₀, a sieve check to 10^10.
What is **not** obtained: an explicit k. The certified point is u₀ = log x₀ = 18 133 579.0743, x₀ ≈ 10^{7 875 313};
for p̄_k ≤ 10^10 (k ≤ 37 919 477) inequality (3) holds without exception (§4).
**q = 11: not settled** (§6): best certified value (under GRH_11) 2.68 < 4.

By-products (§7), not in the brief:
* For most of the next prime moduli (3) fails at tiny, explicit k with no hypothesis at all, e.g.
  **q = 17, a = 1, k = 5** (N̄_5 = 103·137·239·307·409); also q = 23, 29, 37, 41, 43, 47, 59, 61, 67, 71, 73, 79, 83.
  No failure in sieve range for q = 11, 13, 19 (≤ 10^10) and q = 31, 53, 89 (≤ 10^9).
  **This is not a counterexample to anything in the paper**: its Theorem 1.22 is stated only for q ≤ 10 and
  q = 12, 14. It shows that the literal q = 17 analogue of that criterion has a false right-hand side (§7.1).
* q = 17 and q = 19: infinitely many failures too, under GRH_q with certified margins 3.9 and 4.9, and hence
  unconditionally by the argument of Corollary B (§7.2).
* The set {q ≤ 1000 : F_q < 2R_{q,1}} is {1, …, 10, 12, 14, 15, 16, 18, 20, 21, 24, 28, 30, 36, 40, 42, 48, 60}: thirteen
  moduli beyond the paper's list satisfy the hypothesis of its Proposition 1.20.

## 1. Conventions, re-derived from the source

All formulas below were checked against the PDF (page images), not against the triage paraphrase.

* S_{q,a}: integers all of whose prime factors are ≡ a (mod q). N̄_k = N_{q,a}(k) = p̄_1 ⋯ p̄_k, p̄_i the i-th prime
  ≡ a (mod q). θ(x; q, a) = Σ_{p ≤ x, p ≡ a} log p.
* Theorem 1.8: ∏_{p ≤ x, p ≡ a (q)} (1 − 1/p) ~ C(q, a) / (log x)^{1/φ(q)}, with
  C(q, a)^{φ(q)} = e^{−γ} ∏_p (1 − 1/p)^{α(p; q, a)}, α = φ(q) − 1 if p ≡ a (mod q), α = −1 otherwise.
* Inequality (3): n / ( φ(n) (log(φ(q) log n))^{1/φ(q)} ) > 1 / C(q, a).
* f(x; q, a) := (log(φ(q) θ(x; q, a)))^{1/φ(q)} / C(q, a) · ∏_{p ≤ x, p ≡ a} (1 − 1/p)   (p. 3–4). For
  x ∈ [p̄_k, p̄_{k+1}): f(x; q, a) = (log(φ(q) log N̄_k))^{1/φ(q)} φ(N̄_k) / (C(q, a) N̄_k), so
  **(3) holds for N̄_k ⇔ f(x; q, a) < 1 ⇔ log f(x; q, a) < 0** for one (equivalently every) x in [p̄_k, p̄_{k+1}).
  Formula (4): log f = Σ_{p ≤ x} log(1 − 1/p) + (1/φ(q)) log log(φ(q) θ(x; q, a)) − log C(q, a).
* 𝓛(s; q, a) = Σ_χ χ̄(a) L′/L(s, χ). SH_{q,a} (Hypothesis 1.12): its singularities in 0 < Re s < 1 lie on
  Re s = 1/2. GRH_q ⇒ SH_{q,a}; SH_{q,1} ⇔ GRH_q.
* Ind_q(a) = least n > 1 such that a is an n-th power mod q (= 2 iff a is a square); R_{q,a} = #{b mod q :
  b^{Ind_q(a)} ≡ a}. For q = 13 and a square a: m = 2, **R_{13,a} = 2**, threshold 2 R = 4.
* Z(χ) = {ρ : L(ρ, χ) = 0, Re ρ ≥ 0, ρ ≠ 0}; χ′ = the primitive character inducing χ. For q = 13 the eleven
  non-principal characters are primitive; the principal one is induced by the trivial character, so Z(χ₀′) is the
  set of non-trivial zeros of ζ (this is the role of the zeta zeros). The zero at s = 0 of an even L(s, χ) and
  the zeros of 1 − 13^{−s} on Re s = 0 do not enter the sum over Z(χ′).
* **Theorem 1.18** (SH_{q,a}; m = Ind_q(a)), formula (6):

      log f(x; q, a) = (1/φ(q)) [ (1/(√x log x)) Σ_{χ mod q} χ̄(a) Σ_{ρ ∈ Z(χ′)} x^{i Im ρ} / (ρ(ρ − 1))
                                 − m R_{q,a} / ((m − 1) x^{(m−1)/m} log x) ]  +  O( 1/(√x log² x) ).

  The coefficient is the **conjugate** χ̄(a) (the overline is lost by `pdftotext`), and the denominator is
  ρ(ρ − 1), not ρ(1 − ρ).
* **Corollary 1.19** (a a square, SH_{q,a}): there is x₀ with log f(x; q, a) < 0 for x > x₀ if and only if
  lim sup_{x→∞} Σ_χ χ̄(a) Σ_{ρ∈Z(χ′)} x^{i Im ρ}/(ρ(ρ−1)) < 2 R_{q,a}   (8).
  Remark (auditor): from Theorem 1.18 one gets "(8) ⇒ x₀ exists" and "lim sup > 2R ⇒ no x₀"; in the boundary case
  lim sup = 2R exactly the O(1/log x) term decides, so the "only if" direction of the corollary is not justified
  there. This is irrelevant below: only "lim sup > 2R ⇒ log f > 0 for arbitrarily large x" is used, with a strict
  margin, and it is proved directly from Theorem 1.18 in §2.
* F_q := Σ_χ Σ_{ρ∈Z(χ′)} 1/(ρ(1 − ρ))   (9); Table 1: F_11 = 4.26098, F_13 = 6.45484.

**Signs on the critical line.** For ρ = 1/2 + iγ: ρ(ρ − 1) = −(1/4 + γ²). Put, for u = log x,

    S_a(u) := Σ_{χ mod q} χ̄(a) Σ_{ρ = 1/2 + iγ ∈ Z(χ′)} e^{iγu} / (1/4 + γ²).

Under GRH_q the sum (7) equals −S_a(log x), and for m = 2 Theorem 1.18 reads

    φ(q) √x log x · log f(x; q, a) = −S_a(log x) − 2 R_{q,a} + O(1/log x).                                  (★)

So **(3) fails near x when S_a(log x) < −2R_{q,a}**, and |S_a| ≤ F_q. The zeros of L(s, χ̄) are the complex
conjugates of the zeros of L(s, χ). Hence S_a is real: a zero 1/2 + iγ of L(s, χ) and the zero 1/2 − iγ of
L(s, χ̄) contribute together 2 cos(γu − arg χ(a)) / (1/4 + γ²) (for real χ and for ζ: the two zeros ±γ of the
same function, same formula). For a = 1: S_1(u) = Σ cos(γu)/(1/4 + γ²) over all zeros of ζ_{Q(ζ_q)}.

Sanity of the sign: for q = 1 this is Nicolas's formula, log f = −(2 + Σ_γ x^{iγ}/(1/4+γ²))/(√x log x) + …,
negative because F_1 = 0.046 < 2. I also re-derived (★) independently of the statement of Theorem 1.18 from
its ingredients: Prop. 3.1 (log f = K + O(1/x), K = ∫_x^∞ (θ(t;q,a) − t/φ(q)) g(t) dt, g = (1 + log t)/(t² log² t));
ψ − θ = (R_{q,a}/φ(q)) √t + O(√t / log t) for m = 2 (Prop. 3.3); F_s(x) := ∫_x^∞ t^s g(t) dt
= x^{s−1}/((1−s) log x) + O(…/log² x) (Lemma 3.4, Nicolas); the explicit formula ψ(t; q, a) − t/φ(q)
= −(1/φ(q)) Σ_χ χ̄(a) Σ_ρ t^ρ/ρ + … (their Thm 6.1 is McCurley's integrated form). These give
J = ∫_x^∞ (ψ − t/φ) g = (1/(φ(q) √x log x)) Σ χ̄(a) Σ x^{iγ}/(ρ(ρ−1)) + … and the term −(R/φ)·2/(√x log x): (★).
The numerical test of §4 confirms sign, normalisation and the constant 4 against actual primes.

Minor erratum noticed in the source (irrelevant for what follows): the captions of Figures 2 and 3 quote
"C(5,1)^{−1} ≈ 1.2252" and "C(5,3)^{−1} ≈ 0.8060"; from Theorem 1.8 (and from primes up to 10^7,
`sieve.c`) one gets C(5,1) = 1.225238…, C(5,3) = 0.805951…, so those numbers are C, not C^{−1}.

## 2. Theorem A and its proof

Notation: q = 13, a = 1, φ(q) = 12, R = 2, S = S_1.

**Lemma 1 (closed form of F_q; Montgomery–Vaughan Cor. 10.14, 10.18 = the paper's (49), (50)).**
Σ_ρ 1/(ρ(1−ρ)) = 2 + γ − log 4π over the non-trivial zeros of ζ, and for a primitive χ mod q, of parity
κ ∈ {0, 1}: Σ_{ρ∈Z(χ)} 1/(ρ(1−ρ)) = log(q/π) + 2 Re (L′/L)(1, χ) + ψ((1+κ)/2), ψ = Γ′/Γ
(ψ(1) = −γ, ψ(1/2) = −γ − 2 log 2). The series converge absolutely. Evaluated in Arb (`closed_form.py`;
L(1, χ) and L′(1, χ) from the Stieltjes constants of the Hurwitz zeta function):

    F_13 = 6.454845631181664360577825…        F_11 = 4.260981950643087115571535…

(the paper: 6.45484, 4.26098; also F_3, F_5, F_7, F_8, F_9, F_12 of Table 1 are reproduced, `fq_table_1000.txt`).
Under GRH_13 every term is 1/(1/4 + γ²) > 0 and they sum to exactly F_13.

**Lemma 2 (certificate; `certify.py`, output `certificate_q13_a1.json`).** Assume GRH_13 and let
u₀ = 18133579.0743. Then S(u₀) ≤ −4.27166084 < −4.
Proof. The file `zeros_q13.txt` (sha256 `98115fad53b7d7f88e62129fd9b38f6bc97635ad01343ad28f17479b4503d88f`,
pinned in `SHA256SUMS`)
lists 4208 windows (d − 10^{−30}, d + 10^{−30}), each labelled by a function (ζ or χ_13(n, ·)) and certified to
contain a zero 1/2 + iγ, γ > 0, of that function; windows with the same label are disjoint (§3.1). Each window
gives two distinct zeros (the listed one and the conjugate zero of the conjugate function), with joint weight
2w, w = 1/(1/4 + γ²), and joint contribution 2w cos(γu₀) to S(u₀). Ball arithmetic (256 bits) over the windows:

    S_T(u₀) := Σ_windows 2w cos(γu₀) = −4.34330566187385 ± 10^{−22},
    W       := Σ_windows 2w           =  6.38320081487161 ± 10^{−24}.

All remaining zeros have positive weights summing to exactly F_13 − W = 0.07164481631005 (Lemma 1, GRH_13), and
|cos| ≤ 1. Hence S(u₀) ≤ S_T(u₀) + (F_13 − W) = −4.27166084556380 ± 10^{−22}. ∎
(Only "each window contains at least one zero" is used. The completeness of the list below height 400, which is
also certified, is not needed for this inequality. The ± 10^{−22} is the ball radius forced by the deliberately
loose windows of half-width 10^{−30} times u₀; the listed decimals are in fact within about 10^{−42} of the zeros,
by the auditor's Krawczyk re-certification, `audit/logs/recert_q13.log`. An 80-digit mpmath recomputation by
the auditor gives S_T(u₀) = −4.34330566187385447981405661174 and the same bound.)

**Lemma 3 (recurrence; this is the almost periodicity of S, proved directly).** For every ε > 0 and every U there
is τ > U with sup_{u∈ℝ} |S(u + τ) − S(u)| ≤ ε.
Proof. Σ_ρ 1/(1/4+γ²) = F_13 < ∞, so there is a finite set of zeros, with ordinates γ_1, …, γ_N, outside of
which the total weight is < ε/4. Dirichlet's theorem on simultaneous approximation (the tool of the paper's
Remark 1.23(i)): for reals α_1, …, α_N and every integer Q ≥ 1 there is an integer 1 ≤ n ≤ Q^N with
‖nα_j‖ < 1/Q for all j (pigeonhole: two of the Q^N + 1 points ({mα_1}, …, {mα_N}), 0 ≤ m ≤ Q^N, lie in the same
box of side 1/Q; n = the difference of the two m). Take t₀ > U, α_j = γ_j t₀/2π, τ = n t₀ ≥ t₀ > U. Then
|e^{iγ_j τ} − 1| ≤ 2π‖γ_j τ/2π‖ ≤ 2π/Q, so
|S(u + τ) − S(u)| ≤ (2π/Q) F_13 + 2 · ε/4 ≤ ε once Q ≥ 4π F_13/ε. ∎

**Proof of Theorem A.** GRH_13 implies SH_{13,1}, so Theorem 1.18 holds: there are C₀, X₀ with
| 12 √x log x · log f(x; 13, 1) + S(log x) + 4 | ≤ C₀ / log x for x ≥ X₀ (this is (★); the O-constant depends only
on q, a). Put δ = 0.2716 < −4 − S(u₀) (Lemma 2) and fix 0 < ε < δ. For every U, Lemma 3 gives u = u₀ + τ > U with
S(u) ≤ S(u₀) + ε, hence, for x = e^u,

    12 √x log x · log f(x; 13, 1) ≥ −S(u₀) − ε − 4 − C₀/u ≥ δ − ε − C₀/u,

which is positive as soon as u > C₀/(δ − ε). So log f(x; 13, 1) > 0 for arbitrarily large x, and
lim sup 12 √x log x · log f ≥ δ (let ε → 0). For such x let k = π(x; 13, 1) ≥ 1 (x ≥ 53), so x ∈ [p̄_k, p̄_{k+1}) and
f(x; 13, 1) > 1 says exactly that N̄_k/(φ(N̄_k)(log(12 log N̄_k))^{1/12}) < 1/C(13, 1). Since x → ∞ along the
sequence, k = π(x; 13, 1) → ∞, so infinitely many distinct k occur. ∎

Remarks.
1. *Relation between x and k.* No choice of x is needed: f(·; q, a) is constant on [p̄_k, p̄_{k+1}), so failure
   at any x ≥ p̄_1 is failure at k = π(x; q, a).
2. *The error term.* Only its form O(1/(√x log² x)), uniform in x for fixed q, a, is used. Because the implied
   constant is not explicit in the paper, Theorem A does **not** assert that (3) fails exactly at x₀ = e^{u₀}, only
   at x = e^{u₀+τ} for all sufficiently large almost periods τ. (At u₀ ≈ 1.8·10^7 the term C₀/u₀ is below the margin
   unless C₀ > 4.9·10^6, so failure at π(x₀; 13, 1) itself is overwhelmingly likely; proving it needs an explicit
   Theorem 1.18, i.e. explicit GRH-conditional bounds for θ(x; 13, b). Not done.)
3. *Density.* By Bohr's theorem the ε-almost periods of S are relatively dense, so the set of u = log x with
   log f(e^u; 13, 1) > 0 is relatively dense in [U₁, ∞) for some U₁: failures recur with bounded gaps in log x.
4. Theorem A gives lim sup of (7) ≥ 4.2716, i.e. ≥ 0.6618 F_13. It does **not** answer Question 1.24 (is the
   lim sup equal to F_q?).

**Proof of Corollary B.** If GRH_13 holds, Theorem A. If not, SH_{13,1} is false (SH_{q,1} ⇔ GRH_q), so
Θ > 1/2, and for 1 − Θ < b < 1/2 their Theorem 1.16(a) gives log f(x; 13, 1) = Ω_+(x^{−b}) provided 𝓛(s; 13, 1) has
no singularity on the real segment (1 − b, 1). Its singularities in (0, 1) are real zeros of ζ (none) or of
L(s, χ), χ mod 13 non-principal: none. Three sources for this real-zero input: (i) the certified count [N] of
§3.1 (no zero off the critical line and none at s = 1/2 below height 400); (ii) the auditor's independent Arb
check that L(σ, χ) ≠ 0 for σ ∈ [1/2, 1] for all 11 characters (`audit/real_zeros.py`,
`audit/logs/real_zeros_q13.log`; (0, 1/2) follows by the functional equation); (iii) Platt, via their
Corollary 1.13. ∎
I read the proof of Theorem 1.16(a) (their §4–5: Landau's oscillation theorem applied to ∫ J(x) x^{−s} dx) and
found no problem, but did not re-derive it line by line; the auditor checked it line by line (AUDIT §5) and found
it sound. Corollary B depends on it. "Unconditional" means: no unproved hypothesis; it does rest on that
published theorem.

## 3. Computations for Lemma 2

### 3.1 Certified zeros (`lzeros_cert.py`, `zeta_zeros_cert.py`, `assemble.py`, `verify_zeros.py`)

Adapted from `monotonicity/nonmonotone/zeros_cert.py` (odd real characters) to all primitive characters,
complex ones included. Characters are Conrey characters χ_13(n, ·) as in Arb/python-flint 0.9.0;
χ_13(n^{−1} mod 13, ·) is the conjugate (checked on all residues).

* **[Z] localisation.** For each character the Hardy function Z(t, χ) (real for real t, |Z| = |L(1/2 + it, χ)|)
  has certified opposite signs at d ± 10^{−30} (Arb, 224 bits). Positive ordinates only; the zeros of L(s, χ) with
  negative ordinate are the conjugates of the zeros of L(s, χ̄) with positive ordinate. For non-real χ the zero
  set is not symmetric: e.g. χ_13(6, ·) has its first zero at γ = 0.88396…, and L(1/2 − 0.88396 i, χ_13(6, ·)) ≠ 0,
  while χ_13(11, ·) = conj has the zero at −0.88396 i (checked with Arb and independently with mpmath).
* **[N] completeness (not needed for Lemma 2, needed for Corollary B).** For each conjugate pair {χ, χ̄} at a
  common height T: #{zeros of Λ(s, χ) in (−1, 2) × (−T, T)} = (D(χ) + D(χ̄))/π, D = change of arg Λ along
  2 → 2 + iT → 1/2 + iT, using Λ(1 − s̄, χ) = ε(χ) conj Λ(s, χ) (left half = right half) and
  conj Λ(s̄, χ) = Λ(s, χ̄) (lower half for χ = upper half for χ̄). D is enclosed with ball arithmetic (Γ-factor by
  lgamma; Re L > 0 on Re s = 2; the horizontal segment covered by balls on which L ≠ 0). The count equals the
  number of windows of χ and χ̄ at every block top (20 to 28 block tops per pair up to T = 400). Hence exactly one
  simple zero per window, no zero off the line, none at s = 1/2, below T = 400. No grid refinement was ever needed.
* ζ: windows from `acb.zeta_zeros`, each re-certified by the sign change of Z; count N(400) = 202 from
  `arb.zeta_nzeros`.
* Result: 202 (ζ) + 4006 (characters; 364 or 365 per character) = 4208 windows, T = 400.
  `verify_zeros.py 13 --count` re-certifies [Z] for all windows and [N] at T = 400 from the assembled file alone:
  PASS (`logs/verify_q13.log`).
* **Pinning** (audit finding G1 of `nonmonotone-audit/AUDIT.md`): `SHA256SUMS` holds
  `98115fad53b7d7f88e62129fd9b38f6bc97635ad01343ad28f17479b4503d88f  zeros_q13.txt` and
  `eea2f78092d7a0e9e5503c43722d7a77b4849ad58f50d8ab8657a33afb1d93c7  zeros_q11.txt` (and q = 17, 19, §7.2). How a certificate is bound to its zero file: `afcommon.load_zeros(q)` is the only
  reader of zero files; it computes the sha256 of the file on disk and raises unless it equals the value pinned
  in `SHA256SUMS`; `certify.py` then writes that sha256 into the certificate as `"zeros_sha256"`. So a
  certificate cannot be produced from a stale, truncated or edited zero file, and each JSON names the exact file
  it was computed from. `check_bindings.py` re-checks, for all 11 certificate files, that `zeros_sha256` equals
  the pinned hash and the hash of the file now on disk: PASS, and a modified copy of `zeros_q13.txt` is refused
  (`logs/check_bindings.log`). The auditor re-certified the files with exactly these hashes.

### 3.2 F_13 two ways (`f_two_ways.py`)

    closed form                         F_13 = 6.45484563118166
    certified partial sum (|γ| < 400)   W    = 6.38320081487161
    exact tail  F_13 − W                     = 0.07164481631
    density estimate of the tail             = 0.071665   ( (1/πT)[11(log(13T/2π)+1) + (log(T/2π)+1)], T = 400 )

The same agreement holds character by character (each of the eleven tails is 0.00614, the ζ tail 0.00410), which
checks formula (50), the parity bookkeeping, and the zero lists against each other. An mpmath recomputation
(`crosscheck_mpmath.py`, no Arb in the numerics; `logs/crosscheck_mpmath_q13.log`) gives F_13 = 6.4548456311816643606
and S_T(u₀) = −4.343305661873854479814, the same digits as Arb.

### 3.3 Weight distribution and the search (`scan.c`, exploratory, double precision)

2w = 2/(1/4+γ²) for the lowest ordinates of the union: γ = 0.88396 (χ_13(6,·), order 12, odd): **1.939**
(30 % of F_13); then 2.19555: 0.394; 2.27313: 0.369; 2.34547: 0.348; 3.11934 (the quadratic character): 0.200.
Cumulative: 10 zeros 3.92, 20 zeros 4.55, 40 zeros 5.08, 100 zeros 5.61, 788 zeros (γ < 100) 6.22. Required:
−S ≥ 4, i.e. 62 % of F_13.

* u ≈ π/γ_1 alone is not enough: S(3.97) ≈ −2.92 (x = e^{3.97} = 53, the first prime ≡ 1 mod 13).
* A scan over the minima of the dominant cosine (stage 1: 100 heaviest zeros; stage 2: all 4208) reaches −4 quickly.
  No lattice reduction was needed. Successive record lows of S_T (selection): −3.39 (u = 81.5), −3.61 (3401.15),
  −3.70 (4780.54), −3.79 (19109.78), −3.92 (19195.27), −4.028 (436462.76),
  −4.089 (1239915.02), −4.123 (4951517.75), −4.231 (12397185.45), **−4.3433 (18133579.07)**; nothing deeper up to 10^8
  (next best −4.297 at 46529079.52). 99 local minima with S_T < −4.0717 below 10^8 (`scan_out/`).
* Certified points for a = 1 (S_T + tail < −4): u = 1239915.0184 (bound −4.01741, the smallest found),
  u = 18133579.0743 (bound −4.27166, used in Lemma 2). The first candidate, u = 436462.7634 (x ≈ 10^{189553}),
  has S_T = −4.0281 but S_T + tail = −3.956: **not certified** with T = 400 (it would need the tail below 0.028,
  T ≈ 1300).
* Heuristic size (`li_estimate.py`, linear independence + Chernoff bound, not a proof): frequency of
  {S < −4} about 10^{−6}, consistent with the first occurrence near u ≈ 4·10^5.

## 4. Sanity check with actual primes (`sieve.c`, `mertens_const.py`, `compare_primes.py`)

log f(x; 13, 1) was computed from all primes ≡ 1 (mod 13) up to **10^10** (segmented sieve, long double with
compensated summation), with log C(13, 1) = 0.10035428120829340689… (C(13,1) = 1.105562…, 1/C = 0.904516907920896…)
from the L-function expression of Languasco–Zaccagnini type (`mertens_const.py`, Arb).

* **No failure of (3)** for p̄_k ≤ 10^10, i.e. for all k ≤ 37 919 477 (`sieve_out/records_q13_a1_1e10.txt`, nfail = 0).
  With G(x) := 12 √x log x · log f(x; 13, 1): max over k > 10 of G(p̄_k) = −1.3394 at p̄_k = 56 729 011 (k = 281 613);
  by (★) this corresponds to S ≈ −2.66, the scale of the minima of S in that range of u (§3.3).
* **Formula (★) against primes** (samples x = ⌊e^u⌋, u = 9.0, 9.1, …, 23.0; zeros below 400):

      range        G_actual range     G − [−S − 4]: mean / rms / max     G − refined: mean / rms / max     corr (leading / refined)
      u ≥ 12       [−7.57, −1.37]     −0.056 / 0.136 / 0.361             −0.003 / 0.033 / 0.101            0.9975 / 0.99985
      u ≥ 18       [−6.44, −1.41]     +0.034 / 0.062 / 0.091             −0.001 / 0.006 / 0.012            0.9996 / 0.99999
      u ≥ 20       [−6.44, −3.54]     +0.004 / 0.046 / 0.073             +0.001 / 0.006 / 0.012            0.9975 / 0.99996

  "refined" = the same main terms without expanding F_s(x) = x^{s−1}/u + s E_1((1−s)u) in powers of 1/log x, plus
  the terms x^{1/k}, k = 3, 5, 6 of ψ − θ (see the docstring). The leading formula agrees to O(1/log x) and the
  refined one to 0.006 while G varies by several units, e.g. u = 22.1 (x = 3 961 941 421): G = −6.4372,
  −S − 4 = −6.3726, refined −6.4640; u = 23.0: −3.5417, −3.5425, −3.5661. A wrong sign, a missing conjugate or a
  threshold other than 4 would show as an O(1) discrepancy.
* e^{u₀} is out of reach: u₀ = 18 133 579.0743 means x₀ ≈ 10^{7 875 313.3}, k ≈ x₀/(12 u₀) ≈ 10^{7 875 305}, and
  N̄_k ≈ exp(x₀/12). Even the first candidate has x ≈ 10^{189 553}. No explicit counterexample k for q = 13 can
  be exhibited, and none exists below 10^10.

## 5. Other square classes mod 13 (Theorem A′)

**Theorem A′ (GRH_13).** For every square class a ∈ {1, 3, 4, 9, 10, 12} mod 13, inequality (3) with q = 13 fails
for infinitely many k (N̄_k = product of the first k primes ≡ a mod 13), and
lim sup_{x→∞} 12 √x log x · log f(x; 13, a) ≥ margin(a) > 0.

Proof: as for Theorem A, with S_a in place of S (Theorem 1.18 under SH_{13,a}, implied by GRH_13; m = 2,
R_{13,a} = 2; Lemma 3 holds verbatim for S_a since |χ̄(a)| = 1). Certificates (`certify.py 13 a u₀`,
`certificate_q13_a<a>.json`; a window (χ, γ) contributes 2w cos(γu₀ − arg χ(a))):

    a     u₀               S_T(u₀)          bound S_T + tail    margin
    1     18133579.0743    −4.34330566      −4.27166085         0.2717
    3      6533826.3615    −4.23351550      −4.16187068         0.1619
    4     18673543.4487    −4.24466188      −4.17301707         0.1730
    9     14434147.2037    −4.32380278      −4.25215796         0.2522
    10     6560193.3352    −4.31884465      −4.24719984         0.2472
    12     4029543.3045    −4.25945688      −4.18781207         0.1878

(tail = F_13 − W = 0.07164482 in every case). The scans for a ≠ 1 covered u ≤ 2·10^7 (a = 3: ≤ 1.9·10^7,
a = 10: ≤ 7.4·10^6) and were stopped once a comfortable margin was found; the u₀ are not optimal. For a ≠ 1 there
is no analogue of Corollary B here: if GRH_13 failed but SH_{13,a} held, the present certificate would not apply.

Check of the coefficient χ̄(a) against primes, a = 3 (`sieve_out/compare_q13_a3.txt`; primes ≡ 3 mod 13 up to 10^9,
log C(13, 3) = −0.3119342579864…): with the phases as above, G − refined has rms 0.009 (u ≥ 18; corr 0.9999); with
the conjugate phases (those of a = 9) the discrepancy is O(1) (mean −2.3, rms 2.4, corr 0.50). No failure of (3)
for a = 3 and p̄_k ≤ 10^9 (k ≤ 4 236 872; max G = −0.603). The classes a = 4, 9, 10, 12 were not sieved. (Since the
set of square classes is closed under inversion and all six are certified, Theorem A′ as a whole would survive
even an interchange χ(a) ↔ χ̄(a).)

## 6. q = 11

F_11 = 4.26098195 and 2R = 4: one needs S_1(u) < −4, i.e. 93.9 % of full alignment, with total slack
F_11 − 4 = 0.261 for the misalignment of **all** zeros plus the tail.

* Zeros: `zeros_q11.txt` (sha256 `eea2f78092d7a0e9e5503c43722d7a77b4849ad58f50d8ab8657a33afb1d93c7`), 3385 windows,
  complete below T = 400, same certificates ([Z], [N]; `logs/verify_q11.log`: PASS). W = 4.20282824,
  tail F_11 − W = 0.05815371 (density estimate 0.058185).
* Weight distribution: γ_1 = 1.23119 carries 1.133 (27 %); 10 zeros: 2.51; 100 zeros: 3.68; 500 zeros: 4.04;
  1000 zeros: 4.12. So S < −4 requires about the first 500–1000 ordinates to be aligned simultaneously.
* Scan of u ≤ 5.96·10^7 (same method): best S_T = −2.7394 at u = 16133306.6866. **Certified:
  S_1(16133306.6866) ≤ −2.68128** (`certificate_q11_a1.json`; like Lemma 2 this bound assumes GRH_11, i.e. RH and
  GRH for the ten non-principal L(s, χ) mod 11), i.e. **under GRH_11** the lim sup of (7) is ≥ 2.6813 for q = 11.
  **This does not exceed the threshold 4.** Status of q = 11: **open**.
* Heuristic (`li_estimate.py`, linear independence, Chernoff rate; zeros below 400 only): {S_1 < −4} has frequency
  of order 10^{−260} or smaller (10^{−390} if the tail above 400 is taken at its worst), so a witness u would have
  several hundred digits and its certificate would need simultaneous Diophantine approximation in dimension of
  order 10^3 and zeros to several hundred digits. A scan cannot find it. By the same heuristic the frequency of
  {S_1 < −3.3} is about 10^{−23}, so an LLL with u up to 10^25 (the limit of the present 10^{−30} windows) could
  reach about −3.3: still below the threshold. Not attempted.
* If (under GRH_11) the positive ordinates are linearly independent over ℚ, then by Kronecker's theorem
  inf S_1 = −F_11 < −4, so (3) fails infinitely often for q = 11 too. (sup S_1 = F_11 holds trivially, at u = 0,
  and says nothing; an earlier version of this note wrote "sup S_1 = F_11 > 4" here, which was a misstatement.)
  If one wanted to prove the opposite (that (3) holds from some point on), one would need inf S_1 > −4, hence a
  linear relation among the ordinates. Neither is within reach of this computation.
* Sieve: for q = 11, a = 1, inequality (3) holds for all p̄_k ≤ 10^10 (k ≤ 45 504 543), `sieve_out/sweep_a1.txt`.

## 7. Other moduli

These are by-products; they were not in the brief but bear directly on "which q admit a Nicolas-type criterion".

### 7.1 Prime moduli 17 ≤ q ≤ 89, a = 1: explicit failures at small k (no hypothesis)

`sieve.c` for every prime 11 ≤ q ≤ 89 (x ≤ 10^9; 10^10 for q = 11, 13, 19), `sieve_out/sweep_a1.txt`
(long double; C(q,1) from `mertens_const.py`):

    q     first k with (3) false (p̄_k)     number of failing k with p̄_k ≤ X
    11    none up to 10^10                  0
    13    none up to 10^10                  0
    17    k = 5   (409)                     175
    19    none up to 10^10                  0
    23    k = 5   (599)                     40 052
    29    k = 28  (6323)                    6 345
    31    none up to 10^9                   0
    37    k = 282 (107 671)                 31 146
    41    k = 124 (50 021)                  77 984
    43    k = 4   (1033)                    525 948
    47    k = 2   (659)                     139 037
    53    none up to 10^9                   0
    59    k = 2   (827)                     391 678
    61    k = 7   (3539)                    24
    67    k = 5   (3083)                    693 716
    71    k = 2   (853)                     350
    73    k = 1   (293)                     56 633
    79    k = 3   (2213)                    638 372
    83    k = 1   (167)                     3
    89    none up to 10^9                   0

Re-checked in Arb ball arithmetic (`explicit_failure.py`, `logs/explicit_failures.log`): (3) FAILS for
(q, k) = (17, 5), (23, 5), (29, 28), (43, 4), (47, 2). Example: q = 17, N̄_5 = 103·137·239·307·409 = 423464858827:

    N̄_5/(φ(N̄_5)(log(16 log N̄_5))^{1/16}) = 0.9179425463  <  1/C(17,1) = 0.9181134793.

The auditor obtained the same two numbers with an independent constant (log C(17,1) = 0.0854342802 ± 2·10^{−7},
own formula and own sieve): log(LHS/RHS) = −1.86·10^{−4}, about 900 times the error bar; also failures at k = 6, 7,
and at (q, k) = (23, 5), (43, 4) (`audit/logs/explicit17.log`, `audit/logs/explicit_23_43.log`).

**What this does and does not show.** It is **not a counterexample to anything in the paper.** Theorem 1.22
("GRH_q ⇔ (3) for all k ≥ 1") is stated only for q ≤ 10 and q = 12, 14; for q > 14 the paper claims only that (3)
holds for infinitely many k (Cor. 1.15) and the non-square results (Cor. 1.17), and it notes F_p > 4 for
7 < p ≤ 149. What the computation shows is that the literal q = 17 analogue of the criterion,
"GRH_17 ⇔ (3) for all k ≥ 1", has a false right-hand side, so that analogue could hold only if GRH_17 were false.
The same for the other q in the table with a failure. Such finite failures do not decide the asymptotic question
(Corollary 1.19: is log f < 0 from some point on?), which is what Theorem A answers for q = 13.

Novelty, stated with its limits: a failure of (3) for a = 1 is **not found in the arXiv version (which plots
log f only for q ≤ 10, and whose only displayed failures are for non-square a mod 7), in Francis's 2017
Lethbridge thesis, or in the citing literature searched.** The thesis and literature search are the auditor's
(AUDIT §8: thesis text-searched; Crossref, OpenAlex and Semantic Scholar list one citing work,
Aoudjit–Berkane–Dusart 2021 on Robin's criterion, unrelated; web searches on q = 11, 13, 17 and Question 1.24
found nothing). The full journal text was not accessible, so that residual risk remains.

### 7.2 q = 17 and q = 19: infinitely many failures (GRH_q), same method

Zeros below T = 100 (`zeros_q17.txt`, sha256 `217bf7f6dc3021ba806240ad94174a3469138e109b0fff4ef6a04fda74217562`,
1301 windows; `zeros_q19.txt`, sha256 `221c699e97c522d94c73015da7b36a471532a36520f3efd1dacaf22f390cfe39`, 1476
windows; ζ below 400; `verify_zeros.py --count`: PASS for both, `logs/verify_q17.log`, `logs/verify_q19.log`),
F_17 = 13.0206685, F_19 = 18.1646123, R_{q,1} = 2:

    q     u₀            S_T(u₀)        tail       bound S_T + tail     margin
    17    11761.4153    −8.21558101    0.31912875  −7.89645226         3.896
    19    165.7829      −9.26218041    0.36747455  −8.89470586         4.895

Hence, by the proof of Theorem A: **under GRH_17 (resp. GRH_19), (3) with a = 1 fails for infinitely many k.**

**Unconditional version (observation of the auditor, AUDIT §9 MINOR 1).** The argument of Corollary B applies
verbatim: if GRH_q fails then SH_{q,1} fails and Akbary–Francis Theorem 1.16(a) gives log f(x; q, 1) = Ω_+(x^{−b}),
provided no L(s, χ), χ mod q, has a zero in (0, 1). Hence **(3) with a = 1 fails for infinitely many k for q = 17
and for q = 19, with no unproved hypothesis.** Sources of the real-zero input, stated precisely:
* q = 17: (i) the auditor's Arb ball check that L(σ, χ) ≠ 0 on [1/2, 1] for all 15 non-principal characters
  mod 17 (`audit/real_zeros.py`, `audit/logs/real_zeros_q17.log`: "ALL: no real zero in [1/2,1]"); (ii) my count
  [N] at T = 100 (`logs/verify_q17.log`: no zero off the critical line and none at s = 1/2 below height 100);
  (iii) Platt, via the paper's Corollary 1.13 (q ≤ 400 000).
* q = 19: the auditor did not run his check and cites Platt via Corollary 1.13; in addition my count [N] at T = 100
  (`logs/verify_q19.log`) certifies the same for the 17 non-principal characters mod 19. Note the near miss:
  L(1/2, χ_19(13, ·)) = 0.00654 − 0.01799 i, small but certified non-zero (the zero is at height 0.0189564).
For q = 19 nothing fails below 10^10 because χ_19(13, ·) (order 18, odd) has a zero at height γ = 0.0189564
(L(1/2, χ) = 0.00654 − 0.01799 i): the term 7.99 cos(γu) keeps S_1 > 0 until u ≈ π/γ; the first candidate failure
is at u ≈ 165.78, x ≈ 10^{72} (heuristic for the location; the certificate proves only the lim sup statement).
q = 31, 53, 89 (no failure up to 10^9; 4/F_q ≤ 0.11) were not treated; the same pipeline applies.

### 7.3 The set {q : F_q < 2R_{q,1}} for q ≤ 1000 (`fq_table.py`, `fq_table_1000.txt`)

The paper computes F_q for q ≤ 14 (and q = p, 2p) and asks for "the finite set of integers q for which
F_q < 2R_{q,1}" (p. 9). With the closed forms (49), (50) summed over the primitive characters of all conductors
d | q (Arb, 100 bits; Table 1 of the paper is reproduced to all printed digits):

    F_q < 2R_{q,1}  for q ≤ 1000  exactly when
    q ∈ {1, …, 10, 12, 14, 15, 16, 18, 20, 21, 24, 28, 30, 36, 40, 42, 48, 60}        (25 moduli).

New relative to the paper's list: q = 15, 16, 18, 20, 21, 24, 28, 30, 36, 40, 42, 48, 60 (e.g. F_60 = 8.458 < 16,
F_48 = 11.039 < 16, F_36 = 6.258 < 8). For these, Corollary 1.19 gives at once, under GRH_q: (3) fails for at most
finitely many k, for every square class (the analogue of their Prop. 1.21); an "all k" statement as in
Theorem 1.22 would need their §8 computation, not done here. Closest misses: q = 120 (F = 32.473 vs 32),
q = 11, 22 (4.261 vs 4), 56 (22.06 vs 16), 84 (22.70 vs 16), 72 (23.67 vs 16), 13, 26 (6.455 vs 4). For
500 < q ≤ 1000 the smallest ratio F_q/2R_{q,1} is 7.25 (q = 840). This is a table, not a determination: that no
q > 1000 qualifies needs an explicit lower bound for F_q (their (11) is not explicit).

## 8. Trusted inputs and caveats

**Hypotheses.** Theorem A, A′: RH and GRH for the 11 non-principal L(s, χ) mod 13, nothing else (no linear
independence, no simplicity). §7.2, certified lim sup bounds: the same mod 17, mod 19; §6, the bound 2.6813: the
same mod 11. Corollary B, the unconditional statements of §7.2, and §7.1: no unproved hypothesis (the first two
use Akbary–Francis Theorem 1.16(a) and the absence of real zeros).

**On the announced proof.** A proof of RH and of GRH for Dirichlet L-functions was announced on 2026-11-20
(described to us as Lean-checked and unrefereed). We have not examined it, nothing in this note uses it, and no
statement here should be read as unconditional because of it. If it is accepted, Theorem A and A′ lose their
hypothesis; Corollary B does not need it in any case.

**Trusted inputs.**
1. Akbary–Francis, Theorem 1.18 (arXiv v1; proof read in full: Props 3.1, 3.3, 3.5, Lemma 3.4, Thm 6.1, Lemma 6.2,
   Cor. 6.3, Lemma 6.4; I found no error; leading term re-derived in §1; tested against primes in §4), and for
   Corollary B their Theorem 1.16(a) (proof read, not re-derived). Through them: McCurley's explicit formula for
   ∫ψ(t; q, a) dt (their Thm 6.1), Nicolas's Lemma 2.2 (their Lemma 3.4), the PNT in progressions.
2. The Hadamard-product identities of Lemma 1 (Montgomery–Vaughan Cor. 10.14, 10.18; quoted from the paper's
   (49), (50), not from the book; verified numerically per character against the zero sums, §3.2, and by mpmath).
3. Arb / python-flint 0.9.0: `acb_dirichlet_l`, `acb_dirichlet_hardy_z` (including that Z(t, χ) is real with
   |Z| = |L(1/2 + it, χ)| for complex primitive χ), `lgamma`, Stieltjes constants, `zeta_zeros`, `zeta_nzeros`,
   elementary functions. The proof of Lemma 2 uses only: each window contains a zero; the closed form; ball
   arithmetic for 2w cos(γu₀).
4. Classical facts used by the count [N] (only needed for Corollary B and for "exactly one simple zero per
   window"): functional equation, zeros in the critical strip.

**Caveats.**
* No explicit k for q = 13. The certified u₀ corresponds to x₀ ≈ 10^{7 875 313}; (3) holds for all p̄_k ≤ 10^10.
  Theorem A does not say that failure occurs at x₀ itself (non-explicit O-constant; Remark 2 of §2).
* "First failure near x ≈ 10^{189 553}" (u ≈ 436 462.76) is a guess from a non-rigorous scan and is not certified.
  The scans (`scan.c`) are double-precision and non-exhaustive; they only supply candidates, which Arb certifies.
* Question 1.24 (lim sup = F_q?) is not answered; only lim sup ≥ 0.662 F_13 (and ≥ 0.61 F_17, ≥ 0.49 F_19, ≥ 0.63 F_11).
* q = 11 is open (§6). Work order 2 of the triage (explicit exception lists for square classes, q ≤ 14) was not
  started; §7.3 enlarges its scope (13 further moduli).
* The journal version (Math. Comp. 89 (2020), no. 322, 993–1026, DOI 10.1090/mcom/3463) was not read: its full
  text was not accessible (HTTP 403). Statement numbers refer to arXiv v1. Novelty claims are therefore worded as
  "not found in the arXiv version, Francis's 2017 thesis, or the citing literature searched" (search by the
  auditor, AUDIT §8); this covers the q = 13 results, the explicit failures for a = 1, and the 13 extra moduli.
  Before any public claim, the journal PDF should be checked for anything on q = 11, 13, 17.
* The q ≥ 17 failures are not counterexamples to the paper (§7.1).
* `mertens_const.py` checks the branch of log L(1, χ) only at sample points; for a = 1 only real parts enter, and
  the constants are confirmed by the sieve (log f → 0), so §4 and §7.1 do not depend on it.
* The sieve uses long double arithmetic with compensated sums, not interval arithmetic; the five failures listed in
  §7.1 as re-checked are in Arb, the others (and "no failure up to 10^10") are floating point with margins far
  above rounding error but not certified. The sweep over q stops at 89 (the C(q,1) routine is slow for larger q).
* Files in this directory were committed by another process ("Work in progress: intermediate agent output"),
  not by this agent.

## 9. Files and replay

Scripts (python-flint 0.9.0, numpy/scipy/mpmath for the non-rigorous parts; `gcc -O2 -o sieve sieve.c -lm`,
`gcc -O3 -o scan scan.c -lm`; all runs with `nice -n 10`, at most 2 processes):

    lzeros_cert.py q H          certified zeros of all L(s,χ) mod prime q below H  → state/, zeros_q<q>.cert.log
    zeta_zeros_cert.py 400      ζ zeros → state/zeros_zeta.json
    assemble.py q               → zeros_q<q>.txt; its sha256 goes into SHA256SUMS
    verify_zeros.py q --count   re-certify [Z], [N] from the assembled, hash-pinned file
    closed_form.py 13 11        F_q (Lemma 1); f_two_ways.py q: against the zero sum (common T: q = 13, 11)
    make_scan_input.py, scan.c, refine_u.py     exploratory search for u₀ (double precision)
    certify.py q a u₀           THE certificate (Lemma 2) → certificate_q<q>_a<a>.json
    crosscheck_mpmath.py 13 1 18133579.0743     independent recomputation (labels, F_13, S_T(u₀))
    mertens_const.py q a; sieve.c; compare_primes.py     primes vs explicit formula (§4, §5), sweep (§7.1)
    explicit_failure.py q a k   Arb check of (3) at an explicit primorial
    fq_table.py 1000            §7.3;   li_estimate.py q: heuristic frequencies
    check_bindings.py           every certificate_*.json ↔ pinned sha256 of its zero file (added after the audit)

Replay of the proof of Theorem A from the pinned data (about 1 minute):

    python3 verify_zeros.py 13 --count      # [Z] 4208 windows, [N] at T = 400: PASS
    python3 certify.py 13 1 18133579.0743   # CERTIFIED True, margin 0.27166084556

Data: `zeros_q{11,13,17,19}.txt` (+ `.cert.log`), `SHA256SUMS`, `certificate_*.json`
(`certificate_q13_a1_small.json`: u = 1239915.0184; `…candidate436462_NOT_certified.json`: the first candidate;
`certificate_q11_a1.json`: certified = false, the best value for q = 11), `closed_form_q*.json`,
`sieve_out/` (samples and comparisons for q = 13, a = 1 and a = 3; `sweep_a1.txt`; large record files were trimmed;
the compiled `sieve` and `scan` binaries were removed),
`scan_out/` (scan inputs and candidate lists), `fq_table_1000.txt`, `logs/`, `source/` (the PDF and its text).
`explore/` holds the first uncertified experiments (superseded). `AUDIT.md` and `audit/` are the independent
audit (its code, logs, and `audit/SHA256SUMS_audited`, which pins the files as audited; of those only `RESULT.md`
has changed since).

## 10. Review history

**2026-11-24, independent audit (`AUDIT.md`, code and logs in `audit/`).** Verdict: no rigor break. Theorem A,
Theorem A′, Corollary B and the explicit q = 17 failure were re-verified with different code: every window of the
four zero files re-certified by a complex Krawczyk test without `hardy_z` or `dirichlet_char` (0 failures;
ordinates good to about 10^{−42}); F_13, W, the tail and S_T(u₀) recomputed in mpmath at 80 digits (identical
digits; margin 0.27166084556380); the six certificates of Theorem A′ and those for q = 17, 19, 11 reproduced; the
sign, normalisation and the constant 4 of (★) re-tested on an independent sieve to 10^9 (wrong sign: rms 2.60,
wrong threshold: rms 3.93, correct: rms 0.075 for u ≥ 18); Theorem 1.16(a) checked line by line; real zeros
excluded independently for q = 13 and q = 17; C(17,1) and the k = 5 failure recomputed; seven entries of the F_q
table spot-checked. The auditor read the TeX source of arXiv v1, so the conjugation bar and the signs in §1 are
confirmed from the LaTeX.

Changes made to this note in response (no mathematical statement was weakened; nothing was recomputed except the
binding check):

* **G1 (novelty and citation).** Journal version now cited in the header and in §8 (Math. Comp. 89 (2020),
  no. 322, 993–1026, DOI 10.1090/mcom/3463), with the statement that its full text was not accessible to us and
  that numbering refers to arXiv v1. Novelty reworded as "not found in the arXiv version, Francis's 2017 thesis,
  or the citing literature searched", with the auditor's search recorded (§7.1, §8). The sentences "Journal
  version: not checked" and "I have not searched the later literature" are gone.
* **G2 (statement of Theorem A).** The clause saying that RH and GRH mod 13 are covered by the announced proof was
  removed from the statement of Theorem A. Theorem A assumes RH and GRH mod 13, full stop. One remark about the
  announcement remains, in §8 ("On the announced proof"), saying that it is not used.
* **MINOR 1.** §7.2 now states the unconditional "infinitely many k" for q = 17 and q = 19 by the Corollary B
  argument, with the source of the real-zero input given precisely for each q (auditor's check on [1/2, 1] for
  q = 17; Platt via Cor. 1.13 for q = 19; plus this note's own [N] counts at T = 100 for both). §0 and §8 updated.
* **MINOR 2.** §6: "sup S_1 = F_11 > 4" corrected to "inf S_1 = −F_11 < −4" (Kronecker, under linear independence);
  the certified bound "lim sup ≥ 2.6813" for q = 11 now carries its hypothesis GRH_11.
* **MINOR 3.** §1: remark that the "only if" direction of their Corollary 1.19 is not justified when the lim sup
  equals 2R exactly; irrelevant here because of the strict margin, and §2 proves the direction used directly from
  Theorem 1.18.
* **MINOR 4.** Corollary B reworded ("unconditional, using Akbary–Francis Theorem 1.16(a)"; "no criterion of the
  form GRH_13 ⇔ (3) for all k ≥ k₀ can hold unless GRH_13 is false"). The proof paragraph lists the three sources
  of the real-zero input.
* **MINOR 5.** Lemma 2: the 10^{−30} windows are deliberately loose; the ordinates are known to about 10^{−42}.
* **MINOR 6, 7.** No change needed: completeness [N] is not used by any stated theorem except through the
  real-zero input (now triply sourced); the floating-point status of the §7.1 sweep was already stated. The
  auditor's independent re-checks of (17, 5), (23, 5), (43, 4) are recorded in §7.1.
* **q = 17, k = 5.** §0 and §7.1 now say explicitly that this is not a counterexample to anything in the paper
  (Theorem 1.22 is stated only for q ≤ 10 and q = 12, 14) and record the auditor's numbers
  (0.9179425463 against 1/C(17,1) = 0.9181134793).
* **Hash binding.** Confirmed and documented in §3.1: `afcommon.load_zeros` refuses any zero file whose sha256
  differs from `SHA256SUMS`, and `certify.py` writes the sha256 into each certificate. New script
  `check_bindings.py` (not part of any proof) confirms this for all 11 certificate files and includes a negative
  control (`logs/check_bindings.log`). No audited script, zero file or certificate was modified:
  `sha256sum -c audit/SHA256SUMS_audited` reports every file unchanged except `RESULT.md`.

Open after the audit: the journal full text (G1 residual); q = 11; Question 1.24; an explicit version of
Theorem 1.18 (failure at x₀ itself); Work order 2.
