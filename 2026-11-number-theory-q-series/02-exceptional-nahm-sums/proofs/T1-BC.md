# The generalized Nahm sums of the pairs (T1,C_r) and (T1,B_r): Bailey chains with half steps

> **Referee corrections (2026-11-26, `REFEREE-T1-ABC.md`; they supersede any contrary statement below).**
> No mathematical error was found, but **Theorem B is not new**: (a) its identity (0.2) for all r is verbatim
> Warnaar, math/0110307, the unnumbered identity after (5.6) with k = r+1 (reproduced as (2.3) in
> Cherednik–Feigin arXiv:1209.1978); (b) for every r, both parity parts follow from Du–Ji–Shen–Xu
> arXiv:2602.14483 (Theorems 1.1/1.2, 1.5): for even r by fusion, and for odd r — contrary to the sentence
> "I have no direct transformation" below, which is false — by one self-fusion of the short node (Lemma S of
> `KRTW-k2.md`) followed by defusions (checked exactly up to rank 9 by the referee). Corollary B is their
> j = 0 case. What remains ours: the identification of these known theta quotients with the Dynkin pair
> (T1,B_r) of Sun–Wang's Table 1 (listed there with CFT "?"), i.e. Conjecture 1.1 for (T1,B_r) follows from
> known theorems; and a uniform half-step Bailey proof. The Bailey pair called "Slater's G(2)" below is G(2).
> Corollary C (linear-term companions of (T1,C_r) as differences of NS characters of SM(4r+6,4)) is correct
> (referee: r ≤ 5, all linear terms) and plausibly unstated, but it is a short theta dissection of two quoted
> theorems; for r = 1, i = 1 it needs Euler's identity (Theorem C′ is quoted only for r ≥ 2). In §4 (i even)
> the remark "(k → −k−1 in T_2)" should be deleted.


Written 2026-11-26. Exact certificate and numerical checks: `code/sw_certify_T1.py`, log
`certificates/sw_certify_T1.txt` (parts (a)–(d)); see §6.

**What is new and what is not (read this first).**

* (T1,C_r): the sum = product identity is **known** for every r (Theorem C below is Warnaar's identity, quoted by
  Sun–Wang themselves as their (2.8); the companions with linear terms are Warnaar's unpublished identity, proved
  in Wang–Wang arXiv:2402.06253, Theorem 4.1). So Sun–Wang's (3.6) was never open as a statement about q-series:
  it is (2.8) rewritten. What is done here: a short self-contained proof of the B = 0 case (§3), and the explicit
  dictionary to the characters of SM(4r+6,4) for **all r and all r+1 linear terms** (Corollary C, an elementary
  dissection of a theta series), which contains (3.6)–(3.13) and explains why (3.6) prints the same formula twice.
* (T1,B_r): Theorem B gives a closed formula for every r ≥ 1. Sun–Wang's Table 1 lists (T1,B_r) with CFT "?" and
  no reference; of this family only r = 2 (= (T1,C_2)) and r = 3 (Wang–Wang arXiv:2407.21725, Mizuno's Example 7)
  are mentioned there as known. **However** the two theta-quotients of Theorem B coincide with the right-hand
  sides of Du–Ji–Shen–Xu, arXiv:2602.14483, Theorem 1.5 (1) and (3) (identities (1.7), (1.9), proved there for
  other matrices of rank about r/2), and for even r the two sums are related by the fusion lemma (§5). So the
  *modularity* of the (T1,B_r) sum is new at most in the sense that this link was not noticed (Sun–Wang do not
  cite Du–Ji–Shen–Xu); the proof given here is independent of theirs and uniform in r.

## 0. Notation and statements

(a;q)_n, (a;q)_∞, (q)_n = (q;q)_n as usual. All identities are identities of formal power series in q^{1/4}
(resp. q^{1/2}) with integer coefficients; every multiple series below is a summable family (the exponent of
q grows quadratically). N denotes the set of non-negative integers.

**(T1,C_r).** D(C_r) = diag(1,…,1,2), A(T1,C_r)D = C(C_r)⁻¹D = (min(i,j))_{i,j≤r} (checked in `sw_certify_T1.py`).
For r ≥ 1 and 1 ≤ i ≤ r+1 put N_j = n_j + … + n_r and

    S_r^{(i)}(q) = Σ_{n∈N^r}  q^{½(N_1²+…+N_r²) + N_i + N_{i+2} + N_{i+4} + …} / ( (q)_{n_1}⋯(q)_{n_{r−1}} (q²;q²)_{n_r} )

(the linear term is the sum of the N_j with j ≡ i mod 2, i ≤ j ≤ r; it is empty for i = r+1). S_r := S_r^{(r+1)}
is the generalized Nahm sum of the pair (T1,C_r) with B = 0.

**Theorem C** (Warnaar; see §4 for attribution). For all r ≥ 1,

    S_r(q) = ( (−q^{1/2};q)_∞ / (q;q)_∞ ) · Σ_{n∈Z} (−1)^n q^{((2r+3)n² − n)/4}.                                  (0.1)

(By the Jacobi triple product the theta series is (q^{(r+1)/2}, q^{(r+2)/2}, q^{r+3/2}; q^{r+3/2})_∞, and
(−q^{1/2};q)_∞/(q)_∞ = 1/((−q;q)_∞(q^{1/2};q^{1/2})_∞); so (0.1) is Sun–Wang's (2.8).)

**Theorem C′** (Warnaar unpublished; Wang–Wang, arXiv:2402.06253, Theorem 4.1 — *quoted, not reproved here*).
For r ≥ 2, 1 ≤ i ≤ r:  S_r^{(i)}(q) = ((−q^{1/2};q)_∞/(q)_∞) Σ_{n∈Z} (−1)^n q^{((2r+3)n² − (2r+3−2i)n)/4}.

**Corollary C** (character form; Sun–Wang (3.6)–(3.13) for all r). Let χ^{NS}_{ρ,s} be the NS characters of the
N = 1 minimal model SM(4,4r+6) (conventions in §4), 1 ≤ i ≤ r+1 and j = 2r+3−2i. Then

    q^{(j² − (2r+3))/(16(2r+3))} · S_r^{(i)}(q) = χ^{NS}_{1,i} − χ^{NS}_{1,4r+6−i}             (i odd),
                                               = χ^{NS}_{1,2r+3−i} − χ^{NS}_{1,2r+3+i}         (i even).

The effective conformal weights (h − h_min) of the two characters are (j²−1)/(16(2r+3)) and
((4r+6+(−1)^{i+1}j)² − 1)/(16(2r+3)). For i = r+1 (B = 0) they are 0 and ⌊r/2⌋+1, and the prefactor is q^{−c/24},
c = 3(r+1)/(2r+3): this is (3.6) of Sun–Wang, **for every r ≥ 1** (the case i = r+1 uses only Theorem C).
In Kac labels: χ^{NS}_{1,r+1} − χ^{NS}_{1,3r+5} for r even, χ^{NS}_{1,r+2} − χ^{NS}_{1,3r+4} for r odd — the two
cases of (3.6) differ in the labels but not in the weights, which is why the printed formula is the same twice.

**(T1,B_r).** D(B_r) = diag(2,…,2,1), and ½nᵀA(T1,B_r)Dn = ¼(M_1²+…+M_r²) with M_i = 2n_i+…+2n_{r−1}+n_r
(i < r), M_r = n_r (checked in `sw_certify_T1.py`). Put

    R_r(q) = Σ_{n∈N^r}  q^{¼(M_1²+…+M_r²)} / ( (q²;q²)_{n_1}⋯(q²;q²)_{n_{r−1}} (q;q)_{n_r} )          (r ≥ 1)

(for r = 1: R_1 = Σ_n q^{n²/4}/(q)_n, the Nahm sum of (T1,A1)); R_2 = S_2.

**Theorem B.** For all r ≥ 1, R_r = R_r^{even} + R_r^{odd}, where R^{even} (R^{odd}) is the part of the sum with n_r
even (odd), and

    R_r^{even}(q) = ( (−q;q²)_∞ / (q²;q²)_∞ ) · Σ_{m∈Z} (−1)^m q^{(4r+6)m² − 2m},                                 (0.2)
    R_r^{odd}(q)  = q^{r/4} ( (−q²;q²)_∞ / (q²;q²)_∞ ) · Σ_{n∈Z} (−1)^{n(n+1)/2} q^{((2r+3)n² + (2r+1)n)/2}          (0.3)
                 = q^{r/4} ( (−q²;q²)_∞ / (q²;q²)_∞ ) · [ (q^{2r+5},q^{6r+7},q^{8r+12};q^{8r+12})_∞ + q (q^{2r+1},q^{6r+11},q^{8r+12};q^{8r+12})_∞ ].

**Corollary B** (Sun–Wang Conjecture 1.1 for (T1,B_r), B = 0). With c = 3(2r−1)/(2r+3),

    q^{−c/24} R_r(q) = (η(2τ)/(η(τ)η(4τ))) · Σ_{m∈Z} (−1)^m q^{(4r+6)(m − 1/(4r+6))²}
                     + (η(4τ)/η(2τ)²) · Σ_{n∈Z} (−1)^{n(n+1)/2} q^{(2r+3)(n + (2r+1)/(4r+6))²/2},

a sum of two quotients of a unary theta series of weight ½ by an eta quotient of weight ½; hence a modular
function (weight 0) on a congruence subgroup. [Exponents: (−q;q²)_∞/(q²;q²)_∞ = q^{1/8}η(2τ)/(η(τ)η(4τ)),
(−q²;q²)_∞/(q²;q²)_∞ = η(4τ)/η(2τ)², and 1/8 − 1/(4r+6) = c/24 = r/4 − (2r+1)²/(8(2r+3)). Modularity of unary theta
series with a congruence condition and a character: Shimura, Ann. of Math. 97 (1973), Prop. 2.1.]

## 1. Bailey pairs and the half step

A pair of sequences (α_n)_{n≥0}, (β_n)_{n≥0} is a *Bailey pair relative to a (in base q)* if for all n ≥ 0

    β_n = Σ_{r=0}^{n} α_r / ( (q;q)_{n−r} (aq;q)_{n+r} ).                                                      (1.1)

We need only a ∈ {1, q}; y denotes the square root (aq)^{1/2} ∈ {q^{1/2}, q}.

**Lemma 1 (half step).** Let (α,β) be a Bailey pair relative to a and y = (aq)^{1/2}. Then

    α′_n = y^n q^{n(n−1)/2} α_n,      β′_n = (1/(−y;q)_n) Σ_{k=0}^{n} (−y;q)_k y^k q^{k(k−1)/2} β_k / (q)_{n−k}

is a Bailey pair relative to a. (y^n q^{n(n−1)/2} = q^{n²/2} for a = 1 and = q^{n(n+1)/2} for a = q.)

This is Bailey's lemma with ρ_1 → ∞, ρ_2 = −(aq)^{1/2}: step (S2) of Bressoud–Ismail–Stanton, *Change of base
in Bailey pairs*, Ramanujan J. 4 (2000) (arXiv:math/9909053, §3). Since the proof is three lines we give it.

*Proof.* Insert (1.1) for β_k and put k = r + i:

    (−y)_n β′_n = Σ_{r=0}^{n} α_r Σ_{i=0}^{n−r} (−y)_{r+i} y^{r+i} q^{C(r+i,2)} / ( (q)_{n−r−i} (q)_i (aq)_{2r+i} ),      C(m,2) = m(m−1)/2.

With (−y)_{r+i} = (−y)_r(−yq^r)_i, (aq)_{2r+i} = (aq)_{2r}(aq^{2r+1})_i, C(r+i,2) = C(r,2)+C(i,2)+ri and
1/(q)_{n−r−i} = (q^{−(n−r)};q)_i (−1)^i q^{(n−r)i−C(i,2)}/(q)_{n−r} the inner sum is

    ( (−y)_r y^r q^{C(r,2)} / ((q)_{n−r}(aq)_{2r}) ) · ₂φ₁( q^{−(n−r)}, −yq^r ; aq^{2r+1} ; q, −yq^n ).

Here −yq^n = (aq^{2r+1}) q^{n−r}/(−yq^r) because y² = aq, so the q-Chu–Vandermonde sum
₂φ₁(q^{−N}, b; c; q, cq^N/b) = (c/b;q)_N/(c;q)_N (Gasper–Rahman, *Basic Hypergeometric Series*, (1.5.2); an identity
of rational functions) applies with c/b = −yq^r:  ₂φ₁ = (−yq^r)_{n−r}/(aq^{2r+1})_{n−r}. Hence the inner sum equals
y^r q^{C(r,2)} (−y)_n / ((q)_{n−r}(aq)_{n+r}), i.e. β′_n = Σ_r α′_r/((q)_{n−r}(aq)_{n+r}). ∎

**Lemma 2 (chain).** Let a ∈ {1,q}, y = (aq)^{1/2}, and let (α,β) be a Bailey pair relative to a such that all α_n
are power series in q^{1/2} (non-negative order; then so are the β_n). Put β̃_n = (−y;q)_n β_n. Then for every s ≥ 1

    Σ_{N_1≥N_2≥…≥N_s≥0}  ( ∏_{j=1}^{s} y^{N_j} q^{C(N_j,2)} ) β̃_{N_s} / ( (q)_{N_1−N_2}⋯(q)_{N_{s−1}−N_s} )
        = ( (−y;q)_∞ / (aq;q)_∞ ) · Σ_{n≥0} y^{sn} q^{s·C(n,2)} α_n.                                                    (1.2)

*Proof.* (i) *Limit formula.* For a Bailey pair (α,β) as in the lemma,

    Σ_{k≥0} y^k q^{C(k,2)} (−y)_k β_k = ( (−y)_∞/(aq)_∞ ) Σ_{r≥0} y^r q^{C(r,2)} α_r.                             (1.3)

Indeed, by Lemma 1, Σ_{k≤n} (−y)_k y^k q^{C(k,2)} β_k/(q)_{n−k} = (−y)_n Σ_{r≤n} y^r q^{C(r,2)} α_r/((q)_{n−r}(aq)_{n+r}); let
n → ∞ coefficientwise (the terms with index k, resp. r, have order ≥ k²/2 − O(k), so only finitely many
contribute to a given power of q, and 1/(q)_{n−k} → 1/(q)_∞, (aq)_{n+r} → (aq)_∞, (−y)_n → (−y)_∞); multiply by (q)_∞.
(ii) Apply Lemma 1 s−1 times. By induction α^{(i)}_n = y^{in} q^{i·C(n,2)} α_n and, because the factor (−y)_n is
carried along, β̃^{(i+1)}_n = Σ_{k≤n} y^k q^{C(k,2)} β̃^{(i)}_k/(q)_{n−k}, so

    β̃^{(s−1)}_{N_1} = Σ_{N_1≥N_2≥…≥N_s≥0} ( ∏_{j=2}^{s} y^{N_j} q^{C(N_j,2)} ) β̃_{N_s} / ( (q)_{N_1−N_2}⋯(q)_{N_{s−1}−N_s} ).

(iii) Apply (1.3) to (α^{(s−1)},β^{(s−1)}). ∎

## 2. Three Bailey pairs

**Lemma 3.** The following are Bailey pairs (base q; q^{1/2} is a formal square root of q):

    (P1)  a = 1:  α_0 = 1,  α_n = (−1)^n q^{n(3n−1)/4} (1 + q^{n/2})  (n ≥ 1);        β_n = 1/( (−q^{1/2};q)_n (q²;q²)_n );
    (P2)  a = 1:  α_0 = 1,  α_{2m} = (−1)^m q^{3m²−m} (1 + q^{2m}),  α_{2m−1} = 0  (m ≥ 1);   β_n = 1/( (q;q²)_n (q;q)_n );
    (P3)  a = q:  α_n = (−1)^{n(n+1)/2} q^{n(3n+1)/4} (1 + q^{n+1/2})/(1 − q);        β_n = 1/( (q^{1/2};q)_{n+1} (q²;q²)_n ).

(P1) and (P2) are G(1) and C(1) of Slater's list (L. J. Slater, *A new proof of Rogers's transformations of
infinite series*, Proc. London Math. Soc. (2) 53 (1951) 460–475); (P3) is Slater's G(2),
α_n = (−1)^n q^{n(3n+1)/4}(1−q^{n+1/2})/(1−q^{1/2}), β_n = 1/((q²;q²)_n(−q^{3/2};q)_n), with q^{1/2} replaced by −q^{1/2} and
divided by 1−q^{1/2}. [I did not have Slater's paper at hand; the statements were verified exactly, as identities
of rational functions of p = q^{1/2}, for n ≤ 14 — check (a) — and they are derived here from one classical
summation.]

*Proof.* We use the terminating very-well-poised ₆φ₅ sum (Gasper–Rahman (II.21)) in the limit of one
parameter → ∞: for a base p and parameters A, c,

    Σ_{r=0}^{n} V_r(A;p) · ( (c;p)_r (p^{−n};p)_r / ((Ap/c;p)_r (Ap^{n+1};p)_r) ) · (−1)^r p^{C(r,2)} (Ap^{n+1}/c)^r = (Ap;p)_n/(Ap/c;p)_n,       (2.1)
    V_r(A;p) := (1 − Ap^{2r}) (A;p)_r / ( (1 − A)(p;p)_r )      (for A = 1: V_0 = 1, V_r = 1 + p^r for r ≥ 1),

an identity of rational functions of A, c, p. [(II.21) reads Σ_r V_r(A;p)(b,c,p^{−n};p)_r/((Ap/b,Ap/c,Ap^{n+1};p)_r)·(Ap^{n+1}/bc)^r =
(Ap,Ap/bc;p)_n/((Ap/b,Ap/c;p)_n); let b → ∞ termwise in the finite sum: (b;p)_r b^{−r} → (−1)^r p^{C(r,2)}.]

*(P1), (P3).* Take base p = q^{1/2} and c = −p^{−n} in (2.1). Then (c;p)_r(p^{−n};p)_r = (q^{−n};q)_r and
(Ap/c;p)_r(Ap^{n+1};p)_r = (A²q^{n+1};q)_r, and for a := A²

    (q^{−n};q)_r/(aq^{n+1};q)_r = (q)_n (aq)_n (−1)^r q^{C(r,2) − nr} / ( (q)_{n−r} (aq)_{n+r} ).                       (2.2)

The powers of p^n cancel ((p^{2n+1})^r q^{−nr} = p^r), and (2.1) becomes

    Σ_{r=0}^{n} [ V_r(A;p) (−A)^r p^{r(3r−1)/2} ] / ( (q)_{n−r}(aq)_{n+r} ) = (Ap;p)_n / ( (q)_n (aq)_n (−Ap^{n+1};p)_n ) = 1/( (q;q)_n (−Ap;p)_{2n} ),

using (aq;q)_n = (Ap;p)_n(−Ap;p)_n. So α_r = V_r(A;p)(−A)^r p^{r(3r−1)/2}, β_n = 1/((q)_n(−Ap;p)_{2n}) is a Bailey pair
relative to a = A². For A = 1: α_r = (−1)^r p^{r(3r−1)/2}(1+p^r) and (−p;p)_{2n} = (−q^{1/2};q)_n(−q;q)_n, which is (P1).
For A = p (a = q): α_r = (−1)^r p^{r(3r+1)/2}(1−p^{2r+1})/(1−p), β_n = 1/((q)_n(−p²;p)_{2n}) (this is G(2)). The Bailey
relation (1.1) for a = q involves p only through α and β and is an identity of rational functions of p; replacing
p by −p gives the pair α_r = (−1)^{r + r(3r+1)/2} p^{r(3r+1)/2}(1+p^{2r+1})/(1+p), β_n = 1/((q)_n (−p²;−p)_{2n}), and
(−p²;−p)_{2n} = ∏_{i<n}(1+p^{2i+2})(1−p^{2i+3}) = (−q;q)_n (q^{3/2};q)_n. Dividing α and β by 1 − p and using
(−1)^{r+r(3r+1)/2} = (−1)^{r(r+1)/2}, (1−p)(q^{3/2};q)_n = (q^{1/2};q)_{n+1}, (q)_n(−q)_n = (q²;q²)_n gives (P3).

*(P2).* Take base p = q², A = 1 in (2.1), with n replaced by N := ⌊n/2⌋ and c := q^{1−n} if n is even, c := q^{−n}
if n is odd. In both cases {c, p^{−N}} = {q^{−n}, q^{1−n}}, {p/c, p^{N+1}} = {q^{n+1}, q^{n+2}} and p^{N+1}/c = q^{2n+1}; since
(q^{−n};q)_{2m} = (q^{−n},q^{1−n};q²)_m, (q^{n+1};q)_{2m} = (q^{n+1},q^{n+2};q²)_m, and by (2.2) with a = 1 and r = 2m
(q^{−n};q)_{2m}/(q^{n+1};q)_{2m} = (q)_n² q^{m(2m−1) − 2mn}/((q)_{n−2m}(q)_{n+2m}), (2.1) reads

    Σ_m V_m(1;q²) (−1)^m q^{m(m−1) + (2n+1)m + m(2m−1) − 2mn} / ( (q)_{n−2m}(q)_{n+2m} ) = (q²;q²)_N / ( (q)_n² (q²/c;q²)_N ).

The exponent is 3m² − m, so the left side is Σ_m α_{2m}/((q)_{n−2m}(q)_{n+2m}) with α as in (P2). The right side is
1/((q)_n(q;q²)_n): for n = 2N, (q^{n+1};q²)_N = (q;q²)_n/(q;q²)_N and (q)_n = (q;q²)_N(q²;q²)_N; for n = 2N+1,
(q^{n+2};q²)_N = (q;q²)_n/(q;q²)_{N+1} and (q)_n = (q;q²)_{N+1}(q²;q²)_N. ∎

## 3. Proof of Theorem C

In the variables N_1 ≥ … ≥ N_r ≥ 0 (n_j = N_j − N_{j+1}, n_r = N_r),
S_r = Σ q^{½(N_1²+…+N_r²)}/((q)_{N_1−N_2}⋯(q)_{N_{r−1}−N_r}(q²;q²)_{N_r}). This is the left side of (1.2) with a = 1, y = q^{1/2},
s = r and β̃_n = 1/(q²;q²)_n, i.e. β_n = 1/((−q^{1/2};q)_n(q²;q²)_n): the pair (P1). Hence

    S_r = ((−q^{1/2};q)_∞/(q)_∞) ( 1 + Σ_{n≥1} (−1)^n q^{rn²/2 + n(3n−1)/4} (1 + q^{n/2}) ) = ((−q^{1/2};q)_∞/(q)_∞) Σ_{n∈Z} (−1)^n q^{((2r+3)n²−n)/4}. ∎

## 4. Characters of SM(4,4r+6): proof of Corollary C; attributions

*Conventions.* For p, p′ with p′−p even and gcd(p,(p′−p)/2) = 1, SM(p,p′) has c = 3/2 − 3(p′−p)²/(pp′) and, for
1 ≤ ρ ≤ p−1, 1 ≤ s ≤ p′−1, ρ+s even (NS sector),

    χ^{NS}_{ρ,s}(q) = q^{(ρp′−sp)²/(8pp′) − 1/16} · ((−q^{1/2};q)_∞/(q)_∞) · N_{ρ,s},
    N_{ρ,s} := Σ_{k∈Z} ( q^{k(kpp′ + ρp′ − sp)/2} − q^{(kp+ρ)(kp′+s)/2} ),

with χ_{ρ,s} = χ_{p−ρ,p′−s} (Goddard–Kent–Olive; Meurman–Rocha-Caridi; this is the convention under which all the
N = 1 formulas of Sun–Wang were found to hold, see `SUN-WANG-AUDIT.md`). h − c/24 = (ρp′−sp)²/(8pp′) − 1/16, so
the effective weight is ((ρp′−sp)² − d_0²)/(8pp′) with d_0 the minimum of |ρp′−sp|. For (p,p′) = (4,4r+6):
(ρp′−sp)/2 = ρ(2r+3) − 2s, pp′/2 = 8r+12, d_0 = 2, and for ρ = 1

    N_{1,s} = Σ_k ( q^{(8r+12)k² + (2r+3−2s)k} − q^{(8r+12)k² + (2r+3+2s)k + s/2} ),     h_eff(1,s) = ((2r+3−2s)² − 1)/(16(2r+3)).

*Proof of Corollary C.* Let j = 2r+3−2i and Θ_j = Σ_n (−1)^n q^{((2r+3)n²−jn)/4}, so that S_r^{(i)} = ((−q^{1/2})_∞/(q)_∞)Θ_j
(Theorem C for i = r+1, Theorem C′ for i ≤ r). Split n = 4k + a, a ∈ {0,1,2,3}; the four pieces are, with E = 8r+12,

    T_0 = Σ_k q^{Ek² − jk},              T_1 = −Σ_k q^{Ek² + (4r+6−j)k + (2r+3−j)/4},
    T_2 = Σ_k q^{Ek² + (E−j)k + 2r+3 − j/2},     T_3 = −Σ_k q^{Ek² + (12r+18−j)k + (18r+27−3j)/4}.

If i is odd: with s = i, 2r+3−2s = j and 2r+3+2s = 4r+6−j, s/2 = (2r+3−j)/4, so (k → −k in the first sum)
T_0 + T_1 = N_{1,i}. With s″ = 4r+6−i = (6r+9+j)/2: replacing k by k−1 in T_3 and k by −k−1 in T_2 gives
T_3 = −q^{(2r+3+j)/4} Σ_k q^{Ek² − (4r+6+j)k} and T_2 = q^{(2r+3+j)/4} Σ_k q^{Ek² + (E+j)k + s″/2}, i.e.
T_2 + T_3 = −q^{(2r+3+j)/4} N_{1,4r+6−i}. If i is even: with s = 2r+3−i = (2r+3+j)/2, T_0 + T_3 = N_{1,2r+3−i} (same shift of
T_3), and with s = 2r+3+i = (6r+9−j)/2, T_1 + T_2 = −q^{(2r+3−j)/4} N_{1,2r+3+i} (k → −k−1 in T_2). In both cases the
labels are NS labels (1 + s even), and the powers of q are exactly the differences of the exponents
(2r+3−2s)²/(16(2r+3)) of the two characters; the common prefactor is q^{j²/(16(2r+3)) − 1/16}. The weights follow from
h_eff(1,s). For i = r+1: j = 1, weights 0 and ((4r+6±1)²−1)/(16(2r+3)) = (r+2)/2 (r even), (r+1)/2 (r odd). ∎

*Attributions for (T1,C_r).* Theorem C is eq. (5.14) of S. O. Warnaar, *The generalized Borwein conjecture. II.
Refined q-trinomial coefficients*, Discrete Math. 272 (2003) 215–258 (= (5.15) in arXiv:math/0110307v2, where it
is obtained as the L, M → ∞ limit of a polynomial identity proved by his q-trinomial transformations; "For k = 1
this is … entry (20) in Slater's list"). It is quoted as (2.8) in Sun–Wang, and as (3.98) in Wang–Wang
arXiv:2402.06253, who also state and prove the companions (their Theorem 4.1: "The following result was
discovered by Warnaar (unpublished)"). Mizuno (arXiv:2305.02267) and Wang–Wang treat r = 2, 3 in full. The proof
in §3 (Bailey chain of half steps on Slater's G(1)) is the natural one and is presumably known to the experts
(Warnaar remarks that his (5.13) "is a refinement of the Bailey pair G(1)"); I did not find it written out.
Corollary C for general r and i is an immediate consequence that I did not find stated (Sun–Wang state i = r+1
"checked up to r = 7" and the cases r = 2, 3).

## 5. Proof of Theorem B; relation to Du–Ji–Shen–Xu

*Proof of Theorem B.* Write N′_j = n_j + … + n_{r−1} (1 ≤ j ≤ r−1), so M_j = 2N′_j + n_r, and Q := q².

*Even part.* n_r = 2μ. Put L_j = N′_j + μ (j < r), L_r = μ: L_1 ≥ … ≥ L_r ≥ 0, M_j = 2L_j for all j,
(q²;q²)_{n_j} = (Q;Q)_{L_j−L_{j+1}}, and

    R^{even}_r = Σ_{L_1≥…≥L_r≥0} Q^{½(L_1²+…+L_r²)} / ( (Q)_{L_1−L_2}⋯(Q)_{L_{r−1}−L_r} · (q;q)_{2L_r} ).

This is the left side of (1.2) in base Q with a = 1, y = Q^{1/2} = q, s = r, β̃_L = 1/(q;q)_{2L}, i.e.
β_L = 1/((−q;q²)_L(q;q)_{2L}) = 1/((q²;q⁴)_L(q²;q²)_L) (use (q;q)_{2L} = (q;q²)_L(q²;q²)_L): the pair (P2) with q replaced by q².
By Lemma 2,

    R^{even}_r = ((−q;q²)_∞/(q²;q²)_∞) ( 1 + Σ_{m≥1} q^{4rm²} (−1)^m q^{6m²−2m}(1+q^{4m}) ) = ((−q;q²)_∞/(q²;q²)_∞) Σ_{m∈Z} (−1)^m q^{(4r+6)m²−2m}.

*Odd part.* n_r = 2μ+1. Put K_j = N′_j + μ (j < r), K_r = μ; then M_j = 2K_j+1, ¼ΣM_j² = r/4 + Σ_j K_j(K_j+1), and

    R^{odd}_r = q^{r/4} Σ_{K_1≥…≥K_r≥0} Q^{Σ_j K_j(K_j+1)/2} / ( (Q)_{K_1−K_2}⋯(Q)_{K_{r−1}−K_r} · (q;q)_{2K_r+1} ),

the left side of (1.2) in base Q with a = Q, y = Q = q², s = r, β̃_K = 1/(q;q)_{2K+1}, i.e.
β_K = 1/((−q²;q²)_K (q;q)_{2K+1}) = 1/((q;q²)_{K+1}(q⁴;q⁴)_K): the pair (P3) with q replaced by q². By Lemma 2 (aQ = q⁴),

    R^{odd}_r = q^{r/4} ((−q²;q²)_∞/(q⁴;q²)_∞) Σ_{n≥0} q^{rn(n+1)} (−1)^{n(n+1)/2} q^{n(3n+1)/2} (1+q^{2n+1})/(1−q²)
              = q^{r/4} ((−q²;q²)_∞/(q²;q²)_∞) Σ_{n≥0} (−1)^{n(n+1)/2} q^{((2r+3)n²+(2r+1)n)/2} (1 + q^{2n+1}).

The term with q^{2n+1} is the term of index −n−1 of the bilateral series in (0.3) ((2r+3)(n+1)² − (2r+1)(n+1) =
(2r+3)n² + (2r+5)n + 2 and (−1)^{(−n−1)(−n)/2} = (−1)^{n(n+1)/2}). Splitting n = 2m and n = 2m−1 gives
Σ_m (−1)^m q^{(4r+6)m²+(2r+1)m} + q Σ_m (−1)^m q^{(4r+6)m²−(2r+5)m}, and the Jacobi triple product gives the product form. ∎

*Remarks.* (1) r = 1: R_1 = Σ q^{n²/4}/(q)_n; (0.2), (0.3) are Slater-type evaluations of Σ q^{m²}/(q)_{2m} and
Σ q^{m²+m}/(q)_{2m+1} (Rogers). r = 2: R_2 = S_2 (B_2 = C_2), and Theorem B is the dissection of (0.1) for r = 2
into integral and half-integral powers. r = 3: (0.2) is eq. (3.67) of Wang–Wang arXiv:2407.21725 (their Example 7
= (T1,B3)), J_2J_{16,36}/(J_1J_4).

(2) *Du–Ji–Shen–Xu.* In J. Q. D. Du, K. Q. Ji, E. Y. Y. Shen, C. X. Y. Xu, arXiv:2602.14483, Theorem 1.5 (1), (3)
(identities (1.7), (1.9)) the same two theta quotients occur, for generalized Nahm sums of rank ρ with matrices
(1.4), (1.5) there: (1.7) at rank ρ has the right side of Theorem B for r = 2ρ−2, and (1.9) (linear term −n_ρ/2)
that for r = 2ρ−3. For even r = 2s this is explained by the fusion lemma in base q² (Lemma F of
`KRTW-k2.md`): their matrix (1.4) of rank s+1, on the species (1, 2, 4, …, 2s−2, r), becomes A(T1,B_{2s})D after the
s−1 fusions "node 2j+1 := node 1 + node 2j" (j = s−1, …, 1) — exact check for s ≤ 8 in the certificate, part (d).
So for even r Theorem B also follows from their (1.7) and the fusion lemma. For odd r the equality of the
(T1,B_{2ρ−3}) sum with their sum (1.9) follows from Theorem B and (1.9) (checked numerically for ρ ≤ 5), but I have
no direct transformation. Their theorem has r+1 companions (1.8), (1.10); by fusion these give companions of the
(T1,B_{2s}) sum, which I did not work out.

(3) The same method gives nothing new for (A1,B_r) and (A1,C_r): (A1,C_r) is Bressoud's identity ((2.5) in
Sun–Wang; full steps q^{n²} on the Bailey pair with β_n = 1/(q²;q²)_n), and (A1,B_r) is
Baker–Kanade–Russell–Sadowski, Algebr. Comb. 6 (2023), Theorem 3.4 ((2.4) in Sun–Wang); both were re-tested
(`SUN-WANG-AUDIT.md`).

## 6. Verification (`certificates/sw_certify_T1.txt`)

(a) Lemma 3: the three Bailey pairs, exactly (polynomial identities in p = q^{1/2} after clearing denominators),
for n ≤ 14. (b) Lemma 1 applied to (P1) and to (P3): the resulting pair satisfies (1.1), as power series mod
p^{400}, n ≤ 6. (c) Theorem C, Theorem C′ and both forms of Corollary C for all 1 ≤ i ≤ r+1, r ≤ 9, mod q^{100}
(r ≤ 6) / q^{60} (r = 7, 8, 9); the labels and effective weights are checked against the list of weights of
SM(4,4r+6). (d) Theorem B: total, even part, odd part and the product form, r ≤ 9, same orders; Du–Ji–Shen–Xu
(1.7), (1.9) as printed and their equality with R_{2ρ−2}, R_{2ρ−3}, ρ ≤ 5, mod q^{60}; the exact fusion chain.
Independent of this file, `certificates/sw_audit.txt` contains (2.8) and (3.6)–(3.13) of Sun–Wang tested in
their printed form (r ≤ 9, up to q^{200}).
