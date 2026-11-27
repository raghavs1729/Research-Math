# R07. The game CRIM on integer partitions: complete solution of the losing positions; three conjectures PROVED

**Source.** *CRIM: A Natural Game on Integer Partitions*, arXiv:2606.16828 (15 June 2026).
CRIM (Column-Row Impartial Merge): a position is a Young diagram; "a move consists of removing an
arbitrary row or column", the remaining pieces being reattached; normal play (the player who cannot
move loses; the empty partition is the terminal position). `𝒫` = losing positions (previous player
wins), `𝒩` = winning positions. Rank of `λ = [λ_1,…,λ_r]`: `|r − λ_1|`.

**Printed statements settled here.**

* Conjecture `conj:main` (the main conjecture; abstract: "We conjecture that every losing partition
  has even rank"): "Every `𝒫`-position is of even rank." Authors' evidence: partitions of order at most
  50 and several special families.
* Conjecture (Section "Understanding operations on partitions"): "A meld of losing thick hooks is
  losing." (meld: `μ ⋄ λ = [μ_1+λ_1,…,μ_1+λ_m, μ_1,…,μ_k]`; thick hook `Γ^{a,b}_{r,c} = [c^a, b^{r−a}]`,
  `r > a ≥ 2`, `c > b ≥ 2`.)
* Conjecture `conj:bloby`: "Let `λ` and `Γ_{r,c}` be losing and suppose `λ* = Γ_{r,c} ◁ λ` is defined.
  Then `λ* ∈ 𝒫`." (`Γ_{r,c} = [c, 1^{r−1}]` hook; `◁` = stack.)
* Question: "Which partitions consisting of four parts are losing?" (answered below.)
* Conjecture `conj:pself`: "A parity self-conjugate partition is in `𝒩` if and only if it is the stack
  of odd self-conjugate hooks and `[1]`." — only the "if" direction is proved here.

## Encoding
The boundary word `w(λ) ∈ {E,N}^*` of `λ = [λ_1 ≥ … ≥ λ_r > 0]` is
`E^{λ_r} N E^{λ_{r−1}−λ_r} N ⋯ E^{λ_1−λ_2} N` (`w([ ]) = ε`). It starts with `E` and ends with `N`
("canonical"), has `r` letters `N` and `λ_1` letters `E`; so rank `≡ |w(λ)| (mod 2)`. Removing the
`i`-th row deletes the corresponding `N`; removing the `j`-th column deletes the `j`-th `E`; afterwards
leading `N`'s (rows of length 0) and trailing `E`'s (empty columns) are discarded. Thus **a move =
delete one letter, then `canon`**, where `canon(v)` strips leading `N`'s and trailing `E`'s, and every
single-letter deletion is a legal move. Meld is concatenation: `w(μ ⋄ λ) = w(μ)w(λ)`.

## Theorem 1 (solution of CRIM)
Let `B` be the deterministic automaton with states `0,…,16`, start state `0`, and the table

| state | a word reaching it | accepting | on `E` | on `N` |
|---|---|---|---|---|
| 0 | `ε` | yes | 1 | 0 |
| 1 | `E` | yes | 2 | 3 |
| 2 | `EE` | yes | 4 | 5 |
| 3 | `EN` | no | 6 | 7 |
| 4 | `EEE` | yes | 2 | 8 |
| 5 | `EEN` | no | 8 | 2 |
| 6 | `ENE` | no | 9 | 10 |
| 7 | `ENN` | no | 8 | 3 |
| 8 | `EEEN` | no | 5 | 5 |
| 9 | `ENEE` | no | 11 | 12 |
| 10 | `ENEN` | yes | 1 | 13 |
| 11 | `ENEEE` | no | 9 | 14 |
| 12 | `ENEEN` | no | 15 | 9 |
| 13 | `ENENN` | no | 15 | 10 |
| 14 | `ENEEEN` | yes | 16 | 12 |
| 15 | `ENEENE` | no | 12 | 12 |
| 16 | `ENEEENE` | yes | 14 | 9 |

A partition `λ` is a losing position of CRIM if and only if `B` accepts `w(λ)`.

*Proof.* Write `B(v)` for "`B` accepts `v`". Three properties of the table (checked mechanically,
`code/crim_automaton.py`, output `cert/crim_automaton_cert.txt`):

* (i) `δ(0,N) = 0` and state 0 is accepting; (ii) for every state `q`, `δ(q,E)` is accepting iff `q` is.
  Hence `B(v) ⇔ B(canon(v))` for every word `v`, and `B(ε)` holds.
* (iv) For every word `w` that starts with `E` and ends with `N`:
  `B(w)  ⇔  for no position i does B accept w with its i-th letter deleted`.
  This is verified for all (infinitely many) such `w` by exhausting the reachable configurations
  `(δ*(0,w), {δ*(0, w with one letter deleted)})`: the configuration of `wx` is determined by that of
  `w` (`q ↦ δ(q,x)`, `T ↦ δ(T,x) ∪ {q}`); 169 configurations are reachable from `w = E`, and the
  equivalence holds at each of them whenever the last letter is `N`.

Now induct on `|w|` for canonical `w`: `ε` is losing and accepted. For `w ≠ ε`, the options of `w` are
`canon(w − letter i)`, shorter canonical words; `w` is losing iff no option is losing iff (induction)
no option is accepted iff (i, ii) no `w − letter i` is accepted iff (iv) `B(w)`. ∎

The automaton was found by computing Nerode classes from game values of words of length `≤ 16`; the
proof does not depend on how it was found.

## Corollaries (each checked on the automaton, exactly)

1. **`conj:main` is true: every `𝒫`-position has even rank.** In the exhaustive configuration search
   the parity of `|w|` is carried along; no accepted canonical word has odd length
   (`cert/crim_automaton_cert.txt`, line (v)).
2. For even `L ≥ 8` exactly `5·2^{L−5}` of the `2^{L−2}` partitions with `rows + columns = L` are
   losing (proportion 5/8); for `L = 4, 6` the numbers are 2 and 11. Proof: the number `c_L` of accepted
   canonical words of length `L` satisfies a linear recurrence of order `≤ 17` (transfer matrix of `B`),
   `5·2^{L−5}[L even]` one of order 2, so their difference satisfies one of order `≤ 19`; it vanishes for
   all `8 ≤ L ≤ 80` (`cert/crim_counts.txt`), hence for all `L ≥ 8`. (Brute force agrees for `L ≤ 28`.)
3. `δ*(q,E^k) = δ*(q,E^{k+2})` and `δ*(q,N^k) = δ*(q,N^{k+2})` for all states `q` and all `k ≥ 2`.
   Hence the status of `λ` depends only on the classes "1", "even", "odd ≥ 3" of the run lengths of
   `w(λ)`, i.e. of `λ_r`, the differences `λ_i − λ_{i+1}` and the multiplicities of the parts.
4. **"A meld of losing thick hooks is losing" is true.** Losing thick hooks are those with `r ≡ c (mod 2)`
   (the source's Theorem `thm:thick-hook`; reproduced by `B`); their words are
   `E^b N^{r−a} E^{c−b} N^a`. Every such word maps the start state to state 2 and state 2 to state 2
   (all exponent classes enumerated, using 3); state 2 is accepting (`code/crim_conj.py`,
   `cert/crim_conj.txt`).
5. **`conj:bloby` is true.** `w(Γ_{r,c} ◁ λ) = E N^a w(λ) E^b N` with `a = r−1−k`, `b = c−1−λ_1`
   (`k` = number of parts of `λ`); the hook is losing iff `r, c ≥ 2` and `r ≡ c (mod 2)` (source, Lemma
   `lem:hook`), and `λ` losing implies `k ≡ λ_1 (mod 2)` (Corollary 1), so the hypothesis is `a ≡ b (mod 2)`.
   For every class of `a`, a product search over all canonical `w` (pairs of states of `B` started at
   `0` and at `δ*(0, E N^a)`) over the **nonempty** canonical `w` shows: whenever `B(w)`,
   `B(E N^a w E^b N)` for every `b ≡ a (mod 2)` (`code/crim_conj.py`). For the empty partition `λ = [ ]`
   the stack is the hook `Γ_{r,c}` itself, which is losing by hypothesis. (The parity reduction uses
   `λ ≠ [ ]`, so that `k, λ_1 ≥ 1` and `r, c ≥ 2` automatically.)
6. **Four parts.** `λ = [λ_1,λ_2,λ_3,λ_4]` (`λ_4 ≥ 1`) is losing iff the quadruple of classes of
   `(λ_4, λ_3−λ_4, λ_2−λ_3, λ_1−λ_2)` (classes `0`; `1`; `e` = even `≥ 2`; `o` = odd `≥ 3`) is one of the 69
   quadruples listed in `cert/crim_conj2.txt` (of 192).
7. **`conj:pself`, "if" direction.** Every stack of odd self-conjugate hooks and `[1]` is a winning
   position: these are the words obtained from `EN` by `w ↦ E N^a w E^a N`, `a` odd; the closure of the
   corresponding state transformers (`a ∈ {1, 3}` suffices by 3) has 3 elements, none accepting from
   the start state (`code/crim_conj2.py`). The converse (every other parity self-conjugate partition
   is losing) is NOT proved; no counterexample among the 2564 parity self-conjugate partitions of
   `n ≤ 52` (authors: `n ≤ 50`).

## Independent checks (not needed for the proof)
* `code/crim_check.py`: the automaton agrees with a brute-force search that removes rows/columns of
  partitions represented as tuples, for all 66 273 partitions of `n ≤ 34`; it reproduces the source's
  Proposition `thm:2parts` (`λ_1 < 80`), Theorem `thm:3parts` (`λ_1 < 50`; the literal statement
  wrongly includes the column `[1,1,1]`), and accepts the 19 "unexplained losing partitions" of the
  source's appendix that were checked.
* `code/crim.c`: brute-force bitset computation over all boundary words of length `≤ 28`
  (`2^{27}` positions): no losing position of odd length; counts `5·2^{L−5}` for even `8 ≤ L ≤ 28`.

## Not settled
The Sprague–Grundy conjectures for rectairs and staircases, `conj:hookswap` (Conway pairs) and the
"only if" direction of `conj:pself`. The normal-play Grundy function and the misère classes do not
show a small Nerode index (186 resp. 317 classes and growing at prefix length 9), so the method does
not apply directly.

## Novelty (search of 26 Nov 2026, `NOVELTY-2.md`, item 1): NOT NEW
* Zijian Zeng, *A Boundary-Word Classification of Losing Positions in CRIM* (dated 29 Aug 2026; public
  in the GitHub archive `ArtificialZeng/AI-Has-Taste` since 16 Sep 2026; not on arXiv; unrefereed)
  proves: with `w = E^{a_1}N^{b_1}⋯E^{a_k}N^{b_k}`, `L` = number of leading pairs with `a_i = 1`, `b_i` odd,
  `R` = number of trailing pairs with `a_i` odd, `b_i = 1`, `M = min(k, L+R)`: `w` is losing iff `|w|` and
  `M` are even. This gives the main conjecture. Cross-check: this criterion and the automaton of
  Theorem 1 agree on all 2 097 151 canonical words of length `≤ 22`
  (`code/crim_zeng_compare.py`, `cert/crim_zeng_compare.txt`).
* Lean proofs of the main, stack (`conj:bloby`), meld and parity self-conjugate conjectures and of a
  misère classification are claimed in open, unreviewed pull requests #73 and #76 of
  `umaia1234/agentic-conjectures` (25 Oct 2026); the four-part classification is in
  `the-omega-institute/trureturing` PR #26523 (5 Nov 2026).
So Theorem 1 and its corollaries are an independent confirmation (by a different, automaton-based
proof), not a first proof. None of the third-party proofs was read here.

## Addendum: the printed Sprague–Grundy conjecture for `R^k_{r,r−1}` has the wrong parity
`code/crim_sgall.c` computes the Grundy values of all partitions with `rows + columns ≤ 30`
(`cert/crim_sg_conjectures.txt`). The source's conjecture "For `r ≥ 7`: `sg(R^k_{r,r−1}) = 3` if `k = r−2`
and `r` odd, `1` otherwise" fails for every `7 ≤ r ≤ 15` at `k = r−2` (values 1, 3, 1, 3, … for
`r = 7, 8, 9, 10, …`): the value is 3 iff `r` is **even**. (`R^{r−2}_{r,r−1}` is the padded staircase
`ps_{r−1}`, and the source's own conjecture on `ps_n` gives 3 for odd `n = r−1`; so this is an internal
inconsistency.) The conjectures on `R^k_{r,r}` (apart from the edge case `r = 1`, already noted by
others), on staircases `s_n` (`n ≤ 15`) and on padded staircases `ps_n` (`n ≤ 14`), and the corrected
parity for `R^k_{r,r−1}` (`r ≤ 15`) survive.
