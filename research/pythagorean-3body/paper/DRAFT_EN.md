# A computer-assisted proof of the final evolution of the Pythagorean three-body problem

> **Preprint (version 1.2, 2026-10-13; version 1: 2026-10-12).** Internally reviewed; not yet externally
> peer reviewed and not formally (machine-)verified. Section 8 lists the remaining assumptions. Released on
> Zenodo (DOI in Section 9); not yet submitted to a journal. **Change in version 1.2:** the classical
> collision facts that version 1 cited for global existence through binary collisions are replaced by a
> self-contained proof for the planar case (Appendix A.2); see Section 8, item 2.
>
> Jeomseon (Independent researcher), jeomseon@proton.me
>
> The use of AI tools is described in the section "Use of AI tools". Every number is traceable to a file in
> the reproduction package (Section 9).

---

## Abstract

In the Pythagorean (Burrau) three-body problem, point masses 3, 4 and 5 start at rest at the vertices of a
3-4-5 right triangle, each opposite the side of its own length. Numerical studies since Szebehely and Peters
(1967) agree that the motion ends with the mass-3 body escaping and the masses 4 and 5 forming a bound
binary. We describe a computer-assisted argument for this final evolution. A validated Taylor-series
integrator with outward-rounded dyadic interval arithmetic (128-bit mantissas, order 20) and a Lohner–QR
parallelepiped representation with a width-weighted pivoting rule encloses the exact solution on
t ∈ [0, 61.0055] in 31,811 accepted steps. From the stored enclosures at t = 59.50, …, 61.00 we check a
sufficient escape condition (★) — Yoshida's (1972) Theorem 1 with d = m₁m₂/K, which is also a special case of
the corrected form of his Corollary that Yoshida himself published (Erratum 1973; Criterion A′, 1974) — using two independently written exact/interval
checkers. Both certify (★) on every stored box from t = 59.73 to t = 61.00 and on none before. Together with a
hand proof of (★) and of the two potential-bound lemmas it needs, this yields that the mass-3 body escapes
with ρ(t) ≥ ρ₀ + C(t − t₀) for some C > 0 and the 4–5 separation stays below 1200/769 for all later times,
where binary collisions (if any) are continued by Levi-Civita regularization. We report failed attempts and
lower precisions (112 and 96 bits do not reach t = 61) and state the limits of the verification: one-ulp
rounding faults in the interval layer are excluded only by code reading, and no formal verification or
external review has yet been done. Global existence through binary collisions is proved in an appendix
from standard ODE theory (planar case). To our knowledge, no
rigorous proof of this final evolution has been published; our literature search (Section 1.2) cannot
exclude one.

<!-- src: results_r4/r3c_final_128w.json (t_final 2047006249/33554432, steps_accepted 31811);
     Y1_FINAL_JUDGMENT_KO.md; yoshida_check/Y1_RESULT_KO.md §1, §6; LITERATURE_RECHECK_KO.md §3, §5;
     server_archive_20261012/r4_out.tgz → r3c_final_112w.json, r3c_final_96w.json -->

---

## 1. Introduction

### 1.1 The problem

Three bodies with masses m₀ = 3, m₁ = 4, m₂ = 5 (gravitational constant G = 1) start at rest at

  q₀ = (1, 3), q₁ = (−2, −1), q₂ = (1, −1),

so that the mutual distances are |q₀q₁| = 5, |q₀q₂| = 4, |q₁q₂| = 3. The centre of mass is at the origin, the
total momentum and angular momentum are zero, and the energy is exactly

  H = −(3·4/5 + 3·5/4 + 4·5/3) = −769/60.

<!-- src: integrator_review/INTEGRATOR_REVIEW_KO.md §2 row m (initial values); YOSHIDA_REDERIVATION_KO.md §2 (energy) -->

The problem is usually attributed to Burrau (1913). Szebehely and Peters (1967) integrated it numerically with Levi-Civita
regularization of close encounters and reported that the bodies 4 and 5 end as a binary while body 3 recedes
indefinitely. <!-- src: SURVEY_KO.md table 1967; PRIOR_ART_RIGOROUS_KO.md table row Szebehely & Peters -->
Yoshida (1972) applied an analytic hyperbolic–elliptic escape criterion at t = 59.72 to numerical data that
he attributes to Yuasa (1971, private communication; Yoshida 1972, p. 406 — no publication of these data was
found in Crossref or a web search). Later studies are numerical: Aarseth et al. (1994) on the chaotic
neighbourhood of the initial condition; Hirayama, Nishikawa and Shiraishi (2006) and Hirayama (2019) with
high-order Taylor series methods, the latter in quadruple precision (title and bibliographic record checked;
text not read); and the arbitrary-precision Bulirsch–Stoer code Brutus. Boekholt and Portegies Zwart (2015)
needed a tolerance of 10⁻¹⁴ and an 88-bit word length for a solution converged to three decimal places and
report dissolution at about t = 60. Portegies Zwart and Boekholt (2018) used a word length of 136 bits and a
tolerance of at most 10⁻²⁴ so that the time-reversed solution returns to the initial condition, and report a
dissolution time t ≃ 63.4 defined as the time at which the linear behaviour of the phase-space separation
overtakes its exponential behaviour (a fitted model); they note that other escape criteria give different
times. This time is therefore not comparable with the first time t ≈ 59.73 at which our sufficient condition
(★) holds. These are converged numerical solutions, not error-bounded ones. <!-- src: LITERATURE_RECHECK_KO.md §5
     table (which had attributed the 136-bit / 1e-24 / 63.4 values to the group's papers jointly); corrected
     2026-10-12 after an external review comment, by reading arXiv:1411.6671 §3.1 and arXiv:1802.00970 via
     https://ar5iv.labs.arxiv.org/html/1411.6671 and /html/1802.00970 -->

### 1.2 Prior rigorous work

Validated (interval) integration with Taylor methods and Lohner's QR wrapping-control is standard
(Lohner 1988, 2001; Nedialkov, Jackson and Corliss 1999; Zgliczyński 2002; the CAPD library, Kapela et
al. 2021), and computer-assisted proofs exist for other three-body questions — e.g. choreographies and the
figure-eight orbit (Kapela and Zgliczyński 2003; Kapela and Simó 2007) and collision/near-collision orbits of the restricted problem (Capiński, Kepley and
Mireles James 2023). None of the components of our integrator is new.

We searched ADS (abstracts and the citation list of Yoshida 1972), scite, the arXiv API, zbMATH Open (via
its public API) and the general web (English, Japanese, Russian), with both supporting and opposing queries,
for a rigorous proof of the Pythagorean final evolution. **We did not find one.** zbMATH Open returned, for
this problem, only numerical studies (Aarseth et al. 1994; Janssens 2011; Rose and Dullin 2013; Hirayama
2019), and for interval methods in the three-body problem only works on other questions. This is a statement
about our search, not a proof of absence: MathSciNet (subscription) was not searched, nor were thesis
repositories, and validated-numerics groups were not contacted. **To our knowledge, no rigorous proof of the final evolution of the Pythagorean problem has been
published.** <!-- src: LITERATURE_RECHECK_KO.md §5, §6, §7; LITERATURE_RECHECK_ZBMATH_KO.md;
     PRIOR_ART_RIGOROUS_KO.md §6 item 1. Janssens 2011 and Rose–Dullin 2013 listed from the zbMATH note only,
     not cited in the reference list. -->

### 1.3 On Yoshida's Corollary

Yoshida's (1972) Corollary (p. 402) states an escape criterion without the hypothesis r(0) ≤ d, but its proof
uses that hypothesis in its first step. **Yoshida himself corrected this, first in an Erratum (PASJ 25,
p. 285, 1973)**, which replaces condition (26) of the Corollary by d ≥ r(0), and again in 1974 (PASJ 26,
p. 369, Note, calling the 1972 conditions (25), (26) "incorrect" and replacing them by ρ(0) ≥ ad ≥ a r(0),
"Criterion A′"). The same erratum also corrects inequality (40) on p. 405 (Theorem 2, which we do not use).
We did not discover this gap; it was corrected by the author. Our project first noticed the gap
independently (a state at t ≈ 17.118 on the Pythagorean trajectory satisfies the uncorrected hypotheses while
violating the conclusion), and only afterwards found the author's corrections (the 1973 erratum was found
last, while checking the references of this paper). Yoshida's own Pythagorean
application (t = 59.72, r = 0.42 < d₁) satisfies the corrected conditions, so his 1972 numerical conclusion is
not affected; it rests, however, on 1971 floating-point data without error bounds. We use only the corrected
form (★) below, which requires r < d₁ explicitly.
<!-- src: LITERATURE_RECHECK_KO.md §2 rows p.402, 1974 p.369, §3; YOSHIDA_REDERIVATION_KO.md §3 (counterexample);
     Erratum: Crossref record https://api.crossref.org/works/10.1093/pasj/25.2.285 and ADS scan
     1973PASJ...25..285Y (read on screen 2026-10-12 via https://ui.adsabs.harvard.edu/scan/manifest/1973PASJ...25..285Y) -->

### 1.4 Contribution and structure

We give (i) an interval enclosure of the exact solution on [0, 61.0055] (Sections 2, 4); (ii) a proof of the
sufficient escape condition (★) with explicit lemmas (Section 3, Appendix A); (iii) an exact check of (★) on
the enclosures by two independent programs (Section 4.2); (iv) a description of how the code was checked
(Section 5) and of what remains unchecked (Section 8). All numbers, code hashes and run records are in a
reproduction package (Section 9, `REPRODUCE.md`).

---

## 2. Method

### 2.1 Interval arithmetic on dyadic numbers

All rigorous quantities are closed intervals [a, b] whose endpoints are dyadic rationals m·2^e with
|m| < 2^p, where p is the working precision (128 bits in the main run). Each arithmetic operation is
evaluated exactly in rational arithmetic (Python integers and `Fraction`) and the result is rounded
**outward**: the lower endpoint down and the upper endpoint up to p bits. Square roots use integer square
roots (`math.isqrt`) rounded outward. A single precision is fixed for the whole run; every state endpoint is
checked to be dyadic with mantissa ≤ p bits at each step (otherwise the step stops with an error). This
"single precision" rule was adopted after the interval layer's independent checker showed that mixed-precision
inputs make `intersect` slightly widen its result (an issue of contract, not of containment; 0 containment
violations). <!-- src: integrator_review/INTEGRATOR_REVIEW_KO.md §2 rows a2, h1, h2, l; rigorous/R0_DECISION_KO.md; rigorous/R0_FINAL_RUN_KO.md -->

### 2.2 Taylor coefficients

The equations of motion in Cartesian coordinates, q̈_a = Σ_{b≠a} m_b (q_b − q_a)/|q_b − q_a|³, are expanded
in Taylor series by automatic differentiation; the factor |r|^(−3) = s^(−3/2) with s = |r|² uses the
power recursion k s₀ w_k = Σ ((α+1)i − k) s_i w_{k−i}, α = −3/2. The order is p = 20. The Taylor
remainder is enclosed by evaluating the (p+1)-th coefficient over an a-priori enclosure of the trajectory.
<!-- src: integrator_review/INTEGRATOR_REVIEW_KO.md §2 row m, row i; rigorous/r3c_run.py header J1 -->

### 2.3 One step: a-priori enclosure, Lohner–QR update

The state set is a parallelepiped X_k = { c + B r : r ∈ R }, with centre c (exact dyadic vector), frame B
(exact dyadic 12×12 matrix) and interval vector R. One step of size h does:

1. **A-priori (Picard) enclosure.** Find a box P with X̂ + [0, h]·F(P) ⊆ P, where X̂ is the interval hull of the
   enclosure box of X_k and c. Successful evaluation of F(P) requires every pairwise squared distance to have
   a positive lower bound, so **P contains no collision configuration**; a Banach fixed-point argument then
   shows every solution starting in X̂ stays in P on [0, h]. If the check fails, h is halved.
2. **Centre image.** T(c) + Z, with T(c) the order-p Taylor polynomial at c and
   Z = h^{p+1} [x_{[p+1]}](P).
3. **Jacobian enclosure.** [J] ∋ Dφ_h(y) for all y ∈ X̂, from the Taylor coefficients of the variational
   equation (forward-mode differentiation in 12 directions over X̂) plus a remainder term using an a-priori
   enclosure W of the variational matrix. (The code review found that this remainder term is redundant for
   the correctness of the final set, because Z is already evaluated over P; it only makes the set slightly
   larger.)
4. **New centre** c′ = the p-bit floor of the midpoint of T(c) + Z.
5. **New frame.** B′ = Q from a column-pivoted Householder QR factorisation (double precision) of
   A = mid([J])·B; Q is converted exactly to a dyadic matrix. **Its inverse is computed exactly** in rational
   arithmetic (Gauss–Jordan on `Fraction`s) and then enclosed outward. B′ need only be invertible for
   correctness; orthogonality only affects tightness.
6. **New label** R′ = [B′⁻¹]([J]B) R + [B′⁻¹](T(c) + Z − c′).

By the mean-value theorem on the convex set X̂ and inclusion monotonicity, φ_h(X_k) ⊆ X_{k+1}. The step also
returns, for any 0 < τ ≤ h, an enclosure of φ_τ(X_k) ("in-step box"). This is the classical Lohner scheme
in the form of Nedialkov–Jackson–Corliss (1999) and the C¹-Lohner literature (e.g. Walawska and Wilczak 2016);
we claim nothing new for it.
<!-- src: integrator_review/INTEGRATOR_REVIEW_KO.md §1, §2 rows a1–a3, b, c1–c4, d1, e, f, g, §3-2, §5;
     rigorous/lohner_step.py lines 660–695 (read for this draft) -->

### 2.4 Width-weighted pivoting

With the textbook pivoting (largest remaining column norm of A) the run failed at the first close encounter
(Section 6). We instead choose, at each QR stage, the remaining column maximising
(remaining column 2-norm) × w(R_j), the column norm weighted by the width of the corresponding label
component. Lohner recommends a *pivoted* QR factorisation for this purpose (as reported by Revol 2023,
citing Lohner 2001). We have not verified whether Lohner's original pivot uses the same width weighting. Since correctness only needs B′ invertible (step 5),
the pivot rule affects tightness only.
<!-- src: rigorous/lohner_step.py _qr_colpivot and lines 674–678; rigorous/R3C_RESULT_KO.md line 40;
     GATE_R3CP_PLAN_KO.md §0 -->

### 2.5 Step-size rule

Tolerance tol = 2^(−p+16)·max(1, ‖c‖∞). Trial step h = 2^⌊log₂ v⌋ with v = 0.9 (tol/m)^{1/20}, where m is the
largest absolute endpoint of the order-20 coefficient at the centre; h ≤ 1/8. A step is rejected if the
largest component width of Z exceeds tol, and retried with half the step, at most 20 times. The run stops if
the enclosure width exceeds 1, if h < 2⁻⁶⁰, at t ≥ T_end, at the wall-clock limit, or on any exception or
containment violation. The step rule affects efficiency only; correctness rests on Section 2.3.
<!-- src: rigorous/r3c_run.py header J1–J4 -->

### 2.6 Implementation

Pure Python 3 (integers, `fractions.Fraction`, `math.isqrt`), no external numerical library in the
integrator. Files `dyad.py` (interval layer), `taylor_step.py` (Taylor coefficients, Picard enclosure),
`lohner_step.py` (one Lohner step), `r3c_run.py` (driver), `r4_run.py` (thin wrapper: sets T_end = 61 and
writes the enclosures at t = 59.50, …, 61.00 as exact fractions). A C++/GMP port of the Taylor core exists
and was found bit-identical to the Python code on 23,452 test cases. **The results reported here come from the
Python code**; the C++ core was used only for the full re-run R4R (Section 4.1), which reproduced them bit for
bit. <!-- src: GATE_R4_PLAN_KO.md §3; FC1_FINAL_JUDGMENT_KO.md -->

---

## 3. Escape criterion

### 3.1 Notation

Escaping body m₀ = 3; binary m₁ = 4 ≤ m₂ = 5; M = m₀ + m₁ + m₂ = 12. λ = m₂/(m₁+m₂) = 5/9,
ν = m₁/(m₁+m₂) = 4/9, μ = m₀(m₁+m₂)/M = 9/4. r̄ = q₂ − q₁, r = |r̄|. G = (m₁q₁ + m₂q₂)/(m₁+m₂),
ρ̄ = q₀ − G, ρ = |ρ̄|, and **ρ̇ = (ρ̄ · ρ̄̇)/ρ** (radial velocity, not speed). K = −H = 769/60,
d₁ = m₁m₂/K = 1200/769, λd₁ = 2000/2307. Define

  W(ρ, d) = λ/(ρ + νd) + ν/(ρ − λd).

<!-- src: YOSHIDA_REDERIVATION_KO.md §2; proof_review/PROOF_REVIEW_KO.md §2 -->

### 3.2 Statement

**Theorem (★).** Let H < 0 and m₁ ≤ m₂. If at some time t₀ the three strict inequalities

  (C1) r < d₁,  (C2) ρ > λ d₁,  (C3) ρ̇ > √(2M · W(ρ, d₁))

hold, then for all t ≥ t₀ (with binary collisions of m₁, m₂, if any, continued by regularization):
r(t) < d₁, ρ̇(t) > C := √(ρ̇₀² − 2M W(ρ₀, d₁)) > 0, and ρ(t) ≥ ρ₀ + C(t − t₀). In particular body m₀ escapes
and the binary stays bounded (hyperbolic–elliptic motion).

(★) is Yoshida's (1972) Theorem 1 (p. 400) with d = d₁, in which case its third condition coincides with
(C3). It is also the case d = d₁ of the corrected Corollary (Erratum 1973: condition (26) replaced by
d ≥ r(0)) and of Yoshida's (1974) Criterion A′; A′ formally allows the boundary
case r(0) = d₁, which cannot occur together with (C2), (C3) (Appendix A, Remark A.4.2). The 1974 proof of A′
was not read; we use only our own proof. <!-- src: LITERATURE_RECHECK_KO.md §2 rows p.400, p.403, 1974 p.369;
proof_review/PROOF_STAR_FINAL_EN.md A.0, A.4.2; ADS scan 1973PASJ...25..285Y -->

### 3.3 Lemmas

For λ ≥ ν > 0, λ + ν = 1, ρ > 0, r ≥ 0, c = cos θ ∈ [−1, 1], let
r₁ = √(ρ² + ν²r² − 2νρrc), r₂ = √(ρ² + λ²r² + 2λρrc), P = λ/r₁ + ν/r₂,
F = λ(ρ − νrc)/r₁³ + ν(ρ + λrc)/r₂³ = −∂P/∂ρ, F̂(ρ, d) = λ/(ρ + νd)² + ν/(ρ − λd)².

**Lemma A.** If d ≥ 0, ρ > λd, 0 ≤ r ≤ d, then P(ρ, r, c) ≤ W(ρ, d).
**Lemma B.** Under the same hypotheses F(ρ, r, c) ≤ F̂(ρ, d).

Both are equalities at r = d, c = −1. They correspond to Yoshida's (1972) Lemmas 1 and 2 (n = 1). Our proof
(Appendix A) expands ρP and ρ²F in Legendre polynomials, ρ²F = Σ (n+1) sⁿ βₙ Pₙ(c) with s = r/ρ and
βₙ = λνⁿ + ν(−λ)ⁿ, uses (−1)ⁿβₙ ≥ 0 (this is where λ ≥ ν enters) and |Pₙ| ≤ 1, and then monotonicity in s.
The multipole expansion is a standard technique; we claim no novelty for the lemmas.
<!-- src: proof_review/PROOF_REVIEW_KO.md §4 -->

### 3.4 Proof of (★)

The full proof is given in Appendix A. Outline: a bootstrap argument on the closed set
𝒢 = {t : r ≤ d₁ and ρ̇ ≥ 0 on [t₀, t]}. Lemma B and integration against ρ̇ ≥ 0 give
½ρ̇² − MW(ρ, d₁) ≥ ½C²; Lemma A and the (vector) energy identity give r ≤ m₁m₂/(K + ½μC²) < d₁; continuity
then contradicts sup 𝒢 < T*. On 𝒢 the escaper stays at distance ≥ δ = ρ₀ − λd₁ > 0 from both binary
components. The statement (R) of Appendix A.2 then gives that the solution, continued through binary
collisions by Levi-Civita regularization, exists for all time. (R) is **proved** in Appendix A.2 for the
planar case from standard ODE theory: a separation bound keeps the regularized state in a compact set on
finite time, and each unit of regularized time contributes a positive amount of physical time. (Version 1
cited classical facts here instead; see Section 8, item 2.)
<!-- src: proof_review/PROOF_STAR_FINAL_EN.md A.4; proof_review/s3/S3_REGULARIZATION_PROOF_EN.md (v1.2); proof_review/PROOF_REVIEW_KO.md §5, §8 -->

### 3.5 Why the hypothesis r < d₁ matters (sanity checks)

Dropping any ingredient of (★) gives a criterion that is numerically falsified: with the factor 2 omitted, with
λ and ν interchanged in W, or with (C1) removed, the internal review found states satisfying the modified
condition after which body 3 returns (18, 19 and 1 numerically reliable examples, respectively). With (★)
itself, 800 sampled states (400 random, 400 near the boundary) showed no violation of any conclusion up to
t = 40. These are numerical tests, not proofs; they serve to detect sign or factor errors.
<!-- src: proof_review/PROOF_REVIEW_KO.md §6.3 table and §7 item 5 -->

---

## 4. Results

### 4.1 Enclosure on [0, 61] (run R4-A)

| Quantity | Value |
|---|---|
| Precision, order, pivot | 128 bits, p = 20, width-weighted |
| Final time reached | t = 2047006249/33554432 ≈ 61.005540 (stop reason: target time reached) |
| Accepted steps / attempts / rejected | 31,811 / 41,140 / 9,329 |
| Step sizes | 2⁻²⁵ … 2⁻⁵ |
| Largest component width of final enclosure | 10^−9.545 |
| Wall-clock time | 37,015 s ≈ 10.3 h (one core of a Google Cloud n2-custom-12-12288 instance: 12 vCPU, 12 GiB memory, Debian; Python 3.13.5) |
| Containment cross-checks vs. a 60-digit reference trajectory, grid t = k/100, k = 1…6100 | path boxes P: 6,100 PASS; in-step boxes: 1,138 PASS, 4,955 not decidable (box narrower than 10⁻¹⁵; reference inside); step-end boxes: 7 not decidable (reference inside). **Violations: 0** |

<!-- src: results_r4/r3c_final_128w.json (prec, order, pivot, stop_reason C_T_REACHED, t_final, steps_accepted,
     attempts_total, discarded, h_hist_log2, wmax_final_log10, wall_s, grid_counts, python);
     rigorous/r3c_run.py header J5–J6 (meaning of UNDECIDED); GATE_R3CP_PLAN_KO.md §3 (reference grid 60 digits);
     machine type: coordinator message 2026-10-12 (n2-custom-12-12288, 12 vCPU, 12 GiB, Debian) -->

The reference trajectory is a non-rigorous high-precision Taylor integration; the grid comparison is a
consistency check that can expose bugs, not part of the proof. The closest encounter on the trajectory is
between bodies 4 and 5 at t ≈ 15.830 with separation ≈ 4.16 × 10⁻⁴ (non-rigorous reference value).
<!-- src: YOSHIDA_REDERIVATION_KO.md §4 -->

The run stored, for each of the 151 grid times t_k = k/100, k = 5950, …, 6100, two enclosing boxes (302 boxes)
with exact rational endpoints: a wide path box (component widths 1.1 × 10⁻² to 6.3 × 10⁻²) and a narrow
in-step or step-end box (widths 4.8 × 10⁻¹³ to 3.3 × 10⁻⁹).
<!-- src: GATE_Y1_PLAN_KO.md §1; yoshida_check/Y1_RESULT_KO.md §3 header -->

**Second run (R4-B, 144 bits).** The same code at 144-bit precision also reached t = 61.000086 (52,388 accepted
steps, 55,633 s, final largest width 10^−17.416, 0 violations). This is an independent *run*, not an independent
*implementation*. Its stored boxes were **not** checked against (★). In its grid check, 800 in-step and
step-end boxes at t ≤ 8.13, all narrower than 10^−28.9, did not contain the reference point; these are below
the 10⁻¹⁵ decision threshold and counted as undecidable, not as violations. We attribute this to the limited
accuracy of the reference trajectory itself (its nominal tolerance is 10⁻³⁵ and errors grow through close
encounters); **this explanation has not been checked**.
<!-- src: server_archive_20261012/r4_out.tgz → p3b_r3cp/out_r4/r3c_final_144w.json, r3c_grid_144w.jsonl
     (counts computed for this draft from the grid file: 800 UNDECIDED with ref outside, max wmax_log10 −28.875,
     t range 0.01–8.13); integrator_review/INTEGRATOR_REVIEW_KO.md §1-1 (ref_B_grid tolerance 1e-35) -->

**Bit-for-bit reproduction of the whole computation (R4R).** The 128-bit computation on [0, 61] was re-run
from the start with the Taylor-coefficient core replaced by the C++/GMP implementation of Section 2.6
(checked bit-identical to the Python core in 23,452 cases). All 31,811 steps agree in t, step size, number of
rejected attempts and enclosure width (0 mismatches); the stored-box file and the grid-check file have the
same SHA-256 as in the original run (`f65101c1…`, `fe53ded8…`); stop reason, final time and grid counts are
identical. Wall-clock time was 7,974 s (2.2 h) instead of 37,015 s. **Limitation:** this is a second
implementation of the same algorithm — the step rule, the QR frame update, the exact inverse and the
containment checks are the same Python code — not an independent algorithm.
<!-- src: R4R_FINAL_JUDGMENT_KO.md; GATE_R4R_PLAN_KO.md (sealed, not re-read for this draft); archive
     server_archive_20261012/batch3.tgz (SHA-256 648b6aca…, recomputed locally 2026-10-12) -->

### 4.2 Escape check (Y1)

Two programs written for this purpose, sharing no code with each other or with the integrator, evaluated
(C1)–(C3) on all 302 boxes: `y1_exact` (exact rationals; square roots bracketed by integer square roots to
2⁻²⁰⁰) and `y1_mpiv` (mpmath interval arithmetic, 256 bits). A box passes only if each inequality holds at
**every** point of the box (outer bounds). Exact constants were used (K = 769/60, d₁ = 1200/769).

| Result | Value |
|---|---|
| Boxes passing in both checkers | 256 = all boxes at all 128 grid times t = 59.73, …, 61.00 |
| Boxes passing at t ≤ 59.72 | 0 of 46 |
| Disagreements between checkers | 0 of 302 |
| Guaranteed minimum margins (t ≥ 59.73) | C3: 0.0235 (wide box), 0.0374 (narrow box); C1: 0.46; C2: 0.85 |
| Narrow box, t = 59.73: certified C3 margin | [0.03744808107650, 0.03744808108289] (reference value 0.03744808107969 inside) |
| Self-tests | 253 checks, 0 failures; 8 deliberately wrong checkers (λ↔ν, wrong d₁, wrong M, (C1) dropped, …) all detected |

<!-- src: Y1_FINAL_JUDGMENT_KO.md; yoshida_check/Y1_RESULT_KO.md §1, §2.2, §6, §7 -->

The non-rigorous reference trajectory first satisfies (★) at t ≈ 59.7223; the certified boxes switch from
"fail" (narrow-box C3 margin −0.012143 at t = 59.72) to "pass" (+0.037448 at t = 59.73), consistent with this.
As a further consistency check, the exactly known conserved quantities (momentum 0, angular momentum 0,
energy −769/60) lie inside the corresponding interval evaluations of all 302 boxes.
<!-- src: yoshida_check/Y1_RESULT_KO.md §6, §8 -->

### 4.3 Conclusion

Combining 4.1, 4.2 and (★): at t₀ = 59.73 the exact state satisfies (C1)–(C3), hence for all t ≥ 59.73 the
separation of bodies 4 and 5 stays below 1200/769 ≈ 1.5605 and body 3's distance from the binary's centre of
mass grows at least linearly. On [0, 61.0055] the enclosures also show that no collision occurs. This
conclusion is **conditional on the items in Section 8**.

---

## 5. Verification of the implementation

Implementation, checking and review were done in separate AI-tool sessions with restricted file access,
coordinated by the author (see "Use of AI tools" below). The "independence" in this section therefore means
separately written code under a sealed written contract, run and re-run separately — not separate human
teams.

| Check | What was done | Outcome | Source |
|---|---|---|---|
| Interval layer (R0) | Independent checker: 200,000 random cases at precisions 24/53/128/300, boundary cases, sqrt; checked against exact rational values | Containment violations 0 (S1 477,198 checks, S2 434,956, S3 17,000); checker caught 25/25 planted faults | `rigorous/R0_FINAL_RUN_KO.md` |
| Taylor step (R1) | Independent checker vs. 130-digit references: oscillator, circular Kepler, Pythagorean single steps, path boxes | 1600/1600/9/3200 pass, 0 violations; 19/19 planted faults caught | `rigorous/R1_FINAL_RUN_KO.md` |
| Lohner step, checker v1 | Independent checker on both pivot rules | 0 violations, but 10 cases (S2: 8, X4: 2) were refused because the *checker* supplied non-dyadic start states | `GATE_Y1_PLAN_KO.md` §6; `r3c-check-v2-handover` note |
| Lohner step, checker v2 | Fixed input generation (dyadic grid) | Width-pivot rule (used for R4): **PASS, 0 violations, 0 refused** (S1 160, S2 2,356, C 4,098, X4 732 cases). Colnorm rule: 0 violations, unresolved (52 cases refused/not run because boxes grew to contain a collision configuration) | `rigorous/R3C_CHECK_V2_FINAL_KO.md` |
| Code review | Line-by-line review of all five files against a stated contract, plus attack tests and mutation tests against independent mpmath references: 806 one-step cases (oscillator 480, Kepler 216, three-body 60, checker-gap class 50), 47 unit tests, 4 chains (incl. the Pythagorean problem with the R4 settings, t = 0 → 2.0005, 1,309 steps through the first close encounter, 219,912 set-membership checks) | **0 violations**; no path found by which an output box could miss the true solution; 3 issues unrelated to rigour (one robustness, two documentation). Exact inverse step confirmed. Of 12 mutants, 10 were detected; not detected: moving the exact inverse 1 ulp inward, and an exact-midpoint centre (which affects only tightness) | `integrator_review/INTEGRATOR_REVIEW_KO.md` §0, §3, §4-1, §4-2, §6 |
| Determinism | Server run vs. earlier local diagnostic run | First 4,494 steps identical in time, step size, rejections and exact width | `GATE_R4_PLAN_KO.md` §2 |
| Full re-run (R4R) | t = 0 → 61 re-run with the C++ Taylor core | 31,811/31,811 steps identical; box and grid files SHA-256 identical; 2.2 h vs 10.3 h | `R4R_FINAL_JUDGMENT_KO.md` |
| C++ port (FC1) | Bit-identity check of a GMP port of the Taylor core | 23,452 cases identical; used only for the R4R re-run, not for the proof record | `FC1_FINAL_JUDGMENT_KO.md` |
| (★) proof | Independent line-by-line review; new analytic proof of Lemmas A, B; sympy (33/33 identities, 4/4 negative controls); adversarial search of the lemmas (40,000 points, 50 digits); dynamic falsification (800 states) | Logic holds; three presentation fixes required (vector form of the energy identity, citation of (R), explicit definition of the bootstrap set) | `proof_review/PROOF_REVIEW_KO.md` |
| Global existence (R), v1.2 | sympy check of the 10 identities used in Appendix A.2 (conservation laws of the regularized system, recovery of Newton's equations and of (E)); 4 deliberately wrong variants; non-rigorous comparison of the regularized and the Cartesian integration of the Pythagorean problem to t = 30 | 10/10 identities hold; 4/4 wrong variants rejected; positions agree to 7.5 × 10⁻⁸ through a 4–5 approach of 0.0022 | `proof_review/s3/` |
| Literature | Direct reading of Yoshida 1972 (pp. 391–407), 1973 Erratum (p. 285) and 1974 (pp. 367–369) scans; zbMATH Open API search | Corollary gap confirmed and found corrected by the author (1973, 1974) | `LITERATURE_RECHECK_KO.md`, `LITERATURE_RECHECK_ZBMATH_KO.md`, this paper |

<!-- integrator_review §4-1/§4-2 (addendum, file dated 2026-10-12 17:17) now contains the test summary used above. -->

---

## 6. What did not work

We report these because they bound the method.

| Attempt | Outcome | Source |
|---|---|---|
| Axis-aligned boxes (R2), 128 and 256 bits | Enclosure width exceeded 1 at t = 6.968 and 7.055 | `rigorous/R2_RESULT_KO.md`; r2-run note |
| Lohner–QR with textbook column-norm pivoting (R3c), 128 bits | Width exceeded 1 at t = 1.8767, just before the first close encounter (t ≈ 1.879, separation ≈ 0.0097). The pivot order changed at every step during the encounter. Changing only the pivot rule to the width-weighted one passed the same encounter at width ≈ 10⁻²⁹ | `rigorous/R3C_RESULT_KO.md` §9-1, §9-3 |
| Width pivot, 112 bits (R4-C) | Width exceeded 1 at t ≈ 54.498 (15,764 steps, 21,613 s); did not reach the escape window | `r3c_final_112w.json` in `r4_out.tgz` |
| Width pivot, 96 bits (R4-D) | Stopped at t ≈ 27.789: the a-priori box could not exclude a collision configuration (5,102 steps, 9,110 s) | `r3c_final_96w.json` in `r4_out.tgz` |

None of these runs produced a containment violation in its grid checks; they failed by losing tightness,
not by producing wrong boxes. <!-- src: grid_counts in r3c_final_112w.json, r3c_final_96w.json (VIOLATION 0) -->
So 128 bits is, for this step rule and order, close to the minimum that reaches t = 61 (two data points
below; not a proof of minimality). The width growth over [0, 61] is roughly exponential and dominated by
close encounters; we have not quantified the contribution of each encounter.

---

## 7. Discussion

- **Cost.** About 10 h of single-core pure-Python time for the main run, under one second for the escape
  check. The method is entirely standard; the problem-specific choices are the width-weighted pivot, the
  precision, and stopping at t = 61 where the (★) margin is ≈ 0.93 (non-rigorous) instead of the first time it
  holds. <!-- src: YOSHIDA_REDERIVATION_KO.md §4 margins; ROADMAP_AMENDMENT_4_KO.md E1 -->
- **Why no regularization.** The closest approach (≈ 4 × 10⁻⁴) was passed in Cartesian coordinates; the
  a-priori enclosure certifies the absence of collisions on [0, 61.0055]. Regularization enters only through
  the theorem (★) for t > 61.
- **Comparison with Yoshida (1972).** His application used data with ρ, r differing by 1–6 % from our
  reference at t = 59.72; our certified escape time 59.73 is consistent with his 59.72.
  <!-- src: YOSHIDA_REDERIVATION_KO.md §4; LITERATURE_RECHECK_KO.md §2 row p.406 -->

---

## 8. Limitations (what this paper does not establish)

1. **One-ulp rounding faults.** The independent checkers detect containment errors larger than their
   resolution; a fault that rounds 1 ulp inward (2⁻¹²⁸ relative) is excluded only by code reading. Python
   integer, `Fraction`, `>>`, `math.isqrt` and `Fraction(float)` semantics are taken from the language
   specification. <!-- src: integrator_review §3-4, §6 -->
2. **Global existence (R) — proved in version 1.2, not reviewed.** Version 1 cited the classical collision
   facts (Painlevé, Sundman, Levi-Civita, Kustaanheimo–Stiefel, Siegel–Moser) for the statement that binary
   collisions of bodies 4 and 5 can be continued and do not accumulate in finite time, and a verbatim
   statement for L = 0 under the separation bound (SEP) was not found. Version 1.2 replaces this by a
   self-contained proof (Appendix A.2) that uses only the existence, uniqueness and extension theorems for
   ODEs. It covers the planar problem, which is the case of the application; the spatial case is not
   treated. This proof is internal and has not been externally reviewed.
   <!-- src: proof_review/s3/S3_REGULARIZATION_PROOF_EN.md; proof_review/s3/s3_algebra_check.py -->
3. **No formal verification.** Neither the integrator nor the proof of (★) has been machine-checked.
4. **No external review.** All reviews so far were internal to the project.
5. **Literature.** Absence of a prior rigorous proof is a result of a bounded search (Section 1.2);
   **MathSciNet was not searched** (subscription), and zbMATH Open only through its public API. Lemma proofs
   in Yoshida (1972) were not audited step by step; the 1974 paper was read only on pp. 367–369. Hirayama
   (2019) and Lohner's thesis (1988) were identified from bibliographic records; their texts were not read.
6. **Reference trajectory.** The 144-bit grid anomaly in Section 4.1 is unexplained beyond a plausible
   inference.
7. **One integrator algorithm.** The 128- and 144-bit runs share all code. The R4R re-run reproduces the
   128-bit computation bit for bit with a second implementation of the Taylor core, but the rest of the step
   (step rule, QR frame, exact inverse, containment checks) is the same code. Independence lies in the
   checkers and the review, not in a second integrator algorithm.

---

## 9. Reproducibility

All code, plans (sealed by SHA-256 before each run), run records and the server archive are listed with
hashes in `paper/REPRODUCE.md`. The main run needs about 10 h on one core; the escape check takes seconds.
The package is released on Zenodo: https://doi.org/10.5281/zenodo.23317644 (version 1); all versions are
listed under the concept DOI https://doi.org/10.5281/zenodo.23317643. Version 1.2 adds `proof_review/s3/`
(the proof note for Appendix A.2 and its checks).

The released files are redacted copies. The originals are kept unchanged, and all seal hashes refer to them.
In the copies, only the operating-system user name, the cloud instance name and project ID, and any e-mail
address are replaced, by exact string replacement, with fixed placeholders (e.g. `/home/USER`). Numbers,
fractions and computed values are never touched. A public table `REDACTION_MAP.json` lists, for every file, the
SHA-256 of the original and of the copy and the number of replacements per rule. The script
`redact_verify.py` (standard library only) checks that applying the inverse replacement to each copy
reproduces the original byte for byte. It also checks that no placeholder occurred in the originals, which
makes the inverse unambiguous. The output of this check is released with the package. Outside readers can
verify the copy hashes against the table; the byte-level comparison requires the originals and is performed by
the author. <!-- src: author decision 2026-10-12; paper/REPRODUCE.md §8 -->

---

## Use of AI tools

- The design decisions and all final judgments (what to compute, acceptance criteria, which results to accept,
  what to claim) were made by the author.
- Code writing, the independent checkers, the code and proof reviews, the literature searches and the writing
  of the documents (including this paper) were carried out by an AI tool (Anthropic Claude), in multiple
  sessions kept separate from one another.
- Verification structure: each computational step followed a written plan sealed by SHA-256 before it was run;
  implementation and checking were assigned to separate sessions (a checker session did not read the code it
  checked, and vice versa); and a separate coordinating session, directed by the author, independently re-ran key steps (e.g. the checker v2
  final run, the FC1 checks and the Y1 escape check on freshly uploaded files).
- No AI tool is listed as an author. Responsibility for the content rests with the author.
<!-- Policy note (not for the paper text): COPE's position statement (first published 2023-02-13) says AI tools
     cannot be authors and their use should be disclosed. Seen only via search-result summaries; the official
     page https://publicationethics.org/cope-position-statements/ai-author returned HTTP 403 to automated access.
     Venue-specific rules (AI disclosure, pen names) are to be checked when a journal is chosen. -->
<!-- src: coordinator message 2026-10-12 (author's decision); Y1_FINAL_JUDGMENT_KO.md, FC1_FINAL_JUDGMENT_KO.md,
     rigorous/R3C_CHECK_V2_FINAL_KO.md (coordinator re-runs) -->

## Acknowledgements

Computations were run on Google Cloud (free trial credits).

---

## References

Checked 2026-10-12 against Crossref (`https://api.crossref.org/works/<DOI>`), zbMATH Open
(`https://api.zbmath.org/v1/document/_search`), the IPSJ digital library or the ADS scan, as noted per entry.
"Checked" refers to the bibliographic record; the extent to which the text was read is noted separately.
Appendix A has its own reference list with bracketed keys.

1. Aarseth, S. J., Anosova, J. P., Orlov, V. V., Szebehely, V. G. (1994). Global chaoticity in the Pythagorean
   three-body problem. *Celest. Mech. Dyn. Astron.* 58(1), 1–16. doi:10.1007/BF00692114. (Crossref; abstract
   read.)
2. Boekholt, T., Portegies Zwart, S. (2015). On the reliability of N-body simulations. *Comput. Astrophys.
   Cosmol.* 2, 2. doi:10.1186/s40668-014-0005-3; arXiv:1411.6671. (Crossref; pp. 1–5 read.)
3. Burrau, C. (1913). Numerische Berechnung eines Spezialfalles des Dreikörperproblems. *Astron. Nachr.*
   195(6), 113–118. doi:10.1002/asna.19131950602. (Crossref; text not read.)
4. Capiński, M. J., Kepley, S., Mireles James, J. D. (2023). Computer assisted proofs for transverse collision
   and near collision orbits in the restricted three body problem. *J. Differential Equations* 366, 132–191.
   doi:10.1016/j.jde.2023.03.053; arXiv:2205.03922. (Crossref; partly read.)
5. Hirayama, H., Nishikawa, S., Shiraishi, F. (2006). Performance of numerical computation of the
   differential equations by Taylor series [in Japanese: Taylor展開法による常微分方程式の数値計算の性能].
   IPSJ SIG Technical Report 2006-HPC-108, 61–66 (2006-10-06). (IPSJ digital library record
   https://ipsj.ixsq.nii.ac.jp/records/28869; abstract only.)
6. Hirayama, H. (2019). Numerical calculation by quadruple precision higher order Taylor series method of the
   Pythagorean problem of three bodies. In *Integral Methods in Science and Engineering*, 173–183.
   doi:10.1007/978-3-030-16077-7_14. (Crossref; zbMATH; text not read. Volume editors not recorded here.)
7. Kapela, T., Mrozek, M., Wilczak, D., Zgliczyński, P. (2021). CAPD::DynSys: a flexible C++ toolbox for
   rigorous numerical analysis of dynamical systems. *Commun. Nonlinear Sci. Numer. Simul.* 101, 105578.
   doi:10.1016/j.cnsns.2020.105578; arXiv:2010.07097. (Crossref; partly read.)
8. Kapela, T., Simó, C. (2007). Computer assisted proofs for nonsymmetric planar choreographies and for
   stability of the Eight. *Nonlinearity* 20(5), 1241–1255. doi:10.1088/0951-7715/20/5/010. (Crossref; not
   read.)
9. Kapela, T., Zgliczyński, P. (2003). The existence of simple choreographies for the N-body problem — a
   computer-assisted proof. *Nonlinearity* 16(6), 1899–1918. doi:10.1088/0951-7715/16/6/302. (Crossref;
   abstract.)
10. Lohner, R. (1988). Einschließung der Lösung gewöhnlicher Anfangs- und Randwertaufgaben und Anwendungen.
    Dissertation, Universität Karlsruhe, 122 pp. (zbMATH 0663.65074; text not read.)
11. Lohner, R. J. (2001). On the ubiquity of the wrapping effect in the computation of error bounds. In
    *Perspectives on Enclosure Methods*, 201–216. doi:10.1007/978-3-7091-6282-8_12. (Crossref; not read —
    cited for the pivoted-QR recommendation as reported by ref. 13.)
12. Nedialkov, N. S., Jackson, K. R., Corliss, G. F. (1999). Validated solutions of initial value problems for
    ordinary differential equations. *Appl. Math. Comput.* 105(1), 21–68. doi:10.1016/S0096-3003(98)10083-8.
    (Crossref lists only the first author; full author list from zbMATH 0934.65073. Web summary only.)
13. Revol, N. (2023). Affine iterations and wrapping effect. *Acta Cybernetica* 26(1), 129–147.
    doi:10.14232/actacyb.295251. (Crossref; excerpts read via scite.)
14. Szebehely, V., Peters, C. F. (1967). Complete solution of a general problem of three bodies. *Astron. J.*
    72, 876. doi:10.1086/110355. (Crossref; ADS scan OCR read.)
15. Walawska, I., Wilczak, D. (2016). An implicit algorithm for validated enclosures of the solutions to
    variational equations for ODEs. *Appl. Math. Comput.* 291, 303–322. doi:10.1016/j.amc.2016.07.005;
    arXiv:1509.07388. (Crossref; abstract and excerpts.)
16. Yoshida, J. (1972). Improved criteria for hyperbolic-elliptic motion in the general three-body problem.
    *Publ. Astron. Soc. Japan* 24(3), 391–408. doi:10.1093/pasj/24.3.391. (Scan read, pp. 391–407. Crossref
    gives the issue date as 1969-12-25; the journal page and the 1973 erratum give 1972.)
17. Yoshida, J. (1973). Erratum: Improved criteria for hyperbolic-elliptic motion in the general three-body
    problem. *Publ. Astron. Soc. Japan* 25(2), 285. doi:10.1093/pasj/25.2.285. (Crossref record has no author;
    the scan names the corrected paper "by Junzo Yoshida". Full page read on the ADS scan, 2026-10-12.)
18. Yoshida, J. (1974). Improved criteria for hyperbolic-elliptic motion in the general three-body problem. II.
    *Publ. Astron. Soc. Japan* 26(3), 367–377. doi:10.1093/pasj/26.3.367. (Crossref; scan read pp. 367–369.)
19. Yuasa (1971). Private communication, as cited in Yoshida (1972), p. 406. (No publication of these data
    found in Crossref or a web search on 2026-10-12.)
20. Zgliczyński, P. (2002). C¹ Lohner algorithm. *Found. Comput. Math.* 2(4), 429–465.
    doi:10.1007/s102080010025. (Crossref; not read.)
21. Collision and regularization literature: see the reference list of Appendix A ([P1897], [Su1912],
    [LC1920], [KS1965], [SM1971]). Since version 1.2 these are cited for attribution and context only; the
    proof of (R) uses the ODE theorems of [Ha]. Bibliographic records checked 2026-10-12: Painlevé 1897
    (JFM 28.0262.01); Sundman, Acta Math. 36, 105–179 (JFM 43.0826.01 dates it 1912, Crossref
    doi:10.1007/BF02422379 dates it 1913); Levi-Civita, Acta Math. 42 (1920) 99–144, doi:10.1007/BF02404404;
    Kustaanheimo–Stiefel, J. Reine Angew. Math. 218 (1965) 204–219, doi:10.1515/crll.1965.218.204;
    Siegel–Moser 1971 (zbMATH 0312.70017).

22. Portegies Zwart, S. F., Boekholt, T. C. N. (2018). Numerical verification of the microscopic time
    reversibility of Newton's equations of motion: fighting exponential divergence. *Commun. Nonlinear Sci.
    Numer. Simul.* 61, 160–166. doi:10.1016/j.cnsns.2018.02.002; arXiv:1802.00970. (Crossref; text read via
    ar5iv for the Pythagorean values.)

---

## Appendix A. The escape criterion (★): statement and self-contained proof

<!-- Appendix A is inserted verbatim (heading levels only changed) from proof_review/PROOF_STAR_FINAL_EN.md, SHA-256 4dfaea735cf3a90fa07d610be0c46c57dbc95150400630feb215acfac411ce4f (2026-10-12). Version 1.2 (2026-10-13) revises the last paragraph of A.0, A.2, A.4.5, A.6 and the reference list, from proof_review/s3/S3_REGULARIZATION_PROOF_EN.md; PROOF_STAR_FINAL_EN.md itself is unchanged. -->

*Draft for a paper appendix, prepared 2026-10-12. It has not yet been reviewed externally or verified formally (e.g. in a proof assistant).*

### A.0 Scope and attribution

The criterion proved here is the case $d = m_1m_2/K$ of **Theorem 1 of Yoshida (1972)** [Y72, p. 400]. It is also the case $d = m_1m_2/K$ of **criterion A′ of Yoshida (1974)** [Y74, p. 369]; see Remark A.4.2.

The corollary to Theorem 1 in [Y72, p. 402] does not assume $r(0)\le d$. **Yoshida himself corrected this, first in an erratum [Y73] and then in [Y74, p. 369, Note].** The erratum [Y73] replaces condition (26) of the corollary by $d\ge r(0)$; it also corrects inequality (40) on p. 405 of [Y72]. The later Note [Y74, p. 369] calls conditions (25)–(26) of [Y72] "incorrect" and replaces them by $\rho(0)\ge ad\ge a\,r(0)$. The statement below assumes $r(t_0)<d_1$ explicitly, so that issue does not affect it.

Lemmas A and B below are the case $n=1$ of Lemma 1 of [Y72, p. 393] and Lemma 2 of [Y72, p. 396]. Yoshida proves them by induction on $n$ and by analysing the $\theta$-derivative, respectively. We give different, self-contained proofs based on the Legendre (multipole) expansion. This is a standard technique, and we do not claim originality for it.

Everything below is proved in full for the planar problem, using standard ODE theory and elementary analysis. Global existence through binary collisions is proved in A.2 (new in version 1.2; version 1 cited classical collision facts there). The spatial case of A.2 is not treated; A.1, A.3 and A.4 hold in $\mathbb R^2$ and $\mathbb R^3$.

---

### A.1 Setting and notation

Three point masses $m_0,m_1,m_2>0$ move in $\mathbb R^2$ or $\mathbb R^3$ under Newtonian gravity with $G=1$. We work in the barycentric frame, where the total momentum is zero. The body $m_0$ is the candidate escaper. The other two are labelled so that
$$m_1\le m_2 .$$
Define
$$
M=m_0+m_1+m_2,\qquad \lambda=\frac{m_2}{m_1+m_2},\qquad \nu=\frac{m_1}{m_1+m_2},\qquad
m=\frac{m_1m_2}{m_1+m_2},\qquad \mu=\frac{m_0(m_1+m_2)}{M}.
$$
Hence $\lambda\ge\nu>0$, $\lambda+\nu=1$, $\mu M\lambda=m_0m_2$ and $\mu M\nu=m_0m_1$.

**Jacobi vectors.**
- $\mathbf r=\mathbf x_2-\mathbf x_1$.
- $\boldsymbol\rho=\mathbf x_0-\mathbf g$, where $\mathbf g=(m_1\mathbf x_1+m_2\mathbf x_2)/(m_1+m_2)$.
- Then $\mathbf x_2-\mathbf g=\nu\mathbf r$ and $\mathbf x_1-\mathbf g=-\lambda\mathbf r$. Thus the lighter body $m_1$ is the one farther from $\mathbf g$.
- Write
$$
r=|\mathbf r|,\quad \rho=|\boldsymbol\rho|,\quad
r_1=|\mathbf x_0-\mathbf x_2|=|\boldsymbol\rho-\nu\mathbf r|,\quad
r_2=|\mathbf x_0-\mathbf x_1|=|\boldsymbol\rho+\lambda\mathbf r| .
$$
- For $r>0$ let $c=\cos\theta=\boldsymbol\rho\cdot\mathbf r/(\rho r)$. Then
$$
r_1^2=\rho^2+\nu^2r^2-2\nu\rho r c,\qquad r_2^2=\rho^2+\lambda^2r^2+2\lambda\rho r c .
$$

**Radial velocity.** $\dot\rho:=\dfrac{d}{dt}|\boldsymbol\rho|=\dfrac{\boldsymbol\rho\cdot\dot{\boldsymbol\rho}}{\rho}$. This is the radial component, *not* the speed $|\dot{\boldsymbol\rho}|\ge|\dot\rho|$.

**Energy.** In the barycentric frame, and with **vector** velocities,
$$
H=\tfrac12 m|\dot{\mathbf r}|^2+\tfrac12\mu|\dot{\boldsymbol\rho}|^2-\frac{m_1m_2}{r}-\mu M\Big(\frac{\lambda}{r_1}+\frac{\nu}{r_2}\Big). \tag{E}
$$
We assume $H<0$ and put $K=-H>0$ and $d_1=m_1m_2/K$.

**Equation of motion of $\boldsymbol\rho$.**
$$
\ddot{\boldsymbol\rho}=-M\Big(\frac{\lambda}{r_1^3}+\frac{\nu}{r_2^3}\Big)\boldsymbol\rho+\lambda\nu M\Big(\frac1{r_1^3}-\frac1{r_2^3}\Big)\mathbf r=:\Phi(\boldsymbol\rho,\mathbf r). \tag{EM}
$$

**Radial acceleration.** Where $\rho>0$, Lagrange's identity gives
$$
\ddot\rho=\frac{\boldsymbol\rho\cdot\ddot{\boldsymbol\rho}}{\rho}+\frac{|\boldsymbol\rho|^2|\dot{\boldsymbol\rho}|^2-(\boldsymbol\rho\cdot\dot{\boldsymbol\rho})^2}{\rho^3}\;\ge\;\frac{\boldsymbol\rho\cdot\ddot{\boldsymbol\rho}}{\rho}. \tag{RA}
$$
The dropped term equals $|\boldsymbol\rho\times\dot{\boldsymbol\rho}|^2/\rho^3\ge0$. Nothing is assumed about the angular momentum; in particular $L=0$ is allowed.

**Auxiliary functions.** For $\rho>0$, $r\ge0$ and $c\in[-1,1]$ with $r_1,r_2>0$, set
$$
P(\rho,r,c)=\frac{\lambda}{r_1}+\frac{\nu}{r_2},\qquad
F(\rho,r,c)=\frac{\lambda(\rho-\nu rc)}{r_1^3}+\frac{\nu(\rho+\lambda rc)}{r_2^3}=-\partial_\rho P .
$$
- The potential term in (E) is $\mu M P$.
- By (EM), $\boldsymbol\rho\cdot\ddot{\boldsymbol\rho}/\rho=-M\,F(\rho,r,c)$.
- For $r=0$ these reduce to $P=1/\rho$ and $F=1/\rho^2$ whatever the value of $c$.

For $d\ge 0$ and $\rho>\lambda d$, set
$$
W(\rho,d)=\frac{\lambda}{\rho+\nu d}+\frac{\nu}{\rho-\lambda d},\qquad
\hat F(\rho,d)=\frac{\lambda}{(\rho+\nu d)^2}+\frac{\nu}{(\rho-\lambda d)^2}=-\partial_\rho W(\rho,d).
$$
Then $W(\rho,d)=P(\rho,d,-1)$ and $\hat F(\rho,d)=F(\rho,d,-1)$. Geometrically, this is the collinear configuration in which the lighter component $m_1$, at distance $\lambda d$ from $\mathbf g$, lies between $\mathbf g$ and $m_0$.

---

### A.2 Continuation through binary collisions and global existence (planar case)

*Revised in version 1.2. Version 1 cited classical facts here (Painlevé, Sundman, Levi-Civita, Kustaanheimo–Stiefel, Siegel–Moser). This section proves what the proof of (★) needs. It uses only the existence, uniqueness and extension theorems for ODEs with a locally Lipschitz right-hand side [Ha, Ch. II]. The motion is planar. Plane vectors are identified with complex numbers, and $\langle a,b\rangle=\operatorname{Re}(\bar ab)$ is the dot product. The spatial case would need the Kustaanheimo–Stiefel analogue [KS1965]. It is not treated, and the application in A.5 does not need it.*

**The regularized system.** Let $k=m_1+m_2$. The relative motion satisfies $\ddot{\mathbf r}=-k\mathbf r/r^3+\mathbf f$, where
$$
\mathbf f(\boldsymbol\rho,\mathbf r)=m_0\Big[\frac{\boldsymbol\rho-\nu\mathbf r}{r_1^3}-\frac{\boldsymbol\rho+\lambda\mathbf r}{r_2^3}\Big].
$$
Following Levi-Civita [LC1920], use the state $Y=(t,u,w,h,\boldsymbol\rho,\boldsymbol\sigma)\in\mathbb R\times\mathbb C^2\times\mathbb R\times\mathbb C^2$ and write $z=u^2$ in place of $\mathbf r$. Consider, on the open set
$$
D=\{Y:\ |\boldsymbol\rho-\nu z|>0,\ |\boldsymbol\rho+\lambda z|>0\},
$$
the system (with $'=d/ds$)
$$
\begin{aligned} t'&=|u|^2, & u'&=w, & w'&=\tfrac h2u+\tfrac12|u|^2\,\bar u\,\mathbf f(\boldsymbol\rho,z),\\ h'&=\langle 2uw,\mathbf f(\boldsymbol\rho,z)\rangle, & \boldsymbol\rho'&=|u|^2\boldsymbol\sigma, & \boldsymbol\sigma'&=|u|^2\,\Phi(\boldsymbol\rho,z). \end{aligned} \tag{LC}
$$
Its right-hand side is real-analytic on $D$; the only denominators are $r_1^3$ and $r_2^3$. A binary collision $u=0$ is an ordinary point of (LC).

**Lemma A.2.1 (identities).** Along any solution of (LC) in $D$:
1. $I:=2|w|^2-h|u|^2$ is constant.
2. $\mathcal H:=m h+\tfrac12\mu|\boldsymbol\sigma|^2-\mu M P(\boldsymbol\rho,z)$ is constant.

Suppose in addition that $I=k$ and $u\ne0$ on an $s$-interval. Let dots denote $d/dt=|u|^{-2}\,d/ds$. Then on that interval:

3. $\dot z=2w/\bar u$ and $\ddot z=-kz/|z|^3+\mathbf f$;
4. $h=\tfrac12|\dot z|^2-k/|z|$;
5. $\dot{\boldsymbol\rho}=\boldsymbol\sigma$ and $\ddot{\boldsymbol\rho}=\Phi$;
6. $\mathcal H$ is the energy (E).

The Jacobi forms of Newton's equations used here ($\ddot{\mathbf r}=-k\mathbf r/r^3+\mathbf f$ and (EM)) follow from the Cartesian equations.

*Proof.* All items are polynomial or rational identities.

Item 3 by hand:
- $\dot z=2uw/|u|^2=2w/\bar u$.
- $\frac{d}{ds}\dot z=2w'/\bar u-2|w|^2/\bar u^2$.
- Using $2|w|^2=k+h|u|^2$,
$$\ddot z=|u|^{-2}\big(hu/\bar u+|u|^2\mathbf f-(k+h|u|^2)/\bar u^2\big)=\mathbf f-k/(|u|^2\bar u^2)=\mathbf f-kz/|z|^3 .$$

Item 1:
- $\frac{d}{ds}2|w|^2=2h\langle w,u\rangle+2|u|^2\langle w,\bar u\mathbf f\rangle$.
- $\frac{d}{ds}(h|u|^2)=|u|^2\langle 2uw,\mathbf f\rangle+2h\langle u,w\rangle$.
- These are equal, since $\langle w,\bar u\mathbf f\rangle=\langle uw,\mathbf f\rangle$.

All six items, and the Jacobi forms of Newton's equations, were also checked by computer algebra (A.6). $\square$

**Definition (continued solution).** Let the classical solution be given near $t_0$, with $r(t_0)>0$. Choose $u_0$ with $u_0^2=\mathbf r(t_0)$. Put $w_0=\bar u_0\dot{\mathbf r}(t_0)/2$ and $h_0=\tfrac12|\dot{\mathbf r}(t_0)|^2-k/r(t_0)$; then $I=k$.

The *continued solution* is the maximal solution $Y(s)$, $s\in[0,S^*)$, of (LC) in $D$ with $Y(0)=(t_0,u_0,w_0,h_0,\boldsymbol\rho(t_0),\dot{\boldsymbol\rho}(t_0))$. It is read in physical time through $t(s)$. Let
$$T^*:=\lim_{s\to S^*}t(s)\in(t_0,\infty].$$
Then $[t_0,T^*)$ is its maximal interval of existence.

The choice $-u_0$ gives $(t,-u,-w,h,\boldsymbol\rho,\boldsymbol\sigma)$, which is the same motion. Physically, the solution passes a collision as an elastic bounce.

**Lemma A.2.2 (properties).**
- (i) $t(s)$ is $C^1$ and strictly increasing. The zeros of $u$, which are the binary collisions, are isolated in $s$. In particular, collision times do not accumulate in $[t_0,T^*)$.
- (ii) Let $s(t)$ be the inverse function. Then $\boldsymbol\rho(t):=\boldsymbol\rho(s(t))$ is $C^2$ on $[t_0,T^*)$, with $\dot{\boldsymbol\rho}=\boldsymbol\sigma\circ s$ and $\ddot{\boldsymbol\rho}=\Phi(\boldsymbol\rho,\mathbf r)$, also at collision times. The function $\mathbf r(t)=u(s(t))^2$ is continuous.
- (iii) Between collisions, $(\mathbf r,\boldsymbol\rho)$ solves Newton's equations, and (E) holds with the constant value $H$.

*Proof.*
- (i) $t'=|u|^2\ge0$. At a zero of $u$, item 1 of Lemma A.2.1 gives $|u'|^2=k/2>0$. So the zeros of $u$ are isolated, and $t'>0$ except on a discrete set.
- (ii) The substitution rule for the increasing $C^1$ map $t(s)$ gives
$$\boldsymbol\rho(s(t))-\boldsymbol\rho(t_0)=\int_{t_0}^{t}\boldsymbol\sigma(s(\tau))\,d\tau ,$$
with a continuous integrand. The same argument applied to $\boldsymbol\sigma$ gives $\dot{\boldsymbol\sigma}=\Phi\circ s$.
- (iii) This is items 3–6 of Lemma A.2.1, together with item 2. $\square$

**Proposition A.2.3 (statement (R)).** If $\min(r_1,r_2)\ge\delta>0$ along the continued solution, then $T^*=\infty$.

*Proof.* Suppose $T^*<\infty$, and put $\Delta=T^*-t_0$. The following bounds hold for all $s\in[0,S^*)$.

- **(B1) Forces.** $|\mathbf f|\le m_0(r_1^{-2}+r_2^{-2})\le 2m_0/\delta^2=:F$, $P\le1/\delta$ and $|\Phi|\le M/\delta^2$. For the last bound: $\Phi=\tfrac Mk\ddot{\mathbf x}_0$ is an algebraic identity in $(\boldsymbol\rho,z)$, and $|\ddot{\mathbf x}_0|\le m_1/r_2^2+m_2/r_1^2\le k/\delta^2$.
- **(B2) $\boldsymbol\sigma,\boldsymbol\rho$.** $|\boldsymbol\sigma(s)-\boldsymbol\sigma(0)|\le\int|u|^2|\Phi|\,ds\le M\Delta/\delta^2$. So $|\boldsymbol\sigma|\le\Sigma:=|\dot{\boldsymbol\rho}(t_0)|+M\Delta/\delta^2$, and $|\boldsymbol\rho-\boldsymbol\rho(t_0)|\le\Sigma\Delta$.
- **(B3) $h$.** By item 2 of Lemma A.2.1, $mh=H-\tfrac12\mu|\boldsymbol\sigma|^2+\mu MP$. So $|h|\le A:=(|H|+\tfrac12\mu\Sigma^2+\mu M/\delta)/m$.
- **(B4) $u,w$.** Write $R=|u|^2$, which is $r$.
  - Where $R\ge1$ there is no collision. There $|dR/dt|\le|\dot z|=\sqrt{2h+2k/R}\le V:=\sqrt{2A+2k}$.
  - Fix $s$. Look back from $s$ to the last point at which $R\le R_0:=\max(1,r(t_0))$, or to $s=0$ if there is none. In between, $R>R_0\ge1$, so the speed bound applies and $R(s)\le R_0+V\Delta=:U^2$.
  - By item 1 of Lemma A.2.1, $|w|^2=\tfrac12(k+h|u|^2)\le\tfrac12(k+AU^2)=:W_b^2$.

*Case $S^*<\infty$.*
- By (B2)–(B4), $Y(s)$ stays in a closed and bounded set $\mathcal K$ with $r_1,r_2\ge\delta$.
- $\mathcal K$ is compact and contained in $D$.
- A maximal solution with $S^*<\infty$ must leave every compact subset of $D$ [Ha, Ch. II]. Contradiction.

*Case $S^*=\infty$.* Then $\int_0^\infty|u|^2ds=\Delta<\infty$; we show that this integral diverges.
- Choose $\eta\in(0,1]$ with $A\eta^2\le k/2$ and $2\eta+\varepsilon<\sqrt k/2$, where $\varepsilon:=\tfrac A2\eta+\tfrac F2\eta^3$.
- *Every interval $[a,a+1]$ contains a point with $|u|>\eta$.* Suppose instead that $|u|\le\eta$ on all of $[a,a+1]$. Then:
  - $|w'|\le\varepsilon$ there, so $|w-w(a)|\le\varepsilon$;
  - $|w(a)|^2\ge\tfrac12(k-A\eta^2)\ge k/4$;
  - hence $|u(a+1)-u(a)|\ge|w(a)|-\varepsilon\ge\sqrt k/2-\varepsilon>2\eta$.

  This contradicts $|u(a)|,|u(a+1)|\le\eta$.
- For $n\ge1$, take such a point $s^*\in[n,n+1]$.
  - Since $|u'|\le W_b$, we have $|u|\ge\eta/2$ on $[s^*-\ell,s^*+\ell]\subset[n-1,n+2]$, where $\ell=\min(\eta/(2W_b),1)$.
  - Hence $\int_{n-1}^{n+2}|u|^2ds\ge\eta^2\ell/2>0$.
- Summing over $n=1,4,7,\dots$ gives $\int_0^\infty|u|^2ds=\infty$. Contradiction. $\square$

**How A.4 uses this section.**
- Step 1 and Step 3 use the continuity of $r$ and $\dot\rho$, which follows from Lemma A.2.2(ii).
- Step 2 uses $\boldsymbol\rho\in C^2$ with $\ddot{\boldsymbol\rho}=\Phi$, and (E) at times with $r>0$. These are Lemma A.2.2(ii) and (iii).
- Step 4 uses Proposition A.2.3, which is the statement

> **(R)** If $\inf_{[t_0,T^*)}\min(r_1,r_2)>0$, then $T^*=\infty$.

**Regularity used below.** Wherever $\rho>0$, $\rho=|\boldsymbol\rho|$ is $C^2$ and (RA) holds. This follows from Lemma A.2.2(ii).

---

### A.3 Two inequalities

> **Lemma A.** Let $\lambda\ge\nu>0$, $\lambda+\nu=1$, $d\ge0$ and $\rho>\lambda d$. Then for all $0\le r\le d$ and all $c\in[-1,1]$,
> $$P(\rho,r,c)\le W(\rho,d).$$
>
> **Lemma B.** Under the same hypotheses,
> $$F(\rho,r,c)\le \hat F(\rho,d).$$

Both inequalities become equalities at $(r,c)=(d,-1)$, so they are sharp.

Combining Lemma B with (EM) gives the following. Along a solution, if $\rho>\lambda d$ and $r\le d$, then
$$
\frac{\boldsymbol\rho\cdot\ddot{\boldsymbol\rho}}{\rho}=-M F\ \ge\ -M\hat F(\rho,d)=M\,\partial_\rho W(\rho,d). \tag{B′}
$$

**Facts on Legendre polynomials $P_n$.**
- (L1) Generating function [DLMF, Eq. 18.12.11]: $(1-2xz+z^2)^{-1/2}=\sum_{n\ge0}P_n(x)z^n$ for $x\in[-1,1]$ and $|z|<1$.
- (L2) Laplace's integral [DLMF, Eq. 18.10.5]: $P_n(\cos\vartheta)=\frac1\pi\int_0^\pi(\cos\vartheta+i\sin\vartheta\cos\phi)^n\,d\phi$. Since $|\cos\vartheta+i\sin\vartheta\cos\phi|^2=\cos^2\vartheta+\sin^2\vartheta\cos^2\phi\le1$, we get $|P_n(x)|\le1$ on $[-1,1]$.
- (L3) In (L1), replacing $(x,z)$ by $(-x,-z)$ gives $P_n(-x)=(-1)^nP_n(x)$. Setting $x=1$ gives $\sum P_n(1)z^n=(1-z)^{-1}$, so $P_n(1)=1$. Hence $P_n(-1)=(-1)^n$.

**Proof of Lemmas A and B.**

*Step 1 (normalization).*
- Put $s=r/\rho$, $a=\nu s$ and $b=\lambda s$. Since $\lambda r\le\lambda d<\rho$, we have $0\le a\le b<1$.
- Now $r_1=\rho\sqrt{1-2ac+a^2}$ and $\rho-\nu rc=\rho(1-ac)$.
- Similarly $r_2=\rho\sqrt{1-2b(-c)+b^2}$ and $\rho+\lambda rc=\rho(1-b(-c))$.
- Therefore
$$
\rho P=\lambda G(a,c)+\nu G(b,-c),\qquad \rho^2F=\lambda g(a,c)+\nu g(b,-c),
$$
where
$$
G(a,x)=(1-2ax+a^2)^{-1/2},\qquad g(a,x)=(1-ax)(1-2ax+a^2)^{-3/2}.
$$

*Step 2 (expansions).*
- By (L1), $G(a,x)=\sum_{n}a^nP_n(x)$.
- One checks directly that $g=G+a\,\partial_aG$. Hence $g(a,x)=\sum_n(n+1)a^nP_n(x)$.
- For $0\le a<1$ and $|x|\le1$ both series converge absolutely, since $\sum(n+1)a^n|P_n(x)|\le(1-a)^{-2}$ by (L2). So they may be added and regrouped term by term.

*Step 3 (combination).* By Steps 1 and 2 together with (L3),
$$
\rho^2F=\sum_{n\ge0}(n+1)s^n\beta_nP_n(c),\qquad \rho P=\sum_{n\ge0}s^n\beta_nP_n(c),\qquad
\beta_n:=\lambda\nu^n+\nu(-\lambda)^n .
$$
- $\beta_0=1$, and $\beta_1=\lambda\nu-\nu\lambda=0$ (no dipole term about $\mathbf g$).
- **For every $n$, $(-1)^n\beta_n\ge0$.**
  - For even $n$, $(-1)^n\beta_n=\lambda\nu^n+\nu\lambda^n>0$.
  - For odd $n$, $(-1)^n\beta_n=\nu\lambda^n-\lambda\nu^n=\lambda\nu(\lambda^{n-1}-\nu^{n-1})\ge0$, because $\lambda\ge\nu$.
  - This is the only place where $m_1\le m_2$ is used.
- Consequently $|\beta_n|=(-1)^n\beta_n$.

*Step 4 (maximum over the angle).* By (L2) and Step 3,
$$
\rho^2F\le\sum_n(n+1)s^n|\beta_n|=\lambda\sum_n(n+1)(-\nu s)^n+\nu\sum_n(n+1)(\lambda s)^n
=\frac{\lambda}{(1+\nu s)^2}+\frac{\nu}{(1-\lambda s)^2}.
$$
In the same way,
$$
\rho P\le\sum_n s^n|\beta_n|=\frac{\lambda}{1+\nu s}+\frac{\nu}{1-\lambda s}.
$$
The right-hand sides equal $\rho^2F(\rho,r,-1)$ and $\rho P(\rho,r,-1)$, respectively.

*Step 5 (maximum over the separation).*
- Each term $(n+1)s^n|\beta_n|$ and $s^n|\beta_n|$ is nondecreasing in $s\ge0$.
- Since $r\le d$, we have $s\le s_d:=d/\rho<1/\lambda$, so all series still converge at $s_d$.
- Therefore
$$
\rho^2F\le\frac{\lambda}{(1+\nu s_d)^2}+\frac{\nu}{(1-\lambda s_d)^2}=\rho^2\hat F(\rho,d),\qquad
\rho P\le\frac{\lambda}{1+\nu s_d}+\frac{\nu}{1-\lambda s_d}=\rho\,W(\rho,d).
$$
Dividing by $\rho^2$ and by $\rho$ gives the lemmas. $\square$

**Remarks on the range.**
- (i) No upper bound on $\rho$ is needed.
- (ii) The bounds blow up as $\rho\to\lambda d^+$, but they hold at every point with $\rho>\lambda d$. No uniformity is used.
- (iii) $r=0$ is included; this covers the instant of a binary collision.
- (iv) Only $c=\cos\theta$ enters, so the lemmas hold in the plane and in space.
- (v) If $\lambda=\nu$, all odd $\beta_n$ vanish and equality also holds at $c=+1$.

**Remark (second proof of Lemma A).** Fix $r$. Then:
- $\partial_cP=\lambda\nu\rho r\,(r_1^{-3}-r_2^{-3})$, and $r_2^2-r_1^2=r\big((\lambda-\nu)r+2\rho c\big)$.
- So $P$ decreases in $c$ up to $c^*=-(\lambda-\nu)r/(2\rho)$ and increases after it. Its maximum over $[-1,1]$ is therefore attained at $c=\pm1$.
- $P(\rho,r,-1)-P(\rho,r,1)=2\lambda\nu r\big[(\rho^2-\lambda^2r^2)^{-1}-(\rho^2-\nu^2r^2)^{-1}\big]\ge0$.
- Finally, $P(\rho,r,-1)=W(\rho,r)$ and $\partial_dW=\lambda\nu\big[(\rho-\lambda d)^{-2}-(\rho+\nu d)^{-2}\big]\ge0$.

---

### A.4 The criterion

> **Theorem (★).** Let $m_1\le m_2$ and $H<0$, in the barycentric frame. Suppose that at some time $t_0$
>
> - **(C1)** $r(t_0)<d_1=m_1m_2/K$,
> - **(C2)** $\rho(t_0)>\lambda d_1$,
> - **(C3)** $\dot\rho(t_0)>\sqrt{2M\,W(\rho(t_0),d_1)}$.
>
> Let $C:=\sqrt{\dot\rho(t_0)^2-2M\,W(\rho(t_0),d_1)}>0$. Then the solution, continued through binary collisions as in A.2, exists for all $t\ge t_0$. Moreover, for all $t\ge t_0$:
> 1. $r(t)\le \dfrac{m_1m_2}{K+\frac12\mu C^2}<d_1$;
> 2. $\tfrac12\dot\rho(t)^2-M\,W(\rho(t),d_1)\ge\tfrac12C^2$ and $\dot\rho(t)>C$;
> 3. $\rho(t)\ge\rho(t_0)+C\,(t-t_0)$.
>
> In particular $\rho\to\infty$ while the binary separation stays bounded; that is, the motion is hyperbolic–elliptic.

**Proof.**

*Step 1 (bootstrap set).* Let
$$
\mathcal G=\big\{t\in[t_0,T^*):\ r(s)\le d_1\ \text{and}\ \dot\rho(s)\ge0\ \text{for all } s\in[t_0,t]\big\},\qquad \tau=\sup\mathcal G .
$$
$\mathcal G$ is an interval containing $t_0$. By (C1) and (C3), $r(t_0)<d_1$ and $\dot\rho(t_0)>0$, and both $r$ and $\dot\rho$ are continuous. Hence $\tau>t_0$.

*Step 2 (estimates on $\mathcal G$).* Fix $t\in\mathcal G$.

(a) Separation.
- $\rho$ is nondecreasing on $[t_0,t]$. So for $s\in[t_0,t]$ we have $\rho(s)\ge\rho(t_0)>\lambda d_1\ge\lambda r(s)$.
- Put $\delta:=\rho(t_0)-\lambda d_1>0$. Then on $[t_0,t]$
$$
r_2\ge\rho-\lambda r\ge\delta,\qquad r_1\ge\rho-\nu r\ge\rho(t_0)-\nu d_1\ge\delta. \tag{SEP}
$$

(b) Radial velocity.
- For $s\in[t_0,t]$, (RA) and (B′) with $d=d_1$ give $\ddot\rho(s)\ge M\,\partial_\rho W(\rho(s),d_1)$.
- Multiplying by $\dot\rho(s)\ge0$ gives $\dot\rho\ddot\rho\ge\frac{d}{ds}\big[M\,W(\rho(s),d_1)\big]$.
- Integrating (this uses $\rho\in C^2$, see A.2):
$$
\tfrac12\dot\rho(t)^2-M\,W(\rho(t),d_1)\ \ge\ \tfrac12\dot\rho(t_0)^2-M\,W(\rho(t_0),d_1)=\tfrac12 C^2 .
$$
- Since $W>0$ and $\dot\rho(t)\ge0$, it follows that $\dot\rho(t)>C$.

(c) Binary separation.
- If $r(t)=0$ there is nothing to prove.
- Otherwise use (E), $\tfrac12 m|\dot{\mathbf r}|^2\ge0$, $|\dot{\boldsymbol\rho}|\ge|\dot\rho|$, Lemma A (valid since $r(t)\le d_1$ and $\rho(t)>\lambda d_1$), and (b):
$$
\frac{m_1m_2}{r(t)}=K+\tfrac12 m|\dot{\mathbf r}|^2+\tfrac12\mu|\dot{\boldsymbol\rho}|^2-\mu M P
\ \ge\ K+\mu\Big(\tfrac12\dot\rho^2-M\,W(\rho,d_1)\Big)\ \ge\ K+\tfrac12\mu C^2 .
$$
- Hence $r(t)\le m_1m_2/(K+\tfrac12\mu C^2)<d_1$.

*Step 3 ($\tau=T^*$).* Suppose $\tau<T^*$.
- By continuity, $r(\tau)\le d_1$ and $\dot\rho(\tau)\ge0$. Hence $\tau\in\mathcal G$.
- By Step 2, $r(\tau)<d_1$ and $\dot\rho(\tau)>C>0$.
- By continuity these strict inequalities persist on $[\tau,\tau+\varepsilon)$ for some $\varepsilon>0$. So $\tau+\varepsilon/2\in\mathcal G$, which contradicts $\tau=\sup\mathcal G$.
- Therefore $\mathcal G=[t_0,T^*)$.

*Step 4 (global existence).* By (SEP), $\min(r_1,r_2)\ge\delta>0$ on $[t_0,T^*)$. So $T^*=\infty$ by (R).

*Step 5 (conclusion).* Items 1 and 2 hold for every $t\ge t_0$ by Steps 2–4. Item 3 follows by integrating $\dot\rho>C$. $\square$

#### Remarks

**A.4.1 (angular momentum).** The only place where information is discarded is (RA), and the discarded term is nonnegative. The theorem therefore applies to $L=0$, which is the case of the Pythagorean problem.

**A.4.2 (relation to Yoshida's statements).** These relations are based on the statements of [Y72], [Y73] and [Y74] as transcribed from the original scans.

*Theorem 1 of [Y72].* Its hypotheses are:
- (13) $\rho(0)\ge ad>a\,r(0)$ for some $a>\lambda$;
- (14) $\dot\rho(0)>\sqrt{2MW(\rho(0),d)}$;
- (15) $\dot\rho(0)^2-2MW(\rho(0),d)>2H/\mu+2m_1m_2/(\mu d)$.

For $d=d_1$:
- The right side of (15) is $0$, so (15) follows from (14).
- (13) is equivalent to (C1) together with (C2): take $a\in(\lambda,\rho(t_0)/d_1]$.

So (★) is Theorem 1 of [Y72] with $d=m_1m_2/K$.

*Criterion A′ of [Y74].* With $d=d_1$ it requires $\rho(0)\ge ad_1\ge a\,r(0)$, which allows the boundary case $r(0)=d_1$. That case cannot occur together with (C2) and (C3): Step 2(c) at $t=t_0$ would force $r(t_0)<d_1$.

*Corollary of [Y72].* It has the gap described in A.0, which was corrected in [Y73] and again in [Y74, p. 369]. We do not use it.

**A.4.3 (ordering of the masses).** If $m_1>m_2$, relabel the binary components. With $\lambda$ and $\nu$ interchanged in $W$, the function becomes smaller, and Lemmas A–B fail, for example at $(r,c)=(d,-1)$. The ordering is essential.

**A.4.4 (checkability).** (C1)–(C3) are strict inequalities in $(r,\rho,\dot\rho)$ with exact constants. It suffices to check them for every state in a rigorous enclosure of the true state at a single time $t_0$. The value of $K$ must be the exact energy of the true solution, not the energy recomputed from an enclosure.

**A.4.5 (non-accumulation of collisions).** Version 1 gave here only a sketch of why binary collisions cannot accumulate under (SEP). In version 1.2 this is proved: collisions are isolated by Lemma A.2.2(i), and the solution exists for all time by Proposition A.2.3.

---

### A.5 Application to the Pythagorean problem

**Initial data** [Bu1913; SP1967]. Masses $3,4,5$ start at rest at $(1,3)$, $(-2,-1)$ and $(1,-1)$, respectively.
- The barycentre is at the origin, and the total momentum and angular momentum vanish.
- The mutual distances are $|3\text{–}4|=5$, $|3\text{–}5|=4$ and $|4\text{–}5|=3$.
- Hence, exactly,
$$
H=-\Big(\frac{3\cdot4}{5}+\frac{3\cdot5}{4}+\frac{4\cdot5}{3}\Big)=-\frac{769}{60},\qquad K=\frac{769}{60}.
$$
- The energy is conserved along the (regularized) solution, so this $K$ is the value to use at any later time.

**Parameters.** The escaper is $m_0=3$, and the binary is $m_1=4\le m_2=5$:
$$
M=12,\quad \lambda=\tfrac59,\quad \nu=\tfrac49,\quad \mu=\tfrac94,\quad m=\tfrac{20}9,\quad
d_1=\frac{m_1m_2}{K}=\frac{1200}{769}\approx1.560468140442,
$$
$$
\lambda d_1=\frac{2000}{2307}\approx0.866926744690,\qquad \nu d_1=\frac{1600}{2307}\approx0.693541395752 .
$$

**Conditions.** Let $\mathbf x_3,\mathbf x_4,\mathbf x_5$ be the positions of the masses $3,4,5$. Then
$$
r=|\mathbf x_5-\mathbf x_4|,\qquad \boldsymbol\rho=\mathbf x_3-\frac{4\mathbf x_4+5\mathbf x_5}{9},\qquad \dot\rho=\frac{\boldsymbol\rho\cdot\dot{\boldsymbol\rho}}{\rho},
$$
and the conditions of Theorem (★) read
$$
\text{(C1)}\ r<\frac{1200}{769},\qquad \text{(C2)}\ \rho>\frac{2000}{2307},\qquad
\text{(C3)}\ \dot\rho>\sqrt{24\,W(\rho)},\quad W(\rho)=\frac{5/9}{\rho+1600/2307}+\frac{4/9}{\rho-2000/2307}.
$$
In $W$, the coefficient $5/9$ goes with the **larger** distance $\rho+\nu d_1$, and $4/9$ with the **smaller** distance $\rho-\lambda d_1$.

**Conclusion.** If (C1)–(C3) hold for the true state at some $t_0$, then for all $t\ge t_0$:
- the $4$–$5$ separation stays below $m_1m_2/(K+\tfrac12\mu C^2)<1200/769$;
- the body of mass $3$ recedes at least linearly: $\rho(t)\ge\rho(t_0)+C(t-t_0)$.

Here $\mu=9/4$ and $C$ is defined from the state at $t_0$.

Verifying (C1)–(C3) for the actual trajectory is a separate, computer-assisted step described in the main text. It is not part of this appendix.

**Historical note.** Yoshida [Y72, p. 406] already applied the corollary to this problem at $t=59.72$, using an earlier non-rigorous numerical solution. At that point $r\approx0.42<d_1$, so the corrected condition was also satisfied.

---

### A.6 Independent cross-checks (not part of the proof)

These checks look for errors in the proof above; they do not replace it.
- **Computer algebra (SymPy).** We checked every algebraic identity used in A.1 and A.3: (E) and (EM) derived from the Cartesian Newton equations, (RA) in 2D and 3D, the expansions of Steps 1–3 up to finite order, the closed forms of Step 4, the second proof of Lemma A, and, in version 1.2, the 10 identities of Lemma A.2.1 (`proof_review/s3/s3_algebra_check.py`). Deliberately wrong variants were rejected.
- **Regularized vs. Cartesian integration (A.2, non-rigorous).** For the Pythagorean data, the system (LC) integrated in $s$ and the Newton equations integrated in $t$ agree at $t=30$ to $7.5\times10^{-8}$ (binary separation) and $2.1\times10^{-8}$ ($\boldsymbol\rho$), through a 4–5 approach of $0.0022$; the invariant $I$ drifts by $2\times10^{-12}$ (`proof_review/s3/s3_numeric_sanity.py`).
- **Adversarial numerical search** (50-digit arithmetic, 40 000 points concentrated near the boundaries $\rho\to\lambda d^+$, $r\in\{0,d\}$, $c=\pm1$, for several $\lambda\in[\tfrac12,1)$). We found no violation of Lemmas A–B. The maximum ratio was $1+O(10^{-50})$, attained at $(r,c)=(d,-1)$.
- **Dynamical falsification test.** We took 800 states with masses $3,4,5$ and various energies satisfying (C1)–(C3), and integrated them to $t\le40$. No conclusion of Theorem (★) was violated. Variants with the factor $2$ in (C3) omitted, with (C1) dropped, or with $\lambda,\nu$ interchanged in $W$ did produce violations, including the escaper returning.

---

### References

- **[Bu1913]** C. Burrau, *Numerische Berechnung eines Spezialfalles des Dreikörperproblems*, Astron. Nachr. **195** (1913) 113–118. doi:10.1002/asna.19131950602.
- **[DLMF]** NIST Digital Library of Mathematical Functions, https://dlmf.nist.gov/, Eqs. 18.10.5 and 18.12.11 (accessed 2026-10-12).
- **[Ha]** P. Hartman, *Ordinary Differential Equations*, 2nd ed., Classics in Applied Mathematics **38**, SIAM, Philadelphia, 2002. Chapter II (existence, uniqueness, extension of solutions). The exact theorem numbers were not checked against the text.
- **[KS1965]** P. Kustaanheimo, E. Stiefel, *Perturbation theory of Kepler motion based on spinor regularization*, J. Reine Angew. Math. **218** (1965) 204–219.
- **[LC1920]** T. Levi-Civita, *Sur la régularisation du problème des trois corps*, Acta Math. **42** (1920) 99–144. doi:10.1007/BF02404404.
- **[P1897]** P. Painlevé, *Leçons sur la théorie analytique des équations différentielles, professées à Stockholm (1895)*, Paris, 1897.
- **[SM1971]** C. L. Siegel, J. K. Moser, *Lectures on Celestial Mechanics*, Grundlehren der mathematischen Wissenschaften **187**, Springer, Berlin, 1971 (reprinted in Classics in Mathematics, 1995).
- **[Su1912]** K. F. Sundman, *Mémoire sur le problème des trois corps*, Acta Math. **36** (1912) 105–179. doi:10.1007/BF02422379.
- **[SP1967]** V. Szebehely, C. F. Peters, *Complete solution of a general problem of three bodies*, Astron. J. **72** (1967) 876. doi:10.1086/110355.
- **[Y72]** J. Yoshida, *Improved criteria for hyperbolic-elliptic motion in the general three-body problem*, Publ. Astron. Soc. Japan **24** (1972) 391–408. doi:10.1093/pasj/24.3.391.
- **[Y73]** J. Yoshida, *Erratum: Improved criteria for hyperbolic-elliptic motion in the general three-body problem*, Publ. Astron. Soc. Japan **25** (1973) 285. doi:10.1093/pasj/25.2.285.
- **[Y74]** J. Yoshida, *Improved criteria for hyperbolic-elliptic motion in the general three-body problem. II*, Publ. Astron. Soc. Japan **26** (1974) 367–377.
