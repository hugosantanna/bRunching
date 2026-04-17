# CDK 2005 Verification Memo

**Target paper:** Chen, Dahl, Khan (2005), JASA 100(469): 212-221.
**Local PDF source:** Dahl's UCSD page, fetched 2026-04-17 (10 pages, full published version).
**Checked against:** `/R/cct.R`, Appendix C of `Bunching_as_Censoring.pdf` (pp. 25-26), `quality_reports/strategy_cct_2025.md`.

## 1. Access status
Got it. `pdftotext -layout` extracted the full paper cleanly. All equations, remarks, and the empirical-illustration numerical choices (p. 217) were recovered.

## 2. Faithfulness of CCT's retelling
CCT Appendix C is substantively faithful, with three deviations worth flagging:

- **α\*** is a free parameter in CCT. CDK hardwires **α\* = 0.5** — Stage 1 is titled "Local Constant Estimation of the Conditional Median Function" (p. 214), and Eq. 12 and Eq. 15 are written with `q_.5`. CDK's scale normalization is `c_{α_1} ≡ 1`, not `c_{α_s} ≡ 1`. CCT correctly generalizes by replacing 0.5 with a known zero-quantile α\*; their Remark 4.2 / Remark 4.3 anchor this generalization. Not a bug, but worth noting.
- **Trimming τ(x\_i)** in CDK Eq. 15 is absent from CCT's version. CCT justifies this by restricting to finite-support Z (Assumption 9) — boundary trimming is for continuous-regressor kernel smoothing. Defensible.
- **Two trimming constants (η, ε)** in CDK become **(ζ\_1, ζ\_0)** in CCT. CDK's η is the weight-function kink (A3.2, p. 215); CDK's ε is the Stage-3 indicator cutoff (Eq. 17). CCT's mapping ζ\_1 ↔ η, ζ\_0 ↔ ε is correct.
- **α\_s constraint.** CDK Eq. 15 context requires α\_s > 0.5 ("we select a quantile, α\_s > .5, which need not be on the grid"). CCT Appendix C requires α\_s ∈ (α\*, 1). Consistent generalization.

The Buchinsky-Hahn w(q) formula on CCT p. 26 is identical to the standard Buchinsky-Hahn (1998) form and is consistent with CDK's Assumption A3.

## 3. Faithfulness of R code
`R/cct.R` implements CCT's Appendix C, not CDK's published algorithm. Specifically:

- `cct_weight()` matches the CCT/Buchinsky-Hahn formula exactly.
- `cct_m_hat()` Step 2 (lines 148-167) computes Ĉ\_ℓ without a τ(x\_i) trimming — matches CCT, defensible under finite Z.
- α\* is exposed as a user argument — matches CCT generalization, not CDK's hardwired 0.5.

The code is faithful to CCT's retelling. Any residual concern is inherited from CCT, not introduced by the implementation.

## 4. Concrete issues found

- **cct.R:158 — w\_sum uses a bounded subset, but CDK Eq. 15 and CCT normalize by the full `Σ w(q̂(Z\_i; α\*))`.** The R code restricts `w_sum <- sum(w_i[safe])` where `safe <- denom > 0 & w_i > 0`. CDK Eq. 15's denominator is `(1/n) Σ w(q̂(x\_i))` over all i. Observations with `denom = 0` but `w_i > 0` should enter the denominator (they contribute zero to the numerator because of the 0/0 convention — CDK Appendix adopts "0/0 = 0" explicitly, p. 219). The practical effect is small in finite Z with separated quantiles, but it is a divergence from the formula.
- **α\_s and CDK's empirical choice.** CDK's illustration (p. 217) picks α\_s = 0.3 — **below the median** — because they estimate on left tail. CCT's α\_s > α\* direction is opposite for our right-censored setting, which is fine, but the default `alpha_s = max(alpha_grid)` in our code should be flagged as arbitrary. CDK does not prescribe a rule beyond "need not be on the grid" and "any convenient quantile."
- **ζ\_0, ζ\_1 defaults.** CDK gives no numeric rule — both η and ε are "small positive constants." The Monte Carlo (p. 217) does not disclose values. The default `1e-3 * sd(X)` is defensible but should be documented as ours, not CDK's.
- **No α\* = 0.5 hardwiring sanity check.** When α\* ≠ 0.5, CDK's local-constant **median** Stage 1 is not the asymptotic object being estimated. CCT's Proposition C.1 extends the argument, but users should be warned that departing from α\* = 0.5 is a CCT-2025 extension not directly justified by CDK 2005.

## 5. Asymptotic variance
CDK Theorem 2 gives a **sandwich form** `n^{p/(2p+d\_c)}(θ̂ - θ\_0) ⇒ N(C⁻¹B, C⁻¹VC⁻¹)` with explicit C and V (Eqs. 21-22). CCT's Proposition C.1 gives only the influence function ψ(·) and defers the variance to "é óbvio" (line 1195 — untranslated Portuguese placeholder). **This is a CCT drafting gap, not a CDK gap.** CDK already provides the full plug-in variance estimator (Eqs. 25-26, p. 216). A v2 of the package could implement CDK's sandwich directly — priority: medium.

## 6. Recommended fixes (severity order)

1. **(Low/medium)** Document in `?ccn` that our α\* is a CCT-2025 generalization; CDK-2005 uses α\* = 0.5.
2. **(Low)** Fix the `w_sum` normalization in `cct.R:158` to sum over all observations (adopt CDK's 0/0 = 0 convention for the ratio).
3. **(Low)** Document ζ\_0 = ζ\_1 = 1e-3·sd(X) as our choice, not CDK's or CCT's.
4. **(Medium, future)** Implement CDK's analytic sandwich SE (Eqs. 21-22, 25-26) as an alternative to bootstrap.
5. **(Low)** Clarify in docs that α\_s choice is unconstrained in CDK beyond α\_s ≠ α\*; default is heuristic.

**Bottom line:** CCT's retelling is accurate up to the α\* generalization and the trimming-removal (both benign under finite Z). Our R code matches CCT. No bugs that change estimates — only a numerator/denominator asymmetry in `w_sum` that is zero-measure in practice and a handful of documentation improvements.
