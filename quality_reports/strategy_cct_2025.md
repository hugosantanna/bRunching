# Strategy Memo: Implementing CCT-2025 (`method = "cct"`) in bRunching

**Paper:** Caetano, Caetano, Tecchio (November 2025), *Correcting Endogeneity of Treatments with Bunching Using Censoring Strategies.*
**Target parameter:** β(0⁺) := lim_{x↓0} E[(Y(x) − Y(0))/x | X = x] (Eq. 2, p. 5)
**Estimator:** Chen, Dahl, Khan (2005) nonparametric location-scale censored regression for m(Z), plugged into an OLS "generated regressor" representation of Eq. (9), p. 7.

This memo is the build spec for a new `method = "cct"` (or `"cct2025"`) branch of `ccn()`. It is distinct from every existing method in bRunching: the existing methods target a pooled slope β; this one targets the boundary slope β(0⁺) and relies on conditional quantiles of X|Z rather than on per-bin parametric Tobits.

---

## 1. Identification recap (so the coder understands what is being estimated)

Under Assumptions 1–3 (p. 4–5) the identification equation is

```
β(0⁺) = dE[Y | X = 0⁺]/dx − P(X = 0) · Δ / ( E[m(Z)] − E[X] )      (Eq. 9, p. 7)
```

where Δ := E[Y | X = 0] − lim_{x↓0} E[Y | X = x] (Eq. 4, p. 6) and m(Z) := E[X* | Z] in the censoring model X* = m(Z) + η, X = max{0, X*} (Eq. 7–8, p. 7).

The method reduces to identifying m(Z) (more precisely, E[m(Z)]) from conditional quantiles of X | Z under a nonparametric location-scale restriction on η | Z (Assumption 5, p. 10).

Under the additional **Assumption 7** (p. 14) — linearity of E[Y | X] on X > 0 — Eq. (9) collapses to an OLS regression of Y on (1, X, X + 1{X=0}·π̂) with π̂ := (E[m(Z)] − E[X]) / P(X = 0). The coefficient on X is β(0⁺). This is the form we will implement (Section 5, p. 13–15).

---

## 2. Algorithm A: Estimate m̂(Z) via Chen-Dahl-Khan (Appendix C, p. 25–26)

**Inputs.**
- `X` (n×1 numeric, treatment, ≥ 0, with a strictly positive mass at 0)
- `Z` (n×1 categorical; Assumption 9 requires finite support — if the covariate is multivalued continuous the user must pre-discretize, e.g. by hierarchical cluster index as in Section 6; see §7 below)
- `alpha_star` ∈ [0, 1): the known quantile with H⁻¹(α*) = 0 (Assumption 5). The paper typically uses α* = 0.5 (median) or α* = 0.85 (Section 6.2). Default: **α* = 0.5** with automatic fallback if `median(X | Z = z) = 0` for every z (see Remark 4.3 and §7).
- `alpha_grid`: a grid α_1 < … < α_L of quantiles strictly above α* used in the location-scale regression. Paper does not prescribe a specific grid; Chen-Dahl-Khan suggest an equispaced grid. **Default: `alpha_grid = seq(alpha_star + 0.05, 0.95, by = 0.05)`** (TODO: confirm with authors).
- `alpha_s` ∈ (α*, 1): a fixed "scale anchor" quantile used to build Ĉ_ℓ (Eq. in Step 2, p. 26). **Default: `alpha_s = max(alpha_grid)`** — i.e. the top quantile. Paper is silent on the exact choice.
- `zeta_0 > 0`: small trimming constant used to drop small quantile estimates in Eq. (14). **Default: `zeta_0 = 1e-3` × sd(X)** (TODO: paper only says "fix a small constant"; consult CDK 2005 for guidance).
- `zeta_1 > 0`: kink location of the smooth weight w(·). CDK and Buchinsky-Hahn (1998) use `zeta_1 = zeta_0`. **Default: `zeta_1 = zeta_0`**.

**Step 1 — Conditional quantile estimation (Eq. 13, p. 25).**
For each unique value z ∈ supp(Z) and each α ∈ {α*} ∪ alpha_grid, solve

```
q̂(z; α) = argmin_q  Σ_i 1{Z_i = z} · ρ_α(X_i − q)
```

with ρ_α(u) = α|u| + (2α − 1)·u·1{u<0}. Because Z is finite, this is just the (unconditional) sample α-quantile of `X[Z == z]`. Implementation: `quantile(X[Z == z], probs = α, type = 1)` is the check-loss minimizer; `type = 1` or a custom solver (linear program) both work. A linear-programming version via `quantreg::rq(X ~ 1, tau = α, subset = Z == z)` is the cleanest.

**Step 2 — Compute Ĉ_ℓ, ℓ = 1, …, L (p. 26).**
The smooth weight function from Buchinsky-Hahn (1998), used by CDK:

```
w(q) = ( e^{q − 2ζ_1}/(1 + e^{q − 2ζ_1}) − e^{−ζ_1}/(1 + e^{−ζ_1}) )
     × ( (2 + e^{ζ_1} + e^{−ζ_1}) / (e^{ζ_1} − e^{−ζ_1}) ) · 1{ζ_1 < q < 3ζ_1}
     + 1{q > 3ζ_1}
```

(Appendix C, bottom of p. 26.) This is a C¹ approximation to 1{q > 0} that vanishes below ζ_1 and equals 1 above 3ζ_1.

Then, for each ℓ ∈ {1, …, L}:

```
Ĉ_ℓ = ( Σ_i w(q̂(Z_i; α*)) )⁻¹ · Σ_i [ w(q̂(Z_i; α*)) · (q̂(Z_i; α_ℓ) − q̂(Z_i; α*)) / (q̂(Z_i; α_s) − q̂(Z_i; α*)) ]
```

Define Ĉ_ℓ := (1, Ĉ_ℓ)' — a 2-vector.

**Step 3 — Weighted OLS for m̂(z) (Eq. 14, p. 26).**
For each z with q̂(z; α_ℓ) ≥ ζ_0 for at least one ℓ:

```
m̂(z) = κ_1' · ( Σ_ℓ 1{q̂(z; α_ℓ) ≥ ζ_0} · Ĉ_ℓ Ĉ_ℓ' )⁻¹ · ( Σ_ℓ 1{q̂(z; α_ℓ) ≥ ζ_0} · Ĉ_ℓ q̂(z; α_ℓ) )
```

with κ_1 = (1, 0)'. Equivalent R code: stack rows ℓ with valid q̂(z; α_ℓ) ≥ ζ_0 into a design matrix D = [Ĉ_ℓ'] and response y = q̂(z; α_ℓ); run `lm(y ~ 0 + D)`; return first coefficient as m̂(z).

**Step 4 — Return.** A vector `m_hat` of length n with `m_hat[i] = m̂(Z_i)` (or a per-z table joined back on Z).

**Short-circuit (Remark 4.2, p. 11).** If `min_z q̂(z; α*) > 0` for all z (i.e. all α*-quantiles are strictly positive), the coder should skip the location-scale regression and set

```
m̂(z) = q̂(z; α*).
```

Expose this via an argument `locscale = c("auto", "on", "off")`, default `"auto"` (switch off when all α*-quantiles are positive; see §6.1 p. 16 — TV application uses this relaxation).

**Short-circuit (Remark 4.3, p. 12).** If there exists z₁ with `P(X = 0 | Z = z₁) = 0`, then m̂(z₁) = mean(X | Z = z₁) and the H⁻¹(α*) = 0 restriction is unnecessary. This is a useful diagnostic but not strictly needed in v1; flag as a future enhancement.

---

## 3. Algorithm B: Compute β̂(0⁺) (Section 5, p. 13–15)

**Semiparametric form (recommended default, Assumption 7).**

Inputs: `Y`, `X`, `m_hat` (from Algorithm A), optional `W` (controls, Remark 5.2).

1. Compute `pi_hat = sum(m_hat - X) / sum(X == 0)`.
2. Build the generated regressor `Xtilde = X + (X == 0) * pi_hat`.
3. Fit OLS: `lm(Y ~ 1 + X + Xtilde)` — or with controls, `lm(Y ~ 1 + X + Xtilde + W)`.
4. Return the coefficient on `X` as β̂(0⁺).

**Interpretation of the three coefficients.** Per Proposition 5.1 (p. 15):
- intercept = γ_0
- coefficient on X = γ_1 − Δ/π (this is β̂(0⁺))
- coefficient on Xtilde = Δ/π

So the OLS trick cleanly separates the linear-slope part from the bias correction without ever explicitly estimating Δ or dE[Y | X = 0⁺]/dx.

**Nonparametric fallback.** If the user sets `linear = FALSE`:
- Estimate dE[Y | X = 0⁺]/dx with a boundary local-linear regression of Y on X for X > 0 at x = 0. Use `locpol::locLinSmootherC` or `KernSmooth::locpoly` with `deriv = 1` and a bandwidth from `nprobust::lprobust` (default bw rule: MSE-optimal, available via `nprobust::lpbwselect`). Take the right-limit fit at x = 0.
- Estimate `lim_{x↓0} E[Y | X = x]` with the same local-linear estimator at `deriv = 0`.
- Compute `Delta_hat = mean(Y[X == 0]) - (local-linear intercept at 0+)`.
- Plug into Eq. (9): β̂(0⁺) = (deriv estimate) − mean(X == 0) · Delta_hat / (mean(m_hat) − mean(X)).

This fallback is NOT the paper's recommended path — the paper explicitly argues the linear form is what practitioners should use — but exposing it is useful for robustness and matches the discussion on p. 13–14. TODO: decide with author whether to ship in v1 or defer.

---

## 4. Inference: pairs bootstrap (Section 6, Tables 1–2, p. 17, 19)

The paper reports "**Bootstrapped standard errors (500 iterations)**" for both applications. No other inference method is used in the empirical sections. Plug-in analytic SEs are derived in Proposition 5.1 and Proposition C.1 but the paper defers implementation ("Variance plug-in é óbvio a forma, condições seguem de e.g. Newey and McFadden 1994" — p. 27, note this is an author-side TODO in the draft).

**Scheme.** Nonparametric pairs bootstrap.

```
for b in 1..B:
  idx_b <- sample(1:n, n, replace = TRUE)
  re-run Algorithm A on (X[idx_b], Z[idx_b])          # re-estimate m_hat_b
  re-run Algorithm B on (Y[idx_b], X[idx_b], m_hat_b) # re-estimate beta_hat_b
  store beta_hat_b, and E[m(Z)]_b = mean(m_hat_b), E[X]_b, pi_hat_b
se <- sd(beta_hat_b over b)
```

**Defaults:** `B = 500` (paper), cluster = none, report percentile + Normal CI. Expose `B` and `cluster` as arguments. **Do not** cache `m_hat` across bootstrap iterations — every iteration must re-estimate.

**Caveat.** Bootstrap validity for censored quantile regression of the CDK form is not trivial (Chernozhukov-Hong 2002 discuss alternatives). The paper sidesteps this by simply reporting bootstrap SEs; we follow suit but add a note in the documentation. TODO: consider subsampling (Politis-Romano-Wolf) as a robustness option.

**Analytical SE (deferred).** Proposition 5.1 (p. 15) gives the sandwich form:

```
√n(γ̂ − γ) →_d N(0, Σ + (Δ/π)² V_π Γ Γ' − (Λ Γ' + Γ Λ') Δ/π)
```

with Σ, V_π, Γ, Λ defined on p. 15. Proposition C.1 (p. 27) gives the ψ(·) influence function for m̂. These are implementable but nontrivial — punt to v2.

---

## 5. Monte Carlo for unit testing

The paper does **not** contain an explicit Monte Carlo section with a named DGP. This is a gap — the coder will need to construct a test DGP consistent with Assumptions 1, 3, 4, 5 and 7. Recommended setup:

```
n      = 5000
Z      ~ sample(1:K, n, replace = TRUE), K = 20
mu(z)  = -0.5 + 0.15 * z                              # so E[X* | Z=z] = mu(z)
sigma(z) = 1                                          # homoskedastic location-scale, H = Normal
eta    ~ Normal(0, sigma(z))
X_star = mu(Z) + eta
X      = pmax(X_star, 0)                              # censoring at 0
# Linear outcome with known beta(0+)
beta_true  = 1.5
delta_true = 0.8                                      # selection slope
Y0_star = 2 + delta_true * X_star + rnorm(n, sd = 1)
Y       = Y0_star + beta_true * X * (X > 0)           # beta(0+) = beta_true
```

Verify:
- Bunching rate ≈ mean(X == 0) around 10–30%.
- α* = 0.5 gives H⁻¹(0.5) = 0 by symmetry of Normal.
- `ccn(Y, X, Z, method = "cct")` should recover β̂(0⁺) ≈ 1.5 with SE ≈ O(1/√n).

Save this as `tests/testthat/test-cct-mc.R` with a tolerance of 3·SE.

**Note.** Paper may add a Monte Carlo appendix in the revision. Watch for updates and sync the test DGP. TODO.

---

## 6. Empirical application targets (Section 6, p. 15–19)

The paper uses these as "replication targets" for the coder.

**Application 1 — TV hours → non-cognitive skills (Section 6.1, p. 15–17).**
- Data: CDS-PSID 1997, 2002, 2007 waves (same as Caetano et al. 2023c).
- X = weekly TV hours; bunching at 0 ≈ 5%.
- Z options: (a) household income percentile, (b) 50 clusters, (c) 100 clusters (hierarchical clustering on child + family + environment covariates with Gower distance + Ward linkage; see footnote 13).
- Relaxation: because median(X | Z = z) > 0 for all z, the paper **drops the location-scale assumption** (Remark 4.2). So `m̂(z) = median(X | Z = z)` here — the `locscale = "off"` branch.
- Reported results (Table 1, p. 17):

  | Z                | E[X*]  | β̂         |
  |------------------|--------|------------|
  | Household Income | 12.24  | −0.8169    |
  |                  | (0.19) | (0.3029)   |
  | 50 clusters      | 12.05  | −0.7720    |
  |                  | (0.20) | (0.2789)   |
  | 100 clusters     | 12.09  | −0.7786    |
  |                  | (0.20) | (0.3030)   |

**Application 2 — working hours around retirement → overall health (Section 6.2, p. 17–19).**
- Bunching > 50% (heavy). Because of this, α* = 0.5 fails (`median(X | Z) = 0` for many z). Paper uses **α* = 0.85**, with α_1 = 0.85, α_2 = 0.90.
- Z = 50 or 100 clusters on covariates listed in footnote 16.
- Reported results (Table 2, p. 19):

  | Z            | E[X]  | E[X*]   | β̂(0⁺)  |
  |--------------|-------|---------|---------|
  | 50 clusters  | 11.61 | −22.72  | 0.1158  |
  |              | (0.08)| (4.87)  | (0.0138)|
  | 100 clusters | 11.61 | −22.60  | 0.1158  |
  |              | (0.08)| (4.87)  | (0.0138)|

**Replication test:** If the CDS-PSID and HRS-like data are accessible (they are not in the bRunching repo), we can target these numbers in `tests/testthat/test-cct-replication.R`. Otherwise treat them as documentation-only benchmarks.

---

## 7. Assumptions, failure modes, and contrast with `"symmetric"`

| Assumption | What it says | Fails when… | How it compares to `"symmetric"` |
|------------|--------------|-------------|-----------------------------------|
| A1 Local support (p. 4) | X has density on [0, ε₀] | Point mass above 0 or gap after 0 | Same requirement |
| A2 Local differentiability (p. 4) | x ↦ E[Y(x) − Y(0) \| X = x] is right-differentiable at 0 | Kinked or discontinuous effect profile at margin | `"symmetric"` needs global constancy — much stronger |
| A3 Local linear selection (p. 5) | E[Y(0) \| X*] linear on {X* ≤ 0} | Nonlinear selection at margin | `"symmetric"` replaces this with symmetry of the η density — a distributional restriction rather than a conditional-mean one |
| A4 Positive known quantile (p. 10) | P(q(Z; α*) > 0) > 0 | All observations bunch regardless of Z — no identifying variation | `"symmetric"` does not need Z at all |
| A5 Location-scale in η \| Z (p. 10) | F_{η\|Z}(v) = H(v / σ(Z)) with H⁻¹(α*) = 0 | Different covariates have differently-shaped conditional distributions (e.g. bimodal for some z and unimodal for others) | `"symmetric"` requires symmetric density of X*; A5 is strictly weaker but still restrictive |
| A7 Linear E[Y \| X > 0] (p. 14) | E[Y \| X] linear on X > 0 | Curved dose-response (paper itself opts out in §6.2 — see p. 18) | No analog: `"symmetric"` fits a linear model globally |

**When `"cct"` is preferred over `"symmetric"`:** when the density of X* is clearly asymmetric but its conditional-on-Z shape looks similar across z (Figure 1, p. 11). Typical case: right-skewed treatments like TV hours, working hours, cigarettes.

**When `"symmetric"` is preferred:** when no credible covariate Z exists that the researcher can defend as generating the same-shape conditional distribution, but the unconditional X* distribution looks symmetric.

**When both fail:** highly skewed and no valid Z. User should fall back to `"tobit"` / `"het_tobit"` with its distributional assumption.

---

## 8. User-facing API

Recommended addition to `ccn()`:

```r
ccn(outcome, treatment, instrument = NULL, data,
    method = c("symmetric", "uniform", "tobit", "het_tobit", "naive", "cct"),
    # CCT-specific:
    alpha_star  = 0.5,
    alpha_grid  = NULL,      # default seq(alpha_star + 0.05, 0.95, by = 0.05)
    alpha_s     = NULL,      # default max(alpha_grid)
    locscale    = c("auto", "on", "off"),
    zeta_0      = NULL,      # default 1e-3 * sd(X)
    zeta_1      = NULL,      # default zeta_0
    linear      = TRUE,      # Assumption 7; FALSE → boundary local-linear
    controls    = NULL,
    boot_B      = 500,
    ...)
```

Note: when `method = "cct"`, the `instrument` argument is interpreted as Z — the covariate that indexes conditional quantiles — not as an instrument in the IV sense (see Remark 3.1, p. 8, which explicitly disclaims the causal interpretation of Z). The doc must make this clear.

`ccn(...)$beta` returns a single scalar β̂(0⁺) (not a vector as in existing methods). `ccn(...)$m_hat` stores the per-z m̂ estimates.

---

## 9. Ambiguities flagged for the coder (TODOs)

1. **Quantile grid size and spacing.** Paper does not prescribe `alpha_grid`. CDK-2005 originals use a small number of grid points. Default `seq(α* + 0.05, 0.95, 0.05)` is a guess — benchmark against CDK replication code if available.
2. **`alpha_s` choice.** Paper uses `alpha_s ∈ (α*, 1)` but never states how to pick it. Guess: top grid point. Expose as argument.
3. **`zeta_0`, `zeta_1` choices.** Paper says "fix a small constant"; CDK give the smooth-weight formula but no numeric default. Pick a scale-invariant default (`1e-3 · sd(X)`), expose as argument, document the sensitivity.
4. **Nonparametric derivative fallback (Algorithm B, `linear = FALSE`).** Bandwidth selection is not discussed. Use `nprobust::lpbwselect` MSE-optimal default; warn that boundary bias correction is active.
5. **Multi-dimensional W (controls).** Remark 5.2 says "include W linearly in the OLS." Straightforward for continuous W; for factor W the `1(X > 0) · W` block needs to be linearly independent of (1, X, Xtilde) — warn if collinear.
6. **Monte Carlo in paper.** None explicit; we build our own. Cross-check with the authors on whether a Monte Carlo will be added before CRAN submission.
7. **Analytic SE.** Proposition 5.1 gives the sandwich; Proposition C.1 gives the ψ. Not implementing in v1 — defer.
8. **Overidentification.** Remark 4.1 (p. 11) discusses using multiple z₀ or multiple (α_1, α_2) pairs. The CDK estimator already averages across many α_ℓ via Eq. (14); but an outer average across z₀ values is not implemented in our spec. Punt unless performance demands.
9. **Remark 4.3 case (`P(X = 0 | Z = z₁) = 0`).** Useful practical diagnostic — when any z-cell has zero censoring, we can skip location-scale and use sample mean in that cell. Build as a fast-path if we have time.
10. **Author-side draft artifacts.** The PDF has red in-line TODOs from the authors (p. 12: "??We illustrate this case in an application in Section XXXX/or cite paper"; p. 15: "Projection as in Newey's 1994 ECMA"; p. 27: Portuguese note about the variance plug-in). These will likely change in the next revision. **Sync with authors before CRAN release.**

---

## 10. Deliverables for the coder

- `R/cct.R` — entry points: `cct_m_hat()`, `cct_beta()`, `cct_bootstrap()`.
- `R/ccn.R` — add `method = "cct"` branch that dispatches to the above.
- `R/censored_exp.R` — no change (that module computes pooled censored expectations; CCT does not use them).
- `man/` — new Rd for `method = "cct"` documenting the above API.
- `tests/testthat/test-cct-mc.R` — the Monte Carlo from §5, tol = 3·SE.
- `tests/testthat/test-cct-edge.R` — zero-bunching cell (Remark 4.3), all-quantiles-positive short-circuit (Remark 4.2), collinear W check.
- `vignettes/cct.Rmd` — worked example mirroring Section 6.1 as far as synthetic data allows.

---

## 11. Referee-style concerns we should pre-empt

1. "Why Chen-Dahl-Khan rather than Lewbel-Linton (2002) or Powell (1984)?" — Paper explicitly argues CDK's location-scale restriction fits the empirical applications and avoids parametric distributional assumptions (Section 1, p. 2; Section 4, p. 9). We follow the paper. If a user wants alternatives, that is future work — document it.
2. "Bootstrap validity for censored quantile regression." — Chernozhukov-Hong (2002) is the standard reference. Our SEs are paper-consistent but not formally justified. Note in the documentation.
3. "Sensitivity of β̂(0⁺) to α*, α_grid, ζ." — Build a sensitivity function `cct_sensitivity()` (v2) that sweeps these; document in the vignette.

End of memo.
