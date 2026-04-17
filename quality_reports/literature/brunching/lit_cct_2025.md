# Literature Notes — CCT 2025 Implementation

**Target project:** bRunching (R package implementing CCT 2025 extension of CCN 2023)
**Date:** 2026-04-17
**Scope:** Focused pointers for implementation. Not a full lit review.

---

## Task 1 — Chen, Dahl, Khan (2005)

### Exact citation

Chen, Songnian, Gordon B. Dahl, and Shakeeb Khan. 2005. "Nonparametric Identification and Estimation of a Censored Location-Scale Regression Model." *Journal of the American Statistical Association* 100 (469): 212–221.

- **DOI:** 10.1198/016214504000001835
- **JASA landing page:** https://www.tandfonline.com/doi/abs/10.1198/016214504000001835
- **Open-access PDF (author's site, UCSD):** https://econweb.ucsd.edu/~gdahl/papers/censored-regression.pdf

### Estimator summary

CDK assume the latent outcome satisfies a **nonparametric location-scale model**:

  Y* = m(X) + σ(X)·ε,

with ε independent of X, ε normalized to location zero and scale one, and observed Y = max(Y*, c) for a known (fixed) censoring point c (WLOG c = 0). The conditional distribution of Y* | X therefore belongs to a nonparametric location-scale family whose shape H (the CDF of ε) is unknown but fixed across X.

**Identification lever.** When m(X) ≥ c with positive probability (i.e., the location function crosses the censoring point on a set of positive measure), m(·) is identified on the entire support of X — including the region where m(X) < c and conditional censoring probability is 1.

**Estimating equation.** The estimator is built from **conditional quantile regressions at higher quantiles** (quantiles above the censoring threshold, where they coincide with uncensored conditional quantiles of Y*). Given any two such quantile levels α₁, α₂, the location-scale structure yields a linear system that, combined with one boundary condition (H⁻¹(α*) = 0 for a known α*), recovers m(X) at every X — including in the interior censoring region. CDK implement this with **local polynomial** quantile estimates, achieve the optimal nonparametric rate for m(·), and derive a limiting normal distribution. Using more than two quantiles produces overidentification and efficiency gains.

The estimator is **semi-nonparametric**: nonparametric in both the shape H and the scale σ(X), with no symmetry or parametric distributional assumption on the error.

### How CCT 2025 uses it

In CCT's Proposition 4.1 (pp. 11 of the draft), they adapt CDK to the bunching setting: the "censoring point" is the bunching threshold (zero), X* is the latent selection variable, and E[m(Z)] is identified by a linear combination of two conditional quantiles of X | Z at levels α₁(Z), α₂(Z) above the known bunching quantile α*. Appendix D.1 of CCT reproduces the CDK identification argument adapted to their assumption set.

### Existing code / implementations

No official public R, Stata, or Python implementation of CDK (2005) appears to exist. Searches of CRAN, GitHub, SSC, and the authors' pages (Gordon Dahl at UCSD, Shakeeb Khan at Boston College, Songnian Chen at HKUST) return no replication repo for the JASA paper itself.

Partial / tangentially-related code:
- `quantreg` (Koenker) — R package for conditional quantile regression, the workhorse for any reimplementation
- `np` / `npqreg` — Hayfield & Racine, nonparametric quantile regression in R
- Downstream method (Heuchenne & Van Keilegom 2009, AISM, "Estimation in nonparametric location-scale regression models with censored data," https://link.springer.com/article/10.1007/s10463-009-0219-3) — also no public code

Implication for bRunching: a CDK estimator must be coded from scratch, likely layered on `quantreg::rq` (or local-linear quantile via `np`) plus the linear-system inversion in Proposition 4.1.

---

## Task 2 — CCT 2025 replication code

### Paper identity (pulled from the local PDF)

- **Title (actual):** "Correcting Endogeneity of Treatments with Bunching Using Censoring Strategies"
- **Authors:** Carolina Caetano (Univ. of Georgia), Gregorio Caetano (Univ. of Georgia), Otávio Tecchio (MIT)
- **Date on draft:** November 2025
- **Local copy:** `/Users/hsantanna/Library/Mobile Documents/com~apple~CloudDocs/work/research/brunching/Bunching_as_Censoring.pdf`

### Code availability — status

**No public replication code found as of 2026-04-17.** Specifically checked:

- **Gregorio Caetano's GitHub** — https://github.com/GregorioCaetano — 4 public repos (`Bunching`, `Robust-RDD-with-Covariates`, `Multivariate-RDD`, `DummyTest`), none for the 2025 CCT paper. Repos for the 2023 CCN JBES paper and the 2021 Dummy Test FEDS paper exist but are for predecessor papers.
- **Carolina Caetano's website** — carolinacaetano.net — SSL cert mismatch blocked WebFetch; direct browser visit recommended, but Google/Bing searches do not surface a CCT 2025 listing there.
- **Otávio Tecchio's MIT page** — https://economics.mit.edu/people/phd-students/otavio-tecchio — profile only, no papers/code listed (2nd-year PhD student).
- **Blueprint Labs profile** — https://blueprintlabs.mit.edu/team/otavio-tecchio/ — profile only.
- **SSRN / NBER / FEDS searches** — no posting under the current title, the alternate "Bunching as Censoring" title, or any CCT combination. The 2020 CCN FEDS paper (2020-080) and 2021 Dummy-Test FEDS paper (2021-068) exist; no CCT 2025 analog yet.
- **TandF / JBES** — not yet submitted or not yet in the online-first queue.

### Precursor code worth cannibalizing

The 2023 CCN paper (Caetano, Caetano, Nielsen, JBES 42(3): 851-863) is the direct methodological ancestor of CCT 2025 and **has full Stata replication code**:

- **GitHub:** https://github.com/GregorioCaetano/Bunching
- **JBES:** https://www.tandfonline.com/doi/abs/10.1080/07350015.2023.2252471
- **FEDS working paper:** https://www.federalreserve.gov/econres/feds/files/2020080pap.pdf
- **SSRN:** https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3699644

This repo is what bRunching's current implementation is derived from. The CCT 2025 paper reuses the same CDS-PSID TV-watching application (Section 6.1 of the draft, same data/sample as CCN 2023), so the data-handling side of the Stata code maps cleanly; only the CDK-based first step (quantile-based m(Z) estimator) is new in CCT 2025 and has no public code.

### Recommended next step for the bRunching package

Ask Tecchio/Caetano directly (hsantanna88@gmail.com has direct contact with Mazetto → Caetano lab). A draft-stage Stata or Matlab prototype almost certainly exists internally — it would be faster than rebuilding CDK from the JASA paper alone.

---

## Sources

- [Chen, Dahl, Khan (2005) — JASA landing](https://www.tandfonline.com/doi/abs/10.1198/016214504000001835)
- [Chen, Dahl, Khan (2005) — OA PDF, Dahl UCSD site](https://econweb.ucsd.edu/~gdahl/papers/censored-regression.pdf)
- [CCN 2023 JBES](https://www.tandfonline.com/doi/abs/10.1080/07350015.2023.2252471)
- [CCN 2020 FEDS working paper](https://www.federalreserve.gov/econres/feds/files/2020080pap.pdf)
- [CCN Stata repo (ancestor of bRunching)](https://github.com/GregorioCaetano/Bunching)
- [Gregorio Caetano GitHub profile](https://github.com/GregorioCaetano)
- [Otávio Tecchio MIT profile](https://economics.mit.edu/people/phd-students/otavio-tecchio)
- [Heuchenne & Van Keilegom (2009), downstream nonparametric location-scale with censoring](https://link.springer.com/article/10.1007/s10463-009-0219-3)
