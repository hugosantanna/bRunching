{smcl}
{* *! version 0.1.0  2026-04-17}{...}
{vieweralsosee "[R] regress" "help regress"}{...}
{vieweralsosee "[R] tobit" "help tobit"}{...}
{viewerjumpto "Syntax" "brunching##syntax"}{...}
{viewerjumpto "Description" "brunching##description"}{...}
{viewerjumpto "Options" "brunching##options"}{...}
{viewerjumpto "Remarks" "brunching##remarks"}{...}
{viewerjumpto "Examples" "brunching##examples"}{...}
{viewerjumpto "Stored results" "brunching##results"}{...}
{viewerjumpto "References" "brunching##references"}{...}
{title:Title}

{phang}
{bf:brunching} {hline 2} Bunching-as-censoring estimator (CCN 2024 and CCT 2025)

{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:brunching}
{depvar} {it:indepvar}
{ifin}
{cmd:,} {opt method(name)}
[{it:options}]

{synoptset 26 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Main}
{synopt:{opt method(name)}}one of {opt naive}, {opt uniform}, {opt tobit},
 {opt het_tobit}, {opt symmetric}, or {opt cct}.{p_end}
{synopt:{opt z(varname)}}conditioning variable (required for all methods other than {opt naive}; interpreted as Z in the CCT 2025 sense when {opt method(cct)}).{p_end}
{synopt:{opt controls(varlist)}}control variables included additively in the
outcome regression.{p_end}
{synopt:{opt nbins(#)}}number of equal-frequency bins for discretizing {opt z()} in the CCN methods; default is 10.{p_end}
{synopt:{opt swap(name)}}fallback when a bin has >=50% bunching in
{opt method(symmetric)}; one of {opt tobit} (default) or {opt het_tobit}.{p_end}

{syntab:CCT-specific}
{synopt:{opt alpha_star(#)}}known quantile where the conditional quantile of the latent variable equals 0; default 0.5.{p_end}
{synopt:{opt alpha_grid(numlist)}}grid of quantile levels strictly in (alpha_star, 1); default is 8 equally spaced points in (alpha_star, 0.99].{p_end}
{synopt:{opt alpha_s(#)}}scale anchor in (alpha_star, 1); default is max(alpha_grid).{p_end}
{synopt:{opt zeta_0(#)}}trimming constant; default 1e-3 * sd(X).{p_end}
{synopt:{opt zeta_1(#)}}kink location of the Buchinsky-Hahn weight; default zeta_0.{p_end}
{synopt:{opt locscale(name)}}Remark 4.2 short-circuit: {opt auto} (default), {opt on}, or {opt off}.{p_end}
{synopt:{opt boot(#)}}pairs-bootstrap replications; default 0 (no bootstrap). Use Stata's {cmd:bootstrap} prefix for more flexible inference.{p_end}
{synopt:{opt seed(#)}}random seed for the internal bootstrap.{p_end}
{synoptline}


{marker description}{...}
{title:Description}

{pstd}
{cmd:brunching} estimates the causal effect of a treatment variable censored
at zero on an outcome, correcting for the bunching bias that arises when
ordinary least squares treats the bunch mass at zero as a uniform
behavioral response.  This is a Stata port of the R package {bf:bRunching}.

{pstd}
Two estimands are supported:

{phang2}
{bf:Global slope} (CCN 2024).  {cmd:method(naive)}, {cmd:method(uniform)},
{cmd:method(tobit)}, {cmd:method(het_tobit)}, {cmd:method(symmetric)}
construct a censored-expectation correction term from a conditioning
variable {opt z()} and run a corrected OLS regression of the outcome on
the treatment, the correction term, and bin fixed effects.  The
correction term identifies the global slope of Y on the latent
treatment X*.

{phang2}
{bf:Local boundary slope beta(0+)} (CCT 2025).  {cmd:method(cct)} runs
the Caetano-Caetano-Tecchio (2025) estimator: a nonparametric
location-scale regression of X on Z via Chen-Dahl-Khan (2005),
followed by the Proposition 5.1 linear form for {it:beta(0+)} :=
lim_{x|0} E[(Y(x) - Y(0))/x | X = x].

{pstd}
The canonical JBES replication (Caetano, Caetano, Nielsen 2024) uses
"reg Y X i.Z_K ..., noconstant" which omits the base level and forces
the first-bin intercept to zero.  This Stata port uses the equivalent
parameterization with K-1 bin dummies and an intercept, matching the R
package's saturated no-constant form numerically (identical beta and
delta to 10+ decimals).


{marker options}{...}
{title:Options}

{phang}
{opt method(name)} selects the estimator.  For the CCN 2024 methods,
{opt naive} is uncorrected OLS with bin fixed effects; {opt uniform}
imposes a uniform symmetry on the censored mass; {opt tobit} runs a
pooled Tobit with bin dummies to recover the censored expectation;
{opt het_tobit} runs a separate Tobit per bin; {opt symmetric} is the
paper's headline method using tail symmetry with a {opt swap()}
fallback in bins with >=50% bunching.

{phang}
{opt z(varname)} names the conditioning variable.  When {opt method(cct)}
is chosen, {it:varname} is interpreted as Z in Assumption 5 of
Caetano-Caetano-Tecchio (2025) - NOT as an IV-style instrument.

{phang}
{opt nbins(#)} controls the discretization of Z for the CCN methods.
When Z has fewer unique values than {opt nbins()}, it is used as-is.
Ignored for {opt method(cct)}.

{phang}
{opt alpha_star(#)} is the known quantile at which the conditional
quantile of the latent X* equals zero.  Assumption 4 requires
alpha_star > max_z P(X = 0 | Z = z).  For heavy bunching set
alpha_star higher (e.g., 0.85 in the CCT 2025 Section 6.2 application
where overall bunching is around 50%).

{phang}
{opt locscale(auto|on|off)} controls the Remark 4.2 short-circuit.
When "off" or "auto" with all per-cell alpha-star quantiles strictly
positive, the location-scale regression is skipped in favor of
m_hat(z) = q_hat(z; alpha_star).


{marker remarks}{...}
{title:Remarks}

{pstd}
{ul:Standard errors.}  For the CCN methods, {cmd:brunching} returns the
OLS standard errors from the corrected regression.  These ignore the
sampling uncertainty in the censored-expectation estimation; for valid
inference, use Stata's {cmd:bootstrap} prefix:

{phang2}
{cmd:. bootstrap beta=_b[x] delta=e(delta), reps(500): brunching y x, method(symmetric) z(z)}

{pstd}
For {opt method(cct)}, a built-in pairs bootstrap is available via
{opt boot(#)}; results are returned in {cmd:e(boot_se)}, {cmd:e(boot_ci_lo)},
{cmd:e(boot_ci_hi)}.


{marker examples}{...}
{title:Examples}

{pstd}Simulate a DGP with bunching:{p_end}

{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. set obs 5000}{p_end}
{phang2}{cmd:. set seed 42}{p_end}
{phang2}{cmd:. gen z = runiform() * 10}{p_end}
{phang2}{cmd:. gen x_star = -0.5 + 0.3 * z + rnormal()}{p_end}
{phang2}{cmd:. gen x = max(x_star, 0)}{p_end}
{phang2}{cmd:. gen y = 1 + 0.8 * x_star + rnormal()}{p_end}

{pstd}Naive (biased):{p_end}
{phang2}{cmd:. brunching y x, method(naive) z(z)}{p_end}

{pstd}Symmetry correction:{p_end}
{phang2}{cmd:. brunching y x, method(symmetric) z(z) nbins(10)}{p_end}

{pstd}CCT 2025 local boundary slope:{p_end}
{phang2}{cmd:. brunching y x, method(cct) z(z) alpha_star(0.5) boot(200) seed(42)}{p_end}

{pstd}Bootstrapping the CCN global slope:{p_end}
{phang2}{cmd:. bootstrap _b[x], reps(500): brunching y x, method(symmetric) z(z)}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:brunching} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(n_censored)}}number of censored observations{p_end}
{synopt:{cmd:e(n_bins)}}number of bins used (CCN) / unique Z values (CCT){p_end}
{synopt:{cmd:e(n_switched)}}number of bins where symmetry fallback was used{p_end}
{synopt:{cmd:e(cens_exp_mean)}}mean censored expectation (CCN, non-naive){p_end}
{synopt:{cmd:e(delta)}}coefficient on correction term (CCN, non-naive){p_end}
{synopt:{cmd:e(delta_se)}}standard error of delta (CCN, non-naive){p_end}
{synopt:{cmd:e(pi_hat)}}CCT 2025 pi_hat{p_end}
{synopt:{cmd:e(delta_over_pi)}}CCT 2025 coefficient on Xtilde{p_end}
{synopt:{cmd:e(boot_B)}}CCT bootstrap successful reps{p_end}
{synopt:{cmd:e(boot_se)}}CCT bootstrap SE of beta(0+){p_end}
{synopt:{cmd:e(boot_ci_lo) / e(boot_ci_hi)}}CCT 95% percentile CI{p_end}

{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:brunching}{p_end}
{synopt:{cmd:e(method)}}the selected method{p_end}
{synopt:{cmd:e(estimand)}}"beta_global" (CCN) or "beta(0+)" (CCT){p_end}

{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector from the final regression{p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix{p_end}


{marker references}{...}
{title:References}

{phang}
Caetano, C., Caetano, G., and Nielsen, E.  2024.
"Correcting for Endogeneity Biases Using Bunching."
{it:Journal of Business and Economic Statistics}.

{phang}
Caetano, C., Caetano, G., and Tecchio, O.  2025.
"Correcting Endogeneity of Treatments with Bunching Using Censoring Strategies."
Working paper.

{phang}
Chen, S., Dahl, G. B., and Khan, S.  2005.
"Nonparametric Identification and Estimation of a Censored Location-Scale Regression Model."
{it:Journal of the American Statistical Association} 100(469): 212-221.


{title:Authors}

{phang}
Hugo Sant'Anna, University of Alabama at Birmingham.{p_end}
{phang}
Debora Mazetto, University of Alabama at Birmingham.{p_end}


{title:Also see}

{psee}
R package: {browse "https://cran.r-project.org/package=bRunching":bRunching}
