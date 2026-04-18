# brunching (Stata)

Stata port of the R package `bRunching`.  Implements the CCN 2024 (Caetano-
Caetano-Nielsen, JBES) censored-expectation estimators and the CCT 2025
(Caetano-Caetano-Tecchio) local boundary-slope estimator for treatment
variables censored at zero.

## Installation

Copy the contents of this folder into a directory on your ado path, e.g.:

```stata
sysdir               // find PERSONAL or PLUS
cp brunching.ado _brunching_ccn.ado _brunching_cct.ado brunching.sthlp \
   ~/Library/Application\ Support/Stata/ado/plus/b/
```

Or add this folder to the ado path for the current session:

```stata
adopath ++ "/path/to/bRunching/stata"
```

## Usage

```
brunching depvar indepvar [if] [in], method(name) z(varname) [options]
```

Methods:
- `naive`     - uncorrected OLS with bin fixed effects
- `uniform`   - uniform tail-symmetry correction
- `tobit`     - pooled Tobit with bin dummies
- `het_tobit` - per-bin Tobit
- `symmetric` - tail-symmetry with Tobit fallback (CCN 2024 headline method)
- `cct`       - Caetano-Caetano-Tecchio 2025 local boundary slope beta(0+)

See `help brunching` for the full syntax and option list.

## Tests

Run from this folder:

```bash
/Applications/StataNow/StataMP.app/Contents/MacOS/stata-mp -b do tests/test_basic.do
/Applications/StataNow/StataMP.app/Contents/MacOS/stata-mp -b do tests/test_r_vs_stata.do
```

`test_basic.do` exercises each method on a simulated DGP.
`test_r_vs_stata.do` loads `SampleDataWithClusters.dta` (the JBES
replication package's sample) and compares the Stata estimates to the
R package on the same data; all five CCN methods match to 1e-9.

## Parameterization fix

The canonical JBES `ccn_dist.ado` uses `reg Y X i.Z_K ..., noconstant`
which omits the base level and forces the first-bin intercept to zero.
This Stata port uses the equivalent form with K-1 bin dummies plus an
intercept, which matches the R package's saturated no-constant form
numerically (identical beta and delta to 10+ decimals).

## Files

| File | Purpose |
|------|---------|
| `brunching.ado` | Main entry point; parses syntax and dispatches |
| `brunching.sthlp` | Help file |
| `brunching.pkg` | SSC-style package manifest |
| `_brunching_ccn.ado` | CCN 2024 worker (5 methods) |
| `_brunching_cct.ado` | CCT 2025 worker (Mata-backed CDK) |
| `tests/test_basic.do` | Smoke tests on simulated DGP |
| `tests/test_r_vs_stata.do` | Cross-check against R package |

## Dependencies

None.  Pure base Stata + Mata (Stata 14+).
