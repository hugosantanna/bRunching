*! test_r_vs_stata.do - cross-check brunching against the R package
*!
*! Loads the JBES sample data (SampleDataWithClusters.dta) and runs each
*! method through brunching.  Compares to the R reference values baked
*! in below (see /tmp/ccn_repl/r_results.csv line "R_with_controls").

clear all
set more off
version 14.0

adopath ++ "`c(pwd)'/.."

// R reference values (with controls; Z_5 as categorical; nbins = 5)
local r_naive      =  0.0637229954609542
local r_uniform    = -0.0970782468969029
local r_tobit      = -0.0607978427742782
local r_het_tobit  = -0.0839295253260858
local r_symmetric  = -0.115455922165742

local r_delta_uniform    = 0.131998762484478
local r_delta_tobit      = 0.0962169117871671
local r_delta_het_tobit  = 0.113833824244464
local r_delta_symmetric  = 0.14163686197398

local tol = 1e-4  // absolute tolerance, 4+ decimals

local dta "/tmp/ccn_repl/SampleDataWithClusters.dta"
capture confirm file "`dta'"
if _rc {
    di as err "Reference data `dta' not found."
    di as err "Generate via: run_comparison.do in /tmp/ccn_repl"
    exit 601
}

// --------------------------------------------------------------------------
// helper macros to compare scalar results
// --------------------------------------------------------------------------

local fail_count = 0

capture program drop _check
program define _check
    args name obs exp tol
    local diff = abs(`obs' - `exp')
    if `diff' > `tol' {
        di as err "  FAIL: `name'  got=`obs'  want=`exp'  diff=`diff'"
        exit 9
    }
    di as res "  PASS: `name'  got=" %12.8f `obs' "  want=" %12.8f `exp' "  diff=" %9.2e `diff'
end

// --------------------------------------------------------------------------
// Test naive
// --------------------------------------------------------------------------
di as txt _n _dup(78) "=" _n "Cross-check: method = naive" _n _dup(78) "="
use "`dta'", clear
local controls ChildMale ChildWhite ChildBlack ChildHispanic
brunching Y X, method(naive) z(Z_5) controls(`controls') nbins(5)
_check naive `=_b[X]' `r_naive' `tol'

// --------------------------------------------------------------------------
// Test uniform
// --------------------------------------------------------------------------
di as txt _n _dup(78) "=" _n "Cross-check: method = uniform" _n _dup(78) "="
use "`dta'", clear
brunching Y X, method(uniform) z(Z_5) controls(`controls') nbins(5)
_check uniform_beta  `=_b[X]'       `r_uniform'       `tol'
_check uniform_delta `=e(delta)'    `r_delta_uniform' `tol'

// --------------------------------------------------------------------------
// Test tobit (pooled)
// --------------------------------------------------------------------------
di as txt _n _dup(78) "=" _n "Cross-check: method = tobit" _n _dup(78) "="
use "`dta'", clear
brunching Y X, method(tobit) z(Z_5) controls(`controls') nbins(5)
_check tobit_beta  `=_b[X]'    `r_tobit'       `tol'
_check tobit_delta `=e(delta)' `r_delta_tobit' `tol'

// --------------------------------------------------------------------------
// Test het_tobit (per-bin)
// --------------------------------------------------------------------------
di as txt _n _dup(78) "=" _n "Cross-check: method = het_tobit" _n _dup(78) "="
use "`dta'", clear
brunching Y X, method(het_tobit) z(Z_5) controls(`controls') nbins(5)
_check het_tobit_beta  `=_b[X]'    `r_het_tobit'       `tol'
_check het_tobit_delta `=e(delta)' `r_delta_het_tobit' `tol'

// --------------------------------------------------------------------------
// Test symmetric (headline method)
// --------------------------------------------------------------------------
di as txt _n _dup(78) "=" _n "Cross-check: method = symmetric" _n _dup(78) "="
use "`dta'", clear
brunching Y X, method(symmetric) z(Z_5) controls(`controls') nbins(5) swap(tobit)
_check symmetric_beta  `=_b[X]'    `r_symmetric'       `tol'
_check symmetric_delta `=e(delta)' `r_delta_symmetric' `tol'

// --------------------------------------------------------------------------
// Summary
// --------------------------------------------------------------------------
di _n as txt _dup(78) "=" _n "All cross-checks PASSED at tol = `tol'" _n _dup(78) "="
exit 0
