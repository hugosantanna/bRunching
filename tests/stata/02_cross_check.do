*! Cross-check step 2: run ccn_hom_dm on the same data and dump results.
*! Expects test_data.csv in the same folder; writes stata_results.csv.

clear all
set more off

* --- locate files relative to this do-file -----------------------------------
local here "`c(pwd)'"
di "Working dir: `here'"

* load the ado from the repo root (two levels up from tests/stata)
run "../../ccn_hom_dm.ado"

* --- load data ---------------------------------------------------------------
import delimited "test_data.csv", clear varnames(1)

* --- build bin variables expected by ccn_hom_dm ------------------------------
local k = 10

* equal-frequency bins on the instrument z (matches R's quantile-based cut)
xtile X_`k' = z, nquantiles(`k')
gen Z_`k'   = X_`k'
tab X_`k', gen(dumX_`k'_)

* --- run each method and collect results ------------------------------------
tempname memhold
tempfile results
postfile `memhold' str12 method double(beta delta) using `results'

* positional args: model k inv_var dep_var se_type controls swap clusterdelta
ccn_hom_dm naive       `k' x y robust "" tobit 1
post `memhold' ("naive")     (r(BETA)) (.)

ccn_hom_dm het_uniform `k' x y robust "" tobit 1
post `memhold' ("uniform")   (r(BETA)) (r(DELTA))

ccn_hom_dm het_tobit   `k' x y robust "" tobit 1
post `memhold' ("het_tobit") (r(BETA)) (r(DELTA))

ccn_hom_dm symmetric   `k' x y robust "" tobit 1
post `memhold' ("symmetric") (r(BETA)) (r(DELTA))

postclose `memhold'

use `results', clear
list, noobs clean
export delimited "stata_results.csv", replace

di "Stata results written to stata_results.csv"
