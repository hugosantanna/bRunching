*! Monte Carlo in Stata: shift bin-1 Y intercept by alpha_1 in {-2, 0, 2}.
*! 50 reps per alpha; runs each method via ccn_hom_dm.ado.

clear all
set more off

run "../../ccn_hom_dm.ado"

local K = 10
local B = 50
local N = 5000

tempname mh
postfile `mh' double alpha_1 int rep ///
    double(naive uniform het_tobit symmetric) ///
    using "mc_stata_results.dta", replace

local alpha_idx = 0
foreach a of numlist -2 0 2 {
    local ++alpha_idx
    di as txt _n "alpha_1 = `a'"
    forvalues r = 1/`B' {
        clear
        set obs `N'
        set seed `=`r' + 1000 * `alpha_idx''
        gen z      = runiform() * 10
        xtile bin  = z, nquantiles(`K')
        gen x_star = -1 + 0.3 * z + rnormal()
        gen x      = max(x_star, 0)
        * DGP matches package unit test: true beta = 0.8, true delta = 0.5
        gen y      = 1 + `a' * (bin == 1) + 0.8 * x + 0.5 * x_star + rnormal()

        * prep variables that ccn_hom_dm expects (do NOT create cens_ind —
        * ado creates it internally and would collide)
        xtile X_`K' = z, nquantiles(`K')
        gen Z_`K'   = X_`K'
        tab X_`K', gen(dumX_`K'_)

        local naive = .
        local uni   = .
        local het   = .
        local sym   = .

        qui capture ccn_hom_dm naive       `K' x y robust "" tobit 1
        if _rc == 0 local naive = r(BETA)
        qui capture ccn_hom_dm het_uniform `K' x y robust "" tobit 1
        if _rc == 0 local uni = r(BETA)
        qui capture ccn_hom_dm het_tobit   `K' x y robust "" tobit 1
        if _rc == 0 local het = r(BETA)
        qui capture ccn_hom_dm symmetric   `K' x y robust "" tobit 1
        if _rc == 0 local sym = r(BETA)

        post `mh' (`a') (`r') (`naive') (`uni') (`het') (`sym')
        if mod(`r', 10) == 0 di as txt "  rep `r' of `B' done"
    }
}
postclose `mh'

use "mc_stata_results.dta", clear
export delimited "mc_stata_results.csv", replace

di as txt _n "=== Stata Monte Carlo summary ==="
foreach a of numlist -2 0 2 {
    di as txt _n "--- alpha_1 = `a' ---"
    tabstat naive uniform het_tobit symmetric if alpha_1 == `a', ///
        stats(mean sd) format(%9.4f)
}
