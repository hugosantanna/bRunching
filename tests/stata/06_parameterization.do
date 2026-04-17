*! Check whether Stata <-> R disagreement comes from i.Z_k (K-1 dummies)
*! vs ibn.Z_k (all K dummies) under `noconstant`.

clear all
set more off

import delimited "test_data.csv", clear varnames(1)

local k = 10
xtile X_`k' = z, nquantiles(`k')
gen Z_`k' = X_`k'
tab X_`k', gen(dumX_`k'_)

gen cens_ind = (x == 0)
rename x I
rename y S

* build reg_term the same way ccn_hom_dm does for het_uniform
gen imr_term = .
levelsof X_`k', local(levels)
foreach i of local levels {
    qui sum I if I>0 & X_`k' == `i'
    local pos_mean = r(mean)
    qui sum cens_ind if X_`k' == `i'
    local bunching = r(mean)
    qui replace imr_term = -`pos_mean'*(`bunching'/(1-`bunching')) if X_`k' == `i'
}
gen reg_term = imr_term*cens_ind + I

di as txt _n "=== Stata current: i.Z_k (K-1 dummies, base omitted) + noconstant ==="
qui reg S I i.Z_`k' reg_term, noconstant
di "beta     = " _b[I]
di "reg_term = " _b[reg_term]

di as txt _n "=== Stata with ibn.Z_k (all K dummies) + noconstant ==="
qui reg S I ibn.Z_`k' reg_term, noconstant
di "beta     = " _b[I]
di "reg_term = " _b[reg_term]

di as txt _n "=== Stata with dumX_k_* (all K dummies) + noconstant ==="
qui reg S I dumX_`k'_* reg_term, noconstant
di "beta     = " _b[I]
di "reg_term = " _b[reg_term]
