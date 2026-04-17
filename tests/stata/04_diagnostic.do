*! Diagnostic: dump bin assignments, cens_exp, reg_term from Stata's het_uniform
*! so we can compare to R's cens_exp_uniform on the identical rows.

clear all
set more off

run "../../ccn_hom_dm.ado"

import delimited "test_data.csv", clear varnames(1)

local k = 10
xtile X_`k' = z, nquantiles(`k')
gen Z_`k'   = X_`k'
tab X_`k', gen(dumX_`k'_)

gen cens_ind = (x == 0)
rename x I
rename y S

* replicate the het_uniform bin loop EXACTLY as in ccn_hom_dm.ado
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

* per-bin summary: bin id, pos_mean, bunching, cens_exp
preserve
    collapse (mean) pos_mean = I if I > 0, by(X_`k')
    tempfile posmean
    save `posmean'
restore
preserve
    collapse (mean) bunching = cens_ind (first) cens_exp = imr_term, by(X_`k')
    merge 1:1 X_`k' using `posmean', nogen
    gen recomputed = -pos_mean * (bunching / (1 - bunching))
    list, noobs clean
    export delimited "stata_bin_stats.csv", replace
restore

* full row-level dump for reg_term
keep X_`k' I S cens_ind imr_term reg_term
rename X_`k' bin
rename imr_term cens_exp
export delimited "stata_rows.csv", replace
