/* -----------------------------------------------------------------------------
PROJECT: "Is Video Watching Bad for Kids? The Effect of Video Watching on 
         Children's Skills"
AUTHORS: Carolina Caetano, Gregorio Caetano, Debora Mazetto, and Meghan M. Skira
DO-FILE AUTHOR: Debora Mazetto (dmazetto@tntech.edu)
DATE: 2025/09/15
PURPOSE: create discontinuity plots
NOTES: objetive is to create the graphs with the bunching point, a local linear 
       approximation for the rest of the sample, and a discontinuity test.
----------------------------------------------------------------------------- */


********************************************************************************
***                          INITIALIZING COMMANDS                           ***
********************************************************************************
local RHS = "$RHS"
local xtitle = "$xtitle"
local LHS = "$LHS"
local ytitle = "$ytitle"
local figname = "$figname"
local ylabel = "$ylabel"


********************************************************************************
***                            DISCONTINUITY PLOT                            ***
********************************************************************************
use "$modpath\DataEstimation.dta", clear

*We do the exercise up to 8h/day:
keep if `RHS' <= 8
replace `LHS' = . if `LHS' == -9

*Create residualized LHS variable if necessary:
if "$residualized" == "Yes" {
	reg `LHS' `RHS' $controls i.Z_$K if `RHS' > 0
	
	predict double Res_`LHS', residuals
	
	replace Res_`LHS' = Res_`LHS' + _b[`RHS']*`RHS'
	
	label var `LHS' "`ytitle'"
	label var Res_`LHS' "Residual `ytitle'"
	
	local varlist = "`LHS' Res_`LHS'"
}
else {
	display "No residualized variable"
	
	label var `LHS' "`ytitle'"
	local varlist = "`LHS'"
}

*Generate the cut-off for bunching point:
egen x = cut(`RHS'), at(0 0.001(0.1)8.001)

egen tag = tag(x)

label var x "`xtitle'"

tempfile temp1
save `temp1', replace

*Create the graph:
loc i = 1
foreach var of varlist `varlist' {
	use `temp1', clear
	
	*Calculate the CI for the cut-off variable:
	egen mean`var' = mean(`var'), by(x)
	egen sd`var' = sd(`var'), by(x)
	egen n`var' = count(`var'), by(x)
	gen ci`var'_ub = mean`var' + (1.96*sd`var'/sqrt(n`var'))
	gen ci`var'_lb = mean`var' - (1.96*sd`var'/sqrt(n`var'))
	
	*Calculate the CI at the bunching point:
	sum mean`var' if x == 0
	scalar actual_mean = r(mean)
	scalar df = r(N) - 1
	
	gen se`var' = sd`var'/sqrt(n`var')
	sum se`var' if x == 0
	scalar actual_se = r(mean)
	
	*Estimate the local linear for RHS > 0:
	gen zeros = 0
	
	lpoly `var' `RHS' if `RHS' > 0, se(lpoly_se) at(zeros) generate(lpoly_mean) nograph degree(1) bwidth($bandwidth)
	
	scalar lpoly_mean = lpoly_mean
	scalar lpoly_se = lpoly_se
	scalar diff_mean = lpoly_mean - actual_mean
	scalar diff_se = sqrt((lpoly_se^2) + (actual_se^2))
	scalar diff_t = diff_mean/diff_se
	scalar diff_t = abs(diff_t)
	local pvalue = ttail(df,diff_t)
	
	gen pvalue = round(`pvalue', 0.001)
	tostring pvalue, replace force
	replace pvalue = "0" + pvalue
	replace pvalue = substr(pvalue,1,4)
	local pvalue = pvalue
	drop pvalue
	if `pvalue' == 0 {
		local pvalue = "0.000"
	}
	
	*Estimate modified local linear RHS > 0:
	lpoly `var' `RHS' if `RHS' ~= 0, gen(xlp`var' sdlp`var') nograph noscatter se(selp`var') degree(1) bwidth($bandwidth)
	
	tempfile temp2
	save `temp2', replace
	
	*Estimate the local linear at RHS == 0:
	use `temp1', clear
	
	gen zeros = 0 if `RHS' ~= 0 //in 1
	
	lpoly `var' `RHS' if `RHS' ~= 0, gen(x0`var' sd0`var') nograph noscatter se(se0`var') degree(1) bwidth($bandwidth) at(zeros)
	
	keep x0`var' sd0`var' se0`var'
	keep if x0`var' == 0
	append using `temp2'
	replace xlp`var' = x0`var' if _n == 1
	replace sdlp`var' = sd0`var' if _n == 1
	replace selp`var' = se0`var' if _n == 1
	
	reg `var' `RHS' if `RHS' ~= 0
	predict sdlin`var'
	predict selin`var', stdp
	
	gen upper_ci = sdlp`var' + 1.96*selp`var'
	gen lower_ci = sdlp`var' - 1.96*selp`var'
	
	*Make graph:
	twoway (scatter mean`var' x if tag == 1 & x == 0, mcolor(black) msymbol(oh)) ///
		(rcap ci`var'_ub ci`var'_lb x if tag == 1 & x == 0, lcolor(black)) ///
		(rarea lower_ci upper_ci xlp`var' if xlp`var' >= 0 & xlp`var' <= 8, color(gs10%80)) ///
		(line sdlp`var' xlp`var' if xlp`var' >= 0 & xlp`var' <= 8, lcolor(black)), ///
		scheme(s1mono) legend(off) ytitle("`: variable label `var''", size(medium)) ///
		xtitle("`xtitle'", size(medium)) xlabel(0(1)8) `ylabel' ///
		title("P-value of Discontinuity: `pvalue'", size(medlarge))
	
	if `i' == 1 {
		graph export "$figpath\`figname'.png", replace
	}
	else {
		graph export "$figpath\`figname'_Res.png", replace
	}
	
	loc i = `i' + 1
	
	*Clean data set:
	drop zeros lpoly*
	scalar drop _all
	drop mean`var' sd`var' n`var' ci`var'_ub ci`var'_lb se`var' xlp`var' sdlp`var' selp`var' x0`var' sd0`var' se0`var' sdlin`var' selin`var' upper_ci lower_ci
}
drop x tag
