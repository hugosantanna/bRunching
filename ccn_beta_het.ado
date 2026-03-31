capture program drop ccn_beta_het
program define ccn_beta_het, rclass
	local model_type = "`1'"
	local k = "`2'"
	local inv_var = "`3'"
	local dep_var = "`4'"
	local se_type = "`5'"
	local controls = "`6'"	
	local swap_type = "`7'"
	local clusterdelta = "`8'"
	local beta_controls = "`9'"
	local delta_controls = "`10'"
	
	preserve
	
	capture qui gen I = `inv_var'
	capture qui gen S = `dep_var'
	qui gen cens_ind = (I==0)
	
	sum S I `beta_controls' `delta_controls'
	foreach var in `beta_controls' {
		gen I_`var' = I*`var'
		sum `var' I_`var'
	}
	
	if ("`model_type'" == "naiveN") {
		qui sum I
		
		reg S I I_*
		return scalar ESAMP = e(N)
		return scalar DELTA = .	
		return scalar DELTA_SE = .
		return scalar BETA0 = _b[I]
		return scalar SE0 = _se[I]
		
		loc i = 1
		foreach var in `beta_controls' {
			return scalar BETA`i' = _b[I_`var']
			return scalar SE`i' = _se[I_`var']
			loc i = `i' + 1
		}
		di "model: `model_type'; inputs used: K=`k'; inv_var=`inv_var'; dep_var=`dep_var'"
	}
	
	if ("`model_type'" == "naive") {
		qui sum I
		local samp = r(N)
		
		reg S I I_* dumX_`k'_* `controls'
		return scalar ESAMP = e(N)
		return scalar DELTA = .	
		return scalar DELTA_SE = .
		return scalar BETA0 = _b[I]
		return scalar SE0 = _se[I]
		
		loc i = 1
		foreach var in `beta_controls' {
			return scalar BETA`i' = _b[I_`var']
			return scalar SE`i' = _se[I_`var']
			loc i = `i' + 1
		}
		di "model: `model_type'; inputs used: K=`k'; inv_var=`inv_var'; dep_var=`dep_var'"
	}
	
	if ("`model_type'"=="het_uniform") {
		qui gen imr_term = .
		levelsof X_`k', local(levels)
		foreach i of local levels {
			qui sum I if I>0 & X_`k' == `i'
			local pos_mean = r(mean)
			qui sum cens_ind if X_`k' == `i'
			local bunching = r(mean)
			qui replace imr_term = -`pos_mean'*(`bunching'/(1-`bunching')) if X_`k' == `i'
		}

		qui gen reg_term = imr_term*cens_ind + I
		sum imr_term
		return scalar CENS_EXPEC = `r(mean)'
		
		if "`delta_controls'" != "" {
			foreach var in `delta_controls' {
				gen reg_term_`var' = reg_term*`var'
			}
			
			reg S I I_* dumX_`k'_* `controls' reg_term reg_term_*, noconstant
			return scalar ESAMP = e(N)
			return scalar DELTA0 = _b[reg_term]	
			return scalar DELTA_SE0 = _se[reg_term]
			return scalar BETA0 = _b[I]
			return scalar SE0 = _se[I]
			
			loc i = 1
			foreach var in `beta_controls' {
				return scalar BETA`i' = _b[I_`var']
				return scalar SE`i' = _se[I_`var']
				loc i = `i' + 1
			}
			
			loc i = 1
			foreach var in `delta_controls' {
				return scalar DELTA`i' = _b[reg_term_`var']
				return scalar DELTA_SE`i' = _se[reg_term_`var']
				loc i = `i' + 1
			}
		}
		if "`delta_controls'" == "" {
			reg S I I_* dumX_`k'_* `controls' reg_term, noconstant
			return scalar ESAMP = e(N)
			return scalar DELTA0 = _b[reg_term]	
			return scalar DELTA_SE0 = _se[reg_term]
			return scalar BETA0 = _b[I]
			return scalar SE0 = _se[I]
			
			loc i = 1
			foreach var in `beta_controls' {
				return scalar BETA`i' = _b[I_`var']
				return scalar SE`i' = _se[I_`var']
				loc i = `i' + 1
			}
		}
	}
	
	if ("`model_type'" == "het_tobit") {
		*Run Tobit
		qui gen imr_term = .
		levelsof X_`k', local(levels)
		foreach i of local levels {
			qui capture tobit I if X_`k' == `i', ll(0)
			if _rc==0 {
				qui matrix coeff_mat_`i' = e(b)
				qui matrix var_mat_`i' = e(V)
				qui predict imr_temp, e(.,0)
				qui replace imr_term = imr_temp if X_`k' == `i'
				qui drop imr_temp
				local sigma_`i' = coeff_mat_`i'[1,2]^0.5
			}
			else {
				drop if X_`k'==`i'
			}
		}
		*get an imr_term (cens_exp) for each keep
		levelsof X_`k', local(levels)
		foreach i of local levels {
			qui sum imr_term if X_`k' == `i'
			local imr_term_`i' = r(mean)
		}
		qui gen reg_term = imr_term*cens_ind + I
		sum imr_term
		return scalar CENS_EXPEC = `r(mean)'
		
		if "`delta_controls'" != "" {
			foreach var in `delta_controls' {
				gen reg_term_`var' = reg_term*`var'
			}
			
			********************************************************************************
			*Corrected Models
			reg S I I_* dumX_`k'_* `controls' reg_term reg_term_*, noconstant
			return scalar ESAMP = e(N)
			return scalar DELTA0 = _b[reg_term]	
			return scalar DELTA_SE0 = _se[reg_term]
			return scalar BETA0 = _b[I]
			return scalar SE0 = _se[I]
			
			loc i = 1
			foreach var in `beta_controls' {
				return scalar BETA`i' = _b[I_`var']
				return scalar SE`i' = _se[I_`var']
				loc i = `i' + 1
			}
			
			loc i = 1
			foreach var in `delta_controls' {
				return scalar DELTA`i' = _b[reg_term_`var']
				return scalar DELTA_SE`i' = _se[reg_term_`var']
				loc i = `i' + 1
			}
		}
		if "`delta_controls'" == "" {
			********************************************************************************
			*Corrected Models
			reg S I I_* dumX_`k'_* `controls' reg_term, noconstant
			return scalar ESAMP = e(N)
			return scalar DELTA0 = _b[reg_term]	
			return scalar DELTA_SE0 = _se[reg_term]
			return scalar BETA0 = _b[I]
			return scalar SE0 = _se[I]
			
			loc i = 1
			foreach var in `beta_controls' {
				return scalar BETA`i' = _b[I_`var']
				return scalar SE`i' = _se[I_`var']
				loc i = `i' + 1
			}
		}
	}
	
	if ("`model_type'" == "symmetric") {
		*generate a variable that tells whether cluster is censored more than 50%
		levelsof X_`k', local(levels)
		foreach i of local levels {
			local switch_`i' = 0
			qui sum cens_ind if X_`k' == `i', d
			if (r(p50) == 1) {
				local switch_`i' = 1
			}
		}
		*Get the Censored Expectations
		*Step 1 : get the censorted expectations assuming symmetry.
		*These will not be correct within a cluster if the censoring is more than 50%
		gen cens_exp = .
		foreach i of local levels {
			*Step 1 -- find F_{X|Z=z)(0) -- this is just the proportion of observations with Z=z who have X_i = 0
			qui sum cens_ind if X_`k' == `i'
			*qui gen hatF_0_`i' = r(mean)
			qui local hatF_0_`i' = r(mean)
			qui local op_hatF_0_`i' = 1- `hatF_0_`i''

			*Step 2: find the 1-hatF_0_i percentile amon X (I) such that Z = i
			qui cumul I if X_`k' == `i', gen(Icumul_`i') equal
			qui sum Icumul_`i' if X_`k' == `i'
			local Icumul_min = r(min)
			if `op_hatF_0_`i'' > `Icumul_min' { 
				qui sum I if X_`k' == `i' & Icumul_`i' <= `op_hatF_0_`i''
				local step2_`i' = r(max)
			}
			if `op_hatF_0_`i'' <= `Icumul_min' {
				local step2_`i' = 0
			}
			*Step 3: find the average value of I conditional on X >= step2 and Z = i
			qui sum I if I >= `step2_`i'' & X_`k' == `i'
			local step3_`i' = r(mean)

			*Step 4: Calculate the key censored expectation, conditional on Z = i
			local cens_expec_`i' = `step2_`i''- `step3_`i''
			qui replace cens_exp = `cens_expec_`i'' if X_`k' == `i'
		}
		*also want switch codes for cases where the below qreg does not work.
		foreach i of local levels {
			qui capture qreg I if X_`k' == `i', q(`op_hatF_0_`i'')
			if (_rc != 0) {
				local switch_`i' = 1
				local qreg_switch_`i' = 1
			}
		}
		*get total number of switches
		local tot_switch = 0
		foreach i of local levels {
			local tot_switch = `tot_switch' + `switch_`i''
		}
		return scalar CLUS_SWITCH = `tot_switch'
		if (`tot_switch'>0) {
			if ("`swap_type'" == "tobit") {
				*Now use standard Tobit
				qui tobit I dumX_`k'_*, ll(0) noconstant
				matrix coeff_mat_tob = e(b)
				matrix var_mat_tob = e(V)
				qui predict imr_term_tob, e(.,0)
				local sigma = coeff_mat_tob[1,`k'+2]^0.5

				*get an imr_term (cens_exp) for each keep
				foreach i of local levels {
					qui sum imr_term if X_`k' == `i'
					local imr_term_`i' = r(mean)
					*swap in if the switch term is on
					if (`switch_`i'' == 1) {
						local cens_exp_`i' = `imr_term_`i''
						replace cens_exp = `imr_term_`i'' if X_`k' == `i'
					}
				}
			}
			if ("`swap_type'" == "het_tobit") {
				*Run Tobit
				qui gen imr_term = .
				foreach i of local levels {
					if (`switch_`i'' == 1) {
						qui capture tobit I if X_`k' == `i', ll(0)
						if _rc==0 {
							qui matrix coeff_mat_`i' = e(b)
							qui matrix var_mat_`i' = e(V)
							qui predict imr_temp, e(.,0)
							qui replace imr_term = imr_temp if X_`k' == `i'
							qui drop imr_temp
							local sigma_`i' = coeff_mat_`i'[1,2]^0.5
						}
						else {
							drop if X_`k'==`i'
						}
						*get an imr_term (cens_exp) for each keep
						qui sum imr_term if X_`k' == `i'
						local imr_term_`i' = r(mean)
					}
				}
				levelsof X_`k', local(levels) // do it again because some clusters have few observations and het_tobit does not work.
			}
		}
		qui gen reg_term = cens_exp*cens_ind + I
		sum cens_exp
		return scalar CENS_EXPEC = `r(mean)'
		
		if "`delta_controls'" != "" {
			foreach var in `delta_controls' {
				gen reg_term_`var' = reg_term*`var'
			}
			
			********************************************************************************
			*Corrected Models
			reg S I I_* dumX_`k'_* `controls' reg_term reg_term_*, noconstant
			return scalar ESAMP = e(N)
			return scalar DELTA0 = _b[reg_term]	
			return scalar DELTA_SE0 = _se[reg_term]
			return scalar BETA0 = _b[I]
			return scalar SE0 = _se[I]
			
			loc i = 1
			foreach var in `beta_controls' {
				return scalar BETA`i' = _b[I_`var']
				return scalar SE`i' = _se[I_`var']
				loc i = `i' + 1
			}
			
			loc i = 1
			foreach var in `delta_controls' {
				return scalar DELTA`i' = _b[reg_term_`var']
				return scalar DELTA_SE`i' = _se[reg_term_`var']
				loc i = `i' + 1
			}
		}
		if "`delta_controls'" == "" {
			********************************************************************************
			*Corrected Models
			reg S I I_* dumX_`k'_* `controls' reg_term, noconstant
			return scalar ESAMP = e(N)
			return scalar DELTA0 = _b[reg_term]	
			return scalar DELTA_SE0 = _se[reg_term]
			return scalar BETA0 = _b[I]
			return scalar SE0 = _se[I]
			
			loc i = 1
			foreach var in `beta_controls' {
				return scalar BETA`i' = _b[I_`var']
				return scalar SE`i' = _se[I_`var']
				loc i = `i' + 1
			}
		}
	}
	
	restore
end
