capture program drop ccn_hom_dm
program define ccn_hom_dm, rclass
	local model_type = "`1'"
	local k = "`2'"
	local inv_var = "`3'"
	local dep_var = "`4'"
	local se_type = "`5'"
	local controls = "`6'"
	local swap_type = "`7'"
	local clusterdelta = `8'
	preserve
	capture qui rename `inv_var' I
	capture qui rename `dep_var' S
	qui gen cens_ind = (I==0)
	qui sum I
	local samp = r(N)

	if ("`model_type'" == "naiveN") {
		qui sum I
		local samp = r(N)
		qui reg S I
		return scalar DELTA = .	
		return scalar DELTA_SE = .
		return scalar BETA = _b[I]
		return scalar SE = _se[I]
		local beta = _b[I]
		di "model: `model_type'; inputs used: K=`k'; inv_var=`inv_var'; dep_var=`dep_var'; se_type=`se_type'; sample=`e(N)'; beta=`beta'"
	}
	
	if ("`model_type'" == "naive") {
		qui sum I
		local samp = r(N)
		qui reg S I dumX_`k'_* `controls'
		return scalar ESAMP = e(N)
		return scalar SAMP = `samp'		
		return scalar DELTA = .	
		return scalar DELTA_SE = .
		return scalar BETA = _b[I]
		return scalar SE = _se[I]
		local beta = _b[I]
		di "model: `model_type'; inputs used: K=`k'; inv_var=`inv_var'; dep_var=`dep_var'; se_type=`se_type'; sample=`e(N)'; beta=`beta'"
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

		reg S I dumX_`k'_* `controls' reg_term, noconstant
		return scalar ESAMP = e(N)
		return scalar SAMP = `samp'		
		return scalar DELTA = _b[reg_term]	
		return scalar DELTA_SE = _se[reg_term]
		return scalar BETA = _b[I]
		return scalar SE = _se[I]
		local beta = _b[I]
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

		********************************************************************************
		*Corrected Models
		reg S I i.Z_`k' `controls' reg_term		
		return scalar ESAMP = e(N)
		return scalar SAMP = `samp'		
		return scalar DELTA = _b[reg_term]	
		return scalar DELTA_SE = _se[reg_term]
		return scalar BETA = _b[I]
		return scalar SE = _se[I]
		local beta = _b[I]
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
	}
	********************************************************************************
	*Corrected Models
	if ("`model_type'" == "het_uniform") | ("`model_type'" == "het_tobit") | ("`model_type'" == "symmetric") {
		if (`clusterdelta' != 1) {
			qui reg S I dumX_`k'_* `controls' i.Z_`clusterdelta'#c.reg_term, noconstant
			return scalar ESAMP = e(N)
			return scalar SAMP = `samp'		
			return scalar BETA = _b[I]
			
			if `k' == 1 {
				return scalar DELTA = _b[1.Z_`clusterdelta'#c.reg_term]
			}
			else {
				return scalar DELTA = .	
			}
			
			local klist = ""
			levelsof Z_`clusterdelta', local(levelsZ)
			foreach kk of local levelsZ {
				capture return scalar DELTA_`kk'=_b[`kk'.Z_`clusterdelta'#c.reg_term]
				if _rc ~ =0 {
					return scalar DELTA_`kk' = .
				}
			}
			
			if "$BOOT" == "1" { 
				local dlist = ""	
				foreach kk of local levelsZ {
					capture display _b[`kk'.Z_`clusterdelta'#c.reg_term]
					if _rc ~= 0 {
						local dlist = "`dlist'"
					}
					else {
						local dlist = "`dlist' (_b[`kk'.Z_`clusterdelta'#c.reg_term]=${d`kk'_`model_type'})"
					}
				}
				test `dlist'			
			}
			else {
				testparm i.Z_`clusterdelta'*#c.reg_term
			}
			return scalar F = `r(F)'
			return scalar df = `r(df)'
			return scalar df_r = `r(df_r)'
		}
		
		else { // No cluster delta
			reg S I i.Z_`k' `controls' reg_term, noconstant
			return scalar ESAMP = e(N)
			return scalar SAMP = `samp'	
			return scalar BETA = _b[I]
			return scalar SE = _se[I]
			return scalar DELTA = _b[reg_term]
			return scalar DELTA_SE = _se[reg_term]
			return scalar F = .
			return scalar df = .
			return scalar df_r = .
		}
		local beta = _b[I]
		di "model: `model_type'; inputs used: K=`k'; inv_var=`inv_var'; dep_var=`dep_var'; se_type=`se_type'; sample=`e(N)'; beta=`beta'"
	}
	restore
end
