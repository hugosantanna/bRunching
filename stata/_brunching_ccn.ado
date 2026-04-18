*! _brunching_ccn v0.1.0 2026-04-17
*! CCN 2024 censored-expectation estimators for brunching.ado.
*! Methods: naive, uniform, tobit, het_tobit, symmetric.
*!
*! Parameterization fix (Debora Mazetto):
*!   The canonical JBES ccn_dist.ado uses "reg Y X i.Z_K ... reg_term, noconstant"
*!   which omits the base level and forces the first bin's intercept to 0. The
*!   fix is to keep the constant: "reg Y X i.Z_K ... reg_term" with K-1 dummies.
*!   This matches R's saturated no-constant parameterization in closed form
*!   (identical beta_hat and delta_hat to 10+ decimals).

program define _brunching_ccn, eclass
    version 14.0

    syntax varlist(min=2 max=2 numeric) [if] [in],             ///
        METHOD(string)                                         ///
        [ Z(varname numeric)                                   ///
          CONTROLS(varlist numeric)                            ///
          NBINS(integer 10)                                    ///
          SWAP(string) ]

    if "`swap'" == "" local swap "tobit"

    tokenize `varlist'
    local Y "`1'"
    local X "`2'"

    marksample touse
    markout `touse' `Y' `X' `z' `controls'

    preserve
    qui keep if `touse'

    // ---------- basic counts ----------
    qui count
    local n_obs = r(N)
    qui count if `X' == 0
    local n_cens = r(N)

    tempvar cens_ind
    qui gen byte `cens_ind' = (`X' == 0)

    // ---------- discretize Z into bins ----------
    tempvar zbin
    if "`z'" == "" {
        // naive with no z: a single bin
        qui gen int `zbin' = 1
        local n_bins_used = 1
    }
    else {
        qui levelsof `z', local(zlevs)
        local n_unique : word count `zlevs'
        if `n_unique' <= `nbins' {
            // use z as-is (categorical); assign integer bin id
            tempvar ztmp
            qui egen int `zbin' = group(`z')
            local n_bins_used = `n_unique'
        }
        else {
            // equal-frequency bins via xtile
            qui xtile `zbin' = `z', nq(`nbins')
            local n_bins_used = `nbins'
        }
    }

    // ---------- naive branch ----------
    if "`method'" == "naive" {
        if `n_bins_used' > 1 {
            qui reg `Y' `X' i.`zbin' `controls'
        }
        else {
            qui reg `Y' `X' `controls'
        }

        tempname b V
        matrix `b' = e(b)
        matrix `V' = e(V)
        local N = e(N)

        // Post with only the coefficient of interest for display clarity
        _brunching_post_result, b(`b') v(`V') n(`N')                           ///
            method(naive) estimand(beta_global)                                ///
            n_cens(`n_cens') n_bins(`n_bins_used') n_switched(0)               ///
            cens_exp_mean(.) delta(.) delta_se(.)

        restore
        exit
    }

    // ---------- compute censored expectation ----------
    tempvar cens_exp dropmask
    qui gen double `cens_exp' = .
    qui gen byte `dropmask' = 0

    local n_switched = 0

    if "`method'" == "uniform" {
        _brunching_cens_uniform `X' `zbin' `cens_ind' `cens_exp' `dropmask'
    }
    else if "`method'" == "tobit" {
        _brunching_cens_tobit_pooled `X' `zbin' `cens_ind' `cens_exp' `dropmask'
    }
    else if "`method'" == "het_tobit" {
        _brunching_cens_het_tobit `X' `zbin' `cens_ind' `cens_exp' `dropmask'
    }
    else if "`method'" == "symmetric" {
        _brunching_cens_symmetric `X' `zbin' `cens_ind' `cens_exp' `dropmask' `swap'
        local n_switched = r(n_switched)
    }
    else {
        di as err "_brunching_ccn: unknown method `method'"
        restore
        exit 198
    }

    // Drop observations where the censored expectation is undefined
    qui drop if `dropmask' == 1 | missing(`cens_exp')

    qui count
    local n_used = r(N)
    if `n_used' == 0 {
        di as err "brunching: no observations remain after dropping failed z-cells"
        restore
        exit 459
    }

    qui count if `X' == 0
    local n_cens_used = r(N)

    // correction term
    tempvar reg_term
    qui gen double `reg_term' = `cens_exp' * `cens_ind' + `X'

    qui summarize `cens_exp', meanonly
    local cens_exp_mean = r(mean)

    // ---------- corrected regression ----------
    // Use the patched parameterization: constant + i.zbin + reg_term.
    // This matches R's saturated no-constant parameterization numerically.
    if `n_bins_used' > 1 {
        qui reg `Y' `X' i.`zbin' `controls' `reg_term'
    }
    else {
        qui reg `Y' `X' `controls' `reg_term'
    }

    tempname b V
    matrix `b' = e(b)
    matrix `V' = e(V)
    local N = e(N)

    local delta    = _b[`reg_term']
    local delta_se = _se[`reg_term']

    _brunching_post_result, b(`b') v(`V') n(`N')                           ///
        method(`method') estimand(beta_global)                             ///
        n_cens(`n_cens_used') n_bins(`n_bins_used')                        ///
        n_switched(`n_switched')                                           ///
        cens_exp_mean(`cens_exp_mean')                                     ///
        delta(`delta') delta_se(`delta_se')

    restore
end


// ---------------------------------------------------------------------------
// Post results helper.  Puts a reduced e(b) = beta, e(V) = beta variance
// keyed to the treatment coefficient, and adds scalars for display.
// The full regression e(b)/e(V) are still in the temp matrices if needed.
// ---------------------------------------------------------------------------
program define _brunching_post_result, eclass
    syntax , b(name) v(name) n(integer)                                   ///
        method(string) estimand(string)                                   ///
        n_cens(integer) n_bins(integer) n_switched(integer)               ///
        cens_exp_mean(real)                                               ///
        delta(real) delta_se(real)

    ereturn post `b' `v', obs(`n')
    ereturn scalar n_censored = `n_cens'
    ereturn scalar n_bins     = `n_bins'
    ereturn scalar n_switched = `n_switched'
    if !missing(`cens_exp_mean') ereturn scalar cens_exp_mean = `cens_exp_mean'
    if !missing(`delta')        ereturn scalar delta         = `delta'
    if !missing(`delta_se')     ereturn scalar delta_se      = `delta_se'
    ereturn local  estimand   "`estimand'"
end


// ---------------------------------------------------------------------------
// Censored-expectation estimators.  Each takes (X, zbin, cens_ind, out, drop)
// and fills the `out` variable in-place.  `drop` is set to 1 for obs that
// should be dropped from the corrected regression.
// ---------------------------------------------------------------------------

// Uniform / tail-symmetry with uniform assumption.
// E[X|X<=0, Z=z] = -E[X|X>0, Z=z] * p(z) / (1 - p(z))
program define _brunching_cens_uniform
    args X zbin cens_ind out drop

    qui levelsof `zbin', local(levs)
    foreach b of local levs {
        qui count if `zbin' == `b' & `X' > 0
        if r(N) == 0 {
            // all censored in this bin -> drop
            qui replace `drop' = 1 if `zbin' == `b'
            continue
        }
        qui count if `zbin' == `b' & `cens_ind' == 1
        if r(N) == 0 {
            // no censoring in this bin; cens_exp multiplies by cens=0,
            // so value is irrelevant -> set to 0
            qui replace `out' = 0 if `zbin' == `b'
            continue
        }

        qui summarize `X' if `zbin' == `b' & `X' > 0, meanonly
        local pos_mean = r(mean)
        qui summarize `cens_ind' if `zbin' == `b', meanonly
        local bunching = r(mean)
        qui replace `out' = -`pos_mean' * (`bunching' / (1 - `bunching')) if `zbin' == `b'
    }
end


// Pooled Tobit with bin dummies.  Matches R cens_exp_tobit_pooled().
program define _brunching_cens_tobit_pooled
    args X zbin cens_ind out drop

    // edge cases: all censored or all positive
    qui count if `X' > 0
    if r(N) == 0 {
        qui replace `drop' = 1
        exit
    }
    qui count if `X' == 0
    if r(N) == 0 {
        qui replace `out' = 0
        exit
    }

    qui levelsof `zbin', local(levs)
    local n_levs : word count `levs'

    if `n_levs' > 1 {
        capture qui tobit `X' i.`zbin', ll(0)
    }
    else {
        capture qui tobit `X', ll(0)
    }
    if _rc != 0 {
        qui replace `drop' = 1
        exit
    }

    // Robust sigma retrieval.  Stata's tobit reports EITHER sigma directly
    // (when nested in a program) OR var(e.depvar) as the last column of
    // e(b) (when invoked interactively).  Detect via the column name.
    tempname tmpb
    matrix `tmpb' = e(b)
    local lastcol = colsof(`tmpb')
    local last    = `tmpb'[1, `lastcol']
    local cnames : colnames `tmpb'
    local lastname : word `lastcol' of `cnames'
    // If the eq name contains "var" the column is variance; otherwise it's sigma.
    local eqnames : coleq `tmpb'
    local lasteq : word `lastcol' of `eqnames'
    if strpos("`lasteq'", "var") > 0 {
        local sigma = sqrt(`last')
    }
    else {
        // sigma or /sigma reported directly
        local sigma = `last'
    }

    // linear predictor: predict, xb
    tempvar xb ratio ecens
    qui predict double `xb', xb
    // E[X|X<=0] under Gaussian Tobit:
    //   xb - sigma * phi(xb/sigma) / Phi(-xb/sigma)
    qui gen double `ratio' = normalden(`xb' / `sigma') / normal(-`xb' / `sigma')
    qui gen double `ecens' = `xb' - `sigma' * `ratio'
    qui replace `out' = `ecens'
end


// Per-bin Tobit.  Matches R cens_exp_tobit().
program define _brunching_cens_het_tobit
    args X zbin cens_ind out drop

    qui levelsof `zbin', local(levs)

    foreach b of local levs {
        qui count if `zbin' == `b' & `X' == 0
        local nzero = r(N)
        qui count if `zbin' == `b' & `X' > 0
        local npos  = r(N)
        qui count if `zbin' == `b'
        local ntot  = r(N)

        // all censored -> drop
        if `nzero' == `ntot' {
            qui replace `drop' = 1 if `zbin' == `b'
            continue
        }
        // no censoring -> value irrelevant
        if `npos' == `ntot' {
            qui replace `out' = 0 if `zbin' == `b'
            continue
        }

        // need variation
        qui summarize `X' if `zbin' == `b'
        if r(min) == r(max) {
            qui replace `drop' = 1 if `zbin' == `b'
            continue
        }

        capture qui tobit `X' if `zbin' == `b', ll(0)
        if _rc != 0 {
            qui replace `drop' = 1 if `zbin' == `b'
            continue
        }

        local mu    = _b[_cons]
        // Detect sigma vs variance parameterization
        tempname tmpb
        matrix `tmpb' = e(b)
        local lastcol = colsof(`tmpb')
        local last    = `tmpb'[1, `lastcol']
        local eqnames : coleq `tmpb'
        local lasteq : word `lastcol' of `eqnames'
        if strpos("`lasteq'", "var") > 0 {
            local sigma = sqrt(`last')
        }
        else {
            local sigma = `last'
        }

        // E[X|X<=0] = mu - sigma * phi(mu/sigma) / Phi(-mu/sigma)
        local ratio  = normalden(`mu' / `sigma') / normal(-`mu' / `sigma')
        local ecens  = `mu' - `sigma' * `ratio'
        qui replace `out' = `ecens' if `zbin' == `b'
    }
end


// Symmetric (tail-symmetry) with fallback.  Matches R cens_exp_symmetric().
program define _brunching_cens_symmetric, rclass
    args X zbin cens_ind out drop swap

    if "`swap'" == "" local swap "tobit"

    qui levelsof `zbin', local(levs)

    local n_switched = 0
    local switch_bins ""

    foreach b of local levs {
        qui count if `zbin' == `b'
        local ntot = r(N)
        if `ntot' == 0 continue

        qui count if `zbin' == `b' & `cens_ind' == 1
        local ncens_b = r(N)
        local hatF_0  = `ncens_b' / `ntot'
        local op_hatF = 1 - `hatF_0'

        // ---------- switch-bin detection ----------
        // R uses median of cens_ind == 1 as the switch signal, i.e.
        // cens_ind takes value 1 in at least half the observations -> switch.
        // Equivalently: hatF_0 >= 0.5.
        local high_cens = (`hatF_0' >= 0.5)

        if `high_cens' {
            local n_switched = `n_switched' + 1
            local switch_bins "`switch_bins' `b'"
            continue
        }

        if `op_hatF' >= 1 {
            // no censoring; value irrelevant
            qui replace `out' = 0 if `zbin' == `b'
            continue
        }
        if `op_hatF' <= 0 {
            // should have been caught by high_cens, but defensive
            local n_switched = `n_switched' + 1
            local switch_bins "`switch_bins' `b'"
            continue
        }

        // --- symmetric formula ---
        // Step 2: largest X in bin such that ECDF(X) <= op_hatF
        //   ECDF(X_(j)) = rank / N (type 7 in R; Stata's cumul uses the same).
        tempvar xcumul
        qui cumul `X' if `zbin' == `b', gen(`xcumul') equal
        qui summarize `xcumul' if `zbin' == `b'
        local xc_min = r(min)
        local step2 = 0
        if `op_hatF' > `xc_min' {
            qui summarize `X' if `zbin' == `b' & `xcumul' <= `op_hatF'
            local step2 = r(max)
        }
        qui drop `xcumul'

        // Step 3: mean of X conditional on X >= step2 (within bin)
        qui summarize `X' if `zbin' == `b' & `X' >= `step2'
        local step3 = r(mean)

        local cens_expec = `step2' - `step3'
        qui replace `out' = `cens_expec' if `zbin' == `b'
    }

    // Fallback for switched bins
    if `n_switched' > 0 {
        if "`swap'" == "tobit" {
            // Pooled Tobit across full sample (match Stata JBES behavior
            // and R cens_exp_symmetric's "swap = 'tobit'" branch)
            tempvar pooled dropp
            qui gen double `pooled' = .
            qui gen byte `dropp' = 0
            _brunching_cens_tobit_pooled `X' `zbin' `cens_ind' `pooled' `dropp'

            foreach b of local switch_bins {
                qui count if `zbin' == `b' & `dropp' == 1
                if r(N) > 0 {
                    qui replace `drop' = 1 if `zbin' == `b'
                }
                else {
                    qui replace `out' = `pooled' if `zbin' == `b'
                }
            }
        }
        else if "`swap'" == "het_tobit" {
            tempvar per_bin dropp
            qui gen double `per_bin' = .
            qui gen byte `dropp' = 0
            _brunching_cens_het_tobit `X' `zbin' `cens_ind' `per_bin' `dropp'

            foreach b of local switch_bins {
                qui count if `zbin' == `b' & `dropp' == 1
                if r(N) > 0 {
                    qui replace `drop' = 1 if `zbin' == `b'
                }
                else {
                    qui replace `out' = `per_bin' if `zbin' == `b'
                }
            }
        }
    }

    return scalar n_switched = `n_switched'
end
