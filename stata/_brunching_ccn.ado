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
          SWAP(string)                                         ///
          BOOT(integer 0)                                      ///
          SEED(integer -1) ]

    if "`swap'" == "" local swap "tobit"

    tokenize `varlist'
    local Y "`1'"
    local X "`2'"

    marksample touse
    markout `touse' `Y' `X' `z' `controls'

    preserve
    qui keep if `touse'

    // ---------- snapshot the touse sample for the bootstrap ----------
    // Saved before any temp-variable creation so resamples start from a
    // clean dataset containing only Y, X, Z, and controls. Each bootstrap
    // iteration will reload from this snapshot, draw with replacement,
    // and re-run the entire bin/cens_exp/regression pipeline.
    tempfile boot_snap
    if `boot' > 0 {
        qui save "`boot_snap'", replace
    }

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

    // ---------- pairs bootstrap (optional) ----------
    // Implemented as a Stata-level loop because the Tobit methods need
    // Stata's `tobit` command (no Mata implementation). The main fit
    // remains on the original e(), bit-identical to what it was before.
    //
    // We resample from the FULL touse sample (the snapshot saved before
    // dropping failed-bin observations) so that each bootstrap iteration
    // re-runs the entire pipeline: fresh bins, fresh censored expectation,
    // fresh corrected regression. Iterations that fail (rank-deficient
    // design, Tobit non-convergence) are silently skipped and tracked via
    // B_ok.
    tempname B_beta B_delta boot_se_b boot_se_d boot_ci_b boot_ci_d boot_res
    local B_ok = 0
    if `boot' > 0 {
        if `seed' >= 0 set seed `seed'

        // The outer `preserve` from the main fit is already in effect, and
        // Stata does not allow nested preserves in the same program. Each
        // bootstrap iteration loads the snapshot via `use, clear` (which
        // wipes the data but NOT local macros or matrices) so the main-fit
        // matrices `b'/`V' and the harvested locals (n_cens_used, ...,
        // delta_se) survive the loop and can still be passed downstream.
        //
        // CRITICAL: do NOT use `b' as the loop counter -- it shadows the
        // tempname holding the main-fit coefficient matrix.

        mata: `B_beta'  = J(`boot', 1, .)
        mata: `B_delta' = J(`boot', 1, .)

        forvalues bi = 1 / `boot' {
            qui use "`boot_snap'", clear
            qui bsample

            // Recompute bins on the resample.
            tempvar zbin_b cens_ind_b cens_exp_b dropmask_b reg_term_b
            qui gen byte `cens_ind_b' = (`X' == 0)

            qui levelsof `z', local(zlevs_b)
            local n_unique_b : word count `zlevs_b'
            if `n_unique_b' <= `nbins' {
                qui egen int `zbin_b' = group(`z')
                local n_bins_b = `n_unique_b'
            }
            else {
                capture qui xtile `zbin_b' = `z', nq(`nbins')
                if _rc != 0 {
                    drop _all
                    continue
                }
                local n_bins_b = `nbins'
            }

            qui gen double `cens_exp_b' = .
            qui gen byte `dropmask_b' = 0

            capture {
                if "`method'" == "uniform" {
                    qui _brunching_cens_uniform `X' `zbin_b' `cens_ind_b' `cens_exp_b' `dropmask_b'
                }
                else if "`method'" == "tobit" {
                    qui _brunching_cens_tobit_pooled `X' `zbin_b' `cens_ind_b' `cens_exp_b' `dropmask_b'
                }
                else if "`method'" == "het_tobit" {
                    qui _brunching_cens_het_tobit `X' `zbin_b' `cens_ind_b' `cens_exp_b' `dropmask_b'
                }
                else if "`method'" == "symmetric" {
                    qui _brunching_cens_symmetric `X' `zbin_b' `cens_ind_b' `cens_exp_b' `dropmask_b' `swap'
                }
            }
            if _rc != 0 {
                drop _all
                continue
            }

            qui drop if `dropmask_b' == 1 | missing(`cens_exp_b')
            qui count
            if r(N) < 3 {
                drop _all
                continue
            }

            qui gen double `reg_term_b' = `cens_exp_b' * `cens_ind_b' + `X'

            capture {
                if `n_bins_b' > 1 {
                    qui reg `Y' `X' i.`zbin_b' `controls' `reg_term_b'
                }
                else {
                    qui reg `Y' `X' `controls' `reg_term_b'
                }
            }
            if _rc != 0 {
                drop _all
                continue
            }

            // Harvest beta and delta if both are estimable.
            local bx = .
            local dt = .
            capture local bx = _b[`X']
            local rc1 = _rc
            capture local dt = _b[`reg_term_b']
            local rc2 = _rc
            if `rc1' != 0 | `rc2' != 0 {
                drop _all
                continue
            }
            // missing() works on numerics; locals are strings, so coerce.
            if missing(real("`bx'")) | missing(real("`dt'")) {
                drop _all
                continue
            }

            mata: `B_beta'[`bi']  = `bx'
            mata: `B_delta'[`bi'] = `dt'
            local B_ok = `B_ok' + 1

            drop _all
        }

        // Compute SE / CI in Mata.
        mata: _brunching_ccn_boot_summary(`B_beta',  "`boot_se_b'", "`boot_ci_b'")
        mata: _brunching_ccn_boot_summary(`B_delta', "`boot_se_d'", "`boot_ci_d'")
    }

    _brunching_post_result, b(`b') v(`V') n(`N')                           ///
        method(`method') estimand(beta_global)                             ///
        n_cens(`n_cens_used') n_bins(`n_bins_used')                        ///
        n_switched(`n_switched')                                           ///
        cens_exp_mean(`cens_exp_mean')                                     ///
        delta(`delta') delta_se(`delta_se')

    if `boot' > 0 {
        ereturn scalar boot_B  = `B_ok'
        ereturn scalar boot_B_req = `boot'
        ereturn scalar boot_se_beta  = `boot_se_b'[1, 1]
        ereturn scalar boot_ci_beta_lo = `boot_ci_b'[1, 1]
        ereturn scalar boot_ci_beta_hi = `boot_ci_b'[1, 2]
        ereturn scalar boot_se_delta = `boot_se_d'[1, 1]
        ereturn scalar boot_ci_delta_lo = `boot_ci_d'[1, 1]
        ereturn scalar boot_ci_delta_hi = `boot_ci_d'[1, 2]
        // Free Mata containers
        capture mata: mata drop `B_beta' `B_delta'
    }

    restore
end


// ---------------------------------------------------------------------------
// Mata helper for bootstrap summaries: SE = sd, 95% CI via percentile.
// Writes a 1x1 SE scalar matrix and a 1x2 CI matrix into Stata.
// ---------------------------------------------------------------------------
mata:
mata set matastrict off

void _brunching_ccn_boot_summary(real colvector b, string scalar se_name,
                                  string scalar ci_name)
{
    real colvector ok, sorted
    real scalar n, lo, hi, se_val, h_lo, h_hi, i_lo, i_hi, f_lo, f_hi
    real rowvector ci

    ok = select(b, b :!= .)
    n  = length(ok)

    if (n >= 2) se_val = sqrt(variance(ok))
    else        se_val = .

    if (n >= 10) {
        sorted = sort(ok, 1)
        h_lo = (n - 1) * 0.025 + 1
        h_hi = (n - 1) * 0.975 + 1
        i_lo = floor(h_lo); f_lo = h_lo - i_lo
        i_hi = floor(h_hi); f_hi = h_hi - i_hi
        if (i_lo < 1) i_lo = 1
        if (i_lo >= n) lo = sorted[n]
        else           lo = sorted[i_lo] + f_lo * (sorted[i_lo + 1] - sorted[i_lo])
        if (i_hi < 1) i_hi = 1
        if (i_hi >= n) hi = sorted[n]
        else           hi = sorted[i_hi] + f_hi * (sorted[i_hi + 1] - sorted[i_hi])
    }
    else {
        lo = .; hi = .
    }

    st_matrix(se_name, (se_val))
    ci = (lo, hi)
    st_matrix(ci_name, ci)
}

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
