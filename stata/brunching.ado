*! brunching v0.1.0 2026-04-17
*! Stata port of the R package bRunching (CCN 2024 + CCT 2025 bunching-as-censoring estimators)
*! Authors: Hugo Sant'Anna, Debora Mazetto
*!
*! Entry point. Parses the user syntax, validates options, and dispatches
*! to the method-specific worker (_brunching_ccn for CCN 2024 methods, and
*! _brunching_cct for CCT 2025).
*!
*! Usage:
*!   brunching depvar indepvar [if] [in], method(string) z(varname)           ///
*!       [controls(varlist) nbins(#) boot(#) swap(tobit|het_tobit)            ///
*!        alpha_star(#) alpha_grid(numlist) alpha_s(#)                       ///
*!        zeta_0(#) zeta_1(#) locscale(auto|on|off) seed(#)]
*!
*! depvar  = outcome Y
*! indepvar = treatment X censored at 0
*! z()     = conditioning / instrument variable (required for non-naive
*!           CCN methods and for cct)

program define brunching, eclass
    version 14.0

    syntax varlist(min=2 max=2 numeric) [if] [in],  ///
        METHOD(string)                              ///
        [ Z(varname numeric)                        ///
          CONTROLS(varlist numeric)                 ///
          NBINS(integer 10)                         ///
          SWAP(string)                              ///
          BOOT(integer 0)                           ///
          alpha_star(real 0.5)                      ///
          alpha_grid(numlist)                       ///
          alpha_s(real -1)                          ///
          locscale(string)                          ///
          zeta_0(real -1)                           ///
          zeta_1(real -1)                           ///
          SEED(integer -1) ]

    // ---------- split depvar / indepvar ----------
    tokenize `varlist'
    local outcome   "`1'"
    local treatment "`2'"

    // ---------- method validation ----------
    local method = lower("`method'")
    local ok_methods "naive uniform tobit het_tobit symmetric cct"
    if !`:list method in ok_methods' {
        di as err "brunching: method() must be one of: `ok_methods'"
        exit 198
    }

    // swap default
    if "`swap'" == "" local swap "tobit"
    local swap = lower("`swap'")
    if !inlist("`swap'", "tobit", "het_tobit") {
        di as err "brunching: swap() must be tobit or het_tobit"
        exit 198
    }

    // locscale default
    if "`locscale'" == "" local locscale "auto"
    local locscale = lower("`locscale'")
    if !inlist("`locscale'", "auto", "on", "off") {
        di as err "brunching: locscale() must be auto, on, or off"
        exit 198
    }

    // ---------- variable existence ----------
    confirm numeric variable `outcome'
    confirm numeric variable `treatment'
    if "`z'" != "" confirm numeric variable `z'

    // z() required for all non-naive methods
    if "`method'" != "naive" & "`z'" == "" {
        di as err "brunching: method(`method') requires option z()"
        exit 198
    }

    // ---------- handle if / in ----------
    marksample touse
    markout `touse' `outcome' `treatment' `z' `controls'

    // ---------- dispatch ----------
    if "`method'" == "cct" {
        _brunching_cct `outcome' `treatment' if `touse',       ///
            z(`z') controls(`controls')                        ///
            alpha_star(`alpha_star')                           ///
            alpha_grid(`alpha_grid')                           ///
            alpha_s(`alpha_s')                                 ///
            zeta_0(`zeta_0')                                   ///
            zeta_1(`zeta_1')                                   ///
            locscale(`locscale')                               ///
            boot(`boot')                                       ///
            seed(`seed')
    }
    else {
        _brunching_ccn `outcome' `treatment' if `touse',       ///
            method(`method')                                   ///
            z(`z') controls(`controls')                        ///
            nbins(`nbins') swap(`swap')
    }

    // The worker posts results via ereturn. We add a few top-level scalars.
    ereturn local cmd      "brunching"
    ereturn local cmdline  `"brunching `0'"'
    ereturn local method   "`method'"
    ereturn local outcome  "`outcome'"
    ereturn local treatment "`treatment'"
    ereturn local z        "`z'"
    ereturn local controls "`controls'"

    // ---------- display ----------
    _brunching_display
end

// ---------------------------------------------------------------------------
// Display helper: pretty-prints the ereturn contents.
// ---------------------------------------------------------------------------
program define _brunching_display
    di _n as txt "{hline 78}"
    di as txt "bRunching: bunching-as-censoring estimator"
    di as txt "{hline 78}"
    di as txt "Method:     " as res "`e(method)'"
    di as txt "Estimand:   " as res "`e(estimand)'"
    di as txt "Outcome:    " as res "`e(outcome)'"
    di as txt "Treatment:  " as res "`e(treatment)'"
    if "`e(z)'" != "" {
        di as txt "Z variable: " as res "`e(z)'"
    }
    di as txt "N obs:      " as res e(N)
    di as txt "N censored: " as res e(n_censored)
    if !missing(e(n_bins)) {
        di as txt "N bins:     " as res e(n_bins)
    }
    if !missing(e(n_switched)) {
        di as txt "N switched: " as res e(n_switched)
    }
    di _n as txt "Coefficients:"
    ereturn display
end
