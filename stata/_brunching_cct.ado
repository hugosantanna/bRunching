*! _brunching_cct v0.1.0 2026-04-17
*! CCT 2025 Caetano-Caetano-Tecchio bunching-as-censoring estimator.
*! Implements the Chen-Dahl-Khan (2005) location-scale m(Z) procedure and
*! the Proposition 5.1 linear form for beta(0+).
*!
*! Port of R/cct.R.  Uses Mata for the CDK quantile / weighted-OLS step.

program define _brunching_cct, eclass
    version 14.0

    syntax varlist(min=2 max=2 numeric) [if] [in],            ///
        Z(varname numeric)                                    ///
        [ CONTROLS(varlist numeric)                           ///
          alpha_star(real 0.5)                                ///
          alpha_grid(numlist)                                 ///
          alpha_s(real -1)                                    ///
          zeta_0(real -1)                                     ///
          zeta_1(real -1)                                     ///
          locscale(string)                                    ///
          BOOT(integer 0)                                     ///
          SEED(integer -1) ]

    if "`locscale'" == "" local locscale "auto"

    tokenize `varlist'
    local Y "`1'"
    local X "`2'"

    marksample touse
    markout `touse' `Y' `X' `z' `controls'

    preserve
    qui keep if `touse'

    qui count
    local n_obs = r(N)
    qui count if `X' == 0
    local n_cens = r(N)

    // ---------- validate alpha_star ----------
    if `alpha_star' <= 0 | `alpha_star' >= 1 {
        di as err "brunching: alpha_star must be strictly in (0, 1)"
        restore
        exit 198
    }

    // default alpha_grid: seq(alpha_star, 0.99, length.out = 9)[-1]
    if "`alpha_grid'" == "" {
        local step = (0.99 - `alpha_star') / 8
        local alpha_grid ""
        forvalues k = 1/8 {
            local v = `alpha_star' + `k' * `step'
            local alpha_grid "`alpha_grid' `v'"
        }
    }

    // default alpha_s: max(alpha_grid)
    if `alpha_s' < 0 {
        local alpha_s = 0
        foreach a of numlist `alpha_grid' {
            if `a' > `alpha_s' local alpha_s = `a'
        }
    }

    qui summarize `X'
    local sdX = r(sd)
    if `zeta_0' < 0 {
        local zeta_0 = 1e-3 * `sdX'
        if `zeta_0' <= 0 local zeta_0 = 1e-6
    }
    if `zeta_1' < 0 local zeta_1 = `zeta_0'

    foreach a of numlist `alpha_grid' {
        if `a' <= `alpha_star' | `a' >= 1 {
            di as err "brunching: alpha_grid must lie strictly in (alpha_star, 1); got `a'"
            restore
            exit 198
        }
    }
    if `alpha_s' <= `alpha_star' | `alpha_s' >= 1 {
        di as err "brunching: alpha_s must lie strictly in (alpha_star, 1)"
        restore
        exit 198
    }

    qui levelsof `z', local(zlevs)
    local n_unique_z : word count `zlevs'
    if `n_unique_z' < 3 {
        di as err "brunching: method(cct) requires at least 3 unique values of Z; got `n_unique_z'"
        restore
        exit 459
    }

    // ---------- build m_hat(Z) via CDK ----------
    tempvar zid
    qui egen int `zid' = group(`z')

    tempvar m_hat

    mata: _brunching_cct_mhat(                                ///
        "`X'", "`zid'",                                       ///
        "`alpha_star'", "`alpha_grid'", "`alpha_s'",          ///
        `zeta_0', `zeta_1',                                   ///
        "`locscale'",                                         ///
        "`m_hat'"                                             ///
    )

    // ---------- Algorithm B: beta(0+) via generated regressor ----------
    tempvar diff cens_ind xtilde
    qui gen double `diff' = `m_hat' - `X' if !missing(`m_hat')
    qui gen byte  `cens_ind' = (`X' == 0)

    qui count if `cens_ind' == 1 & !missing(`m_hat')
    local n_cens_used = r(N)
    if `n_cens_used' == 0 {
        di as err "brunching: no censored observations (sum(X == 0) == 0); beta(0+) not identified"
        restore
        exit 459
    }

    qui summarize `diff' if !missing(`m_hat'), meanonly
    local sum_diff = r(sum)
    local pi_hat = `sum_diff' / `n_cens_used'

    qui gen double `xtilde' = `X' + `cens_ind' * `pi_hat'

    qui reg `Y' `X' `xtilde' `controls' if !missing(`m_hat')

    tempname b V
    matrix `b' = e(b)
    matrix `V' = e(V)
    local Nused = e(N)

    local beta_hat      = _b[`X']
    local delta_over_pi = _b[`xtilde']
    local gamma0        = _b[_cons]

    // ---------- bootstrap SE ----------
    local boot_se = .
    local boot_lo = .
    local boot_hi = .
    local B_ok = 0

    if `boot' > 0 {
        if `seed' >= 0 set seed `seed'
        tempname bootres
        mata: _brunching_cct_bootstrap(                            ///
            "`Y'", "`X'", "`zid'",                                 ///
            "`controls'",                                          ///
            "`alpha_star'", "`alpha_grid'", "`alpha_s'",           ///
            `zeta_0', `zeta_1',                                    ///
            "`locscale'",                                          ///
            `boot',                                                ///
            "`bootres'"                                            ///
        )
        // bootres is a 1x4 matrix: [se, lo, hi, B_ok]
        local boot_se = `bootres'[1, 1]
        local boot_lo = `bootres'[1, 2]
        local boot_hi = `bootres'[1, 3]
        local B_ok    = `bootres'[1, 4]
    }

    // ---------- post ----------
    ereturn post `b' `V', obs(`Nused')
    ereturn local  estimand      "beta(0+)"
    ereturn scalar n_censored    = `n_cens'
    ereturn scalar n_bins        = `n_unique_z'
    ereturn scalar n_switched    = 0
    ereturn scalar pi_hat        = `pi_hat'
    ereturn scalar delta_over_pi = `delta_over_pi'
    ereturn scalar gamma0        = `gamma0'
    ereturn scalar alpha_star    = `alpha_star'
    ereturn scalar alpha_s       = `alpha_s'
    ereturn scalar zeta_0        = `zeta_0'
    ereturn scalar zeta_1        = `zeta_1'
    ereturn local  alpha_grid    "`alpha_grid'"
    ereturn local  locscale_used "`locscale'"
    if `boot' > 0 {
        ereturn scalar boot_B     = `B_ok'
        ereturn scalar boot_se    = `boot_se'
        ereturn scalar boot_ci_lo = `boot_lo'
        ereturn scalar boot_ci_hi = `boot_hi'
    }

    restore
end


// ---------------------------------------------------------------------------
// Mata core
// ---------------------------------------------------------------------------
mata:
mata set matastrict off

// Sample quantile, type 1 (inverse of empirical CDF), matching R's
// stats::quantile(x, probs = alpha, type = 1).
real scalar _brunching_quantile_t1(real colvector x, real scalar alpha)
{
    real colvector xs
    real scalar n, h, j
    if (length(x) == 0) return(.)
    xs = sort(x, 1)
    n  = length(xs)
    h  = n * alpha
    if (h <= 0) return(xs[1])
    j = ceil(h)
    if (j > n) j = n
    return(xs[j])
}

// Buchinsky-Hahn smooth weight: C^1 approx to 1{q > 0}
real colvector _brunching_cct_weight(real colvector q, real scalar z1)
{
    real scalar i, qi, num1, scale
    real colvector out
    out = J(rows(q), 1, 0)
    for (i = 1; i <= rows(q); i++) {
        qi = q[i]
        if (qi >= 3 * z1) {
            out[i] = 1
        }
        else if (qi > z1) {
            num1  = exp(qi - 2 * z1) / (1 + exp(qi - 2 * z1))    ///
                    - exp(-z1) / (1 + exp(-z1))
            scale = (2 + exp(z1) + exp(-z1)) / (exp(z1) - exp(-z1))
            out[i] = num1 * scale
        }
    }
    return(out)
}

// ---- helper: find column index in a sorted vector matching a value ----
real scalar _brunching_find_col(real colvector v, real scalar target)
{
    real scalar k, L
    L = length(v)
    for (k = 1; k <= L; k++) {
        if (abs(v[k] - target) < 1e-12) return(k)
    }
    return(1)
}

// ---- core m_hat computation returning a colvector ----
// Returns J(0,1,.) on degenerate input.
real colvector _brunching_cct_mhat_core(real colvector X, real colvector Zid,
                                        real scalar alpha_star,
                                        real colvector alpha_grid,
                                        real scalar alpha_s,
                                        real scalar zeta_0, real scalar zeta_1,
                                        string scalar ls_arg)
{
    real colvector z_levels, alphas_all, q_star, q_s, m_by_z, mhat_obs
    real colvector q_star_i, q_s_i, w_i, denom, safe, C_ell, cols_ell
    real colvector q_ell_i, rr, q_row, keep_ell, z_idx, c_k, yy, bb
    real matrix Q, Q_ell, DD, XtX, inv_mat
    real scalar K, L, k, ell, col_star, col_s, n, w_sum, use_locscale
    real colvector xk

    z_levels = sort(uniqrows(Zid), 1)
    K = length(z_levels)
    L = length(alpha_grid)
    if (K < 2) return(J(0, 1, .))

    alphas_all = sort(uniqrows(alpha_star \ alpha_grid \ alpha_s), 1)

    Q = J(K, length(alphas_all), .)
    for (k = 1; k <= K; k++) {
        xk = select(X, Zid :== z_levels[k])
        if (length(xk) == 0) continue
        for (ell = 1; ell <= length(alphas_all); ell++) {
            Q[k, ell] = _brunching_quantile_t1(xk, alphas_all[ell])
        }
    }

    col_star = _brunching_find_col(alphas_all, alpha_star)
    col_s    = _brunching_find_col(alphas_all, alpha_s)
    cols_ell = J(L, 1, .)
    for (ell = 1; ell <= L; ell++) {
        cols_ell[ell] = _brunching_find_col(alphas_all, alpha_grid[ell])
    }

    q_star = Q[., col_star]
    q_s    = Q[., col_s]
    Q_ell  = J(K, L, .)
    for (ell = 1; ell <= L; ell++) Q_ell[., ell] = Q[., cols_ell[ell]]

    // locscale decision
    if (ls_arg == "off") {
        use_locscale = 0
    }
    else if (ls_arg == "on") {
        use_locscale = 1
    }
    else {
        // auto: use full CDK if any q_star(z) <= 0
        use_locscale = 0
        for (k = 1; k <= K; k++) {
            if (q_star[k] <= 0) {
                use_locscale = 1
                k = K  // break
            }
        }
    }

    // short-circuit path: m_hat(z) = q_star(z)
    if (use_locscale == 0) {
        m_by_z = q_star
    }
    else {
        // Step 2: C_ell weights from observation-level design
        n = length(X)
        z_idx = J(n, 1, 0)
        for (k = 1; k <= K; k++) z_idx = z_idx :+ (Zid :== z_levels[k]) :* k

        q_star_i = q_star[z_idx]
        q_s_i    = q_s[z_idx]
        w_i      = _brunching_cct_weight(q_star_i, zeta_1)
        denom    = q_s_i :- q_star_i
        safe     = (denom :> 0) :* (w_i :> 0)
        w_sum    = sum(w_i :* safe)

        if (w_sum <= 0) {
            // degenerate: fall back to q_star
            m_by_z = q_star
        }
        else {
            // compute C_ell
            C_ell = J(L, 1, 0)
            for (ell = 1; ell <= L; ell++) {
                q_ell_i = Q_ell[z_idx, ell]
                // ratio, replacing unsafe entries with 0
                rr = J(n, 1, 0)
                for (k = 1; k <= n; k++) {
                    if (safe[k] != 0) rr[k] = (q_ell_i[k] - q_star_i[k]) / denom[k]
                }
                C_ell[ell] = sum(w_i :* rr :* safe) / w_sum
            }

            // Step 3: weighted OLS per z-cell
            m_by_z = J(K, 1, .)
            for (k = 1; k <= K; k++) {
                q_row    = Q_ell[k, .]'
                keep_ell = (q_row :>= zeta_0)
                if (sum(keep_ell) < 2) {
                    if (q_star[k] > 0) m_by_z[k] = q_star[k]
                    continue
                }
                c_k = select(C_ell, keep_ell)
                yy  = select(q_row,  keep_ell)
                DD  = J(length(c_k), 2, 1)
                DD[., 2] = c_k
                XtX = DD' * DD
                inv_mat = invsym(XtX)
                if (inv_mat[1,1] == .) {
                    if (q_star[k] > 0) m_by_z[k] = q_star[k]
                    continue
                }
                bb = inv_mat * (DD' * yy)
                m_by_z[k] = bb[1]
            }
        }
    }

    // assemble per-observation m_hat
    n = length(X)
    mhat_obs = J(n, 1, .)
    for (k = 1; k <= K; k++) {
        for (ell = 1; ell <= n; ell++) {
            if (Zid[ell] == z_levels[k]) mhat_obs[ell] = m_by_z[k]
        }
    }
    return(mhat_obs)
}


// ---- Stata-entry version: writes result to new Stata variable ----
void _brunching_cct_mhat(string scalar Xname, string scalar Zidname,
                         string scalar astar_str, string scalar grid_str,
                         string scalar as_str,
                         real scalar zeta_0, real scalar zeta_1,
                         string scalar locscale_arg,
                         string scalar mhatvar)
{
    real colvector X, Zid, alpha_grid, mhat_obs
    real scalar alpha_star, alpha_s, idx

    X   = st_data(., Xname)
    Zid = st_data(., Zidname)
    alpha_star = strtoreal(astar_str)
    alpha_s    = strtoreal(as_str)
    alpha_grid = strtoreal(tokens(grid_str))'

    mhat_obs = _brunching_cct_mhat_core(X, Zid, alpha_star, alpha_grid,
                                         alpha_s, zeta_0, zeta_1, locscale_arg)
    idx = st_addvar("double", mhatvar)
    st_store(., idx, mhat_obs)
}


// ---- bootstrap: returns SE and 95% CI via a 1x4 Stata matrix [se lo hi B_ok] ----
void _brunching_cct_bootstrap(string scalar Yname, string scalar Xname,
                              string scalar Zidname,
                              string scalar Cnames,
                              string scalar astar_str, string scalar grid_str,
                              string scalar as_str,
                              real scalar zeta_0, real scalar zeta_1,
                              string scalar locscale_arg,
                              real scalar B,
                              string scalar resname)
{
    real colvector Y, X, Zid, alpha_grid, beta_b, idx, Yb, Xb, Zb
    real colvector mhat_b, keep, Yk, Xk, mk, censk, Xtilde, bb, beta_ok, sorted
    real matrix C, Cb, D, Ck, XtX, inv_mat
    real scalar alpha_star, alpha_s, n, b, has_C, nc, pi_hat, nk, p, B_ok
    real scalar se_val, lo, hi, h_lo, h_hi, i_lo, i_hi, f_lo, f_hi, nb

    Y   = st_data(., Yname)
    X   = st_data(., Xname)
    Zid = st_data(., Zidname)
    alpha_star = strtoreal(astar_str)
    alpha_s    = strtoreal(as_str)
    alpha_grid = strtoreal(tokens(grid_str))'
    has_C = (Cnames != "")
    if (has_C) C = st_data(., Cnames)

    n = length(Y)
    beta_b = J(B, 1, .)

    for (b = 1; b <= B; b++) {
        idx = ceil(runiform(n, 1) :* n)
        Yb = Y[idx]
        Xb = X[idx]
        Zb = Zid[idx]
        if (has_C) Cb = C[idx, .]

        if (length(uniqrows(Zb)) < 2) continue

        mhat_b = _brunching_cct_mhat_core(Xb, Zb, alpha_star, alpha_grid,
                                           alpha_s, zeta_0, zeta_1, locscale_arg)
        if (length(mhat_b) == 0) continue

        // Algorithm B on resample
        keep = (mhat_b :!= .)
        if (sum(keep) < 3) continue
        Yk    = select(Yb, keep)
        Xk    = select(Xb, keep)
        mk    = select(mhat_b, keep)
        censk = (Xk :== 0)
        nc    = sum(censk)
        if (nc == 0) continue
        pi_hat = sum(mk :- Xk) / nc
        Xtilde = Xk :+ censk :* pi_hat

        nk = length(Yk)
        p  = 3
        if (has_C) p = p + cols(C)
        D = J(nk, p, 1)
        D[., 2] = Xk
        D[., 3] = Xtilde
        if (has_C) {
            Ck = select(Cb, keep)
            D[., 4..p] = Ck
        }
        XtX = D' * D
        inv_mat = invsym(XtX)
        if (inv_mat[1,1] == .) continue
        bb = inv_mat * (D' * Yk)
        beta_b[b] = bb[2]
    }

    // summaries
    B_ok = sum(beta_b :!= .)
    beta_ok = select(beta_b, beta_b :!= .)

    if (B_ok >= 2) se_val = sqrt(variance(beta_ok))
    else            se_val = .

    if (B_ok >= 10) {
        sorted = sort(beta_ok, 1)
        nb = length(sorted)
        h_lo = (nb - 1) * 0.025 + 1
        h_hi = (nb - 1) * 0.975 + 1
        i_lo = floor(h_lo); f_lo = h_lo - i_lo
        i_hi = floor(h_hi); f_hi = h_hi - i_hi
        if (i_lo < 1) i_lo = 1
        if (i_lo >= nb) lo = sorted[nb]
        else            lo = sorted[i_lo] + f_lo * (sorted[i_lo + 1] - sorted[i_lo])
        if (i_hi < 1) i_hi = 1
        if (i_hi >= nb) hi = sorted[nb]
        else            hi = sorted[i_hi] + f_hi * (sorted[i_hi + 1] - sorted[i_hi])
    }
    else {
        lo = .; hi = .
    }

    // Pass results back via a Stata matrix
    real rowvector res
    res = (se_val, lo, hi, B_ok)
    st_matrix(resname, res)
}

end
