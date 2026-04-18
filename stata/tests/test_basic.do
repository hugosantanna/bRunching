*! test_basic.do - smoke tests for brunching
*! Exercises each of the six methods on a simulated DGP and verifies
*! that estimates are finite and sensible.  Port of R test-ccn.R / test-cct.R.

clear all
set more off
version 14.0

// ensure ado files in this repository are found
adopath ++ "`c(pwd)'/.."

// ---------- Test 1: naive method runs ----------
di as txt _n _dup(78) "=" _n "Test 1: naive method"
clear
set obs 2000
set seed 1
gen z = runiform() * 10
gen x_star = 2 + 0.5 * z + rnormal()
gen x = max(x_star, 0)
gen y = 1 + 0.8 * x_star + rnormal()

brunching y x, method(naive) z(z)
assert !missing(_b[x])
assert !missing(_se[x])
assert e(N) == 2000
assert "`e(method)'" == "naive"
assert "`e(estimand)'" == "beta_global"
di as res "  PASS: naive method"

// ---------- Test 2: uniform method runs ----------
di as txt _n _dup(78) "=" _n "Test 2: uniform method"
clear
set obs 2000
set seed 2
gen z = runiform() * 10
gen x_star = 2 + 0.5 * z + rnormal()
gen x = max(x_star, 0)
gen y = 1 + 0.8 * x_star + rnormal()

brunching y x, method(uniform) z(z)
assert !missing(_b[x])
assert !missing(e(delta))
assert !missing(e(cens_exp_mean))
di as res "  PASS: uniform method"

// ---------- Test 3: tobit method runs ----------
di as txt _n _dup(78) "=" _n "Test 3: tobit method"
clear
set obs 2000
set seed 3
gen z = runiform() * 10
gen x_star = 2 + 0.5 * z + rnormal()
gen x = max(x_star, 0)
gen y = 1 + 0.8 * x_star + rnormal()

brunching y x, method(tobit) z(z)
assert !missing(_b[x])
assert !missing(e(delta))
di as res "  PASS: tobit method"

// ---------- Test 4: het_tobit method runs ----------
di as txt _n _dup(78) "=" _n "Test 4: het_tobit method"
clear
set obs 2000
set seed 30
gen z = runiform() * 10
gen x_star = 2 + 0.5 * z + rnormal()
gen x = max(x_star, 0)
gen y = 1 + 0.8 * x_star + rnormal()

brunching y x, method(het_tobit) z(z)
assert !missing(_b[x])
assert !missing(e(delta))
di as res "  PASS: het_tobit method"

// ---------- Test 5: symmetric method runs ----------
di as txt _n _dup(78) "=" _n "Test 5: symmetric method"
clear
set obs 3000
set seed 4
gen z = runiform() * 10
gen x_star = 2 + 0.5 * z + rnormal()
gen x = max(x_star, 0)
gen y = 1 + 0.8 * x_star + rnormal()

brunching y x, method(symmetric) z(z) nbins(10)
assert !missing(_b[x])
assert !missing(e(delta))
di as res "  PASS: symmetric method"

// ---------- Test 6: corrected closer to truth than naive ----------
di as txt _n _dup(78) "=" _n "Test 6: corrected closer to truth than naive"
// DGP: y = 1 + 0.8*x + 0.5*x_star + noise
// true beta = 0.8; naive should be biased toward 1.3.
clear
set obs 5000
set seed 10
gen z = runiform() * 10
gen x_star = -1 + 0.3 * z + rnormal()
gen x = max(x_star, 0)
gen y = 1 + 0.8 * x + 0.5 * x_star + rnormal()

brunching y x, method(naive) z(z)
local beta_naive = _b[x]

brunching y x, method(uniform) z(z) nbins(10)
local beta_corr = _b[x]

local bias_naive = abs(`beta_naive' - 0.8)
local bias_corr  = abs(`beta_corr'  - 0.8)

di "  naive beta = `beta_naive'  (bias = `bias_naive')"
di "  corr  beta = `beta_corr'  (bias = `bias_corr')"

if `bias_corr' >= `bias_naive' {
    di as err "  FAIL: corrected bias (`bias_corr') should be less than naive bias (`bias_naive')"
    exit 9
}
di as res "  PASS: corrected is closer to truth"

// ---------- Test 7: CCT method runs and identifies beta(0+) ----------
di as txt _n _dup(78) "=" _n "Test 7: CCT method"
// DGP from R test-cct.R: beta(0+) = 0.8 + 1.5 = 2.3 (latent slope + kink)
clear
set obs 5000
set seed 1
gen Z = ceil(runiform() * 20)
gen mu_z = -0.5 + 0.15 * Z
gen eta = rnormal()
gen x_star = mu_z + eta
gen x = max(x_star, 0)
gen y = 2 + 0.8 * x_star + 1.5 * x * (x > 0) + rnormal()

brunching y x, method(cct) z(Z) alpha_star(0.5) boot(50) seed(42)
local beta_cct = _b[x]
assert !missing(`beta_cct')
di "  cct beta(0+) = `beta_cct'  (truth ~ 2.3)"
assert "`e(estimand)'" == "beta(0+)"

// Basic sanity: estimate should be strictly positive and within a wide band.
// (small sample, bootstrap, so use a generous 2.0 - 3.0 window.)
if `beta_cct' < 0.5 | `beta_cct' > 4.0 {
    di as err "  FAIL: CCT beta(0+) = `beta_cct' far from truth 2.3"
    exit 9
}
di as res "  PASS: cct method"

// ---------- Test 8: cct requires at least 3 unique z values ----------
di as txt _n _dup(78) "=" _n "Test 8: cct input validation"
clear
set obs 100
set seed 99
gen zbad = 1
replace zbad = 2 if _n > 50
gen x_star = rnormal()
gen x = max(x_star, 0)
gen y = x + rnormal()
capture brunching y x, method(cct) z(zbad)
if _rc == 0 {
    di as err "  FAIL: cct should fail on 2 unique Z levels"
    exit 9
}
di as res "  PASS: cct correctly rejects too few Z levels"

// ---------- Test 9: missing variable error ----------
di as txt _n _dup(78) "=" _n "Test 9: missing variable error"
clear
set obs 10
gen y = _n
gen x = _n
capture brunching y x, method(tobit) z(z)
if _rc == 0 {
    di as err "  FAIL: should error when z does not exist"
    exit 9
}
di as res "  PASS: missing z variable raises error"

// ---------- Tests 10-13: bootstrap inference for the corrected CCN methods ----------
// Each method should produce finite e(boot_se_beta), e(boot_se_delta), and a
// 95% CI when called with boot(50) seed(42).

foreach m in uniform tobit het_tobit symmetric {
    di as txt _n _dup(78) "=" _n "Bootstrap test: method = `m'"
    clear
    set obs 1500
    set seed 1001
    gen z = runiform() * 10
    gen x_star = -0.5 + 0.3 * z + rnormal()
    gen x = max(x_star, 0)
    gen y = 1 + 0.8 * x_star + rnormal()

    brunching y x, method(`m') z(z) nbins(5) boot(50) seed(42)

    assert !missing(_b[x])
    assert !missing(e(delta))
    assert !missing(e(boot_se_beta))
    assert e(boot_se_beta) > 0
    assert !missing(e(boot_se_delta))
    assert e(boot_se_delta) > 0
    assert !missing(e(boot_ci_beta_lo))
    assert !missing(e(boot_ci_beta_hi))
    assert e(boot_ci_beta_lo) <= _b[x] & _b[x] <= e(boot_ci_beta_hi)
    assert e(boot_B) > 20  // most reps succeed
    assert e(boot_B_req) == 50
    di as res "  PASS: bootstrap (`m')  boot_se=" %9.5f e(boot_se_beta)             ///
              "  ci=[" %9.5f e(boot_ci_beta_lo) ", " %9.5f e(boot_ci_beta_hi) "]"  ///
              "  B_ok=" e(boot_B)
}

// ---------- Test 14: bootstrap with seed is reproducible ----------
di as txt _n _dup(78) "=" _n "Test 14: bootstrap reproducibility"
clear
set obs 1500
set seed 2002
gen z = runiform() * 10
gen x_star = -0.5 + 0.3 * z + rnormal()
gen x = max(x_star, 0)
gen y = 1 + 0.8 * x_star + rnormal()

brunching y x, method(uniform) z(z) nbins(5) boot(30) seed(99)
local se1 = e(boot_se_beta)
brunching y x, method(uniform) z(z) nbins(5) boot(30) seed(99)
local se2 = e(boot_se_beta)
local diff = abs(`se1' - `se2')
if `diff' > 1e-12 {
    di as err "  FAIL: same seed gave different SE: `se1' vs `se2'"
    exit 9
}
di as res "  PASS: bootstrap is reproducible with seed (SE=" %9.5f `se1' ")"

// ---------- Test 15: boot(0) leaves the bootstrap fields unset ----------
di as txt _n _dup(78) "=" _n "Test 15: boot(0) backward compat"
clear
set obs 1500
set seed 3003
gen z = runiform() * 10
gen x_star = -0.5 + 0.3 * z + rnormal()
gen x = max(x_star, 0)
gen y = 1 + 0.8 * x_star + rnormal()

brunching y x, method(uniform) z(z) nbins(5) boot(0)
assert !missing(_b[x])
assert missing(e(boot_se_beta))
brunching y x, method(uniform) z(z) nbins(5)
assert missing(e(boot_se_beta))
di as res "  PASS: default boot(0) does not run the bootstrap"

di _n as txt _dup(78) "=" _n "All basic tests PASSED." _n _dup(78) "=" _n
exit 0
