#' Cross-check step 3: compare R and Stata estimates side by side.

out_dir <- "tests/stata"
if (!file.exists(file.path(out_dir, "r_results.csv")))     out_dir <- "."
if (!file.exists(file.path(out_dir, "stata_results.csv"))) out_dir <- "."

r_res  <- read.csv(file.path(out_dir, "r_results.csv"))
st_res <- read.csv(file.path(out_dir, "stata_results.csv"))

names(r_res)[-1]  <- paste0("r_",     names(r_res)[-1])
names(st_res)[-1] <- paste0("stata_", names(st_res)[-1])

cmp <- merge(r_res, st_res, by = "method", sort = FALSE)
cmp <- cmp[match(c("naive","uniform","het_tobit","symmetric"), cmp$method), ]

cmp$beta_diff  <- cmp$r_beta  - cmp$stata_beta
cmp$delta_diff <- cmp$r_delta - cmp$stata_delta

cat("\n=== R vs Stata (true beta = 0.8) ===\n")
print(cmp, row.names = FALSE, digits = 6)

tol <- 1e-4
failed <- abs(cmp$beta_diff) > tol
if (any(failed, na.rm = TRUE)) {
  cat("\nMISMATCH on beta (tolerance ", tol, "):\n", sep = "")
  print(cmp[failed, c("method","r_beta","stata_beta","beta_diff")], row.names = FALSE)
  quit(status = 1)
} else {
  cat("\nPASS: all beta estimates agree within ", tol, ".\n", sep = "")
}
