#!/usr/bin/env bash
# Cross-check runner: R generates data + R estimates -> Stata estimates -> diff.
set -euo pipefail

cd "$(dirname "$0")"

STATA=/Applications/StataNow/StataMP.app/Contents/MacOS/stata-mp

echo "=== Step 1: R (generate data + R estimates) ==="
Rscript 01_gen_and_r.R

echo
echo "=== Step 2: Stata cross-check ==="
"$STATA" -b do 02_cross_check.do
# Stata batch-mode writes to 02_cross_check.log; dump it so errors are visible
if [[ -f 02_cross_check.log ]]; then
  tail -n 50 02_cross_check.log
fi

echo
echo "=== Step 3: compare ==="
Rscript 03_compare.R
