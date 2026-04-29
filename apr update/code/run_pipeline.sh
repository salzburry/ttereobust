#!/usr/bin/env bash
## run_pipeline.sh
##
## End-to-end AIPTW pipeline: simulate -> summarise -> render Quarto.
## Set R, B, WORKERS, METHOD, RESCALE_TIME via environment variables to
## override defaults. Examples:
##
##   ./run_pipeline.sh                           # protocol defaults
##   R=50 B=50 ./run_pipeline.sh                 # smoke-test sized
##   METHOD=both WORKERS=16 ./run_pipeline.sh    # full grid, both arms
##   RESCALE_TIME=1 ./run_pipeline.sh            # coarse grid for slow pods

set -euo pipefail

R_REPS="${R:-1900}"
B_BOOT="${B:-200}"
WORKERS="${WORKERS:-1}"
METHOD="${METHOD:-aiptw_plr}"
RESCALE_TIME="${RESCALE_TIME:-0.25}"

cd "$(dirname "$0")"

echo "[pipeline] R=${R_REPS} B=${B_BOOT} workers=${WORKERS} method=${METHOD} rescale_time=${RESCALE_TIME}"

# 1. Simulate.
Rscript simulate_aiptw.R \
  --R "${R_REPS}" \
  --B "${B_BOOT}" \
  --workers "${WORKERS}" \
  --method "${METHOD}" \
  --rescale-time "${RESCALE_TIME}"

# 2. Aggregate to results/summary.csv.
Rscript summarise_aiptw.R

# 3. Render the Quarto report (HTML).
if command -v quarto >/dev/null 2>&1; then
  quarto render reports/AIPTW_report.qmd
  echo "[pipeline] report: reports/AIPTW_report.html"
else
  echo "[pipeline] quarto not on PATH; skipping report render."
  echo "[pipeline] install quarto and re-run: quarto render reports/AIPTW_report.qmd"
fi

echo "[pipeline] done"
