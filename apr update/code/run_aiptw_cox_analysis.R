## run_aiptw_cox_analysis.R
##
## Thin wrapper that runs simulate_aiptw.R with --method aiptw_cox over
## the full protocol grid. Sister to run_aiptw_plr_analysis.R; uses
## riskRegression::ate() with a Cox PH outcome model. Written to match
## the working group's run_<method>_analysis.R naming convention.
##
## Defaults are protocol-grade. Override on the CLI as needed.
##
## Usage:
##   Rscript run_aiptw_cox_analysis.R
##   Rscript run_aiptw_cox_analysis.R --R 50 --B 50
##   Rscript run_aiptw_cox_analysis.R --workers 16

argv <- commandArgs(trailingOnly = TRUE)

fill_arg <- function(flag, val) {
  if (any(argv == flag)) argv else c(argv, flag, val)
}
argv <- fill_arg("--method", "aiptw_cox")
argv <- fill_arg("--R",      "1900")
argv <- fill_arg("--B",      "200")

old <- commandArgs
commandArgs <- function(trailingOnly = FALSE) argv
on.exit(commandArgs <- old, add = TRUE)

source("simulate_aiptw.R")
main()
