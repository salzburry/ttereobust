## run_aiptw_cox_analysis.R
##
## Thin wrapper that runs simulate_aiptw.R with --method aiptw_cox over
## the full protocol grid. Sister to run_aiptw_plr_analysis.R; uses
## riskRegression::ate() with a Cox PH outcome model and a marginal
## censoring model.
##
## Outputs land under results/raw/aiptw_cox/<scenario>.csv so this
## wrapper can coexist with run_aiptw_plr_analysis.R in the same
## results directory.
##
## Usage:
##   Rscript run_aiptw_cox_analysis.R
##   Rscript run_aiptw_cox_analysis.R --R 50 --B 50
##   Rscript run_aiptw_cox_analysis.R --workers 16

local({
  argv <- commandArgs(trailingOnly = TRUE)
  fill_arg <- function(av, flag, val) {
    if (any(av == flag)) av else c(av, flag, val)
  }
  argv <- fill_arg(argv, "--method", "aiptw_cox")
  argv <- fill_arg(argv, "--R",      "1900")
  argv <- fill_arg(argv, "--B",      "200")

  old <- commandArgs
  commandArgs <<- function(trailingOnly = FALSE) {
    if (trailingOnly) argv else c("RScript", "--", argv)
  }
  on.exit(commandArgs <<- old, add = TRUE)

  source("simulate_aiptw.R")
  main()
})
