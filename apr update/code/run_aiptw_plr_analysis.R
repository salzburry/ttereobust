## run_aiptw_plr_analysis.R
##
## Thin wrapper that runs simulate_aiptw.R with --method aiptw_plr over
## the full protocol grid. Defaults to protocol-grade settings; pass any
## simulate_aiptw.R flag through and it overrides the default.
##
## Outputs land under results/raw/aiptw_plr/<scenario>.csv so this
## wrapper can coexist with run_aiptw_cox_analysis.R in the same
## results directory.
##
## Usage:
##   Rscript run_aiptw_plr_analysis.R                 # full grid, defaults
##   Rscript run_aiptw_plr_analysis.R --R 50 --B 50   # smoke test
##   Rscript run_aiptw_plr_analysis.R --workers 16    # parallel

local({
  argv <- commandArgs(trailingOnly = TRUE)
  fill_arg <- function(av, flag, val) {
    if (any(av == flag)) av else c(av, flag, val)
  }
  argv <- fill_arg(argv, "--method",       "aiptw_plr")
  argv <- fill_arg(argv, "--R",            "1900")
  argv <- fill_arg(argv, "--B",            "200")
  argv <- fill_arg(argv, "--rescale-time", "0.25")

  # Override commandArgs so that simulate_aiptw.R picks up our prefilled
  # flags. on.exit() inside this local() block restores the original.
  old <- commandArgs
  commandArgs <<- function(trailingOnly = FALSE) {
    if (trailingOnly) argv else c("RScript", "--", argv)
  }
  on.exit(commandArgs <<- old, add = TRUE)

  source("simulate_aiptw.R")
  main()
})
