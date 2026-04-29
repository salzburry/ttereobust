## run_aiptw_plr_analysis.R
##
## Thin wrapper that runs simulate_aiptw.R with --method aiptw_plr over
## the full protocol grid using the YAML configs in config/. Written to
## match the working group's run_<method>_analysis.R naming convention so
## the cross-method Comparison.qmd can pick up these outputs the same
## way it does for IPTW / G-comp / TMLE.
##
## Defaults are the protocol-grade settings (R = 1900, B = 200,
## --rescale-time 0.25). Override with --R, --B, etc. as needed.
##
## Usage:
##   Rscript run_aiptw_plr_analysis.R                 # full grid, defaults
##   Rscript run_aiptw_plr_analysis.R --R 50 --B 50   # smoke test
##   Rscript run_aiptw_plr_analysis.R --workers 16    # parallel

argv <- commandArgs(trailingOnly = TRUE)
default_args <- c("--method", "aiptw_plr",
                  "--R", "1900",
                  "--B", "200",
                  "--rescale-time", "0.25")

# Forward the user's args; default_args fills in anything they didn't set.
fill_arg <- function(flag, val) {
  if (any(argv == flag)) argv else c(argv, flag, val)
}
argv <- fill_arg("--method",       "aiptw_plr")
argv <- fill_arg("--R",            "1900")
argv <- fill_arg("--B",            "200")
argv <- fill_arg("--rescale-time", "0.25")

# Stash and re-pass so simulate_aiptw.R's parse_cli_args() picks them up.
old <- commandArgs
commandArgs <- function(trailingOnly = FALSE) argv
on.exit(commandArgs <- old, add = TRUE)

source("simulate_aiptw.R")
main()
