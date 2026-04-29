## utils/cor_matrix.R
## Build the 8x8 covariance matrix for the simulated covariates.
##
## Per Apr 28 methods document Section 6.2:
##   "an 8x8 correlation matrix that specifies none to high levels of
##    correlation to evaluate the impact of collinearity"
##
## Implementation choice: the L1..L6 block is given a compound-symmetric
## correlation rho_L. W and O are kept independent of each other and of the
## L block. Setting rho_L = 0 reproduces the previous diag(8) behaviour.

make_sigma <- function(rho_L = 0,
                        n_L = 6,
                        n_extra = 2) {
  stopifnot(rho_L >= 0, rho_L < 1)
  k <- n_L + n_extra
  sigma <- diag(k)
  if (rho_L > 0) {
    L_block <- matrix(rho_L, n_L, n_L)
    diag(L_block) <- 1
    sigma[seq_len(n_L), seq_len(n_L)] <- L_block
  }
  sigma
}
