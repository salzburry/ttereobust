### parameters for simulation:
this.seed <- 2026
# set covariate params
N.Lcovs.linear = 4
N.Lcovs.sq = 2
mu = rep(0,8) # Lcovs.linear + Lcovs.sq + 2 (W and O)
sigma = diag(8)
 
# exposure model params
alpha.L = c(0.005, 0.001, -0.3, 0.1, 0.2, 0.15) # params for exposure model c(L_1, ..., L_{Lcovs.linear + Lcovs.sq})
alpha.W = 1 # coefficient / log-HRs for outcome model
 
# outcome model params
coeff.A = log(0.75)
coeff.L = log(c(0.95, 1.1, 0.95, 0.9, 1.05, 1.01))
coeff.Lsq = log(c(0.945, 1.05))
coeff.O = log(0.9) # outcome-only confounder
coeff.W = log(0.9) # confounder of BOTH exposure and outcome (matches sim_data.R DGM)
 
# params for generating event and censoring times
gamma.tte = 1
lambda.tte = 0.3
lambda.cens = 0.2
admin.cens = 10
