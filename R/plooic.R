# Calculate parameter-to-data ratio

library(tidyverse)
library(greta)
library(loo)
library(extraDistr)

source(file = "R/1_prep_inputs.R")

# read models
mod_native_ABA <- read_rds(file = "out/mod_native_ABA_greta_LVs.rds")
mod_exotic_ABA <- read_rds(file = "out/mod_exotic_ABA_greta_LVs.rds")
mod_native_WAS <- read_rds(file = "out/mod_native_WAS_greta_LVs.rds")
mod_exotic_WAS <- read_rds(file = "out/mod_exotic_WAS_greta_LVs.rds")
mod_all <- read_rds(file = "out/mod_all_greta_LVs.rds")


# calculate effective number of parameters
# function to retrieve log-likelihood
get_loglik <- function(model, origin, landuse, ...) {
  # subset data
  species <- rownames(Z)[Z[, "Exotic1"] %in% origin]
  plots <- rownames(X)[X[, "Type1"] %in% landuse]
  Y_dat <- Y[plots, species]
  # cleanup
  Y_dat <- Y_dat[rowSums(Y_dat > 0) > 0, colSums(Y_dat > 0) > 0]
  patch <- as.numeric(as.factor(Pi[rownames(Y_dat), "patch"]))
  trial_size <- comb.dat$n_trials[rownames(Y_dat)]

  X_var <- c("size", "dist", "nitrogen", "phosphorous", "potassium")
  X_dat <- X[rownames(Y_dat), X_var]

  Z_var <- c("LPC", "WD", "logSM", "Hmax")
  Z_dat <- Z[colnames(Y_dat), Z_var]

  n_site <- nrow(Y_dat) # number of sites
  n_patch <- length(unique(patch)) # number of patches
  n_sp <- ncol(Y_dat) # number of species
  n_cov <- ncol(X_dat) # number of env. covariates

  # prediction
  kappa <- model$mod$target_greta_arrays$kappa
  Sigma_U <- model$mod$target_greta_arrays$Sigma_U
  mu_beta <- kappa %*% t(cbind(1, Z_dat))
  z_beta <- model$mod$target_greta_arrays$z_beta
  beta <- t(mu_beta) + t(z_beta) %*% Sigma_U
  loadings <- model$mod$target_greta_arrays$loadings
  lvs <- model$mod$target_greta_arrays$lvs
  random <- lvs %*% t(loadings)
  fixed <- cbind(1, X_dat) %*% t(beta)
  linpred <- fixed + random
  mu <- imultilogit(linpred)
  theta <- model$mod$visible_greta_arrays$theta
  conc <- mu * theta

  # log likelihood of abundances
  nsim <- dim(model$draws[[1]])[1]
  n_site <- nrow(Y_dat)
  conc <- calculate(conc, values = model$draws, nsim = nsim)[[1]]
  LL_array <- array(NA, dim = c(nsim, n_site))
  for (k in seq_len(nsim)) {
    LL_array[k, ] <- ddirmnom(
      cbind(Y_dat, trial_size - rowSums(Y_dat)),
      size = trial_size,
      alpha = conc[k, , ],
      log = T
    )
  }
  LL_mat <- t(apply(LL_array, 1, as.numeric))

  return(LL_mat)
}

# function to calculate pLOOIC
get_pLOOIC <- function(loglik) {
  r_eff <- relative_eff(exp(loglik), rep(1, nrow(loglik)))
  loo(loglik, r_eff = r_eff, is_method = "sis")
}

# calculate pLOOIC
pLOOIC_nABA <-
  get_loglik(mod_native_ABA, origin = 1, landuse = 1) |> get_pLOOIC()
pLOOIC_eABA <-
  get_loglik(mod_exotic_ABA, origin = -1, landuse = 1) |> get_pLOOIC()
pLOOIC_nWAS <-
  get_loglik(mod_native_WAS, origin = 1, landuse = -1) |> get_pLOOIC()
pLOOIC_eWAS <-
  get_loglik(mod_exotic_WAS, origin = -1, landuse = -1) |> get_pLOOIC()
pLOOIC_all <-
  get_loglik(mod_all, origin = c(-1, 1), landuse = c(-1, 1)) |> get_pLOOIC()

# calculate proportion of number of effective parameter to total sample size
pLOOIC_nABA$estimates["p_loo", "Estimate"] / prod(dim((mod_native_ABA$y)))
pLOOIC_eABA$estimates["p_loo", "Estimate"] / prod(dim((mod_exotic_ABA$y)))
pLOOIC_nWAS$estimates["p_loo", "Estimate"] / prod(dim((mod_native_WAS$y)))
pLOOIC_eWAS$estimates["p_loo", "Estimate"] / prod(dim((mod_exotic_WAS$y)))
pLOOIC_all$estimates["p_loo", "Estimate"] / prod(dim((mod_all$y)))
