# Measure Kullback-Leibler divergence between train and test datasets

library(greta)
library(tidyverse)
library(kldest)

# read models
mod_native_ABA <- read_rds(file = "out/mod_native_ABA_greta_LVs.rds")
mod_exotic_ABA <- read_rds(file = "out/mod_exotic_ABA_greta_LVs.rds")
mod_native_WAS <- read_rds(file = "out/mod_native_WAS_greta_LVs.rds")
mod_exotic_WAS <- read_rds(file = "out/mod_exotic_WAS_greta_LVs.rds")
mod_all <- read_rds(file = "out/mod_all_greta_LVs.rds")

# read data
load(file = "Data/Clean/comb.dat.Rdata")

# function to calculate KL divergence
KLD <- function(test, train) {
  # global variables
  X_var <- c("size", "dist", "nitrogen", "phosphorous", "potassium")
  Z_var <- c("LPC", "WD", "logSM", "Hmax")

  # test dataset variables
  plot_test <- rownames(test$y)
  species_test <- colnames(test$y)
  X_test <- comb.dat$X[plot_test, X_var, drop = FALSE]
  Z_test <- comb.dat$Z[species_test, Z_var, drop = FALSE]

  # train dataset variables
  plot_train <- rownames(train$y)
  species_train <- colnames(train$y)
  X_train <- comb.dat$X[plot_train, X_var, drop = FALSE]
  Z_train <- comb.dat$Z[species_train, Z_var, drop = FALSE]

  # KL divergence
  KLD_X <- kldest:::kld_est_kde(X_test, X_train)
  KLD_Z <- kldest:::kld_est_kde(Z_test, Z_train)
  return(list(X = KLD_X, Z = KLD_Z))
}

# combine and save outputs
KLDs <- list(
  nABA_to_nABA = KLD(mod_native_ABA, mod_native_ABA),
  eABA_to_nABA = KLD(mod_native_ABA, mod_exotic_ABA),
  nWAS_to_nABA = KLD(mod_native_ABA, mod_native_WAS),
  eWAS_to_nABA = KLD(mod_native_ABA, mod_exotic_WAS),
  all_to_nABA = KLD(mod_native_ABA, mod_all),

  nABA_to_eABA = KLD(mod_exotic_ABA, mod_native_ABA),
  eABA_to_eABA = KLD(mod_exotic_ABA, mod_exotic_ABA),
  nWAS_to_eABA = KLD(mod_exotic_ABA, mod_native_WAS),
  eWAS_to_eABA = KLD(mod_exotic_ABA, mod_exotic_WAS),
  all_to_eABA = KLD(mod_exotic_ABA, mod_all),

  nABA_to_nWAS = KLD(mod_native_WAS, mod_native_ABA),
  eABA_to_nWAS = KLD(mod_native_WAS, mod_exotic_ABA),
  nWAS_to_nWAS = KLD(mod_native_WAS, mod_native_WAS),
  eWAS_to_nWAS = KLD(mod_native_WAS, mod_exotic_WAS),
  all_to_nWAS = KLD(mod_native_WAS, mod_all),

  nABA_to_eWAS = KLD(mod_exotic_WAS, mod_native_ABA),
  eABA_to_eWAS = KLD(mod_exotic_WAS, mod_exotic_ABA),
  nWAS_to_eWAS = KLD(mod_exotic_WAS, mod_native_WAS),
  eWAS_to_eWAS = KLD(mod_exotic_WAS, mod_exotic_WAS),
  all_to_eWAS = KLD(mod_exotic_WAS, mod_all)
)
write_rds(KLDs, file = "out/KLDs.rds")
