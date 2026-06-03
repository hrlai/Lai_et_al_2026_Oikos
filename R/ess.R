# Calculate effective sample size

library(tidyverse)
library(greta)
library(coda)


# read models
mod_native_ABA <- read_rds(file = "out/mod_native_ABA_greta_LVs.rds")
mod_exotic_ABA <- read_rds(file = "out/mod_exotic_ABA_greta_LVs.rds")
mod_native_WAS <- read_rds(file = "out/mod_native_WAS_greta_LVs.rds")
mod_exotic_WAS <- read_rds(file = "out/mod_exotic_WAS_greta_LVs.rds")
mod_all <- read_rds(file = "out/mod_all_greta_LVs.rds")


# calculate effective sample size of the posteriors
ess_native_ABA <- effectiveSize(mod_native_ABA$draws)
ess_exotic_ABA <- effectiveSize(mod_exotic_ABA$draws)
ess_native_WAS <- effectiveSize(mod_native_WAS$draws)
ess_exotic_WAS <- effectiveSize(mod_exotic_WAS$draws)
ess_all <- effectiveSize(mod_all$draws)


# examine ESS
# function to examine fourth-corner coefs
prop_ess_kappa <- function(ess, threshold = 10) {
  proportions(table(ess[str_detect(names(ess), "kappa")] > threshold))
}

prop_ess_kappa(ess_native_ABA, 100)
prop_ess_kappa(ess_exotic_ABA, 100)
prop_ess_kappa(ess_native_WAS, 100)
prop_ess_kappa(ess_exotic_WAS, 100)
prop_ess_kappa(ess_all, 100)
