library(tidyverse)
library(greta)

# source code that wrangle data
source(file = "R/1_prep_inputs.R")
source(file = "R/fit_greta.R")
source(file = "R/randortho.R")
source(file = "R/geweke.R")

# fit greta models to different training sets
mod_native_ABA <-
    fit_greta(
        origin = 1,
        landuse = 1,
        family = "dirichlet_multinomial",
        chains = 1,
        warmup = 3000,
        sampler = hmc(30, 35),
        one_by_one = TRUE
    )
mod_exotic_ABA <-
    fit_greta(
        origin = -1,
        landuse = 1,
        family = "dirichlet_multinomial",
        chains = 1,
        warmup = 3000,
        sampler = hmc(30, 35),
        one_by_one = TRUE
    )
mod_native_WAS <-
    fit_greta(
        origin = 1,
        landuse = -1,
        family = "dirichlet_multinomial",
        chains = 1,
        warmup = 3000,
        sampler = hmc(30, 35),
        one_by_one = TRUE
    )
mod_exotic_WAS <-
    fit_greta(
        origin = -1,
        landuse = -1,
        family = "dirichlet_multinomial",
        chains = 1,
        warmup = 3000,
        sampler = hmc(30, 35),
        one_by_one = TRUE
    )
mod_all <-
    fit_greta(
        origin = c(-1, 1),
        landuse = c(-1, 1),
        family = "dirichlet_multinomial",
        chains = 1,
        warmup = 3000,
        sampler = hmc(30, 35)
    )

# save models
write_rds(mod_native_ABA, file = "out/mod_native_ABA_greta_LVs.rds")
write_rds(mod_exotic_ABA, file = "out/mod_exotic_ABA_greta_LVs.rds")
write_rds(mod_native_WAS, file = "out/mod_native_WAS_greta_LVs.rds")
write_rds(mod_exotic_WAS, file = "out/mod_exotic_WAS_greta_LVs.rds")
write_rds(mod_all, file = "out/mod_all_greta_LVs.rds")
