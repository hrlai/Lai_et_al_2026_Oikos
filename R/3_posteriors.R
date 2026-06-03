library(tidyverse)
library(tidybayes)
library(greta)
library(GGally)

source(file = "R/1_prep_inputs.R")

# read models
mod_native_ABA <- read_rds(file = "out/mod_native_ABA_greta_LVs.rds")
mod_exotic_ABA <- read_rds(file = "out/mod_exotic_ABA_greta_LVs.rds")
mod_native_WAS <- read_rds(file = "out/mod_native_WAS_greta_LVs.rds")
mod_exotic_WAS <- read_rds(file = "out/mod_exotic_WAS_greta_LVs.rds")
mod_all <- read_rds(file = "out/mod_all_greta_LVs.rds")


# Fourth corners ----------------------------------------------------------

env_match <-
  data.frame(
    k = 1:6,
    env = c("beta0", "size", "dist", "nitrogen", "phosphorous", "potassium")
  )
trait_match <-
  data.frame(
    m = 1:5,
    trait = c("T0", "LPC", "WD", "logSM", "Hmax")
  )

summarise_kappa <- function(draws, origin, landuse) {
  draws %>%
    gather_draws(kappa[k, m]) %>%
    median_qi(.width = 0.9) %>%
    left_join(env_match) %>%
    left_join(trait_match) %>%
    mutate(origin = origin, landuse = landuse)
}

kappa <-
  bind_rows(
    summarise_kappa(mod_native_ABA$draws, "Native", "AL"),
    summarise_kappa(mod_native_WAS$draws, "Native", "WW"),
    summarise_kappa(mod_exotic_ABA$draws, "Exotic", "AL"),
    summarise_kappa(mod_exotic_WAS$draws, "Exotic", "WW"),
    summarise_kappa(mod_all$draws, "All", "All")
  ) %>%
  mutate(
    origin = fct_relevel(origin, c("Native", "All", "Exotic")),
    landuse = factor(landuse, levels = c("AL", "All", "WW")),
    env = factor(
      env,
      levels = c(
        "beta0",
        "size",
        "dist",
        "nitrogen",
        "phosphorous",
        "potassium"
      ),
      labels = c(
        expression(beta[0]),
        "Patch~size",
        "Distance~from~old~growth",
        "Soil~nitrogen",
        "Soil~phosphorous",
        "Soil~potassium"
      )
    ),
    trait = factor(
      trait,
      levels = c("T0", "LPC", "WD", "logSM", "Hmax"),
      labels = c("T[0]", "LPC", "WD", "SDM", "H[max]")
    ),
    sig = ifelse(sign(.lower) == sign(.upper), "yes", "no")
  )


# Plot --------------------------------------------------------------------

png("fig/fourth_corners.png", width = 7, height = 5, units = "in", res = 300)
ggplot(kappa) +
  geom_stripped_rows(aes(y = env), even = "lightgrey", odd = "white") +
  geom_vline(xintercept = 0, colour = "darkgrey") +
  geom_errorbar(
    aes(
      y = env,
      xmin = .lower,
      xmax = .upper,
      alpha = sig,
      group = interaction(origin, landuse),
    ),
    linewidth = 1,
    width = 0,
    position = position_dodge(.5, orientation = "y")
  ) +
  geom_point(
    aes(
      x = .value,
      y = env,
      shape = landuse,
      fill = origin,
      alpha = sig
    ),
    size = 2,
    position = position_dodge(.5)
  ) +
  facet_wrap(~trait, nrow = 1, scale = "free_x", labeller = "label_parsed") +
  scale_fill_viridis_d(
    option = "turbo",
    begin = 0.2,
    end = 0.8,
    name = "Origin:"
  ) +
  scale_shape_manual(values = c(24, 21, 25), name = "Land use:") +
  scale_alpha_manual(values = c(0.3, 1)) +
  scale_y_discrete(
    limits = rev,
    labels = scales::label_parse(),
    expand = c(0, 0)
  ) +
  guides(
    fill = guide_legend(override.aes = list(shape = 21)),
    shape = guide_legend(override.aes = list(fill = "white")),
    alpha = "none"
  ) +
  labs(x = "Trait-environment interaction coefficients") +
  theme_classic() +
  theme(
    legend.position = "top",
    legend.justification = 5,
    legend.margin = margin(0, 0, 0, 0),
    plot.title = element_text(size = 12),
    strip.text = element_text(size = 11),
    strip.background = element_blank(),
    axis.text.y = element_text(colour = "black", size = 11),
    axis.text.x = element_text(colour = "black"),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.line.y = element_blank(),
    legend.text = element_text(size = 10)
  )
dev.off()


# function to calculate KL divergence in posterior space
KLD_posterior <- function(mod_test, mod_train, ndraws = 100) {
  param_test <-
    mod_test$draws %>%
    spread_draws(kappa[k, m], ndraws = ndraws) %>%
    select(.draw, k, m, kappa) %>%
    pivot_wider(names_from = c(k, m), values_from = kappa)

  param_train <-
    mod_train$draws %>%
    spread_draws(kappa[k, m], ndraws = ndraws) %>%
    select(.draw, k, m, kappa) %>%
    pivot_wider(names_from = c(k, m), values_from = kappa)

  # KL divergence
  kldest::kld_est_nn(param_test, param_train)
}

KLDs_posterior <- list(
  nABA_to_nABA = KLD_posterior(mod_native_ABA, mod_native_ABA),
  eABA_to_nABA = KLD_posterior(mod_native_ABA, mod_exotic_ABA),
  nWAS_to_nABA = KLD_posterior(mod_native_ABA, mod_native_WAS),
  eWAS_to_nABA = KLD_posterior(mod_native_ABA, mod_exotic_WAS),
  all_to_nABA = KLD_posterior(mod_native_ABA, mod_all),

  nABA_to_eABA = KLD_posterior(mod_exotic_ABA, mod_native_ABA),
  eABA_to_eABA = KLD_posterior(mod_exotic_ABA, mod_exotic_ABA),
  nWAS_to_eABA = KLD_posterior(mod_exotic_ABA, mod_native_WAS),
  eWAS_to_eABA = KLD_posterior(mod_exotic_ABA, mod_exotic_WAS),
  all_to_eABA = KLD_posterior(mod_exotic_ABA, mod_all),

  nABA_to_nWAS = KLD_posterior(mod_native_WAS, mod_native_ABA),
  eABA_to_nWAS = KLD_posterior(mod_native_WAS, mod_exotic_ABA),
  nWAS_to_nWAS = KLD_posterior(mod_native_WAS, mod_native_WAS),
  eWAS_to_nWAS = KLD_posterior(mod_native_WAS, mod_exotic_WAS),
  all_to_nWAS = KLD_posterior(mod_native_WAS, mod_all),

  nABA_to_eWAS = KLD_posterior(mod_exotic_WAS, mod_native_ABA),
  eABA_to_eWAS = KLD_posterior(mod_exotic_WAS, mod_exotic_ABA),
  nWAS_to_eWAS = KLD_posterior(mod_exotic_WAS, mod_native_WAS),
  eWAS_to_eWAS = KLD_posterior(mod_exotic_WAS, mod_exotic_WAS),
  all_to_eWAS = KLD_posterior(mod_exotic_WAS, mod_all)
)
write_rds(KLDs_posterior, file = "out/KLDs_posterior.rds")
