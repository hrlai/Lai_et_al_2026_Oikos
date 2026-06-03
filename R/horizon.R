# Prediction horizon

library(tidyverse)
library(here)
library(ggpubr)
library(viridis)


# cross_vals <- read_rds(file = "out/cross_vals_LVs.rds")
KLDs <- read_rds(file = "out/KLDs.rds")
KLDs_posterior <- read_rds(file = "out/KLDs_posterior.rds")

# read models
mod_native_ABA <- read_rds(file = "out/mod_native_ABA_greta_LVs.rds")
mod_exotic_ABA <- read_rds(file = "out/mod_exotic_ABA_greta_LVs.rds")
mod_native_WAS <- read_rds(file = "out/mod_native_WAS_greta_LVs.rds")
mod_exotic_WAS <- read_rds(file = "out/mod_exotic_WAS_greta_LVs.rds")
mod_all <- read_rds(file = "out/mod_all_greta_LVs.rds")

# summarise loglik to block levels
LL_summary <-
  lapply(cross_vals2, function(x) {
    quantile(rowMeans(x$LL_mat), probs = c(0.05, 0.5, 0.95))
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame %>%
  rownames_to_column("tmp") %>%
  separate(tmp, c("from", "to"), sep = "_to_") %>%
  mutate(in_sample = from == to)

# tidy up KL divergences
KLD_X <-
  unlist(lapply(KLDs, function(x) x$X)) %>%
  data.frame(KLD_X = .) %>%
  rownames_to_column("tmp") %>%
  separate(tmp, c("from", "to"), sep = "_to_")
KLD_Z <-
  unlist(lapply(KLDs, function(x) x$Z)) %>%
  data.frame(KLD_Z = .) %>%
  rownames_to_column("tmp") %>%
  separate(tmp, c("from", "to"), sep = "_to_")
KLD_posterior <-
  unlist(KLDs_posterior) %>%
  data.frame(KLD_posterior = .) %>%
  rownames_to_column("tmp") %>%
  separate(tmp, c("from", "to"), sep = "_to_") %>%
  mutate(KLD_posterior = ifelse(KLD_posterior == -Inf, 0, KLD_posterior))

Dist_LL <-
  LL_summary %>%
  left_join(KLD_X) %>%
  left_join(KLD_Z) %>%
  left_join(KLD_posterior)


# Visualise ---------------------------------------------------------------
Dist_LL_plot <-
  Dist_LL %>%
  filter(to != "all", from != "all") %>%
  mutate(
    from_landuse = str_extract(from, "(ABA|WAS)"),
    from_origin = str_extract(from, "(e|n)")
  )

p_Y_env_KL <-
  ggplot(Dist_LL_plot) +
  geom_hline(yintercept = 0, colour = "grey", linewidth = 0.3) +
  geom_errorbar(
    aes(KLD_X, ymin = `5%`, ymax = `95%`),
    width = 0,
    colour = "darkgrey"
  ) +
  geom_point(
    aes(KLD_X, `50%`, fill = from_origin, shape = from_landuse),
    size = 2
  ) +
  scale_fill_manual(values = c(turbo(1, 1, 0.85), turbo(1, 1, 0.15))) +
  scale_shape_manual(values = c(24, 25)) +
  labs(
    x = "Environment-data divergence",
    y = "Prediction accuracy\n(log likelihood ratio)"
  ) +
  theme_classic() +
  theme(legend.position = "none")
p_Y_trait_KL <-
  ggplot(Dist_LL_plot) +
  geom_hline(yintercept = 0, colour = "grey", linewidth = 0.3) +
  geom_errorbar(
    aes(KLD_Z, ymin = `5%`, ymax = `95%`),
    width = 0,
    colour = "darkgrey"
  ) +
  geom_point(
    aes(KLD_Z, `50%`, fill = from_origin, shape = from_landuse),
    size = 2
  ) +
  scale_fill_manual(values = c(turbo(1, 1, 0.85), turbo(1, 1, 0.15))) +
  scale_shape_manual(values = c(24, 25)) +
  labs(
    x = "Trait-data divergence",
    y = "Prediction accuracy\n(log likelihood ratio)"
  ) +
  theme_classic() +
  theme(legend.position = "none")
p_Y_param_KL <-
  ggplot(Dist_LL_plot) +
  geom_hline(yintercept = 0, colour = "grey", linewidth = 0.3) +
  geom_errorbar(
    aes(KLD_posterior, ymin = `5%`, ymax = `95%`),
    width = 0,
    colour = "darkgrey"
  ) +
  geom_point(
    aes(KLD_posterior, `50%`, fill = from_origin, shape = from_landuse),
    size = 2
  ) +
  scale_fill_manual(values = c(turbo(1, 1, 0.85), turbo(1, 1, 0.15))) +
  scale_shape_manual(values = c(24, 25)) +
  labs(
    x = "Parameter-estimate divergence",
    y = "Prediction accuracy\n(log likelihood ratio)"
  ) +
  theme_classic() +
  theme(legend.position = "none")

png(here("fig/horizon.png"), width = 3.5, height = 8, units = "in", res = 300)
ggarrange(
  p_Y_env_KL,
  p_Y_trait_KL,
  p_Y_param_KL,
  ncol = 1,
  nrow = 3,
  labels = paste0("(", letters[1:3], ")"),
  hjust = 0
)
dev.off()


# Quick regression
horizon_dat <-
  Dist_LL %>%
  filter(to != "all", from != "all") %>%
  mutate(
    accuracy = `50%`,
    KLD_X_s = KLD_X / sd(KLD_X),
    KLD_Z_s = KLD_Z / sd(KLD_Z)
  )
horizon_fit <-
  lm(`50%` ~ 0 + KLD_X_s * KLD_Z_s, data = horizon_dat)
summary(horizon_fit)

horizon_fit_robust <-
  quantreg::rq(`50%` ~ 0 + KLD_X_s * KLD_Z_s, data = horizon_dat, tau = 0.5)
summary(horizon_fit_robust)


# 3D plot
library(scatterplot3d)

with(
  horizon_dat,
  scatterplot3d(
    KLD_X,
    KLD_Z,
    `50%`,
    angle = 10,
    # box = FALSE,
    # highlight.3d = TRUE,
    # type="h"
  )
)


# Predictor space ---------------------------------------------------------

load(here("Data/Clean/comb.dat.Rdata"))

X_pca <-
  prcomp(
    comb.dat$X[, c("size", "dist", "nitrogen", "phosphorous", "potassium")],
    scale. = TRUE
  )
X_pca_coord <-
  predict(X_pca) %>%
  cbind(Type = comb.dat$X[, "Type1"]) %>%
  as.data.frame()
X_pca_arrows <- data.frame(
  x = rep(0, 5),
  y = rep(0, 5),
  xend = X_pca$rotation[, 1],
  yend = X_pca$rotation[, 2]
)
X_arrow_scale <- 4
p_X_pca <-
  ggplot() +
  geom_point(
    data = X_pca_coord,
    aes(PC1, PC2, shape = as.factor(Type), fill = as.factor(Type)),
    size = 2
  ) +
  geom_segment(
    data = X_pca_arrows,
    aes(
      x * X_arrow_scale,
      y * X_arrow_scale,
      xend = xend * X_arrow_scale,
      yend = yend * X_arrow_scale
    ),
    arrow = arrow(length = unit(0.03, "npc"))
  ) +
  geom_text(
    data = X_pca_arrows,
    aes(
      xend * X_arrow_scale,
      yend * X_arrow_scale,
      label = c("Patch size", "Distance", "N", "P", "K")
    ),
    vjust = c(0, 1, 0, 0, 0)
  ) +
  scale_shape_manual(
    values = c(25, 24),
    name = "Land use:",
    labels = c("WW", "AL")
  ) +
  scale_fill_manual(
    values = c("white", "grey"),
    name = "Land use:",
    labels = c("WW", "AL")
  ) +
  labs(x = "Environment PC1 (36%)", y = "Environment PC2 (29%)") +
  theme_classic() +
  theme(
    legend.position = "bottom",
    legend.margin = margin(0, 0, 0, 0),
    legend.title = element_text(size = 10)
  )

Z_pca <-
  prcomp(comb.dat$Z[, c("LPC", "WD", "logSM", "Hmax")], scale. = TRUE)
Z_pca_coord <-
  predict(Z_pca) %>%
  cbind(Exotic = comb.dat$Z[, "Exotic1"]) %>%
  as.data.frame()
Z_pca_arrows <- data.frame(
  x = rep(0, 4),
  y = rep(0, 4),
  xend = Z_pca$rotation[, 1],
  yend = Z_pca$rotation[, 2]
)
Z_arrow_scale <- 3
p_Z_pca <-
  ggplot() +
  geom_point(
    data = Z_pca_coord,
    aes(PC1, PC2, fill = as.factor(Exotic)),
    pch = 21,
    size = 2
  ) +
  geom_segment(
    data = Z_pca_arrows,
    aes(
      x * Z_arrow_scale,
      y * Z_arrow_scale,
      xend = xend * Z_arrow_scale,
      yend = yend * Z_arrow_scale
    ),
    arrow = arrow(length = unit(0.03, "npc"))
  ) +
  geom_text(
    data = Z_pca_arrows,
    aes(
      xend * Z_arrow_scale,
      yend * Z_arrow_scale,
      label = c("LPC", "WD", "SDM", "Hmax")
    ),
    hjust = 1
  ) +
  scale_fill_manual(
    values = viridis::turbo(2, begin = 0.8, end = 0.2),
    name = "Origin:",
    labels = c("Exotic", "Native")
  ) +
  labs(x = "Trait PC1 (40%)", y = "Trait PC2 (23%)") +
  theme_classic() +
  theme(
    legend.position = "bottom",
    legend.margin = margin(0, 0, 0, 0),
    legend.title = element_text(size = 10)
  )


png(
  here("fig/predictor_space.png"),
  width = 7,
  height = 3.5,
  units = "in",
  res = 300
)
ggarrange(
  p_Z_pca,
  p_X_pca,
  ncol = 2,
  nrow = 1,
  heights = c(1, 1.2),
  labels = paste0("(", letters[1:2], ")"),
  hjust = 0
)
dev.off()
