# Cross validations

library(greta)
library(loo)
library(scoringRules)
library(tidyverse)
library(tidybayes)
library(scales)
library(ggpubr)


source(file = "R/1_prep_inputs.R")
source(file = "R/ic.R")

# read models
mod_native_ABA <- readRDS(file = "out/mod_native_ABA_greta_LVs.rds")
mod_exotic_ABA <- readRDS(file = "out/mod_exotic_ABA_greta_LVs.rds")
mod_native_WAS <- readRDS(file = "out/mod_native_WAS_greta_LVs.rds")
mod_exotic_WAS <- readRDS(file = "out/mod_exotic_WAS_greta_LVs.rds")
mod_all <- readRDS(file = "out/mod_all_greta_LVs.rds")


# function to predict onto new data
crossval <- function(model, origin, landuse, ...) {
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
  mu_beta <- kappa %*% t(cbind(1, Z_dat))

  beta <- t(mu_beta)

  fixed <- cbind(1, X_dat) %*% t(beta)
  linpred <- fixed

  mu <- imultilogit(linpred)
  theta <- model$mod$visible_greta_arrays$theta
  conc <- mu * theta
  Y_hat <- dirichlet_multinomial(cbind(trial_size), conc)

  Y_hat_sim <- calculate(Y_hat, values = model$draws, nsim = 1000)
  Y_hat_sim2 <- calculate(Y_hat, values = model$draws, nsim = 1000)

  # relative abundances
  p <- dirichlet(conc)
  # species richness (not run because we don't know how many species in OTHERS)
  SD0 <- rowSums(Y_hat > 0)
  # Inverse Simpson diversity
  SD2 <- 1 / rowSums(p^2)

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

  return(list(
    prediction = Y_hat_sim,
    LL_mat = LL_mat,
    newY = Y_dat,
    newX = X_dat,
    newtraits = Z_dat
  ))
}

# cross validate and calculate prediction accuracy
nABA_to_nABA <- crossval(mod_native_ABA, origin = 1, landuse = 1)
nABA_to_eABA <- crossval(mod_native_ABA, origin = -1, landuse = 1)
nABA_to_nWAS <- crossval(mod_native_ABA, origin = 1, landuse = -1)
nABA_to_eWAS <- crossval(mod_native_ABA, origin = -1, landuse = -1)

eABA_to_nABA <- crossval(mod_exotic_ABA, origin = 1, landuse = 1)
eABA_to_eABA <- crossval(mod_exotic_ABA, origin = -1, landuse = 1)
eABA_to_nWAS <- crossval(mod_exotic_ABA, origin = 1, landuse = -1)
eABA_to_eWAS <- crossval(mod_exotic_ABA, origin = -1, landuse = -1)

nWAS_to_nABA <- crossval(mod_native_WAS, origin = 1, landuse = 1)
nWAS_to_eABA <- crossval(mod_native_WAS, origin = -1, landuse = 1)
nWAS_to_nWAS <- crossval(mod_native_WAS, origin = 1, landuse = -1)
nWAS_to_eWAS <- crossval(mod_native_WAS, origin = -1, landuse = -1)

eWAS_to_nABA <- crossval(mod_exotic_WAS, origin = 1, landuse = 1)
eWAS_to_eABA <- crossval(mod_exotic_WAS, origin = -1, landuse = 1)
eWAS_to_nWAS <- crossval(mod_exotic_WAS, origin = 1, landuse = -1)
eWAS_to_eWAS <- crossval(mod_exotic_WAS, origin = -1, landuse = -1)

all_to_nABA <- crossval(mod_all, origin = 1, landuse = 1)
all_to_eABA <- crossval(mod_all, origin = -1, landuse = 1)
all_to_nWAS <- crossval(mod_all, origin = 1, landuse = -1)
all_to_eWAS <- crossval(mod_all, origin = -1, landuse = -1)

all_to_all <- crossval(mod_all, origin = c(-1, 1), landuse = c(-1, 1))

# combine and save outputs
cross_vals <- list(
  nABA_to_nABA = nABA_to_nABA,
  eABA_to_nABA = eABA_to_nABA,
  nWAS_to_nABA = nWAS_to_nABA,
  eWAS_to_nABA = eWAS_to_nABA,
  all_to_nABA = all_to_nABA,

  nABA_to_eABA = nABA_to_eABA,
  eABA_to_eABA = eABA_to_eABA,
  nWAS_to_eABA = nWAS_to_eABA,
  eWAS_to_eABA = eWAS_to_eABA,
  all_to_eABA = all_to_eABA,

  nABA_to_nWAS = nABA_to_nWAS,
  eABA_to_nWAS = eABA_to_nWAS,
  nWAS_to_nWAS = nWAS_to_nWAS,
  eWAS_to_nWAS = eWAS_to_nWAS,
  all_to_nWAS = all_to_nWAS,

  nABA_to_eWAS = nABA_to_eWAS,
  eABA_to_eWAS = eABA_to_eWAS,
  nWAS_to_eWAS = nWAS_to_eWAS,
  eWAS_to_eWAS = eWAS_to_eWAS,
  all_to_eWAS = all_to_eWAS,

  all_to_all = all_to_all
)
write_rds(cross_vals, file = "out/cross_vals_LVs.rds")

# visualise cross-validation accuracies
cross_vals <- read_rds(file = "out/cross_vals_LVs.rds")

# subtract by all_to_* (log likelihood ratios)
cross_vals2 <- cross_vals
cross_vals2$nABA_to_nABA$LL_mat <-
  cross_vals2$nABA_to_nABA$LL_mat - cross_vals2$all_to_nABA$LL_mat
cross_vals2$eABA_to_nABA$LL_mat <-
  cross_vals2$eABA_to_nABA$LL_mat - cross_vals2$all_to_nABA$LL_mat
cross_vals2$nWAS_to_nABA$LL_mat <-
  cross_vals2$nWAS_to_nABA$LL_mat - cross_vals2$all_to_nABA$LL_mat
cross_vals2$eWAS_to_nABA$LL_mat <-
  cross_vals2$eWAS_to_nABA$LL_mat - cross_vals2$all_to_nABA$LL_mat
cross_vals2$nABA_to_eABA$LL_mat <-
  cross_vals2$nABA_to_eABA$LL_mat - cross_vals2$all_to_eABA$LL_mat
cross_vals2$eABA_to_eABA$LL_mat <-
  cross_vals2$eABA_to_eABA$LL_mat - cross_vals2$all_to_eABA$LL_mat
cross_vals2$nWAS_to_eABA$LL_mat <-
  cross_vals2$nWAS_to_eABA$LL_mat - cross_vals2$all_to_eABA$LL_mat
cross_vals2$eWAS_to_eABA$LL_mat <-
  cross_vals2$eWAS_to_eABA$LL_mat - cross_vals2$all_to_eABA$LL_mat
cross_vals2$nABA_to_nWAS$LL_mat <-
  cross_vals2$nABA_to_nWAS$LL_mat - cross_vals2$all_to_nWAS$LL_mat
cross_vals2$eABA_to_nWAS$LL_mat <-
  cross_vals2$eABA_to_nWAS$LL_mat - cross_vals2$all_to_nWAS$LL_mat
cross_vals2$nWAS_to_nWAS$LL_mat <-
  cross_vals2$nWAS_to_nWAS$LL_mat - cross_vals2$all_to_nWAS$LL_mat
cross_vals2$eWAS_to_nWAS$LL_mat <-
  cross_vals2$eWAS_to_nWAS$LL_mat - cross_vals2$all_to_nWAS$LL_mat
cross_vals2$nABA_to_eWAS$LL_mat <-
  cross_vals2$nABA_to_eWAS$LL_mat - cross_vals2$all_to_eWAS$LL_mat
cross_vals2$eABA_to_eWAS$LL_mat <-
  cross_vals2$eABA_to_eWAS$LL_mat - cross_vals2$all_to_eWAS$LL_mat
cross_vals2$nWAS_to_eWAS$LL_mat <-
  cross_vals2$nWAS_to_eWAS$LL_mat - cross_vals2$all_to_eWAS$LL_mat
cross_vals2$eWAS_to_eWAS$LL_mat <-
  cross_vals2$eWAS_to_eWAS$LL_mat - cross_vals2$all_to_eWAS$LL_mat

LL_summary <-
  lapply(cross_vals2, function(x) {
    quantile(rowMeans(x$LL_mat), probs = c(0.05, 0.5, 0.95))
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame %>%
  rownames_to_column("tmp") %>%
  separate(tmp, c("from", "to"), sep = "_to_") %>%
  mutate(
    in_sample = from == to,
    from = factor(
      from,
      levels = c("nABA", "nWAS", "eABA", "eWAS", "all"),
      labels = c("Native\nAL", "Native\nWW", "Exotic\nAL", "Exotic\nWW", "All")
    ),
    to = factor(
      to,
      levels = c("nABA", "nWAS", "eABA", "eWAS", "all"),
      labels = c("Native\nAL", "Native\nWW", "Exotic\nAL", "Exotic\nWW", "All")
    )
  )

CRPS_summary <-
  lapply(cross_vals, function(x) median(x$CRPS, na.rm = TRUE)) %>%
  do.call(rbind, .) %>%
  as.data.frame %>%
  rownames_to_column("tmp") %>%
  separate(tmp, c("from", "to"), sep = "_to_") %>%
  mutate(
    in_sample = from == to,
    from = factor(
      from,
      levels = c("nABA", "nWAS", "eABA", "eWAS", "all"),
      labels = c(
        "Native\nABA",
        "Native\nWAS",
        "Exotic\nABA",
        "Exotic\nWAS",
        "All"
      )
    ),
    to = factor(
      to,
      levels = c("nABA", "nWAS", "eABA", "eWAS", "all"),
      labels = c(
        "Native\nABA",
        "Native\nWAS",
        "Exotic\nABA",
        "Exotic\nWAS",
        "All"
      )
    )
  )

SD0_CRPS_summary <-
  lapply(cross_vals, function(x) median(x$SD0_CRPS, na.rm = TRUE)) %>%
  do.call(rbind, .) %>%
  as.data.frame %>%
  rownames_to_column("tmp") %>%
  separate(tmp, c("from", "to"), sep = "_to_") %>%
  mutate(
    in_sample = from == to,
    from = factor(
      from,
      levels = c("nABA", "nWAS", "eABA", "eWAS", "all"),
      labels = c(
        "Native\nABA",
        "Native\nWAS",
        "Exotic\nABA",
        "Exotic\nWAS",
        "All"
      )
    ),
    to = factor(
      to,
      levels = c("nABA", "nWAS", "eABA", "eWAS", "all"),
      labels = c(
        "Native\nABA",
        "Native\nWAS",
        "Exotic\nABA",
        "Exotic\nWAS",
        "All"
      )
    )
  )

SD2_CRPS_summary <-
  lapply(cross_vals, function(x) median(x$SD2_CRPS, na.rm = TRUE)) %>%
  do.call(rbind, .) %>%
  as.data.frame %>%
  rownames_to_column("tmp") %>%
  separate(tmp, c("from", "to"), sep = "_to_") %>%
  mutate(
    in_sample = from == to,
    from = factor(
      from,
      levels = c("nABA", "nWAS", "eABA", "eWAS", "all"),
      labels = c(
        "Native\nABA",
        "Native\nWAS",
        "Exotic\nABA",
        "Exotic\nWAS",
        "All"
      )
    ),
    to = factor(
      to,
      levels = c("nABA", "nWAS", "eABA", "eWAS", "all"),
      labels = c(
        "Native\nABA",
        "Native\nWAS",
        "Exotic\nABA",
        "Exotic\nWAS",
        "All"
      )
    )
  )

# Plot --------------------------------------------------------------------

plot_transfer_matrix <- function(data, legend.title, col.direction) {
  # reference log-likelihood (in sample of the "all" model)
  midpoint <- data$`50%`[data$to == "All"]
  plot_data <-
    data %>%
    filter(from != "All")
  # plot
  p_out <-
    ggplot(plot_data) +
    geom_tile(
      aes(
        x = to,
        y = from,
        fill = `50%`,
        colour = in_sample
      ),
      linewidth = 1,
      width = 0.95,
      height = 0.95
    ) +
    geom_text(aes(x = to, y = from, label = sprintf("%.1f", round(`50%`, 1)))) +
    scale_x_discrete(position = "top") +
    scale_y_discrete(limits = rev) +
    scale_colour_manual(values = c(NA, "black"), na.value = NA) +
    scale_fill_viridis_c(
      option = "turbo",
      direction = col.direction,
      begin = 0.1,
      end = 0.9,
      transform = scales::transform_modulus(p = 0.5),
      name = legend.title
    ) +
    coord_equal() +
    guides(colour = "none") +
    labs(x = "To", y = "Transfer from") +
    theme_minimal() +
    theme(
      axis.ticks = element_blank(),
      axis.text = element_text(size = 10, colour = "black"),
      axis.title = element_text(size = 12),
      plot.margin = margin(0.1, 0.1, 0.1, 0.1),
      legend.margin = margin(0, 0, 0, 0),
      legend.position = "bottom",
      panel.grid = element_blank()
    )
  p_out
}

png("fig/transfer_matrix.png", width = 4, height = 4, units = "in", res = 300)
plot_transfer_matrix(
  LL_summary,
  "Prediction accuracy\n(log likelihood ratio)",
  col.direction = -1
)
dev.off()
