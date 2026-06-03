# Script to fit greta model to training dataset
# origin: 1 (native) or -1 (exotic)
# landuse: 1 (ABA) or -1 (WAS)
# family: "negative_binomial", "zero_inflated_negative_binomial" or "dirichlet_multinomial"

fit_greta <- function(origin, landuse, family, ...) {
    # subset data
    species <- rownames(Z)[Z[, "Exotic1"] %in% origin]
    plots <- rownames(X)[X[, "Type1"] %in% landuse]

    Y_dat <- Y[plots, species]
    # cleanup
    Y_dat <- Y_dat[rowSums(Y_dat > 0) > 0, colSums(Y_dat > 0) > 0]
    patch <- as.numeric(as.factor(Pi[rownames(Y_dat), "patch"]))

    X_var <- c("size", "dist", "nitrogen", "phosphorous", "potassium")
    X_dat <- X[rownames(Y_dat), X_var]

    Z_var <- c("LPC", "WD", "logSM", "Hmax")
    Z_dat <- Z[colnames(Y_dat), Z_var]

    n_site <- nrow(Y_dat) # number of sites
    n_patch <- length(unique(patch)) # number of patches
    n_sp <- ncol(Y_dat) # number of species
    n_cov <- ncol(X_dat) # number of env. covariates
    n_lv <- 2 # number of LVs
    n_trait <- ncol(Z_dat) # number of traits

    trial_size <- n_trials[rownames(Y_dat)]

    # species parameters
    # sub-model parameters have normal prior distributions, including intercepts
    z_kappa <- normal(0, 1, dim = c(n_cov + 1, n_trait + 1))
    sd_kappa <- 10
    kappa <- sd_kappa * z_kappa

    # parameters of the base model are a function of the parameters of the
    # sub-model
    sd_beta_unif <- uniform(0, pi / 2, dim = n_cov + 1)
    sd_beta <- 2.5 * tan(sd_beta_unif)

    # prior on the correlation between the varying coefficient
    Omega <- lkj_correlation(3, n_cov + 1)
    Omega_U <- chol(Omega)
    Sigma_U <- sweep(Omega_U, 2, sd_beta, "*")

    mu_beta <- kappa %*% t(cbind(1, Z_dat))
    z_beta <- normal(0, 1, dim = c(n_cov + 1, n_sp))
    beta <- t(mu_beta) + t(z_beta) %*% Sigma_U

    # loadings
    loadings <- zeros(n_sp, n_lv)
    diag(loadings) <- normal(0, 1, dim = n_lv, truncation = c(0, Inf))
    lower_idx <- lower.tri(loadings)
    loadings[lower_idx] <- normal(0, 1, dim = sum(lower_idx))

    # latent variables
    lvs <- normal(0, 1, dim = c(n_site, n_lv))

    # linear predictor
    fixed <- cbind(1, X_dat) %*% t(beta)
    random <- lvs %*% t(loadings)
    linpred <- fixed + random

    if (family == "negative_binomial") {
        log_phi <- normal(0, 2.5)
        size <- exp(-log_phi)
        log_eta_size <- linpred + log_phi
        prob <- ilogit(-log_eta_size)

        Y_data <- as_data(Y_dat)
        distribution(Y_data) <- negative_binomial(size, prob)

        mod <- model(
            alpha,
            beta,
            alpha_bar,
            alphap,
            alphap_bar,
            kappa,
            lvs,
            loadings,
            log_phi,
            sd_alpha,
            sd_alphap,
            sd_beta,
            Sigma_U
        )
    } else if (family == "zero_inflated_negative_binomial") {
        log_phi <- normal(0, 2.5)
        size <- exp(-log_phi)
        log_eta_size <- linpred + log_phi
        prob <- ilogit(-log_eta_size)

        # probability of excess zeroes
        zi_a <- normal(0, 2.5)
        zi_b <- normal(0, 2.5, truncation = c(0, Inf))
        logit_zi <- zi_a - zi_b * linpred
        zi_prob <- ilogit(logit_zi)

        Y_data <- as_data(Y_dat)
        distribution(Y_data) <-
            zero_inflated_negative_binomial(size, prob, zi_prob)

        mod <- model(
            alpha,
            beta,
            alpha_bar,
            kappa,
            zi_a,
            zi_b,
            log_phi,
            sd_alpha,
            sd_alphap,
            sd_beta,
            Sigma_U
        )
    } else if (family == "dirichlet_multinomial") {
        mu <- imultilogit(linpred)
        theta <- exponential(10)
        conc <- mu * theta

        Y_dat_all <- cbind(Y_dat, Others = trial_size - rowSums(Y_dat))
        Y_dat_all <- as_data(Y_dat_all)

        distribution(Y_dat_all) <- dirichlet_multinomial(
            cbind(trial_size),
            conc
        )

        mod <- model(
            beta,
            kappa,
            lvs,
            loadings,
            theta,
            sd_beta,
            z_beta,
            Sigma_U
        )
    }

    # sampling
    draws <- mcmc(mod, ...)

    return(list(mod = mod, draws = draws, y = Y_dat, trial_size = trial_size))
}
