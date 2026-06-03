library(novelforestSG)
library(here)
library(ggpubr)


# Read species and plot list from Lai et al. 2021 AVS
spp_used <- read.csv(here("Data/Clean/spp.used.csv"))
plot_used <- read.csv(here("Data/Clean/plot.used.csv"))

# Species list with origin status
SpList <-
  read.csv(here("Data/Raw/SL_Use.csv"), header = T) %>%
  # update species origin to Lindsay et al. (2022)
  mutate(
    Exotic = case_when(
      Species == "Mangifera odorata" ~ 1,
      .default = Exotic
    )
  )

# X matrix (environmental data) -------------------------------------------

X.df <-
  read.csv(here("Data/Raw/envi.csv"), header = T) %>%
  # summarise to plot level
  group_by(Type, patch, plot) %>%
  summarise_at(vars(size, dist, nitrogen:potassium), ~ mean(log(.))) %>%
  ungroup %>%
  as.data.frame

X.df.s <- X.df %>%
  mutate_at(
    vars(size, dist, nitrogen, phosphorous, potassium),
    ~ as.numeric(scale(.))
  )
X <- model.matrix(
  ~ 1 + Type * (size + dist + nitrogen + phosphorous + potassium),
  contrasts.arg = list(Type = "contr.sum"),
  data = X.df.s
)[, -1]
rownames(X) <- X.df[, 3]


# Y matrix (community data) -----------------------------------------------

Y.df <-
  # read.csv(here("Data/Raw/AlexYee_dbh.csv")) %>%
  novelforest_data$trees %>%
  mutate(species = word(species, 1, 2)) %>%
  mutate(dbh_2011 = as.numeric(dbh_2011)) %>%
  filter(plot != "IK1N", stem == 1, dbh_2011 >= 5) %>%
  rename(Plot = plot, Species = species) %>%
  group_by(Plot, Species) %>%
  summarise(Count = n())

n_trials <- with(Y.df, tapply(Count, Plot, sum))


# Traits ------------------------------------------------------------------
trait <-
  read.csv(here("Data/Raw/trait.csv")) %>%
  rename(Species = X) %>%
  filter(Species %in% unique(Y.df$Species)) %>%
  select(Species, logSLA, Th, LDMC, WD, logSM, Hmax) %>%
  left_join(select(SpList, Species, Exotic)) %>%
  # convert Exotic status to sum contrast
  mutate(Exotic = ifelse(Exotic == 0, -1, 1))

trait %>% select(Species, Exotic) %>% group_by(Exotic) %>% count

# obtain species with complete trait data
trait <- trait[complete.cases(trait), ]

# reduce leaf trait dimension
# Mandai trait data and build PCA model
leaf_trait <- trait[, c("logSLA", "LDMC", "Th")]
leaf_pca <- prcomp(leaf_trait, scale. = TRUE)
lpc1 <- predict(leaf_pca)[, 1]

# compile traits
Z.df <- as.data.frame(cbind(Species = trait[, 1], LPC = lpc1, trait[, 5:8]))

# scale traits and model matrix
Z.df.s <- Z.df %>%
  mutate_at(vars(LPC, WD, logSM, Hmax), ~ as.numeric(scale(.))) %>%
  mutate(Exotic = as.factor(Exotic))
Z <- model.matrix(
  ~ 1 + Exotic * (LPC + WD + logSM + Hmax),
  contrasts.arg = list(Exotic = "contr.sum"),
  data = Z.df.s
)[, -1]
rownames(Z) <- Z.df.s$Species

# remove species without trait in community data
Y.df <-
  Y.df %>%
  filter(Species %in% rownames(Z))

Y <- with(Y.df, tapply(Count, list(Plot, Species), sum, default = 0))

# then filter the trait matrix to only use species in Y
Z <- Z[colnames(Y), ]


# Random (nested) effects -------------------------------------------------

# plots within patch
Pi <-
  select(X.df, plot, patch) %>%
  filter(plot %in% rownames(Y)) %>%
  droplevels
rownames(Pi) <- Pi$plot
Pi$patch <- as.numeric(as.factor(Pi$patch))
Pi$plot <- as.numeric(as.factor(Pi$plot))
Pi <- as.matrix(Pi)

n_trials <- n_trials[rownames(Y)]


# Consolidate all into a list object --------------------------------------
## make sure X and Y have the same row number
X <- X[match(rownames(Y), rownames(X)), ]
identical(nrow(X), nrow(Y))
identical(rownames(Z), colnames(Y))
which(rowSums(Y) == 0)

Pi <- Pi[match(rownames(Y), rownames(Pi)), , drop = FALSE]

comb.dat <- list(X = X, Y = Y, Z = Z, Pi = Pi, n_trials = n_trials)

save(comb.dat, file = here("Data/Clean/comb.dat.Rdata"))


# Supplementary table -----------------------------------------------------

# Species information
tab_SI <-
  Z %>%
  as.data.frame() %>%
  rownames_to_column("Species") %>%
  select(Species, Exotic1) %>%
  mutate(
    Origin = ifelse(Exotic1 == 1, "Native", "Exotic"),
    .keep = "unused"
  ) %>%
  # join abundance by land use type
  left_join(
    reshape2::melt(Y) %>%
      rename(plot = Var1, Species = Var2) %>%
      left_join(
        X.df %>%
          mutate(Type = ifelse(Type == "ABA", "AL", "WW")) %>%
          select(Type, plot)
      ) %>%
      group_by(Species, Type) %>%
      summarise(N = sum(value)) %>%
      ungroup %>%
      pivot_wider(names_from = Type, values_from = N)
  )

write_csv(tab_SI, "out/tab_SI.csv")


# Trait and environment overlaps

p_freq_env <-
  read.csv(here("Data/Raw/envi.csv"), header = T) %>%
  group_by(Type, patch, plot) %>%
  summarise_at(vars(size, dist, nitrogen:potassium), ~ exp(mean(log(.)))) %>%
  ungroup %>%
  filter(plot %in% rownames(X)) %>%
  mutate(dist = dist / 1000) %>%
  select(Type, nitrogen, phosphorous, potassium, size, dist) %>%
  pivot_longer(cols = -Type, names_to = "env") %>%
  mutate(
    env = case_match(
      env,
      "dist" ~ "Dist. from old-growth (km)",
      "size" ~ "Patch size (sq. km)",
      "nitrogen" ~ "Nitrogen (mg/kg)",
      "phosphorous" ~ "Phosphorous (mg/kg)",
      "potassium" ~ "Potassium (mg/kg)"
    )
  ) %>%
  ggplot() +
  facet_wrap(~env, scales = "free_x", nrow = 1, strip.position = "bottom") +
  geom_freqpoly(aes(value, colour = Type)) +
  labs(y = "Frequency") +
  scale_colour_viridis_d(option = "turbo", begin = 0.2, end = 0.8) +
  scale_y_continuous(labels = scales::label_number(accuracy = 1)) +
  scale_x_log10() +
  theme_bw() +
  coord_cartesian(expand = FALSE) +
  theme(
    legend.position = "top",
    strip.placement = "outside",
    strip.background = element_blank(),
    axis.title.x = element_blank(),
    panel.grid = element_blank()
  )

p_freq_trait <-
  Z.df %>%
  filter(Species %in% rownames(Z)) %>%
  mutate(
    Hmax = Hmax / 100,
    Origin = ifelse(Exotic == 1, "Exotic", "Native")
  ) %>%
  select(Origin, LPC, WD, logSM, Hmax) %>%
  pivot_longer(cols = -Origin, names_to = "trait") %>%
  mutate(
    trait = case_match(
      trait,
      "LPC" ~ "Leaf PC1",
      "WD" ~ "Wood density (mg/mm^3)",
      "logSM" ~ "log Seed dry mass (log mg)",
      "Hmax" ~ "Max. height (m)"
    )
  ) %>%
  ggplot() +
  facet_wrap(~trait, scales = "free_x", nrow = 1, strip.position = "bottom") +
  geom_freqpoly(aes(value, colour = Origin)) +
  labs(y = "Frequency") +
  scale_colour_viridis_d(
    option = "viridis",
    begin = 0.2,
    end = 0.8,
    direction = -1
  ) +
  scale_y_continuous(labels = scales::label_number(accuracy = 1)) +
  theme_bw() +
  coord_cartesian(expand = FALSE) +
  theme(
    legend.position = "top",
    strip.placement = "outside",
    strip.background = element_blank(),
    axis.title.x = element_blank(),
    panel.grid = element_blank()
  )

png(
  here("fig/predictor_distribution.png"),
  width = 7,
  height = 5,
  units = "in",
  res = 300
)
ggarrange(p_freq_env, p_freq_trait, nrow = 2, labels = c("(a)", "(b)"))
dev.off()
