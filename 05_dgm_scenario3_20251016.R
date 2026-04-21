
#================================================================================================
#
# Scenario 3 DGM: Effect modifier is a confounder, correlated with a mediator and a collider
#
#================================================================================================

rm(list = ls())
# Set your working directory to the project folder before running this script.
library(faux)

#===============================================================================
# Settings
#===============================================================================

n <- 100000
R <- 100

mu4 <- rep(2, 4)
sd4 <- rep(1, 4)

mu6 <- rep(2, 6)
sd6 <- rep(1, 6)

#===============================================================================
# Effect modifier is also a confounder
# X = 0.3G + 0.3E + 0.1(G*E) + ex
# Y = E + ey
#===============================================================================

set.seed(123)

all_sim_data <- vector("list", R)

for (i in seq_len(R)) {
  
  dat <- rnorm_multi(
    n        = n,
    mu       = mu4,
    sd       = sd4,
    r        = 0,
    varnames = c("G", "E", "ex", "ey")
  )
  
  G  <- dat$G
  E  <- dat$E
  ex <- dat$ex
  ey <- dat$ey
  
  X <- 0.3 * G + 0.3 * E + 0.1 * (G * E) + ex
  Y <- E + ey
  
  sim_data <- data.frame(G = G, E = E, X = X, Y = Y)
  
  fit <- lm(X ~ G * E, data = sim_data)
  b_GE <- coef(fit)[["G:E"]]
  
  sim_data$X1 <- sim_data$X - b_GE * (sim_data$G * sim_data$E)
  
  all_sim_data[[i]] <- sim_data
}

sim_data_conf <- all_sim_data

#===============================================================================
# Effect modifier is related to the mediator
# X = 0.3G + 0.3E + 0.1(G*E) + U + ex
# M = 0.5X + E + em
# Y = M + U + ey
#===============================================================================

set.seed(123)

all_sim_data <- vector("list", R)

for (i in seq_len(R)) {
  
  sim_data <- rnorm_multi(
    n        = n,
    mu       = mu6,
    sd       = sd6,
    r        = 0,
    varnames = c("G", "E", "U", "ex", "em", "ey")
  )
  
  sim_data$X <- 0.3 * sim_data$G + 0.3 * sim_data$E + 0.1 * (sim_data$G * sim_data$E) + sim_data$U + sim_data$ex
  sim_data$M <- 0.5 * sim_data$X + sim_data$E + sim_data$em
  sim_data$Y <- sim_data$M + sim_data$U + sim_data$ey
  
  fit <- lm(X ~ G * E, data = sim_data)
  b_GE <- coef(fit)[["G:E"]]
  
  sim_data$X1 <- sim_data$X - b_GE * (sim_data$G * sim_data$E)
  
  all_sim_data[[i]] <- sim_data
}

sim_data_med <- all_sim_data

#===============================================================================
# Effect modifier is related to a collider
# X = 0.3G + 0.3E + 0.1(G*E) + U + ex
# Y = U + ey
# C = 0.5X + E + 0.5Y + ec
#===============================================================================

set.seed(123)

all_sim_data <- vector("list", R)

for (i in seq_len(R)) {
  
  sim_data <- rnorm_multi(
    n        = n,
    mu       = mu6,
    sd       = sd6,
    r        = 0,
    varnames = c("G", "E", "U", "ex", "ec", "ey")
  )
  
  sim_data$X <- 0.3 * sim_data$G + 0.3 * sim_data$E + 0.1 * (sim_data$G * sim_data$E) + sim_data$U + sim_data$ex
  sim_data$Y <- sim_data$U + sim_data$ey
  sim_data$C <- 0.5 * sim_data$X + sim_data$E + 0.5 * sim_data$Y + sim_data$ec
  
  fit <- lm(X ~ G * E, data = sim_data)
  b_GE <- coef(fit)[["G:E"]]
  
  sim_data$X1 <- sim_data$X - b_GE * (sim_data$G * sim_data$E)
  
  all_sim_data[[i]] <- sim_data
}

sim_data_col <- all_sim_data

#===============================================================================
# Save data
#===============================================================================

save(
  sim_data_conf,
  sim_data_med,
  sim_data_col,
  file = "rdata_nlmr_gxe_sim_s3_addBias_20251016.RData"
)
