

#=======================================================================================
#
# Scenario 2 DGM: GxE correction under null, linear, and nonlinear X-Y relationship
#
#=======================================================================================

rm(list = ls())
# Set your working directory to the project folder before running this script.
library(faux)


#===============================================================================
# Settings
#===============================================================================

n <- 100000
R <- 100
set.seed(123)

mu_vec <- rep(2, 5)
sd_vec <- rep(1, 5)

#===============================================================================
# Null model
# X = 0.3G + 0.3E + 0.1(G*E) + U + ex
# Y = U + ey
#===============================================================================

all_sim_data <- vector("list", R)

for (i in seq_len(R)) {
  
  dat <- rnorm_multi(
    n        = n,
    mu       = mu_vec,
    sd       = sd_vec,
    r        = 0,
    varnames = c("G", "E", "U", "ex", "ey")
  )
  
  G  <- dat$G
  E  <- dat$E
  U  <- dat$U
  ex <- dat$ex
  ey <- dat$ey
  
  X <- 0.3 * G + 0.3 * E + 0.1 * (G * E) + U + ex
  Y <- U + ey
  
  sim_data <- data.frame(G = G, E = E, U = U, X = X, Y = Y)
  
  fit <- lm(X ~ G * E, data = sim_data)
  b_GE <- coef(fit)[["G:E"]]
  
  sim_data$X1 <- sim_data$X - b_GE * (sim_data$G * sim_data$E)
  
  all_sim_data[[i]] <- sim_data

}

sim_data_null <- all_sim_data

#===============================================================================
# Linear model
# X = 0.3G + 0.3E + 0.1(G*E) + U + ex
# Y = 0.5X + U + ey
#===============================================================================

all_sim_data <- vector("list", R)

for (i in seq_len(R)) {
  
  dat <- rnorm_multi(
    n        = n,
    mu       = mu_vec,
    sd       = sd_vec,
    r        = 0,
    varnames = c("G", "E", "U", "ex", "ey")
  )
  
  G  <- dat$G
  E  <- dat$E
  U  <- dat$U
  ex <- dat$ex
  ey <- dat$ey
  
  X <- 0.3 * G + 0.3 * E + 0.1 * (G * E) + U + ex
  Y <- 0.5 * X + U + ey
  
  sim_data <- data.frame(G = G, E = E, U = U, X = X, Y = Y)
  
  fit <- lm(X ~ G * E, data = sim_data)
  b_GE <- coef(fit)[["G:E"]]
  
  sim_data$X1 <- sim_data$X - b_GE * (sim_data$G * sim_data$E)
  
  all_sim_data[[i]] <- sim_data

}

sim_data_linear <- all_sim_data

#===============================================================================
# Nonlinear model
# X = 0.3G + 0.3E + 0.1(G*E) + U + ex
# Y = 0.1X^2 + U + ey
#===============================================================================

all_sim_data <- vector("list", R)

for (i in seq_len(R)) {
  
  dat <- rnorm_multi(
    n        = n,
    mu       = mu_vec,
    sd       = sd_vec,
    r        = 0,
    varnames = c("G", "E", "U", "ex", "ey")
  )
  
  G  <- dat$G
  E  <- dat$E
  U  <- dat$U
  ex <- dat$ex
  ey <- dat$ey
  
  X <- 0.3 * G + 0.3 * E + 0.1 * (G * E) + U + ex
  Y <- 0.1 * X^2 + U + ey
  
  sim_data <- data.frame(G = G, E = E, U = U, X = X, Y = Y)
  
  fit <- lm(X ~ G * E, data = sim_data)
  b_GE <- coef(fit)[["G:E"]]
  
  sim_data$X1 <- sim_data$X - b_GE * (sim_data$G * sim_data$E)
  
  all_sim_data[[i]] <- sim_data
}

sim_data_nonlinear <- all_sim_data

#===============================================================================
# Save data
#===============================================================================

save(
  sim_data_null,
  sim_data_linear,
  sim_data_nonlinear,
  file = "rdata_nlmr_gxe_sim_s2_strategy_20251016.RData"
)
