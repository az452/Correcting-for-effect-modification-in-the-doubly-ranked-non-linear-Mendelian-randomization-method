
#===============================================================================
#
# simulations in the supplementary materials: DGM
#
#===============================================================================

rm(list = ls())
# Set your working directory to the project folder before running this script.
library(faux)



#===============================================================================
#
# Scenario 2 supplementary condition: threshold effect
#
#===============================================================================

n <- 100000
R <- 100
set.seed(123)

all_sim_data <- vector("list", R)

for (i in seq_len(R)) {
  
  dat <- rnorm_multi(
    n = n,
    mu = c(0, 0, 0, 0, 0),
    sd = c(0.5, 1, 1, 1, 1),
    r = 0,
    varnames = c("G", "E", "U", "ex", "ey")
  )
  
  G  <- dat$G
  E  <- dat$E
  U  <- dat$U
  ex <- dat$ex
  ey <- dat$ey
  
  X <- 0.5 * G + 0.5 * E + 0.1 * (G * E) + U + ex
  
  Y <- ifelse(
    X <= 0,
    U + ey,
    -0.1 * X^2 + U + ey
  )
  
  sim_data <- data.frame(
    G   = G,
    E   = E,
    U   = U,
    X   = X,
    Y   = Y,
    rep = i
  )
  
  fit <- lm(X ~ G * E, data = sim_data)
  b_GE <- coef(fit)[["G:E"]]
  
  sim_data$X1 <- sim_data$X - b_GE * (sim_data$G * sim_data$E)
  
  all_sim_data[[i]] <- sim_data
}

sim_data_thresholdUpdate <- all_sim_data

save(
  sim_data_thresholdUpdate,
  file = "rdata_nlmr_gxe_sim_s4_thresholdUpdate_20260322.RData"
)



#===============================================================================
#
# Scenario 3 supplementary condition: effect modifier is also a collider
#
#===============================================================================

n <- 100000
R <- 100
set.seed(123)

all_sim_data <- vector("list", R)

for (i in seq_len(R)) {
  
  dat <- rnorm_multi(
    n = n,
    mu = rep(2, 7),
    sd = rep(1, 7),
    r = 0,
    varnames = c("G", "U", "U1", "U2", "ex", "ey", "eE")
  )
  
  G  <- dat$G
  U  <- dat$U
  U1 <- dat$U1
  U2 <- dat$U2
  ex <- dat$ex
  ey <- dat$ey
  eE <- dat$eE
  
  E <- 0.1 * U1 + 0.1 * U2 + eE
  X <- 0.3 * G + 0.3 * E + 0.3 * (G * E) + U + 0.3 * U1 + ex
  Y <- U + 0.3 * U2 + ey
  
  sim_data <- data.frame(
    G = G,
    U = U,
    U1 = U1,
    U2 = U2,
    E = E,
    X = X,
    Y = Y,
    rep = i
  )
  
  fit <- lm(X ~ G * E, data = sim_data)
  b_GE <- coef(fit)[["G:E"]]
  
  sim_data$X1 <- sim_data$X - b_GE * (sim_data$G * sim_data$E)
  
  all_sim_data[[i]] <- sim_data
}

sim_data_colliderUpdate <- all_sim_data

save(
  sim_data_colliderUpdate,
  file = "rdata_nlmr_gxe_sim_s4_colliderUpdate_20260330.RData"
)

#===============================================================================
#
# Supplementary scenario 1: imperfectly measured effect modifiers
#
#===============================================================================

n <- 100000
R <- 100
set.seed(123)

rho1 <- 1
rho2 <- 0.8
rho3 <- 0.4
rho4 <- 0

sim_data_imperfectE <- vector("list", R)

for (i in seq_len(R)) {
  
  dat <- rnorm_multi(
    n = n,
    mu = rep(2, 5),
    sd = rep(1, 5),
    r = 0,
    varnames = c("G", "E", "U", "ex", "ey")
  )
  
  G  <- dat$G
  E  <- dat$E
  U  <- dat$U
  ex <- dat$ex
  ey <- dat$ey
  
  X <- 0.3 * G + 0.3 * E + 0.1 * (G * E) + U + ex
  Y <- U + ey
  
  eta1 <- rnorm(n, mean = 0, sd = 1)
  eta2 <- rnorm(n, mean = 0, sd = 1)
  eta3 <- rnorm(n, mean = 0, sd = 1)
  eta4 <- rnorm(n, mean = 0, sd = 1)
  
  E1 <- 2 + rho1 * (E - 2) + sqrt(1 - rho1^2) * eta1
  E2 <- 2 + rho2 * (E - 2) + sqrt(1 - rho2^2) * eta2
  E3 <- 2 + rho3 * (E - 2) + sqrt(1 - rho3^2) * eta3
  E4 <- 2 + rho4 * (E - 2) + sqrt(1 - rho4^2) * eta4
  
  fit1 <- lm(X ~ G * E1)
  fit2 <- lm(X ~ G * E2)
  fit3 <- lm(X ~ G * E3)
  fit4 <- lm(X ~ G * E4)
  
  b_GE1 <- coef(fit1)[["G:E1"]]
  b_GE2 <- coef(fit2)[["G:E2"]]
  b_GE3 <- coef(fit3)[["G:E3"]]
  b_GE4 <- coef(fit4)[["G:E4"]]
  
  X1 <- X - b_GE1 * (G * E1)
  X2 <- X - b_GE2 * (G * E2)
  X3 <- X - b_GE3 * (G * E3)
  X4 <- X - b_GE4 * (G * E4)
  
  sim_data <- data.frame(
    G    = G,
    E    = E,
    U    = U,
    ex   = ex,
    ey   = ey,
    X    = X,
    Y    = Y,
    E1   = E1,
    E2   = E2,
    E3   = E3,
    E4   = E4,
    X1   = X1,
    X2   = X2,
    X3   = X3,
    X4   = X4,
    rho1 = rho1,
    rho2 = rho2,
    rho3 = rho3,
    rho4 = rho4,
    rep  = i
  )
  
  sim_data_imperfectE[[i]] <- sim_data
  
  cat("replication ", i, " of ", R, " is done\n", sep = "")
}

save(
  sim_data_imperfectE,
  file = "rdata_nlmr_gxe_sim_s4_imperfectE_20260320.RData"
)

#===============================================================================
#
# Supplementary scenario 2: log-transformed modifier
#
#===============================================================================

n <- 100000
R <- 100
set.seed(123)

all_sim_data <- vector("list", R)

for (i in seq_len(R)) {
  
  dat <- rnorm_multi(
    n = n,
    mu = rep(2, 5),
    sd = rep(1, 5),
    r = 0,
    varnames = c("G", "E", "U", "ex", "ey")
  )
  
  G  <- dat$G
  E  <- dat$E
  U  <- dat$U
  ex <- dat$ex
  ey <- dat$ey
  
  E1 <- log10(E + 8)
  
  X <- 0.3 * G + 0.3 * E + 0.1 * (G * E) + U + ex
  Y <- U + ey
  
  sim_data <- data.frame(
    G   = G,
    E   = E,
    E1  = E1,
    U   = U,
    X   = X,
    Y   = Y,
    rep = i
  )
  
  fit <- lm(X ~ G * E1, data = sim_data)
  b_GE1 <- coef(fit)[["G:E1"]]
  
  sim_data$X1 <- sim_data$X - b_GE1 * (sim_data$G * sim_data$E1)
  
  all_sim_data[[i]] <- sim_data
}

sim_data_logE <- all_sim_data

save(
  sim_data_logE,
  file = "rdata_nlmr_gxe_sim_s4_logE_20260320.RData"
)

#===============================================================================
#
# Supplementary scenario 2: exp-transformed modifier
#
#===============================================================================

n <- 100000
R <- 100
set.seed(123)

all_sim_data <- vector("list", R)

for (i in seq_len(R)) {
  
  dat <- rnorm_multi(
    n = n,
    mu = rep(2, 5),
    sd = rep(1, 5),
    r = 0,
    varnames = c("G", "E", "U", "ex", "ey")
  )
  
  G  <- dat$G
  E  <- dat$E
  U  <- dat$U
  ex <- dat$ex
  ey <- dat$ey
  
  E1 <- exp(E)
  
  X <- 0.3 * G + 0.3 * E + 0.1 * (G * E) + U + ex
  Y <- U + ey
  
  sim_data <- data.frame(
    G   = G,
    E   = E,
    E1  = E1,
    U   = U,
    X   = X,
    Y   = Y,
    rep = i
  )
  
  fit <- lm(X ~ G * E1, data = sim_data)
  b_GE1 <- coef(fit)[["G:E1"]]
  
  sim_data$X1 <- sim_data$X - b_GE1 * (sim_data$G * sim_data$E1)
  
  all_sim_data[[i]] <- sim_data
}

sim_data_expE <- all_sim_data

save(
  sim_data_expE,
  file = "rdata_nlmr_gxe_sim_s4_expE_20260325.RData"
)
