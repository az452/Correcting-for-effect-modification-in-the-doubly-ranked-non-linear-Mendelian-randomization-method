
#===============================================================================
#
# Scenario 1 DGM: net GxE interaction on GxE-induced bias
#
#===============================================================================

rm(list = ls())
# Set your working directory to the project folder before running this script.
library(faux)

#===============================================================================
# Settings
#===============================================================================

n <- 100000
R <- 100

set.seed(123)

scenario_rhos <- c(
  sim_rho0  = 0,
  sim_rho03 = 0.3,
  sim_rho06 = 0.6
)

#===============================================================================
# Simulate data
#===============================================================================

for (scenario_name in names(scenario_rhos)) {
  
  rho <- scenario_rhos[[scenario_name]]
  sim_list <- vector("list", R)
  
  for (i in 1:R) {
    
    cat("Scenario:", scenario_name, "| Replicate:", i, "of", R, "\n")
    
    # Generate variables
    G  <- rnorm(n, mean = 2, sd = 1)
    U  <- rnorm(n, mean = 2, sd = 1)
    ex <- rnorm(n, mean = 2, sd = 1)
    ey <- rnorm(n, mean = 2, sd = 1)
    
    ee <- rnorm_multi(
      n = n,
      mu = c(2, 2),
      sd = c(1, 1),
      r = rho,
      varnames = c("E1", "E2")
    )
    
    E1 <- ee$E1
    E2 <- ee$E2
    
    # Generate exposure and outcome
    X <- 0.3 * G + 0.3 * E1 + 0.2 * E2 + 0.1 * (G * E1) - 0.1 * (G * E2) + U + ex
    Y <- U + ey
    
    # Store data
    df <- data.frame(G, U, E1, E2, ex, ey, X, Y)
    
    # GxE-corrected exposure
    fit <- lm(X ~ G + E1 + E2 + G:E1 + G:E2, data = df)
    b_GE1 <- coef(fit)["G:E1"]
    b_GE2 <- coef(fit)["G:E2"]
    
    df$X1 <- df$X - b_GE1 * df$G * df$E1 - b_GE2 * df$G * df$E2
    
    sim_list[[i]] <- df
  }
  
  assign(scenario_name, sim_list)
}

#===============================================================================
# Save data
#===============================================================================

save(
  sim_rho0,
  sim_rho03,
  sim_rho06,
  file = "rdata_nlmr_gxe_sim_data_modifier_20251009.RData"
)
