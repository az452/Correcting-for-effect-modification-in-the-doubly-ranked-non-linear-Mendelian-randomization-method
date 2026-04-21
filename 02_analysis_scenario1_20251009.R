
#===============================================================================
#
# Scenario 1 data analysis: net GxE interaction on GxE-induced bias
#
#===============================================================================

rm(list = ls())
# Set your working directory to the project folder before running this script.

library(dplyr)
library(SUMnlmr)
source("helpers/create_nlmr_summary2025.R")

#===============================================================================
# Load data
#===============================================================================

load("rdata_nlmr_gxe_sim_data_modifier_20251009.RData")

scenario_data <- list(
  sim_results_rho0  = sim_rho0,
  sim_results_rho03 = sim_rho03,
  sim_results_rho06 = sim_rho06
)

scenario_rho_labels <- c(
  sim_results_rho0  = "0.0",
  sim_results_rho03 = "0.3",
  sim_results_rho06 = "0.6"
)

#===============================================================================
# Analyse data
#===============================================================================

for (scenario_name in names(scenario_data)) {
  
  sim_list <- scenario_data[[scenario_name]]
  rho_label <- scenario_rho_labels[[scenario_name]]
  R <- length(sim_list)
  
  scenario_results <- vector("list", R)
  
  for (i in 1:R) {
    
    dat <- sim_list[[i]]
    
    cat("Scenario:", rho_label, "| Replicate:", i, "of", R, "\n")
    
    # Linear MR
    sx <- summary(lm(X ~ G, data = dat))
    bx   <- sx$coefficients["G", "Estimate"]
    bxse <- sx$coefficients["G", "Std. Error"]
    
    sy <- summary(lm(Y ~ G, data = dat))
    by   <- sy$coefficients["G", "Estimate"]
    byse <- sy$coefficients["G", "Std. Error"]
    
    m0 <- data.frame(
      bx    = bx,
      bxse  = bxse,
      by    = by,
      byse  = byse,
      xmean = mean(dat$X),
      xmin  = min(dat$X),
      xmax  = max(dat$X)
    ) %>%
      mutate(
        rep     = i,
        lace    = by / bx,
        lace_se = byse / bx,
        method  = "lmr",
        strata  = 0
      )
    
    # Doubly-ranked NLMR using X
    m1 <- create_nlmr_summary2025(
      y = dat$Y,
      x = dat$X,
      g = dat$G,
      covar = NULL,
      family = "gaussian",
      strata_method = "ranked",
      q = 10,
      seed = 123
    )$summary %>%
      mutate(
        rep     = i,
        lace    = by / bx,
        lace_se = byse / bx,
        method  = "rank+X",
        strata  = seq_len(10)
      )
    
    scenario_results[[i]] <- bind_rows(m0, m1)
  }
  
  assign(scenario_name, bind_rows(scenario_results))
}

#===============================================================================
# Save results
#===============================================================================

save(
  sim_results_rho0,
  sim_results_rho03,
  sim_results_rho06,
  file = "rdata_nlmr_gxe_sim_results_modifier_20251009.RData"
)
