

#================================================================================================
#
# Scenario 3 DGM: Effect modifier is a confounder, correlated with a mediator and a collider
#
#================================================================================================


rm(list = ls())
# Set your working directory to the project folder before running this script.

library(dplyr)
library(SUMnlmr)

source("helpers/create_nlmr_summary2025.R")
source("helpers/sim_methods.R")

#===============================================================================
# Load data
#===============================================================================

load("rdata_nlmr_gxe_sim_s3_addBias_20251016.RData")

R <- 100

#===============================================================================
# Effect modifier is a confounder
#===============================================================================

sim_results <- data.frame()

for (i in seq_len(R)) {
  
  sim_data <- sim_data_conf[[i]]
  results <- sim_methods(sim_data, rep_num = i)
  sim_results <- rbind(sim_results, results)
  
  cat("replication ", i, " of ", R, " is done\n", sep = "")
}

sim_conf <- sim_results

#===============================================================================
# Effect modifier is related to the mediator
#===============================================================================

sim_results <- data.frame()

for (i in seq_len(R)) {
  
  sim_data <- sim_data_med[[i]]
  results <- sim_methods(sim_data, rep_num = i)
  sim_results <- rbind(sim_results, results)
  
  cat("replication ", i, " of ", R, " is done\n", sep = "")
}

sim_med <- sim_results

#===============================================================================
# Effect modifier is related to a collider
#===============================================================================

sim_results <- data.frame()

for (i in seq_len(R)) {
  
  sim_data <- sim_data_col[[i]]
  results <- sim_methods(sim_data, rep_num = i)
  sim_results <- rbind(sim_results, results)
  
  cat("replication ", i, " of ", R, " is done\n", sep = "")
}

sim_col <- sim_results

#===============================================================================
# Save results
#===============================================================================

save(
  sim_conf,
  sim_med,
  sim_col,
  file = "rdata_nlmr_gxe_sim_s3_addBias_results_20251016.RData"
)
