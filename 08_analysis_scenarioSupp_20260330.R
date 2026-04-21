

#===============================================================================
#
# Overview
#
# This script analyses supplementary simulation scenarios and saves the results.
#
# It includes:
#   1) data analysis for threshold-effect simulation
#   2) data analysis for collider-modifier simulation
#   3) data analysis for imperfectly measured modifiers
#   4) data analysis for log-transformed modifier
#   5) data analysis for exp-transformed modifier
#   6) compute target effects for the threshold-effect simulation
#
# Output files:
#   - rdata_nlmr_gxe_sim_s4_thresholdUpdate_results_20260322.RData
#   - rdata_nlmr_gxe_sim_s4_colliderUpdate_results_20260330.RData
#   - rdata_nlmr_gxe_sim_s4_imperfectE_results_20260320.RData
#   - rdata_nlmr_gxe_sim_s4_logE_results_20260320.RData
#   - rdata_nlmr_gxe_sim_s4_expE_results_20260325.RData
#   - rdata_nlmr_gxe_sim_s4_thresholdUpdate_target_20260322.RData
#
#===============================================================================



#===============================================================================
#
# simulations in the supplementary materials: data analysis
#
#===============================================================================

rm(list = ls())
# Set your working directory to the project folder before running this script.

library(dplyr)
library(SUMnlmr)

source("helpers/create_nlmr_summary2025.R")
source("helpers/sim_methods.R")



#===============================================================================
#
# Scenario 2 supplementary condition: threshold effect
#
#===============================================================================

load("rdata_nlmr_gxe_sim_s4_thresholdUpdate_20260322.RData")

R <- 100
sim_results_thresholdUpdate <- data.frame()

for (i in seq_len(R)) {
  
  sim_data <- sim_data_thresholdUpdate[[i]]
  results <- sim_methods(sim_data, rep_num = i)
  sim_results_thresholdUpdate <- rbind(sim_results_thresholdUpdate, results)
  
  cat("replication ", i, " of ", R, " is done\n", sep = "")
}

save(
  sim_results_thresholdUpdate,
  file = "rdata_nlmr_gxe_sim_s4_thresholdUpdate_results_20260322.RData"
)



#===============================================================================
#
# Scenario 3 supplementary condition: effect modifier is also a collider
#
#===============================================================================

load("rdata_nlmr_gxe_sim_s4_colliderUpdate_20260330.RData")

R <- 100
sim_results_colliderUpdate <- data.frame()

for (i in seq_len(R)) {
  
  sim_data <- sim_data_colliderUpdate[[i]]
  results <- sim_methods(sim_data, rep_num = i)
  sim_results_colliderUpdate <- rbind(sim_results_colliderUpdate, results)
  
  cat("replication ", i, " of ", R, " is done\n", sep = "")
}

save(
  sim_results_colliderUpdate,
  file = "rdata_nlmr_gxe_sim_s4_colliderUpdate_results_20260330.RData"
)



#===============================================================================
#
# Supplementary scenario 1: imperfectly measured effect modifiers
#
#===============================================================================

load("rdata_nlmr_gxe_sim_s4_imperfectE_20260320.RData")

R <- 100
sim_results_imperfectE <- data.frame()

for (i in seq_len(R)) {
  
  sim_data <- sim_data_imperfectE[[i]]
  results <- sim_methods_multiX(sim_data, rep_num = i)
  sim_results_imperfectE <- rbind(sim_results_imperfectE, results)
  
  cat("replication ", i, " of ", R, " is done\n", sep = "")
}

save(
  sim_results_imperfectE,
  file = "rdata_nlmr_gxe_sim_s4_imperfectE_results_20260320.RData"
)


#===============================================================================
#
# Supplementary scenario 2: log-transformed modifier
#
#===============================================================================

load("rdata_nlmr_gxe_sim_s4_logE_20260320.RData")

R <- 100
sim_results_logE <- data.frame()

for (i in seq_len(R)) {
  
  sim_data <- sim_data_logE[[i]]
  results <- sim_methods(sim_data, rep_num = i)
  sim_results_logE <- rbind(sim_results_logE, results)
  
  cat("replication ", i, " of ", R, " is done\n", sep = "")
}

save(
  sim_results_logE,
  file = "rdata_nlmr_gxe_sim_s4_logE_results_20260320.RData"
)


#===============================================================================
#
# Supplementary scenario 2: exp-transformed modifier
#
#===============================================================================

load("rdata_nlmr_gxe_sim_s4_expE_20260325.RData")

R <- 100
sim_results_expE <- data.frame()

for (i in seq_len(R)) {
  
  sim_data <- sim_data_expE[[i]]
  results <- sim_methods(sim_data, rep_num = i)
  sim_results_expE <- rbind(sim_results_expE, results)
  
  cat("replication ", i, " of ", R, " is done\n", sep = "")
}

save(
  sim_results_expE,
  file = "rdata_nlmr_gxe_sim_s4_expE_results_20260325.RData"
)



#===============================================================================
#
# compute target effects for threshold-effect simulation
#
#===============================================================================

rm(list = ls())
# Set your working directory to the project folder before running this script.

library(dplyr)
library(fda)

source("helpers/create_nlmr_summary2025.R")

#===============================================================================
# Load data
#===============================================================================

load("rdata_nlmr_gxe_sim_s4_thresholdUpdate_20260322.RData")

#===============================================================================
# Add ranked strata to each replicate
#   strata_rank0: ranked on X
#   strata_rank1: ranked on X1
#===============================================================================

sim_data_thresholdUpdate <- lapply(seq_along(sim_data_thresholdUpdate), function(i) {
  
  cat("Processing replicate", i, "of", length(sim_data_thresholdUpdate), "\n")
  
  dat <- sim_data_thresholdUpdate[[i]]
  
  dat$strata_rank0 <- create_nlmr_summary2025(
    y = dat$Y,
    x = dat$X,
    g = dat$G,
    covar = NULL,
    family = "gaussian",
    strata_method = "ranked",
    q = 10,
    seed = 123
  )$strata
  
  dat$strata_rank1 <- create_nlmr_summary2025(
    y = dat$Y,
    x = dat$X,
    xs = dat$X1,
    g = dat$G,
    covar = NULL,
    family = "gaussian",
    strata_method = "ranked",
    q = 10,
    seed = 123
  )$strata
  
  dat
})

#===============================================================================
# True derivative
#   0        if x <= 0
#   -0.2 * x if x > 0
#===============================================================================

nonlinear_dif <- function(x) {
  ifelse(x <= 0, 0, -0.2 * x)
}

#===============================================================================
# Simpson integration on equally spaced grid
#===============================================================================

simpson_equal <- function(x, y) {
  n <- length(x)
  
  if (n != length(y)) stop("x and y must have the same length")
  if (n < 2) return(NA_real_)
  
  h <- x[2] - x[1]
  
  if (n %% 2 == 0) {
    n1 <- n - 1
    
    simpson_part <- h / 3 * (
      y[1] + y[n1] +
        4 * sum(y[seq(2, n1 - 1, by = 2)]) +
        2 * sum(y[seq(3, n1 - 2, by = 2)])
    )
    
    trap_part <- (x[n] - x[n1]) * (y[n] + y[n1]) / 2
    return(simpson_part + trap_part)
  }
  
  h / 3 * (
    y[1] + y[n] +
      4 * sum(y[seq(2, n - 1, by = 2)]) +
      2 * sum(y[seq(3, n - 2, by = 2)])
  )
}

#===============================================================================
# Target estimate for one stratum
#===============================================================================

target <- function(selected_dat,
                   exp_col = "X",
                   iv_col = "G",
                   hprime_fun = nonlinear_dif,
                   grid_length = 1000,
                   nbasis = 54,
                   lambda = 0) {
  
  exp <- selected_dat[[exp_col]]
  iv  <- selected_dat[[iv_col]]
  
  ok <- is.finite(exp) & is.finite(iv)
  exp <- exp[ok]
  iv  <- iv[ok]
  
  if (length(exp) < 2) return(NA_real_)
  if (length(unique(iv)) < 2) return(NA_real_)
  
  var_iv <- var(iv)
  if (!is.finite(var_iv) || var_iv <= 0) return(NA_real_)
  
  range_exp <- c(min(exp), max(exp))
  if (!all(is.finite(range_exp)) || diff(range_exp) <= 0) return(NA_real_)
  
  At <- function(t) {
    cov(iv, as.numeric(exp >= t)) / var_iv
  }
  
  time_obs <- seq(range_exp[1], range_exp[2], length.out = grid_length)
  Aobs <- vapply(time_obs, At, numeric(1))
  
  bs_basis <- create.bspline.basis(
    rangeval = c(0, range_exp[2] - range_exp[1]),
    nbasis = nbasis
  )
  
  par_obj <- fdPar(fdobj = bs_basis, Lfdobj = NULL, lambda = lambda)
  
  smoothed_fd <- smooth.basis(
    argvals = time_obs - range_exp[1],
    y = Aobs,
    fdParobj = par_obj
  )$fd
  
  smoothed_At <- as.numeric(
    eval.fd(
      evalarg = seq(0, range_exp[2] - range_exp[1], length.out = grid_length),
      fdobj = smoothed_fd
    )
  )
  
  hdx <- hprime_fun(time_obs)
  intres <- simpson_equal(time_obs, smoothed_At * hdx)
  
  bx_ <- unname(coef(lm(exp ~ iv))[2])
  if (!is.finite(bx_) || abs(bx_) < 1e-12) return(NA_real_)
  
  intres / bx_
}

#===============================================================================
# Compute targets for all replicates and strata
#===============================================================================

target_results <- bind_rows(
  lapply(seq_along(sim_data_thresholdUpdate), function(i) {
    
    dat <- sim_data_thresholdUpdate[[i]]
    
    cat("Replicate", i, "of", length(sim_data_thresholdUpdate), "\n")
    
    rank0_res <- tibble(
      rep = i,
      method = "rank+X",
      strata = 1:10,
      target = sapply(1:10, function(s) {
        cat("  method = rank+X, strata =", s, "\n")
        target(dat %>% filter(strata_rank0 == s))
      })
    )
    
    rank1_res <- tibble(
      rep = i,
      method = "rank+(X-GV)",
      strata = 1:10,
      target = sapply(1:10, function(s) {
        cat("  method = rank+(X-GV), strata =", s, "\n")
        target(dat %>% filter(strata_rank1 == s))
      })
    )
    
    bind_rows(rank0_res, rank1_res)
  })
)

#===============================================================================
# Save results
#===============================================================================

save(
  target_results,
  file = "rdata_nlmr_gxe_sim_s4_thresholdUpdate_target_20260322.RData"
)
