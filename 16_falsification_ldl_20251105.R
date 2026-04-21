


#===============================================================================
#
# Falsification test for LDL
#
#===============================================================================

rm(list = ls())
# Set your working directory to the project folder before running this script.

library(dplyr)
library(SUMnlmr)

source("helpers/sim_nulloutcome.R")
source("helpers/create_nlmr_summary2025.R")

#===============================================================================
# Load and prepare data
#===============================================================================

load("rdata_nlmr_ldl_20250428.RData")

pheno <- pheno %>%
  mutate(sex = as.numeric(sex) - 1)

df <- pheno[complete.cases(pheno), ]

#===============================================================================
# Define variables
#===============================================================================

g <- df$ldlgrs
exp <- df$ldl
expsd <- (exp - mean(exp, na.rm = TRUE)) / sd(exp, na.rm = TRUE)

mod1 <- df$towndepriv
mod2 <- df$age
mod3 <- df$sex
mod4 <- df$fi

#===============================================================================
# Simulation loop
#===============================================================================

set.seed(123)
rep <- 100

for (i in 1:rep) {
  
  dat <- sim_nulloutcome(g = g, exp = expsd)
  
  #------------------------------
  # G×E-corrected exposures
  #------------------------------
  
  model1   <- lm(dat$x ~ g * mod1)
  mm1      <- model.matrix(~ g * mod1)
  coefs1   <- coef(model1)
  int1     <- grep("^g:mod1", colnames(mm1))
  int_eff1 <- as.vector(mm1[, int1, drop = FALSE] %*% coefs1[int1])
  dat$x1   <- dat$x - int_eff1
  
  model2   <- lm(dat$x ~ g * mod2)
  mm2      <- model.matrix(~ g * mod2)
  coefs2   <- coef(model2)
  int2     <- grep("^g:mod2", colnames(mm2))
  int_eff2 <- as.vector(mm2[, int2, drop = FALSE] %*% coefs2[int2])
  dat$x2   <- dat$x - int_eff2
  
  model3   <- lm(dat$x ~ g * mod3)
  mm3      <- model.matrix(~ g * mod3)
  coefs3   <- coef(model3)
  int3     <- grep("^g:mod3", colnames(mm3))
  int_eff3 <- as.vector(mm3[, int3, drop = FALSE] %*% coefs3[int3])
  dat$x3   <- dat$x - int_eff3
  
  model4   <- lm(dat$x ~ g * mod4)
  mm4      <- model.matrix(~ g * mod4)
  coefs4   <- coef(model4)
  int4     <- grep("^g:mod4", colnames(mm4))
  int_eff4 <- as.vector(mm4[, int4, drop = FALSE] %*% coefs4[int4])
  dat$x4   <- dat$x - int_eff4
  
  model_all1   <- lm(dat$x ~ g * mod1 + g * mod2 + g * mod3 + g * mod4)
  mm_all1      <- model.matrix(~ g * mod1 + g * mod2 + g * mod3 + g * mod4)
  coefs_all1   <- coef(model_all1)
  int_all1     <- grep("^g:mod", colnames(mm_all1))
  int_eff_all1 <- as.vector(mm_all1[, int_all1, drop = FALSE] %*% coefs_all1[int_all1])
  dat$x5       <- dat$x - int_eff_all1
  
  model_all2   <- lm(dat$x ~ g * mod1 * mod2 * mod3 * mod4)
  mm_all2      <- model.matrix(~ g * mod1 * mod2 * mod3 * mod4)
  coefs_all2   <- coef(model_all2)
  int_all2     <- grep("^g:", colnames(mm_all2))
  int_eff_all2 <- as.vector(mm_all2[, int_all2, drop = FALSE] %*% coefs_all2[int_all2])
  dat$x6       <- dat$x - int_eff_all2
  
  model_all1_minus_tdi   <- lm(dat$x ~ g * mod2 + g * mod3 + g * mod4)
  mm_all1_minus_tdi      <- model.matrix(~ g * mod2 + g * mod3 + g * mod4)
  coefs_all1_minus_tdi   <- coef(model_all1_minus_tdi)
  int_all1_minus_tdi     <- grep("^g:mod", colnames(mm_all1_minus_tdi))
  int_eff_all1_minus_tdi <- as.vector(
    mm_all1_minus_tdi[, int_all1_minus_tdi, drop = FALSE] %*%
      coefs_all1_minus_tdi[int_all1_minus_tdi]
  )
  dat$x7 <- dat$x - int_eff_all1_minus_tdi
  
  model_all2_minus_tdi   <- lm(dat$x ~ g * mod2 * mod3 * mod4)
  mm_all2_minus_tdi      <- model.matrix(~ g * mod2 * mod3 * mod4)
  coefs_all2_minus_tdi   <- coef(model_all2_minus_tdi)
  int_all2_minus_tdi     <- grep("^g:", colnames(mm_all2_minus_tdi))
  int_eff_all2_minus_tdi <- as.vector(
    mm_all2_minus_tdi[, int_all2_minus_tdi, drop = FALSE] %*%
      coefs_all2_minus_tdi[int_all2_minus_tdi]
  )
  dat$x8 <- dat$x - int_eff_all2_minus_tdi
  
  #------------------------------
  # MR and NLMR
  #------------------------------
  
  model_x <- lm(x ~ g, data = dat)
  summary_x <- summary(model_x)
  bx   <- summary_x$coefficients["g", "Estimate"]
  bxse <- summary_x$coefficients["g", "Std. Error"]
  
  model_y <- lm(y ~ g, data = dat)
  summary_y <- summary(model_y)
  by   <- summary_y$coefficients["g", "Estimate"]
  byse <- summary_y$coefficients["g", "Std. Error"]
  
  m1 <- data.frame(
    exposure  = "ldl",
    outcome   = "sim-ve",
    modifier  = NA,
    method    = "lmr",
    replicate = i,
    strata    = 0,
    bx        = bx,
    bxse      = bxse,
    by        = by,
    byse      = byse,
    xmean     = mean(dat$x),
    xmin      = min(dat$x),
    xmax      = max(dat$x)
  )
  
  m2 <- create_nlmr_summary2025(
    y = dat$y,
    x = dat$x,
    g = dat$g,
    covar = NULL,
    family = "gaussian",
    strata_method = "residual",
    controlsonly = FALSE,
    q = 10
  )$summary %>%
    mutate(
      exposure  = "ldl",
      outcome   = "sim-ve",
      modifier  = NA,
      method    = "residual",
      replicate = i,
      strata    = row_number()
    )
  
  m3 <- create_nlmr_summary2025(
    y = dat$y,
    x = dat$x,
    g = dat$g,
    covar = NULL,
    family = "gaussian",
    strata_method = "ranked",
    q = 10,
    seed = 123
  )$summary %>%
    mutate(
      exposure  = "ldl",
      outcome   = "sim-ve",
      modifier  = NA,
      method    = "rank",
      replicate = i,
      strata    = row_number()
    )
  
  m4 <- create_nlmr_summary2025(
    y = dat$y,
    x = dat$x,
    g = dat$g,
    xs = dat$x1,
    covar = NULL,
    family = "gaussian",
    strata_method = "ranked",
    q = 10,
    seed = 123
  )$summary %>%
    mutate(
      exposure  = "ldl",
      outcome   = "sim-ve",
      modifier  = "tdi",
      method    = "rank(tdi)",
      replicate = i,
      strata    = row_number()
    )
  
  m5 <- create_nlmr_summary2025(
    y = dat$y,
    x = dat$x,
    g = dat$g,
    xs = dat$x2,
    covar = NULL,
    family = "gaussian",
    strata_method = "ranked",
    q = 10,
    seed = 123
  )$summary %>%
    mutate(
      exposure  = "ldl",
      outcome   = "sim-ve",
      modifier  = "age",
      method    = "rank(age)",
      replicate = i,
      strata    = row_number()
    )
  
  m6 <- create_nlmr_summary2025(
    y = dat$y,
    x = dat$x,
    g = dat$g,
    xs = dat$x3,
    covar = NULL,
    family = "gaussian",
    strata_method = "ranked",
    q = 10,
    seed = 123
  )$summary %>%
    mutate(
      exposure  = "ldl",
      outcome   = "sim-ve",
      modifier  = "sex",
      method    = "rank(sex)",
      replicate = i,
      strata    = row_number()
    )
  
  m7 <- create_nlmr_summary2025(
    y = dat$y,
    x = dat$x,
    g = dat$g,
    xs = dat$x4,
    covar = NULL,
    family = "gaussian",
    strata_method = "ranked",
    q = 10,
    seed = 123
  )$summary %>%
    mutate(
      exposure  = "ldl",
      outcome   = "sim-ve",
      modifier  = "fi",
      method    = "rank(fi)",
      replicate = i,
      strata    = row_number()
    )
  
  m8 <- create_nlmr_summary2025(
    y = dat$y,
    x = dat$x,
    g = dat$g,
    xs = dat$x5,
    covar = NULL,
    family = "gaussian",
    strata_method = "ranked",
    q = 10,
    seed = 123
  )$summary %>%
    mutate(
      exposure  = "ldl",
      outcome   = "sim-ve",
      modifier  = "all1",
      method    = "rank(all1)",
      replicate = i,
      strata    = row_number()
    )
  
  m9 <- create_nlmr_summary2025(
    y = dat$y,
    x = dat$x,
    g = dat$g,
    xs = dat$x6,
    covar = NULL,
    family = "gaussian",
    strata_method = "ranked",
    q = 10,
    seed = 123
  )$summary %>%
    mutate(
      exposure  = "ldl",
      outcome   = "sim-ve",
      modifier  = "all2",
      method    = "rank(all2)",
      replicate = i,
      strata    = row_number()
    )
  
  m10 <- create_nlmr_summary2025(
    y = dat$y,
    x = dat$x,
    g = dat$g,
    xs = dat$x7,
    covar = NULL,
    family = "gaussian",
    strata_method = "ranked",
    q = 10,
    seed = 123
  )$summary %>%
    mutate(
      exposure  = "ldl",
      outcome   = "sim-ve",
      modifier  = "all1-tdi",
      method    = "rank(all1-tdi)",
      replicate = i,
      strata    = row_number()
    )
  
  m11 <- create_nlmr_summary2025(
    y = dat$y,
    x = dat$x,
    g = dat$g,
    xs = dat$x8,
    covar = NULL,
    family = "gaussian",
    strata_method = "ranked",
    q = 10,
    seed = 123
  )$summary %>%
    mutate(
      exposure  = "ldl",
      outcome   = "sim-ve",
      modifier  = "all2-tdi",
      method    = "rank(all2-tdi)",
      replicate = i,
      strata    = row_number()
    )
  
  temp <- bind_rows(m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11) %>%
    mutate(
      lace    = by / bx,
      lace_se = byse / bx
    )
  
  if (i == 1) {
    results <- temp
  } else {
    results <- rbind(results, temp)
  }
  
  cat("rep ", i, "/", rep, " is done\n", sep = "")
}

#===============================================================================
# Save results
#===============================================================================

ldl_simNeg <- results

save(
  ldl_simNeg,
  file = "rdata_nlmr_ldl_negControl_results_20251105.RData"
)
