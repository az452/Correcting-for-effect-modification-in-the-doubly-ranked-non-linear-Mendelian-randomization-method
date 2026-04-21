
#===============================================================================
#
# Empirical negative and positive controls for LDL
#
#===============================================================================

rm(list = ls())
# Set your working directory to the project folder before running this script.

library(dplyr)
library(SUMnlmr)

source("helpers/create_nlmr_summary2025.R")

#===============================================================================
# Load and prepare data
#===============================================================================

load("rdata_nlmr_ldl_20250428.RData")

pheno <- pheno %>%
  mutate(sex = as.numeric(sex) - 1)

df_all <- pheno[complete.cases(pheno), ]

#===============================================================================
# Define variables
#===============================================================================

g   <- df_all$ldlgrs
exp <- df_all$ldl

mod_fi  <- df_all$fi
mod_age <- df_all$age
mod_sex <- df_all$sex
mod_tdi <- df_all$towndepriv

cov0_base <- data.frame(
  age    = df_all$age,
  sex    = as.factor(df_all$sex),
  centre = df_all$centre,
  df_all[, c(10:49)]
)

#===============================================================================
# G×E-corrected exposures
#===============================================================================

model_all1 <- lm(exp ~ g * mod_fi + g * mod_age + g * mod_sex)
b1a <- coef(model_all1)["g:mod_fi"]
b2a <- coef(model_all1)["g:mod_age"]
b3a <- coef(model_all1)["g:mod_sex"]

exp4 <- exp - (b1a * g * mod_fi + b2a * g * mod_age + b3a * g * mod_sex)

model_all2 <- lm(exp ~ g * mod_fi * mod_age * mod_sex)
coefs <- coef(model_all2)

b1b   <- coefs["g:mod_fi"]
b2b   <- coefs["g:mod_age"]
b3b   <- coefs["g:mod_sex"]
b12b  <- coefs["g:mod_fi:mod_age"]
b13b  <- coefs["g:mod_fi:mod_sex"]
b23b  <- coefs["g:mod_age:mod_sex"]
b123b <- coefs["g:mod_fi:mod_age:mod_sex"]

exp5 <- exp - (
  b1b   * g * mod_fi +
    b2b   * g * mod_age +
    b3b   * g * mod_sex +
    b12b  * g * mod_fi * mod_age +
    b13b  * g * mod_fi * mod_sex +
    b23b  * g * mod_age * mod_sex +
    b123b * g * mod_fi * mod_age * mod_sex
)

#===============================================================================
# Age as empirical negative control
#===============================================================================

out_age <- df_all$age
cov_age <- model.matrix(~ ., data = dplyr::select(cov0_base, -age))[, -1, drop = FALSE]

dat_age <- data.frame(
  x  = exp,
  g  = g,
  y  = out_age,
  x4 = exp4,
  x5 = exp5
)

model_x <- lm(x ~ g + cov_age, data = dat_age)
summary_x <- summary(model_x)
bx   <- summary_x$coefficients["g", "Estimate"]
bxse <- summary_x$coefficients["g", "Std. Error"]

model_y <- lm(y ~ g + cov_age, data = dat_age)
summary_y <- summary(model_y)
by   <- summary_y$coefficients["g", "Estimate"]
byse <- summary_y$coefficients["g", "Std. Error"]

m1_age <- data.frame(
  exposure = "ldl",
  outcome  = "age",
  modifier = NA,
  method   = "lmr",
  strata   = 0,
  bx       = bx,
  bxse     = bxse,
  by       = by,
  byse     = byse,
  xmean    = mean(dat_age$x),
  xmin     = min(dat_age$x),
  xmax     = max(dat_age$x)
)

m2_age <- create_nlmr_summary2025(
  y             = dat_age$y,
  x             = dat_age$x,
  g             = dat_age$g,
  covar         = cov_age,
  family        = "gaussian",
  strata_method = "residual",
  controlsonly  = FALSE,
  q             = 10
)$summary %>%
  mutate(
    exposure = "ldl",
    outcome  = "age",
    modifier = NA,
    method   = "residual",
    strata   = row_number()
  )

m3_age <- create_nlmr_summary2025(
  y             = dat_age$y,
  x             = dat_age$x,
  g             = dat_age$g,
  covar         = cov_age,
  family        = "gaussian",
  strata_method = "ranked",
  q             = 10,
  seed          = 123
)$summary %>%
  mutate(
    exposure = "ldl",
    outcome  = "age",
    modifier = NA,
    method   = "rank",
    strata   = row_number()
  )

m4_age <- create_nlmr_summary2025(
  y             = dat_age$y,
  x             = dat_age$x,
  g             = dat_age$g,
  xs            = dat_age$x4,
  covar         = cov_age,
  family        = "gaussian",
  strata_method = "ranked",
  q             = 10,
  seed          = 123
)$summary %>%
  mutate(
    exposure = "ldl",
    outcome  = "age",
    modifier = "all1",
    method   = "rank(all1)",
    strata   = row_number()
  )

m5_age <- create_nlmr_summary2025(
  y             = dat_age$y,
  x             = dat_age$x,
  g             = dat_age$g,
  xs            = dat_age$x5,
  covar         = cov_age,
  family        = "gaussian",
  strata_method = "ranked",
  q             = 10,
  seed          = 123
)$summary %>%
  mutate(
    exposure = "ldl",
    outcome  = "age",
    modifier = "all2",
    method   = "rank(all2)",
    strata   = row_number()
  )

ldl_age <- bind_rows(m1_age, m2_age, m3_age, m4_age, m5_age) %>%
  mutate(
    lace    = by / bx,
    lace_se = byse / bx
  )

#===============================================================================
# Sex as empirical negative control
#===============================================================================

out_sex <- df_all$sex
cov_sex <- model.matrix(~ ., data = dplyr::select(cov0_base, -sex))[, -1, drop = FALSE]

dat_sex <- data.frame(
  x  = exp,
  g  = g,
  y  = out_sex,
  x4 = exp4,
  x5 = exp5
)

model_x <- lm(x ~ g + cov_sex, data = dat_sex)
summary_x <- summary(model_x)
bx   <- summary_x$coefficients["g", "Estimate"]
bxse <- summary_x$coefficients["g", "Std. Error"]

model_y <- glm(y ~ g + cov_sex, data = dat_sex, family = binomial)
summary_y <- summary(model_y)
by   <- summary_y$coefficients["g", "Estimate"]
byse <- summary_y$coefficients["g", "Std. Error"]

m1_sex <- data.frame(
  exposure = "ldl",
  outcome  = "sex",
  modifier = NA,
  method   = "lmr",
  strata   = 0,
  bx       = bx,
  bxse     = bxse,
  by       = by,
  byse     = byse,
  xmean    = mean(dat_sex$x),
  xmin     = min(dat_sex$x),
  xmax     = max(dat_sex$x)
)

m2_sex <- create_nlmr_summary2025(
  y             = dat_sex$y,
  x             = dat_sex$x,
  g             = dat_sex$g,
  covar         = cov_sex,
  family        = "binomial",
  strata_method = "residual",
  controlsonly  = FALSE,
  q             = 10
)$summary %>%
  mutate(
    exposure = "ldl",
    outcome  = "sex",
    modifier = NA,
    method   = "residual",
    strata   = row_number()
  )

m3_sex <- create_nlmr_summary2025(
  y             = dat_sex$y,
  x             = dat_sex$x,
  g             = dat_sex$g,
  covar         = cov_sex,
  family        = "binomial",
  strata_method = "ranked",
  q             = 10,
  seed          = 123
)$summary %>%
  mutate(
    exposure = "ldl",
    outcome  = "sex",
    modifier = NA,
    method   = "rank",
    strata   = row_number()
  )

m4_sex <- create_nlmr_summary2025(
  y             = dat_sex$y,
  x             = dat_sex$x,
  g             = dat_sex$g,
  xs            = dat_sex$x4,
  covar         = cov_sex,
  family        = "binomial",
  strata_method = "ranked",
  q             = 10,
  seed          = 123
)$summary %>%
  mutate(
    exposure = "ldl",
    outcome  = "sex",
    modifier = "all1",
    method   = "rank(all1)",
    strata   = row_number()
  )

m5_sex <- create_nlmr_summary2025(
  y             = dat_sex$y,
  x             = dat_sex$x,
  g             = dat_sex$g,
  xs            = dat_sex$x5,
  covar         = cov_sex,
  family        = "binomial",
  strata_method = "ranked",
  q             = 10,
  seed          = 123
)$summary %>%
  mutate(
    exposure = "ldl",
    outcome  = "sex",
    modifier = "all2",
    method   = "rank(all2)",
    strata   = row_number()
  )

ldl_sex <- bind_rows(m1_sex, m2_sex, m3_sex, m4_sex, m5_sex) %>%
  mutate(
    lace    = by / bx,
    lace_se = byse / bx
  )

#===============================================================================
# CHD as positive control
#===============================================================================

out_chd <- df_all$chd
cov_chd <- model.matrix(~ ., data = cov0_base)[, -1, drop = FALSE]

dat_chd <- data.frame(
  x  = exp,
  g  = g,
  y  = out_chd,
  x4 = exp4,
  x5 = exp5
)

model_x <- lm(x ~ g + cov_chd, data = dat_chd)
summary_x <- summary(model_x)
bx   <- summary_x$coefficients["g", "Estimate"]
bxse <- summary_x$coefficients["g", "Std. Error"]

model_y <- glm(y ~ g + cov_chd, data = dat_chd, family = binomial)
summary_y <- summary(model_y)
by   <- summary_y$coefficients["g", "Estimate"]
byse <- summary_y$coefficients["g", "Std. Error"]

m1_chd <- data.frame(
  exposure = "ldl",
  outcome  = "CAD",
  modifier = NA,
  method   = "lmr",
  strata   = 0,
  bx       = bx,
  bxse     = bxse,
  by       = by,
  byse     = byse,
  xmean    = mean(dat_chd$x),
  xmin     = min(dat_chd$x),
  xmax     = max(dat_chd$x)
)

m2_chd <- create_nlmr_summary2025(
  y             = dat_chd$y,
  x             = dat_chd$x,
  g             = dat_chd$g,
  covar         = cov_chd,
  family        = "binomial",
  strata_method = "residual",
  controlsonly  = FALSE,
  q             = 10
)$summary %>%
  mutate(
    exposure = "ldl",
    outcome  = "CAD",
    modifier = NA,
    method   = "residual",
    strata   = row_number()
  )

m3_chd <- create_nlmr_summary2025(
  y             = dat_chd$y,
  x             = dat_chd$x,
  g             = dat_chd$g,
  covar         = cov_chd,
  family        = "binomial",
  strata_method = "ranked",
  q             = 10,
  seed          = 123
)$summary %>%
  mutate(
    exposure = "ldl",
    outcome  = "CAD",
    modifier = NA,
    method   = "rank",
    strata   = row_number()
  )

m4_chd <- create_nlmr_summary2025(
  y             = dat_chd$y,
  x             = dat_chd$x,
  g             = dat_chd$g,
  xs            = dat_chd$x4,
  covar         = cov_chd,
  family        = "binomial",
  strata_method = "ranked",
  q             = 10,
  seed          = 123
)$summary %>%
  mutate(
    exposure = "ldl",
    outcome  = "CAD",
    modifier = "all1",
    method   = "rank(all1)",
    strata   = row_number()
  )

m5_chd <- create_nlmr_summary2025(
  y             = dat_chd$y,
  x             = dat_chd$x,
  g             = dat_chd$g,
  xs            = dat_chd$x5,
  covar         = cov_chd,
  family        = "binomial",
  strata_method = "ranked",
  q             = 10,
  seed          = 123
)$summary %>%
  mutate(
    exposure = "ldl",
    outcome  = "CAD",
    modifier = "all2",
    method   = "rank(all2)",
    strata   = row_number()
  )

ldl_cad <- bind_rows(m1_chd, m2_chd, m3_chd, m4_chd, m5_chd) %>%
  mutate(
    lace    = by / bx,
    lace_se = byse / bx
  )

#===============================================================================
# Save results
#===============================================================================

save(
  ldl_age,
  ldl_sex,
  ldl_cad,
  file = "rdata_nlmr_ldl_EmpNegPosControls_results_20251105.RData"
)
