
#===============================================================================
#
# Empirical negative controls for vitamin D using focused score
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

load("/rds/user/az452/hpc-work/re_analysis/rdata_vitd_crp_20250312.RData")

df0 <- vitd.crp %>%
  left_join(dat[, c("app20175", "towndepriv")], by = "app20175") %>%
  mutate(sex = ifelse(sex == 1, 1, 0))

df <- df0[complete.cases(df0), ]

#===============================================================================
# Define variables
#===============================================================================

g   <- df$vitd21es_grs
exp <- df$vitd

mod1 <- df$towndepriv
mod2 <- df$age
mod3 <- df$sex
mod4 <- df$fi
mod5 <- factor(df$vitd_month)

cov0 <- data.frame(
  age     = df$age,
  sex     = df$sex,
  centre  = as.factor(df$centre2),
  array   = df$snp_array,
  fasting = as.factor(df$fast_q6),
  month   = as.factor(df$vitd_month),
  aliquot = as.factor(df$vitd_ali),
  birthns = as.factor(df$birth_ns),
  birthew = as.factor(df$birth_ew),
  df[, c(15:54)]
)

#===============================================================================
# G×E-corrected exposures
#===============================================================================

model5      <- lm(exp ~ g * mod5)
mm_month    <- model.matrix(~ g * mod5)
coefs_month <- coef(model5)

int_cols_month <- grep("^g:mod5", colnames(mm_month))
int_eff_month  <- as.vector(mm_month[, int_cols_month, drop = FALSE] %*%
                              coefs_month[int_cols_month])

exp1 <- exp - int_eff_month

model_all1 <- lm(exp ~ g * mod1 + g * mod2 + g * mod3 + g * mod4 + g * mod5)
mm_all1    <- model.matrix(~ g * mod1 + g * mod2 + g * mod3 + g * mod4 + g * mod5)
coefs_all1 <- coef(model_all1)

int_cols_all1 <- grep("^g:mod", colnames(mm_all1))
int_eff_all1  <- as.vector(mm_all1[, int_cols_all1, drop = FALSE] %*%
                             coefs_all1[int_cols_all1])

exp2 <- exp - int_eff_all1

#===============================================================================
# Age as empirical negative control
#===============================================================================

out <- df$age
cov <- model.matrix(~ ., data = dplyr::select(cov0, -age))[, -1, drop = FALSE]

dat <- data.frame(
  x  = exp,
  g  = g,
  y  = out,
  x1 = exp1,
  x2 = exp2
)

model_x   <- lm(x ~ g + cov, data = dat)
summary_x <- summary(model_x)
bx   <- summary_x$coefficients["g", "Estimate"]
bxse <- summary_x$coefficients["g", "Std. Error"]

model_y   <- lm(y ~ g + cov, data = dat)
summary_y <- summary(model_y)
by   <- summary_y$coefficients["g", "Estimate"]
byse <- summary_y$coefficients["g", "Std. Error"]

m1 <- data.frame(
  exposure = "vitd",
  outcome  = "age",
  modifier = NA,
  method   = "lmr",
  strata   = 0,
  bx       = bx,
  bxse     = bxse,
  by       = by,
  byse     = byse,
  xmean    = mean(dat$x),
  xmin     = min(dat$x),
  xmax     = max(dat$x)
)

m2 <- create_nlmr_summary2025(
  y             = dat$y,
  x             = dat$x,
  g             = dat$g,
  covar         = cov,
  family        = "gaussian",
  strata_method = "residual",
  controlsonly  = FALSE,
  q             = 10
)$summary %>%
  mutate(
    exposure = "vitd",
    outcome  = "age",
    modifier = NA,
    method   = "residual",
    strata   = row_number()
  )

m3 <- create_nlmr_summary2025(
  y             = dat$y,
  x             = dat$x,
  g             = dat$g,
  covar         = cov,
  family        = "gaussian",
  strata_method = "ranked",
  q             = 10,
  seed          = 123
)$summary %>%
  mutate(
    exposure = "vitd",
    outcome  = "age",
    modifier = NA,
    method   = "rank",
    strata   = row_number()
  )

m4 <- create_nlmr_summary2025(
  y             = dat$y,
  x             = dat$x,
  g             = dat$g,
  xs            = dat$x1,
  covar         = cov,
  family        = "gaussian",
  strata_method = "ranked",
  q             = 10,
  seed          = 123
)$summary %>%
  mutate(
    exposure = "vitd",
    outcome  = "age",
    modifier = "month",
    method   = "rank(month)",
    strata   = row_number()
  )

m5 <- create_nlmr_summary2025(
  y             = dat$y,
  x             = dat$x,
  g             = dat$g,
  xs            = dat$x2,
  covar         = cov,
  family        = "gaussian",
  strata_method = "ranked",
  q             = 10,
  seed          = 123
)$summary %>%
  mutate(
    exposure = "vitd",
    outcome  = "age",
    modifier = "all1",
    method   = "rank(all1)",
    strata   = row_number()
  )

vitd_age <- bind_rows(m1, m2, m3, m4, m5) %>%
  mutate(
    lace    = by / bx,
    lace_se = byse / bx
  )

#===============================================================================
# Sex as empirical negative control
#===============================================================================

out <- df$sex
cov <- model.matrix(~ ., data = dplyr::select(cov0, -sex))[, -1, drop = FALSE]

dat <- data.frame(
  x  = exp,
  g  = g,
  y  = out,
  x1 = exp1,
  x2 = exp2
)

model_x   <- lm(x ~ g + cov, data = dat)
summary_x <- summary(model_x)
bx   <- summary_x$coefficients["g", "Estimate"]
bxse <- summary_x$coefficients["g", "Std. Error"]

model_y   <- glm(y ~ g + cov, data = dat, family = binomial)
summary_y <- summary(model_y)
by   <- summary_y$coefficients["g", "Estimate"]
byse <- summary_y$coefficients["g", "Std. Error"]

m1 <- data.frame(
  exposure = "vitd",
  outcome  = "sex",
  modifier = NA,
  method   = "lmr",
  strata   = 0,
  bx       = bx,
  bxse     = bxse,
  by       = by,
  byse     = byse,
  xmean    = mean(dat$x),
  xmin     = min(dat$x),
  xmax     = max(dat$x)
)

m2 <- create_nlmr_summary2025(
  y             = dat$y,
  x             = dat$x,
  g             = dat$g,
  covar         = cov,
  family        = "binomial",
  strata_method = "residual",
  controlsonly  = FALSE,
  q             = 10
)$summary %>%
  mutate(
    exposure = "vitd",
    outcome  = "sex",
    modifier = NA,
    method   = "residual",
    strata   = row_number()
  )

m3 <- create_nlmr_summary2025(
  y             = dat$y,
  x             = dat$x,
  g             = dat$g,
  covar         = cov,
  family        = "binomial",
  strata_method = "ranked",
  q             = 10,
  seed          = 123
)$summary %>%
  mutate(
    exposure = "vitd",
    outcome  = "sex",
    modifier = NA,
    method   = "rank",
    strata   = row_number()
  )

m4 <- create_nlmr_summary2025(
  y             = dat$y,
  x             = dat$x,
  g             = dat$g,
  xs            = dat$x1,
  covar         = cov,
  family        = "binomial",
  strata_method = "ranked",
  q             = 10,
  seed          = 123
)$summary %>%
  mutate(
    exposure = "vitd",
    outcome  = "sex",
    modifier = "month",
    method   = "rank(month)",
    strata   = row_number()
  )

m5 <- create_nlmr_summary2025(
  y             = dat$y,
  x             = dat$x,
  g             = dat$g,
  xs            = dat$x2,
  covar         = cov,
  family        = "binomial",
  strata_method = "ranked",
  q             = 10,
  seed          = 123
)$summary %>%
  mutate(
    exposure = "vitd",
    outcome  = "sex",
    modifier = "all1",
    method   = "rank(all1)",
    strata   = row_number()
  )

vitd_sex <- bind_rows(m1, m2, m3, m4, m5) %>%
  mutate(
    lace    = by / bx,
    lace_se = byse / bx
  )

#===============================================================================
# Save results
#===============================================================================

save(
  vitd_age,
  vitd_sex,
  file = "rdata_nlmr_vitd_EmpNegPosControls_results_focusedScore_20251205.RData"
)

