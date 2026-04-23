
#===============================================================================
#
# Falsification test for vitamin D using replicated score
#
#===============================================================================

rm(list = ls())
# Set your working directory to the project folder before running this script.

library(dplyr)
library(SUMnlmr)

source("helpers/custom_summary.R")
source("helpers/sim_nulloutcome.R")
source("helpers/create_nlmr_summary2025.R")

#===============================================================================
# Load and prepare data
#===============================================================================

load("rdata_vitd_crp_20250312.RData")

vitd.crp <- vitd.crp %>%
  left_join(dat[, c("app20175", "towndepriv")], by = "app20175") %>%
  mutate(
    sex = ifelse(sex == 1, 1, 0),
    season = case_when(
      vitd_month %in% 1:3   ~ "winter",
      vitd_month %in% 4:6   ~ "spring",
      vitd_month %in% 7:9   ~ "summer",
      vitd_month %in% 10:12 ~ "autumn",
      TRUE                  ~ NA_character_
    ),
    season = factor(season, levels = c("summer", "autumn", "winter", "spring")),
    vitd_month = factor(vitd_month, levels = 1:12)
  )

df0 <- vitd.crp[, c(
  "vitdSun_grs", "vitd", "towndepriv", "age",
  "sex", "fi", "season", "vitd_month"
)]
df <- df0[complete.cases(df0), ]

#===============================================================================
# Define variables
#===============================================================================

g <- df$vitdSun_grs
exp <- df$vitd
expsd <- (exp - mean(exp, na.rm = TRUE)) / sd(exp, na.rm = TRUE)

tdi    <- df$towndepriv
age    <- df$age
sex    <- df$sex
fi     <- df$fi
season <- df$season
month  <- df$vitd_month

#===============================================================================
# Simulation loop
#===============================================================================

set.seed(123)
rep <- 100

for (i in 1:rep) {
  
  #------------------------------
  # Simulate negative-control outcome
  #------------------------------
  
  dat <- sim_nulloutcome(g = g, exp = expsd)
  
  #------------------------------
  # G×E-corrected exposures
  #------------------------------
  
  model_tdi   <- lm(dat$x ~ g * tdi)
  mm_tdi      <- model.matrix(~ g * tdi)
  coefs_tdi   <- coef(model_tdi)
  int_tdi     <- grep("^g:tdi", colnames(mm_tdi))
  int_eff_tdi <- as.vector(mm_tdi[, int_tdi, drop = FALSE] %*% coefs_tdi[int_tdi])
  dat$x_tdi   <- dat$x - int_eff_tdi
  
  model_age   <- lm(dat$x ~ g * age)
  mm_age      <- model.matrix(~ g * age)
  coefs_age   <- coef(model_age)
  int_age     <- grep("^g:age", colnames(mm_age))
  int_eff_age <- as.vector(mm_age[, int_age, drop = FALSE] %*% coefs_age[int_age])
  dat$x_age   <- dat$x - int_eff_age
  
  model_sex   <- lm(dat$x ~ g * sex)
  mm_sex      <- model.matrix(~ g * sex)
  coefs_sex   <- coef(model_sex)
  int_sex     <- grep("^g:sex", colnames(mm_sex))
  int_eff_sex <- as.vector(mm_sex[, int_sex, drop = FALSE] %*% coefs_sex[int_sex])
  dat$x_sex   <- dat$x - int_eff_sex
  
  model_fi    <- lm(dat$x ~ g * fi)
  mm_fi       <- model.matrix(~ g * fi)
  coefs_fi    <- coef(model_fi)
  int_fi      <- grep("^g:fi", colnames(mm_fi))
  int_eff_fi  <- as.vector(mm_fi[, int_fi, drop = FALSE] %*% coefs_fi[int_fi])
  dat$x_fi    <- dat$x - int_eff_fi
  
  model_month   <- lm(dat$x ~ g * month)
  mm_month      <- model.matrix(~ g * month)
  coefs_month   <- coef(model_month)
  int_month     <- grep("^g:month", colnames(mm_month))
  int_eff_month <- as.vector(mm_month[, int_month, drop = FALSE] %*% coefs_month[int_month])
  dat$x_month   <- dat$x - int_eff_month
  
  model_season   <- lm(dat$x ~ g * season)
  mm_season      <- model.matrix(~ g * season)
  coefs_season   <- coef(model_season)
  int_season     <- grep("^g:season", colnames(mm_season))
  int_eff_season <- as.vector(mm_season[, int_season, drop = FALSE] %*% coefs_season[int_season])
  dat$x_season   <- dat$x - int_eff_season
  
  model_all1   <- lm(dat$x ~ g * tdi + g * age + g * sex + g * fi + g * month)
  mm_all1      <- model.matrix(~ g * tdi + g * age + g * sex + g * fi + g * month)
  coefs_all1   <- coef(model_all1)
  int_all1     <- grep("^g:", colnames(mm_all1))
  int_eff_all1 <- as.vector(mm_all1[, int_all1, drop = FALSE] %*% coefs_all1[int_all1])
  dat$x_all1   <- dat$x - int_eff_all1
  
  model_all2   <- lm(dat$x ~ g * tdi * age * sex * fi * month)
  mm_all2      <- model.matrix(~ g * tdi * age * sex * fi * month)
  coefs_all2   <- coef(model_all2)
  int_all2     <- grep("^g:", colnames(mm_all2))
  int_eff_all2 <- as.vector(mm_all2[, int_all2, drop = FALSE] %*% coefs_all2[int_all2])
  dat$x_all2   <- dat$x - int_eff_all2
  
  #------------------------------
  # MR and NLMR
  #------------------------------
  
  model_x   <- lm(x ~ g, data = dat)
  summary_x <- summary(model_x)
  bx   <- summary_x$coefficients["g", "Estimate"]
  bxse <- summary_x$coefficients["g", "Std. Error"]
  
  model_y   <- lm(y ~ g, data = dat)
  summary_y <- summary(model_y)
  by   <- summary_y$coefficients["g", "Estimate"]
  byse <- summary_y$coefficients["g", "Std. Error"]
  
  m1 <- data.frame(
    exposure  = "vitd",
    outcome   = "sim-ve",
    modifier  = NA,
    method    = "lmr",
    replicate = i,
    strata    = 0,
    bx   = bx,
    bxse = bxse,
    by   = by,
    byse = byse,
    xmean = mean(dat$x),
    xmin  = min(dat$x),
    xmax  = max(dat$x)
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
      exposure  = "vitd",
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
      exposure  = "vitd",
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
    xs = dat$x_tdi,
    covar = NULL,
    family = "gaussian",
    strata_method = "ranked",
    q = 10,
    seed = 123
  )$summary %>%
    mutate(
      exposure  = "vitd",
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
    xs = dat$x_age,
    covar = NULL,
    family = "gaussian",
    strata_method = "ranked",
    q = 10,
    seed = 123
  )$summary %>%
    mutate(
      exposure  = "vitd",
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
    xs = dat$x_sex,
    covar = NULL,
    family = "gaussian",
    strata_method = "ranked",
    q = 10,
    seed = 123
  )$summary %>%
    mutate(
      exposure  = "vitd",
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
    xs = dat$x_fi,
    covar = NULL,
    family = "gaussian",
    strata_method = "ranked",
    q = 10,
    seed = 123
  )$summary %>%
    mutate(
      exposure  = "vitd",
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
    xs = dat$x_season,
    covar = NULL,
    family = "gaussian",
    strata_method = "ranked",
    q = 10,
    seed = 123
  )$summary %>%
    mutate(
      exposure  = "vitd",
      outcome   = "sim-ve",
      modifier  = "season",
      method    = "rank(season)",
      replicate = i,
      strata    = row_number()
    )
  
  m9 <- create_nlmr_summary2025(
    y = dat$y,
    x = dat$x,
    g = dat$g,
    xs = dat$x_month,
    covar = NULL,
    family = "gaussian",
    strata_method = "ranked",
    q = 10,
    seed = 123
  )$summary %>%
    mutate(
      exposure  = "vitd",
      outcome   = "sim-ve",
      modifier  = "month",
      method    = "rank(month)",
      replicate = i,
      strata    = row_number()
    )
  
  m10 <- create_nlmr_summary2025(
    y = dat$y,
    x = dat$x,
    g = dat$g,
    xs = dat$x_all1,
    covar = NULL,
    family = "gaussian",
    strata_method = "ranked",
    q = 10,
    seed = 123
  )$summary %>%
    mutate(
      exposure  = "vitd",
      outcome   = "sim-ve",
      modifier  = "all1",
      method    = "rank(all1)",
      replicate = i,
      strata    = row_number()
    )
  
  m11 <- create_nlmr_summary2025(
    y = dat$y,
    x = dat$x,
    g = dat$g,
    xs = dat$x_all2,
    covar = NULL,
    family = "gaussian",
    strata_method = "ranked",
    q = 10,
    seed = 123
  )$summary %>%
    mutate(
      exposure  = "vitd",
      outcome   = "sim-ve",
      modifier  = "all2",
      method    = "rank(all2)",
      replicate = i,
      strata    = row_number()
    )
  
  #------------------------------
  # Combine results
  #------------------------------
  
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

vitd_simNeg <- results

save(
  vitd_simNeg,
  file = "rdata_nlmr_vitd_negControl_results_20251105.RData"
)

