

#===============================================================================
#
# GxE interaction for vitamin D using focused score
#
#===============================================================================

rm(list = ls())
# Set your working directory to the project folder before running this script.

library(dplyr)

#===============================================================================
# Load and prepare data
#===============================================================================

load("/rds/user/az452/hpc-work/re_analysis/rdata_vitd_crp_20250312.RData")

dat2 <- vitd.crp %>%
  left_join(dat[, c("app20175", "towndepriv")], by = "app20175") %>%
  mutate(
    season = case_when(
      vitd_month %in% 1:3   ~ "Winter",
      vitd_month %in% 4:6   ~ "Spring",
      vitd_month %in% 7:9   ~ "Summer",
      vitd_month %in% 10:12 ~ "Autumn",
      TRUE                  ~ NA_character_
    ),
    season = factor(season, levels = c("Winter", "Spring", "Summer", "Autumn"))
  )

#===============================================================================
# G-X by Townsend deprivation index
#===============================================================================

x <- dat2$vitd
g <- scale(dat2$vitd21es_grs)
v <- dat2$towndepriv

c <- data.frame(
  age    = dat2$age,
  sex    = dat2$sex,
  centre = as.factor(dat2$centre2),
  month  = as.factor(dat2$vitd_month),
  dat2[, 15:54]
)

data_for_model <- data.frame(x = x, g = g, v = v, c)
data_for_model <- data_for_model[complete.cases(data_for_model), ]

model <- lm(x ~ g * v + ., data = data_for_model)
coef_table <- summary(model)$coefficients
p_value_gv <- coef_table["g:v", "Pr(>|t|)"]

quartiles <- quantile(data_for_model$v, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)

data_for_model$v_quartile <- cut(
  data_for_model$v,
  breaks = quartiles,
  include.lowest = TRUE,
  labels = c("Q1", "Q2", "Q3", "Q4")
)

covariate_names <- names(c)
formula_str <- paste("x ~ g + v +", paste(covariate_names, collapse = " + "))
formula <- as.formula(formula_str)

strata_results <- NULL

for (quartile in c("Q1", "Q2", "Q3", "Q4")) {
  subset_data <- data_for_model[data_for_model$v_quartile == quartile, ]
  model_q <- lm(formula, data = subset_data)
  
  beta_g <- coef(model_q)["g"]
  se_g <- summary(model_q)$coefficients["g", "Std. Error"]
  
  v_percentiles <- quantile(subset_data$v, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
  
  temp_df <- data.frame(
    strata = quartile,
    b = beta_g,
    se = se_g,
    v25 = v_percentiles[1],
    v50 = v_percentiles[2],
    v75 = v_percentiles[3]
  )
  
  strata_results <- rbind(strata_results, temp_df)
}

vitd.tdi.focus <- list(
  pval_gv = p_value_gv,
  strata_results = strata_results
)

#===============================================================================
# G-X by age
#===============================================================================

x <- dat2$vitd
g <- scale(dat2$vitd21es_grs)
v <- dat2$age

c <- data.frame(
  sex    = dat2$sex,
  centre = as.factor(dat2$centre2),
  month  = as.factor(dat2$vitd_month),
  dat2[, 15:54]
)

data_for_model <- data.frame(x = x, g = g, v = v, c)
data_for_model <- data_for_model[complete.cases(data_for_model), ]

model <- lm(x ~ g * v + ., data = data_for_model)
coef_table <- summary(model)$coefficients
p_value_gv <- coef_table["g:v", "Pr(>|t|)"]

quartiles <- quantile(data_for_model$v, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)

data_for_model$age_quartile <- cut(
  data_for_model$v,
  breaks = quartiles,
  include.lowest = TRUE,
  labels = c("Q1", "Q2", "Q3", "Q4")
)

covariate_names <- names(c)
formula_str <- paste("x ~ g + v +", paste(covariate_names, collapse = " + "))
formula <- as.formula(formula_str)

strata_results <- NULL

for (quartile in c("Q1", "Q2", "Q3", "Q4")) {
  subset_data <- subset(data_for_model, age_quartile == quartile)
  if (nrow(subset_data) == 0L) next
  
  model_q <- lm(formula, data = subset_data)
  beta_g <- coef(model_q)["g"]
  se_g   <- summary(model_q)$coefficients["g", "Std. Error"]
  
  v_percentiles <- quantile(subset_data$v, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
  
  temp_df <- data.frame(
    strata = quartile,
    b      = beta_g,
    se     = se_g,
    v25    = v_percentiles[1],
    v50    = v_percentiles[2],
    v75    = v_percentiles[3]
  )
  
  strata_results <- rbind(strata_results, temp_df)
}

vitd.age.focus <- list(
  pval_gv = p_value_gv,
  strata_results = strata_results
)

#===============================================================================
# G-X by sex
#===============================================================================

x <- dat2$vitd
g <- scale(dat2$vitd21es_grs)
v <- dat2$sex

c <- data.frame(
  age    = dat2$age,
  centre = as.factor(dat2$centre2),
  month  = as.factor(dat2$vitd_month),
  dat2[, 15:54]
)

data_for_model <- data.frame(x = x, g = g, v = v, c)
data_for_model <- data_for_model[complete.cases(data_for_model), ]

model <- lm(x ~ g * v + ., data = data_for_model)
coef_table <- summary(model)$coefficients
p_value_gv <- coef_table["g:v", "Pr(>|t|)"]

sex_levels <- sort(unique(data_for_model$v))

data_for_model$sex_strata <- factor(
  data_for_model$v,
  levels = sex_levels,
  labels = if (length(sex_levels) == 2) c("Male", "Female") else as.character(sex_levels)
)

covariate_names <- names(c)
formula_str <- paste("x ~ g +", paste(covariate_names, collapse = " + "))
formula <- as.formula(formula_str)

strata_results <- NULL

for (stratum in levels(data_for_model$sex_strata)) {
  subset_data <- subset(data_for_model, sex_strata == stratum)
  if (nrow(subset_data) == 0L) next
  
  model_q <- lm(formula, data = subset_data)
  beta_g <- coef(model_q)["g"]
  se_g   <- summary(model_q)$coefficients["g", "Std. Error"]
  
  temp_df <- data.frame(
    strata = stratum,
    b      = beta_g,
    se     = se_g
  )
  
  strata_results <- rbind(strata_results, temp_df)
}

vitd.sex.focus <- list(
  pval_gv = p_value_gv,
  strata_results = strata_results
)

#===============================================================================
# G-X by frailty index
#===============================================================================

x <- dat2$vitd
g <- scale(dat2$vitd21es_grs)
v <- dat2$fi

c <- data.frame(
  age    = dat2$age,
  sex    = dat2$sex,
  centre = as.factor(dat2$centre2),
  month  = as.factor(dat2$vitd_month),
  dat2[, 15:54]
)

data_for_model <- data.frame(x = x, g = g, v = v, c)
data_for_model <- data_for_model[complete.cases(data_for_model), ]

model <- lm(x ~ g * v + ., data = data_for_model)
coef_table <- summary(model)$coefficients
p_value_gv <- coef_table["g:v", "Pr(>|t|)"]

quartiles <- quantile(data_for_model$v, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)

data_for_model$v_quartile <- cut(
  data_for_model$v,
  breaks = quartiles,
  include.lowest = TRUE,
  labels = c("Q1", "Q2", "Q3", "Q4")
)

covariate_names <- names(c)
formula_str <- paste("x ~ g + v +", paste(covariate_names, collapse = " + "))
formula <- as.formula(formula_str)

strata_results <- NULL

for (quartile in c("Q1", "Q2", "Q3", "Q4")) {
  subset_data <- subset(data_for_model, v_quartile == quartile)
  if (nrow(subset_data) == 0L) next
  
  model_q <- lm(formula, data = subset_data)
  beta_g <- coef(model_q)["g"]
  se_g   <- summary(model_q)$coefficients["g", "Std. Error"]
  
  v_percentiles <- quantile(subset_data$v, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
  
  temp_df <- data.frame(
    strata = quartile,
    b      = beta_g,
    se     = se_g,
    v25    = v_percentiles[1],
    v50    = v_percentiles[2],
    v75    = v_percentiles[3]
  )
  
  strata_results <- rbind(strata_results, temp_df)
}

vitd.fi.focus <- list(
  pval_gv = p_value_gv,
  strata_results = strata_results
)

#===============================================================================
# G-X by month of blood collection
#===============================================================================

x <- dat2$vitd
g <- scale(dat2$vitd21es_grs)
v <- factor(dat2$vitd_month)

c <- data.frame(
  age    = dat2$age,
  sex    = dat2$sex,
  centre = factor(dat2$centre2),
  dat2[, 15:54]
)

data_for_model <- data.frame(
  x      = x,
  g      = g,
  v      = v,
  season = dat2$season,
  c
)

data_for_model <- data_for_model[complete.cases(data_for_model), ]

covariate_names <- names(c)
covar_part <- paste(covariate_names, collapse = " + ")

form_int    <- as.formula(paste("x ~ g * v +", covar_part))
form_no_int <- as.formula(paste("x ~ g + v +", covar_part))

model_int    <- lm(form_int, data = data_for_model)
model_no_int <- lm(form_no_int, data = data_for_model)

anova_res  <- anova(model_no_int, model_int, test = "F")
p_value_gv <- tail(anova_res$`Pr(>F)`, 1)

form_strat <- as.formula(paste("x ~ g +", covar_part))

strata_results <- NULL
season_levels <- levels(data_for_model$season)

for (season_level in season_levels) {
  subset_data <- subset(data_for_model, season == season_level)
  if (nrow(subset_data) == 0L) next
  
  model_strat <- lm(form_strat, data = subset_data)
  beta_g <- coef(model_strat)["g"]
  se_g   <- summary(model_strat)$coefficients["g", "Std. Error"]
  
  temp_df <- data.frame(
    strata = season_level,
    b      = beta_g,
    se     = se_g,
    n      = nrow(subset_data)
  )
  
  strata_results <- rbind(strata_results, temp_df)
}

strata_results <- strata_results %>%
  mutate(
    strata = factor(strata, levels = c("Winter", "Spring", "Summer", "Autumn"))
  ) %>%
  arrange(strata)

vitd.month.focus <- list(
  pval_gv = p_value_gv,
  strata_results = strata_results
)

#===============================================================================
# Save results
#===============================================================================

save(
  vitd.tdi.focus,
  vitd.age.focus,
  vitd.sex.focus,
  vitd.fi.focus,
  vitd.month.focus,
  file = "rdata_vitd_interaction_results_focusedScore_20250405_grsSD.RData"
)
