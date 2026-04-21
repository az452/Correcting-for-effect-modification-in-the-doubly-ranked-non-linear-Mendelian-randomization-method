
#===============================================================================
#
# GxE interaction for LDL
#
#===============================================================================

rm(list = ls())
# Set your working directory to the project folder before running this script.

library(dplyr)

#===============================================================================
# Load data
#===============================================================================

load("rdata_nlmr_ldl_20250428.RData")

dat2 <- pheno

#===============================================================================
# G-X by Townsend deprivation index
#===============================================================================

x <- dat2$ldl
g <- dat2$ldlgrs
v <- dat2$towndepriv

c <- data.frame(
  age    = dat2$age,
  sex    = dat2$sex,
  centre = as.factor(dat2$centre),
  dat2[, 10:49]
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

ldl.tdi <- list(
  pval_gv = p_value_gv,
  strata_results = strata_results
)

#===============================================================================
# G-X by age
#===============================================================================

x <- dat2$ldl
g <- dat2$ldlgrs
v <- dat2$age

c <- data.frame(
  sex    = dat2$sex,
  centre = as.factor(dat2$centre),
  dat2[, 10:49]
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
  subset_data <- data_for_model[data_for_model$age_quartile == quartile, ]
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

ldl.age <- list(
  pval_gv = p_value_gv,
  strata_results = strata_results
)

#===============================================================================
# G-X by sex
#===============================================================================

dat2_sex <- pheno %>%
  mutate(sex = ifelse(as.numeric(sex) == 2, 1, 2))

x <- dat2_sex$ldl
g <- dat2_sex$ldlgrs
v <- dat2_sex$sex

c <- data.frame(
  age    = dat2_sex$age,
  centre = as.factor(dat2_sex$centre),
  dat2_sex[, 10:49]
)

data_for_model <- data.frame(x = x, g = g, v = v, c)
data_for_model <- data_for_model[complete.cases(data_for_model), ]

model <- lm(x ~ g * v + ., data = data_for_model)
coef_table <- summary(model)$coefficients
p_value_gv <- coef_table["g:v", "Pr(>|t|)"]

data_for_model$sex_strata <- factor(
  data_for_model$v,
  levels = c(1, 2),
  labels = c("Male", "Female")
)

covariate_names <- names(c)
formula_str <- paste("x ~ g +", paste(covariate_names, collapse = " + "))
formula <- as.formula(formula_str)

strata_results <- NULL

for (stratum in c("Male", "Female")) {
  subset_data <- data_for_model[data_for_model$sex_strata == stratum, ]
  model_q <- lm(formula, data = subset_data)
  
  beta_g <- coef(model_q)["g"]
  se_g <- summary(model_q)$coefficients["g", "Std. Error"]
  
  temp_df <- data.frame(
    strata = stratum,
    b = beta_g,
    se = se_g
  )
  
  strata_results <- rbind(strata_results, temp_df)
}

ldl.sex <- list(
  pval_gv = p_value_gv,
  strata_results = strata_results
)

#===============================================================================
# G-X by frailty index
#===============================================================================

x <- dat2$ldl
g <- dat2$ldlgrs
v <- dat2$fi

c <- data.frame(
  age    = dat2$age,
  sex    = dat2$sex,
  centre = as.factor(dat2$centre),
  dat2[, 10:49]
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

ldl.fi <- list(
  pval_gv = p_value_gv,
  strata_results = strata_results
)

#===============================================================================
# Save results
#===============================================================================

save(
  ldl.tdi,
  ldl.age,
  ldl.sex,
  ldl.fi,
  file = "rdata_ldl_interaction_results_20250428.RData"
)
