

#=====================================================================================================
#
# Figure S10: heterogeneity measured by normalized MSE for LDL-C emperical negative/positive controls
# including Age, Sex, and CAD
#
#=====================================================================================================

rm(list = ls())
# Set your working directory to the project folder before running this script.

library(dplyr)
library(ggplot2)


#========================================
# Load data
#========================================
load("rdata_nlmr_ldl_EmpNegPosControls_results_20251105.RData")

#========================================
# Data processing
#========================================
lmr_ref <- bind_rows(ldl_age, ldl_sex, ldl_cad) %>%
  filter(method == "lmr") %>%
  group_by(exposure, outcome) %>%
  summarise(
    lace_LMR = lace[1],
    .groups = "drop"
  )

dat <- bind_rows(ldl_age, ldl_sex, ldl_cad) %>%
  filter(method %in% c("rank", "rank(all1)", "rank(all2)")) %>%
  mutate(
    method = factor(method, levels = c("rank", "rank(all1)", "rank(all2)")),
    Methods = recode(
      method,
      "rank"       = "Rank",
      "rank(all1)" = "Rank (FO)",
      "rank(all2)" = "Rank (HO)"
    ),
    Methods = factor(Methods, levels = c("Rank", "Rank (FO)", "Rank (HO)"))
  )

mse_norm <- dat %>%
  left_join(lmr_ref, by = c("exposure", "outcome")) %>%
  mutate(
    lace_ref = if_else(outcome == "CAD", lace_LMR, 0),
    sq_dev = (lace - lace_ref)^2
  ) %>%
  group_by(exposure, outcome, Methods) %>%
  summarise(
    mse = mean(sq_dev, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(exposure, outcome) %>%
  mutate(
    mse_rank = mse[Methods == "Rank"][1],
    mse_norm = mse / mse_rank
  ) %>%
  ungroup() %>%
  mutate(
    Exposure = "LDL-C",
    Outcome = dplyr::recode(
      outcome,
      "age" = "Age\n(Negative control)",
      "sex" = "Sex\n(Negative control)",
      "CAD" = "CAD\n(Positive control)"
    ),
    Outcome = factor(
      Outcome,
      levels = c(
        "Age\n(Negative control)",
        "Sex\n(Negative control)",
        "CAD\n(Positive control)"
      )
    )
  )


#========================================
# Plot settings
#========================================
method_colors <- c(
  "LMR"          = "#CD202CFF",
  "Rank"         = "#ffa600",
  "Rank (TDI)"   = "#c7522a",
  "Rank (age)"   = "#a98467",
  "Rank (sex)"   = "#ffa0c5",
  "Rank (FI)"    = "#32dba9",
  "Rank (MET)"   = "#476f95",
  "Rank (month)" = "#E74C3C",
  "Rank (FO)"    = "#893f71",
  "Rank (HO)"    = "#92ba92"
)


#========================================
# Plot
#========================================
plot.MSE_ldl <- ggplot(mse_norm, aes(x = Outcome, y = mse_norm, fill = Methods)) +
  geom_bar(
    stat = "identity",
    width = 0.85,
    position = position_dodge(width = 0.87)
  ) +
  labs(x = NULL, y = "Normalized MSE", title = NULL) +
  scale_fill_manual(values = method_colors, drop = FALSE) +
  theme_minimal() +
  theme(
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 10),
    axis.title.y = element_text(size = 12),
    legend.text  = element_text(size = 12),
    legend.title = element_text(size = 13)
  )

plot.MSE_ldl

#========================================
# Save figure
#========================================
ggsave(
  plot = plot.MSE_ldl,
  filename = "manuscript/figureS10.svg",
  device = "svg",
  scale = 0.9,
  width = 10,
  height = 5,
  units = "in"
)
