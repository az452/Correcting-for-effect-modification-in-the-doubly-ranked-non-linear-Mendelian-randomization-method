
#===============================================================================
#
# Figure S2: MSE for nonlinear X-Y in simulation scenario 2
#
#===============================================================================


rm(list = ls())
# Set your working directory to the project folder before running this script.

library(dplyr)
library(ggplot2)

#========================================
# Load data
#========================================
load("rdata_nlmr_gxe_sim_s2_strategy_results_20251016.RData")
load("rdata_nlmr_gxe_sim_s2_strategy_non-linear_target_20260324.RData")


#========================================
# Prepare LACE results
#========================================
dat_base <- sim_nonlinear %>%
  filter(
    method %in% c("rank+X", "rank+(X-GV)"),
    strata != 0
  ) %>%
  mutate(
    Methods = recode(
      method,
      "rank+X"      = "Rank",
      "rank+(X-GV)" = "Rank (GxE)"
    ),
    Methods = factor(Methods, levels = c("Rank", "Rank (GxE)"))
  )

#========================================
# Prepare targets
#========================================
target_by_stratum <- target_results %>%
  filter(
    method %in% c("rank+X", "rank+(X-GV)"),
    strata != 0
  ) %>%
  mutate(
    Methods = recode(
      method,
      "rank+X"      = "Rank",
      "rank+(X-GV)" = "Rank (GxE)"
    ),
    Methods = factor(Methods, levels = c("Rank", "Rank (GxE)"))
  ) %>%
  group_by(Methods, strata) %>%
  summarise(
    target = median(target, na.rm = TRUE),
    .groups = "drop"
  )

#========================================
# Compute MSE
#========================================
mse_rep <- dat_base %>%
  left_join(target_by_stratum, by = c("Methods", "strata")) %>%
  mutate(
    sq_error = (lace - target)^2
  ) %>%
  group_by(rep, Methods) %>%
  summarise(
    mse = mean(sq_error, na.rm = TRUE),
    .groups = "drop"
  )


#========================================
# Plot settings
#========================================
method_pal <- c(
  "Rank"       = "#E18727FF",
  "Rank (GxE)" = "#E74C3CFF"
)


#========================================
# Plot
#========================================

plot_mse <- ggplot(
  mse_rep,
  aes(x = Methods, y = mse, fill = Methods, colour = Methods)
) +
  geom_boxplot(
    width = 0.45,
    outlier.size = 0.6,
    show.legend = FALSE
  ) +
  stat_summary(
    fun = median,
    geom = "crossbar",
    width = 0.50,
    colour = "white",
    linewidth = 0.3,
    show.legend = FALSE
  ) +
  labs(
    x = NULL,
    y = "MSE"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 14)
  ) +
  scale_fill_manual(values = method_pal, drop = FALSE) +
  scale_color_manual(values = method_pal, drop = FALSE) +
  scale_y_continuous(
    limits = c(0, 0.003),
    breaks = seq(0, 0.003, by = 0.0005)
  )

plot_mse

#========================================
# Save figure
#========================================
ggsave(
  plot = plot_mse,
  filename = "manuscript/figureS2.svg",
  device = "svg",
  scale = 0.9,
  width = 6,
  height = 6,
  units = "in"
)
