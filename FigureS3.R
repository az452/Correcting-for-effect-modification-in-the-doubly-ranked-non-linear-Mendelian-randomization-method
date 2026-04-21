
#======================================================================================
#
# Figure S3: additional nonlinear form of (Threshold effect) for simulation scenario 2
#
#======================================================================================


rm(list = ls())
# Set your working directory to the project folder before running this script.

library(ggplot2)
library(dplyr)

#========================================
# Load data
#========================================
load("rdata_nlmr_gxe_sim_s4_thresholdUpdate_results_20260322.RData")
load("rdata_nlmr_gxe_sim_s4_thresholdUpdate_target_20260322.RData")

#========================================
# Data processing
#========================================
dat <- sim_results_thresholdUpdate %>%
  filter(method != "residual") %>%
  mutate(
    strata_label = factor(strata, levels = 0:10, labels = c("Overall", 1:10)),
    method = factor(method, levels = c("lmr", "rank+X", "rank+(X-GV)")),
    Methods = recode(
      method,
      "lmr"         = "LMR",
      "rank+X"      = "Rank",
      "rank+(X-GV)" = "Rank (GxE)"
    ),
    Methods = factor(Methods, levels = c("LMR", "Rank", "Rank (GxE)"))
  )

laceTarget_rank <- target_results %>%
  filter(method == "rank+X") %>%
  group_by(strata) %>%
  summarise(
    laceTarget = median(target, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    strata_label = factor(strata, levels = 1:10, labels = c(1:10))
  )

laceTarget_gxe <- target_results %>%
  filter(method == "rank+(X-GV)") %>%
  group_by(strata) %>%
  summarise(
    laceTarget = median(target, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    strata_label = factor(strata, levels = 1:10, labels = c(1:10))
  )

#========================================
# Plot settings
#========================================
method_pal <- c(
  "LMR"        = "#0072B5FF",
  "Rank"       = "#E18727FF",
  "Rank (GxE)" = "#E74C3CFF"
)

pos <- position_dodge(width = 0.65)
offset <- 0.18

#========================================
# Plot
#========================================
plot <- ggplot(dat, aes(x = strata_label, y = lace, fill = Methods, colour = Methods)) +
  geom_boxplot(
    width = 0.40,
    position = pos,
    outlier.size = 0.6,
    show.legend = TRUE
  ) +
  stat_summary(
    fun = median,
    geom = "crossbar",
    width = 0.50,
    colour = "white",
    linewidth = 0.3,
    position = pos,
    show.legend = FALSE
  ) +
  geom_point(
    data = laceTarget_rank,
    aes(x = strata_label, y = laceTarget),
    position = position_nudge(x = -offset),
    color = "red",
    shape = 16,
    size = 2,
    inherit.aes = FALSE
  ) +
  geom_point(
    data = laceTarget_gxe,
    aes(x = strata_label, y = laceTarget),
    position = position_nudge(x = offset),
    color = "blue",
    shape = 16,
    size = 2,
    inherit.aes = FALSE
  ) +
  labs(x = "Strata", y = "LACE estimates") +
  theme_minimal() +
  theme(legend.position = "right") +
  scale_fill_manual(values = method_pal, drop = FALSE) +
  scale_color_manual(values = method_pal, drop = FALSE)

plot

#========================================
# Save figure
#========================================
ggsave(
  plot = plot,
  filename = "manuscript/figureS3.svg",
  device = "svg",
  scale = 0.9,
  width = 10,
  height = 5,
  units = "in"
)
