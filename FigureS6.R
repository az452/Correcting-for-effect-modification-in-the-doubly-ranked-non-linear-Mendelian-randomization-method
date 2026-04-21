
#===============================================================================
#
# Figure S6: simulation scenario 2: added condition, E is a collider
#
#===============================================================================

rm(list = ls())
# Set your working directory to the project folder before running this script.

library(ggplot2)
library(dplyr)

#========================================
# Load data
#========================================
load("rdata_nlmr_gxe_sim_s4_colliderUpdate_results_20260330.RData")


#========================================
# Data processing
#========================================
dat <- sim_results_colliderUpdate %>%
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

#========================================
# Plot settings
#========================================
method_pal <- c(
  "LMR"        = "#0072B5FF",
  "Rank"       = "#E18727FF",
  "Rank (GxE)" = "#E74C3CFF"
)

pos <- position_dodge(width = 0.65)

#========================================
# Plot
#========================================
plot <- ggplot(dat, aes(x = strata_label, y = lace, fill = Methods, colour = Methods)) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.5,
    color = "black"
  ) +
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
  filename = "manuscript/figureS6.svg",
  device = "svg",
  scale = 0.9,
  width = 10,
  height = 5,
  units = "in"
)
