
#==========================================================================================
#
# Figure S4: Supplementary simulation scenario 1: imperfectly measured effect modifier
#
#==========================================================================================

rm(list = ls())
# Set your working directory to the project folder before running this script.

library(ggplot2)
library(dplyr)

#========================================
# Load data
#========================================
load("rdata_nlmr_gxe_sim_s4_imperfectE_results_20260320.RData")

#========================================
# Data processing
#========================================
dat <- sim_results_imperfectE %>%
  filter(method != "residual") %>%
  mutate(
    strata_label = factor(strata, levels = 0:10, labels = c("Overall", 1:10)),
    method = factor(
      method,
      levels = c("lmr", "rank+X", "rank+X1", "rank+X2", "rank+X3", "rank+X4")
    ),
    Methods = recode(
      method,
      "lmr"     = "LMR",
      "rank+X"  = "Rank",
      "rank+X1" = "Rank (GxE,1)",
      "rank+X2" = "Rank (GxE,0.8)",
      "rank+X3" = "Rank (GxE,0.4)",
      "rank+X4" = "Rank (GxE,0)"
    ),
    Methods = factor(
      Methods,
      levels = c(
        "LMR",
        "Rank",
        "Rank (GxE,1)",
        "Rank (GxE,0.8)",
        "Rank (GxE,0.4)",
        "Rank (GxE,0)"
      )
    )
  )

#========================================
# Plot settings
#========================================
method_pal <- c(
  "LMR"            = "#0072B5FF",
  "Rank"           = "#E18727FF",
  "Rank (GxE,1)"   = "#E74C3CFF",
  "Rank (GxE,0.8)" = "#EC7063FF",
  "Rank (GxE,0.4)" = "#F1948AFF",
  "Rank (GxE,0)"   = "#F5B7B1FF"
)

pos <- position_dodge(width = 0.65)

#========================================
# Plot
#========================================
plot <- ggplot(dat, aes(x = strata_label, y = lace, fill = Methods, colour = Methods)) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.2,
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
  filename = "manuscript/figureS4.svg",
  device = "svg",
  scale = 0.9,
  width = 16,
  height = 5,
  units = "in"
)
