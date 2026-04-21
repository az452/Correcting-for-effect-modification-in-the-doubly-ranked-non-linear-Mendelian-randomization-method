

#=====================================================================================
#
# Figure S5: supplementary simulation scenario 2: transformation of E (log and exp)
#
#=====================================================================================

rm(list = ls())
# Set your working directory to the project folder before running this script.

library(ggplot2)
library(dplyr)

#========================================
# Load data
#========================================
load("rdata_nlmr_gxe_sim_s2_strategy_results_20251016.RData")
load("rdata_nlmr_gxe_sim_s4_logE_results_20260320.RData")
load("rdata_nlmr_gxe_sim_s4_expE_results_20260325.RData")

#========================================
# Data processing
#========================================
dat1 <- sim_null %>%
  filter(method != "residual")

dat2 <- sim_results_logE %>%
  filter(method == "rank+(X-GV)") %>%
  mutate(method = "rank+log(X-GV)")

dat3 <- sim_results_expE %>%
  filter(method == "rank+(X-GV)") %>%
  mutate(method = "rank+exp(X-GV)")

dat <- bind_rows(dat1, dat2, dat3) %>%
  mutate(
    strata_label = factor(strata, levels = 0:10, labels = c("Overall", 1:10)),
    method = factor(
      method,
      levels = c("lmr", "rank+X", "rank+(X-GV)", "rank+log(X-GV)", "rank+exp(X-GV)")
    ),
    Methods = recode(
      method,
      "lmr"            = "LMR",
      "rank+X"         = "Rank",
      "rank+(X-GV)"    = "Rank (GxE)",
      "rank+log(X-GV)" = "Rank (GxE,log)",
      "rank+exp(X-GV)" = "Rank (GxE,exp)"
    ),
    Methods = factor(
      Methods,
      levels = c("LMR", "Rank", "Rank (GxE)", "Rank (GxE,log)", "Rank (GxE,exp)")
    )
  )

#========================================
# Plot settings
#========================================
method_pal <- c(
  "LMR"            = "#0072B5FF",
  "Rank"           = "#E18727FF",
  "Rank (GxE)"     = "#E74C3CFF",
  "Rank (GxE,log)" = "#EC7063FF",
  "Rank (GxE,exp)" = "#F1948AFF"
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
  filename = "manuscript/figureS5.svg",
  device = "svg",
  scale = 0.9,
  width = 16,
  height = 5,
  units = "in"
)
