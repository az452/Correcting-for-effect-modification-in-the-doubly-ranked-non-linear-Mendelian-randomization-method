#===============================================================================
#
# Figure 3: results from simulation scenario 2
#
#===============================================================================

rm(list = ls())
# Set your working directory to the project folder before running this script.
library(ggplot2)
library(dplyr)
library(ggpubr)



#========================================
# Load data
#========================================
load("rdata_nlmr_gxe_sim_s2_strategy_results_20251016.RData")
load("rdata_nlmr_gxe_sim_s2_strategy_non-linear_target_20260324.RData")



#========================================
# Data processing
#========================================
dat1 <- sim_null %>%
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

dat2 <- sim_linear %>%
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

dat3 <- sim_nonlinear %>%
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
# Plot 1: null
#========================================
plot1 <- ggplot(dat1, aes(x = strata_label, y = lace, fill = Methods, colour = Methods)) +
  geom_hline(
    yintercept = 0,
    linetype   = "dashed",
    linewidth  = 0.5,
    color      = "black"
  ) +
  geom_boxplot(
    width        = 0.40,
    position     = pos,
    outlier.size = 0.6,
    show.legend  = TRUE
  ) +
  stat_summary(
    fun         = median,
    geom        = "crossbar",
    width       = 0.50,
    colour      = "white",
    linewidth   = 0.3,
    position    = pos,
    show.legend = FALSE
  ) +
  labs(x = "Strata", y = "LACE estimates") +
  theme_minimal() +
  theme(legend.position = "right") +
  scale_fill_manual(values = method_pal, drop = FALSE) +
  scale_color_manual(values = method_pal, drop = FALSE)

plot1

#========================================
# Plot 2: linear
#========================================
bx0 <- dat2 %>%
  filter(strata == 0) %>%
  pull(bx) %>%
  median()

plot2 <- ggplot(dat2, aes(x = strata_label, y = lace, fill = Methods, colour = Methods)) +
  geom_hline(
    yintercept = bx0,
    linetype   = "dashed",
    linewidth  = 0.5,
    color      = "black"
  ) +
  geom_boxplot(
    width        = 0.40,
    position     = pos,
    outlier.size = 0.6,
    show.legend  = TRUE
  ) +
  stat_summary(
    fun         = median,
    geom        = "crossbar",
    width       = 0.50,
    colour      = "white",
    linewidth   = 0.3,
    position    = pos,
    show.legend = FALSE
  ) +
  labs(x = "Strata", y = "LACE estimates") +
  theme_minimal() +
  theme(legend.position = "right") +
  scale_fill_manual(values = method_pal, drop = FALSE) +
  scale_color_manual(values = method_pal, drop = FALSE)

plot2

#========================================
# Plot 3: nonlinear
#========================================

# Target estimates for Rank
laceTarget1c <- target_results %>%
  filter(method == "rank+X") %>%
  group_by(strata) %>%
  summarise(
    laceTarget = median(target, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    strata_label = factor(strata, levels = 1:10, labels = c(1:10))
  )

# Target estimates for Rank (GxE)
laceTarget2c <- target_results %>%
  filter(method == "rank+(X-GV)") %>%
  group_by(strata) %>%
  summarise(
    laceTarget = median(target, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    strata_label = factor(strata, levels = 1:10, labels = c(1:10))
  )

offset <- 0.18

plot3 <- ggplot(dat3, aes(x = strata_label, y = lace, fill = Methods, colour = Methods)) +
  geom_boxplot(
    width        = 0.40,
    position     = pos,
    outlier.size = 0.6,
    show.legend  = TRUE
  ) +
  stat_summary(
    fun         = median,
    geom        = "crossbar",
    width       = 0.50,
    colour      = "white",
    linewidth   = 0.3,
    position    = pos,
    show.legend = FALSE
  ) +
  geom_point(
    data = laceTarget1c,
    aes(x = strata_label, y = laceTarget),
    position = position_nudge(x = -offset),
    color = "red",
    shape = 16,
    size = 2,
    inherit.aes = FALSE
  ) +
  geom_point(
    data = laceTarget2c,
    aes(x = strata_label, y = laceTarget),
    position = position_nudge(x =  offset),
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

plot3

#========================================
# Combine plots
#========================================
comb <- ggarrange(
  plot1,
  plot2,
  plot3,
  ncol = 1,
  heights = c(10, 10, 10),
  align = "v",
  labels = c("A", "B", "C"),
  font.label = list(size = 10, face = "bold"),
  common.legend = TRUE,
  legend = "right"
)

comb

#========================================
# Save figure
#========================================
ggsave(
  comb,
  filename = "manuscript/figure3.svg",
  scale = 0.9,
  width = 15,
  height = 10,
  units = "in"
)
