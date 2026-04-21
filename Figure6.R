
#========================================
# Figure 6: falsification tests
#========================================

rm(list = ls())
# Set your working directory to the project folder before running this script.

library(dplyr)
library(ggplot2)
library(ggpubr)

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
  "Rank (met)"   = "#476f95",
  "Rank (month)" = "#E74C3C",
  "Rank (FO)"    = "#893f71",
  "Rank (HO)"    = "#92ba92"
)

pos <- position_dodge(width = 0.55)

#========================================
# Panel A and B: Vitamin D
#========================================

load("rdata_nlmr_vitd_negControl_results_20251105_focusedScore.RData")
vitd_simNeg.focus <- vitd_simNeg

load("rdata_nlmr_vitd_negControl_results_20251105.RData")

vitd <- vitd_simNeg %>%
  filter(
    method %in% c("lmr", "rank", "rank(month)", "rank(all1)")
  ) %>%
  mutate(
    strata_label = factor(strata, levels = 0:10, labels = c("Overall", 1:10)),
    method = factor(method, levels = c("lmr", "rank", "rank(month)", "rank(all1)")),
    Methods = recode(
      method,
      "lmr"         = "LMR",
      "rank"        = "Rank",
      "rank(month)" = "Rank (month)",
      "rank(all1)"  = "Rank (FO)"
    ),
    Methods = factor(Methods, levels = c("LMR", "Rank", "Rank (month)", "Rank (FO)"))
  )

vitd.focus <- vitd_simNeg.focus %>%
  filter(
    method %in% c("lmr", "rank", "rank(month)", "rank(all1)")
  ) %>%
  mutate(
    strata_label = factor(strata, levels = 0:10, labels = c("Overall", 1:10)),
    method = factor(method, levels = c("lmr", "rank", "rank(month)", "rank(all1)")),
    Methods = recode(
      method,
      "lmr"         = "LMR",
      "rank"        = "Rank",
      "rank(month)" = "Rank (month)",
      "rank(all1)"  = "Rank (FO)"
    ),
    Methods = factor(Methods, levels = c("LMR", "Rank", "Rank (month)", "Rank (FO)"))
  )

vitd.plot <- ggplot(
  vitd,
  aes(x = strata_label, y = lace, fill = Methods, colour = Methods)
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.4,
    color = "black"
  ) +
  geom_boxplot(
    width = 0.4,
    position = pos,
    outlier.size = 0.6,
    show.legend = TRUE
  ) +
  stat_summary(
    fun = median,
    geom = "crossbar",
    width = 0.5,
    colour = "white",
    linewidth = 0.3,
    position = pos,
    show.legend = FALSE
  ) +
  labs(x = "Strata", y = "Sim -ve") +
  theme_minimal() +
  theme(legend.position = "right") +
  scale_fill_manual(values = method_colors) +
  scale_color_manual(values = method_colors)

vitd.focus.plot <- ggplot(
  vitd.focus,
  aes(x = strata_label, y = lace, fill = Methods, colour = Methods)
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.4,
    color = "black"
  ) +
  geom_boxplot(
    width = 0.4,
    position = pos,
    outlier.size = 0.6,
    show.legend = TRUE
  ) +
  stat_summary(
    fun = median,
    geom = "crossbar",
    width = 0.5,
    colour = "white",
    linewidth = 0.3,
    position = pos,
    show.legend = FALSE
  ) +
  labs(x = "Strata", y = "Sim -ve") +
  theme_minimal() +
  theme(legend.position = "right") +
  scale_fill_manual(values = method_colors) +
  scale_color_manual(values = method_colors)

comb_vitd <- ggarrange(
  vitd.plot + theme(legend.position = "right"),
  vitd.focus.plot + theme(legend.position = "none"),
  ncol = 1,
  heights = c(10, 10),
  align = "v",
  labels = c("A", "B"),
  font.label = list(size = 10, face = "bold"),
  common.legend = TRUE,
  legend = "right"
)

#========================================
# Panel C: BMI
#========================================

load("rdata_nlmr_bmi_negControl_results_20251105.RData")

bmi <- bmi_simNeg %>%
  filter(
    method %in% c("lmr", "rank", "rank(fi)", "rank(all1)")
  ) %>%
  mutate(
    strata_label = factor(strata, levels = 0:10, labels = c("Overall", 1:10)),
    method = factor(method, levels = c("lmr", "rank", "rank(fi)", "rank(all1)")),
    Methods = recode(
      method,
      "lmr"        = "LMR",
      "rank"       = "Rank",
      "rank(fi)"   = "Rank (FI)",
      "rank(all1)" = "Rank (FO)"
    ),
    Methods = factor(Methods, levels = c("LMR", "Rank", "Rank (FI)", "Rank (FO)"))
  )

bmi.plot <- ggplot(
  bmi,
  aes(x = strata_label, y = lace, fill = Methods, colour = Methods)
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.4,
    color = "black"
  ) +
  geom_boxplot(
    width = 0.4,
    position = pos,
    outlier.size = 0.6,
    show.legend = TRUE
  ) +
  stat_summary(
    fun = median,
    geom = "crossbar",
    width = 0.5,
    colour = "white",
    linewidth = 0.3,
    position = pos,
    show.legend = FALSE
  ) +
  labs(x = "Strata", y = "Sim -ve") +
  theme_minimal() +
  theme(legend.position = "right") +
  scale_fill_manual(values = method_colors) +
  scale_color_manual(values = method_colors)

#========================================
# Panel D: LDL
#========================================

load("rdata_nlmr_ldl_negControl_results_20251105.RData")

ldl <- ldl_simNeg %>%
  filter(
    method %in% c("lmr", "rank", "rank(all1)", "rank(all2)")
  ) %>%
  mutate(
    strata_label = factor(strata, levels = 0:10, labels = c("Overall", 1:10)),
    method = factor(method, levels = c("lmr", "rank", "rank(all1)", "rank(all2)")),
    Methods = recode(
      method,
      "lmr"        = "LMR",
      "rank"       = "Rank",
      "rank(all1)" = "Rank (FO)",
      "rank(all2)" = "Rank (HO)"
    ),
    Methods = factor(Methods, levels = c("LMR", "Rank", "Rank (FO)", "Rank (HO)"))
  )

ldl.plot <- ggplot(
  ldl,
  aes(x = strata_label, y = lace, fill = Methods, colour = Methods)
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.4,
    color = "black"
  ) +
  geom_boxplot(
    width = 0.4,
    position = pos,
    outlier.size = 0.6,
    show.legend = TRUE
  ) +
  stat_summary(
    fun = median,
    geom = "crossbar",
    width = 0.5,
    colour = "white",
    linewidth = 0.3,
    position = pos,
    show.legend = FALSE
  ) +
  labs(x = "Strata", y = "Sim -ve") +
  theme_minimal() +
  theme(legend.position = "right") +
  scale_fill_manual(values = method_colors) +
  scale_color_manual(values = method_colors)

#========================================
# Combine plots
#========================================

bmi_ldl_comb <- ggarrange(
  bmi.plot,
  ldl.plot,
  ncol = 1,
  heights = c(10, 10),
  align = "v",
  labels = c("C", "D"),
  font.label = list(size = 10, face = "bold"),
  common.legend = FALSE
)

plot_combined <- ggarrange(
  comb_vitd,
  bmi_ldl_comb,
  ncol = 1,
  widths = c(1, 1),
  align = "h"
)

plot_combined



#========================================
# Save figure
#========================================
ggsave(
  plot = plot_combined,
  filename = "manuscript/figure6.svg",
  device = "svg",
  scale = 0.9,
  width = 16,
  height = 14,
  units = "in"
)
