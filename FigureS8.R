

#=====================================================================================
#
# Figure S8: results from VitD empirical negative controls using replicated score
#
#=====================================================================================

rm(list = ls())
# Set your working directory to the project folder before running this script.

library(dplyr)
library(ggplot2)
library(ggpubr)


#========================================
# Load data
#========================================
load("rdata_nlmr_vitd_EmpNegPosControls_results_20251205.RData")


#========================================
# Data processing
#========================================

# Age
dat.age <- vitd_age %>%
  filter(method %in% c("lmr", "rank", "rank(month)", "rank(all1)")) %>%
  mutate(
    strata_lab = factor(strata, levels = 0:10, labels = c("Overall", 1:10)),
    bxlci    = bx   - qnorm(0.975) * bxse,
    bxuci    = bx   + qnorm(0.975) * bxse,
    bylci    = by   - qnorm(0.975) * byse,
    byuci    = by   + qnorm(0.975) * byse,
    lace_lci = lace - qnorm(0.975) * lace_se,
    lace_uci = lace + qnorm(0.975) * lace_se,
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

# Sex
dat.sex <- vitd_sex %>%
  filter(method %in% c("lmr", "rank", "rank(month)", "rank(all1)")) %>%
  mutate(
    strata_lab = factor(strata, levels = 0:10, labels = c("Overall", 1:10)),
    bxlci    = bx   - qnorm(0.975) * bxse,
    bxuci    = bx   + qnorm(0.975) * bxse,
    bylci    = by   - qnorm(0.975) * byse,
    byuci    = by   + qnorm(0.975) * byse,
    lace_lci = lace - qnorm(0.975) * lace_se,
    lace_uci = lace + qnorm(0.975) * lace_se,
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

pos <- position_dodge(width = 0.5)


#========================================
# Plot A: Age
#========================================
lci <- dat.age %>%
  filter(method == "lmr") %>%
  pull(lace_lci)

uci <- dat.age %>%
  filter(method == "lmr") %>%
  pull(lace_uci)

plot.age <- ggplot(dat.age, aes(x = strata_lab, y = lace, color = Methods)) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.7,
    color = "black"
  ) +
  geom_point(size = 2.0, position = pos) +
  geom_linerange(
    aes(ymin = lace_lci, ymax = lace_uci),
    linewidth = 0.8,
    position = pos,
    show.legend = FALSE
  ) +
  annotate(
    "rect",
    xmin = -Inf, xmax = Inf,
    ymin = lci, ymax = uci,
    fill = "grey60", alpha = 0.2
  ) +
  labs(x = "Strata", y = "Age\n(Negative control)") +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.text  = element_text(size = 13),
    legend.title = element_text(size = 14),
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  ) +
  scale_color_manual(values = method_colors)

plot.age

#========================================
# Plot B: Sex
#========================================
or_lci <- dat.sex %>%
  filter(method == "lmr") %>%
  pull(lace_lci) %>%
  exp()

or_uci <- dat.sex %>%
  filter(method == "lmr") %>%
  pull(lace_uci) %>%
  exp()

plot.sex <- ggplot(dat.sex, aes(x = strata_lab, y = exp(lace), color = Methods)) +
  geom_hline(
    yintercept = exp(0),
    linetype = "dashed",
    linewidth = 0.7,
    color = "black"
  ) +
  geom_point(size = 2.0, position = pos) +
  geom_linerange(
    aes(ymin = exp(lace_lci), ymax = exp(lace_uci)),
    linewidth = 0.8,
    position = pos
  ) +
  annotate(
    "rect",
    xmin = -Inf, xmax = Inf,
    ymin = or_lci, ymax = or_uci,
    fill = "grey60", alpha = 0.2
  ) +
  coord_trans(y = "log") +
  labs(x = "Strata", y = "Odds ratio for being a male\n(Negative control)") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  ) +
  scale_color_manual(values = method_colors)

plot.sex

#========================================
# Combine plots
#========================================
comb <- ggarrange(
  plot.age,
  plot.sex,
  ncol = 1,
  heights = c(10, 10),
  align = "v",
  labels = c("A", "B"),
  font.label = list(size = 14, face = "bold"),
  common.legend = TRUE,
  legend = "right"
)

comb



#========================================
# Save figure
#========================================
ggsave(
  plot = comb,
  filename = "manuscript/figureS8.svg",
  device = "svg",
  scale = 0.9,
  width = 16,
  height = 9,
  units = "in"
)
