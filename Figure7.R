#===============================================================================
#
# Figure 7: normalised MSE from the falsification test
#
#===============================================================================

rm(list = ls())
# Set your working directory to the project folder before running this script.

library(dplyr)
library(ggplot2)

#========================================
# Load data
#========================================

load("rdata_nlmr_vitd_negControl_results_20251105_focusedScore.RData")
vitd_simNeg.focus <- vitd_simNeg %>%
  mutate(exposure = "vitd[focus]")

load("rdata_nlmr_vitd_negControl_results_20251105.RData")
load("rdata_nlmr_bmi_negControl_results_20251206.RData")
load("rdata_nlmr_ldl_negControl_results_20251105.RData")

#========================================
# Data processing
#========================================

vitd <- vitd_simNeg %>%
  filter(method != "residual",
         method != "rank(season)") %>%
  mutate(
    strata_label = factor(strata, levels = 0:10, labels = c("Overall", 1:10)),
    method = factor(method, levels = c(
      "lmr", "rank", "rank(tdi)", "rank(age)", "rank(sex)",
      "rank(fi)", "rank(month)", "rank(all1)", "rank(all2)"
    )),
    Methods = recode(
      method,
      "lmr"         = "LMR",
      "rank"        = "Rank",
      "rank(tdi)"   = "Rank (TDI)",
      "rank(age)"   = "Rank (age)",
      "rank(sex)"   = "Rank (sex)",
      "rank(fi)"    = "Rank (FI)",
      "rank(month)" = "Rank (month)",
      "rank(all1)"  = "Rank (FO)",
      "rank(all2)"  = "Rank (HO)"
    ),
    Methods = factor(
      Methods,
      levels = c("LMR", "Rank", "Rank (TDI)", "Rank (age)", "Rank (sex)",
                 "Rank (FI)", "Rank (month)", "Rank (FO)", "Rank (HO)")
    )
  )

vitd.focus <- vitd_simNeg.focus %>%
  filter(method != "residual",
         method != "rank(season)") %>%
  mutate(
    strata_label = factor(strata, levels = 0:10, labels = c("Overall", 1:10)),
    method = factor(method, levels = c(
      "lmr", "rank", "rank(tdi)", "rank(age)", "rank(sex)",
      "rank(fi)", "rank(month)", "rank(all1)", "rank(all2)"
    )),
    Methods = recode(
      method,
      "lmr"         = "LMR",
      "rank"        = "Rank",
      "rank(tdi)"   = "Rank (TDI)",
      "rank(age)"   = "Rank (age)",
      "rank(sex)"   = "Rank (sex)",
      "rank(fi)"    = "Rank (FI)",
      "rank(month)" = "Rank (month)",
      "rank(all1)"  = "Rank (FO)",
      "rank(all2)"  = "Rank (HO)"
    ),
    Methods = factor(
      Methods,
      levels = c("LMR", "Rank", "Rank (TDI)", "Rank (age)", "Rank (sex)",
                 "Rank (FI)", "Rank (month)", "Rank (FO)", "Rank (HO)")
    )
  )

bmi <- bmi_simNeg %>%
  filter(method != "residual") %>%
  mutate(
    strata_label = factor(strata, levels = 0:10, labels = c("Overall", 1:10)),
    method = factor(method, levels = c(
      "lmr", "rank", "rank(tdi)", "rank(age)", "rank(sex)",
      "rank(fi)", "rank(met)", "rank(all1)", "rank(all2)"
    )),
    Methods = recode(
      method,
      "lmr"        = "LMR",
      "rank"       = "Rank",
      "rank(tdi)"  = "Rank (TDI)",
      "rank(age)"  = "Rank (age)",
      "rank(sex)"  = "Rank (sex)",
      "rank(fi)"   = "Rank (FI)",
      "rank(met)"  = "Rank (MET)",
      "rank(all1)" = "Rank (FO)",
      "rank(all2)" = "Rank (HO)"
    ),
    Methods = factor(
      Methods,
      levels = c("LMR", "Rank", "Rank (TDI)", "Rank (age)", "Rank (sex)",
                 "Rank (FI)", "Rank (MET)", "Rank (FO)", "Rank (HO)")
    )
  )

ldl <- ldl_simNeg %>%
  filter(method != "residual",
         method != "rank(all1-tdi)",
         method != "rank(all2-tdi)") %>%
  mutate(
    strata_label = factor(strata, levels = 0:10, labels = c("Overall", 1:10)),
    method = factor(method, levels = c(
      "lmr", "rank", "rank(tdi)", "rank(age)", "rank(sex)",
      "rank(fi)", "rank(all1)", "rank(all2)"
    )),
    Methods = recode(
      method,
      "lmr"        = "LMR",
      "rank"       = "Rank",
      "rank(tdi)"  = "Rank (TDI)",
      "rank(age)"  = "Rank (age)",
      "rank(sex)"  = "Rank (sex)",
      "rank(fi)"   = "Rank (FI)",
      "rank(all1)" = "Rank (FO)",
      "rank(all2)" = "Rank (HO)"
    ),
    Methods = factor(
      Methods,
      levels = c("LMR", "Rank", "Rank (TDI)", "Rank (age)", "Rank (sex)",
                 "Rank (FI)", "Rank (FO)", "Rank (HO)")
    )
  )

#========================================
# MSE by replicate and method
#========================================

vitd.mse <- vitd_simNeg %>%
  filter(method != "residual",
         method != "lmr",
         method != "rank(season)") %>%
  group_by(exposure, outcome, replicate, method) %>%
  summarise(
    mse = mean(lace^2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    method = factor(method, levels = c(
      "lmr", "rank", "rank(tdi)", "rank(age)", "rank(sex)",
      "rank(fi)", "rank(month)", "rank(all1)", "rank(all2)"
    )),
    Methods = recode(
      method,
      "lmr"         = "LMR",
      "rank"        = "Rank",
      "rank(tdi)"   = "Rank (TDI)",
      "rank(age)"   = "Rank (age)",
      "rank(sex)"   = "Rank (sex)",
      "rank(fi)"    = "Rank (FI)",
      "rank(month)" = "Rank (month)",
      "rank(all1)"  = "Rank (FO)",
      "rank(all2)"  = "Rank (HO)"
    ),
    Methods = factor(
      Methods,
      levels = c("LMR", "Rank", "Rank (TDI)", "Rank (age)", "Rank (sex)",
                 "Rank (FI)", "Rank (month)", "Rank (FO)", "Rank (HO)")
    )
  )

vitd.focus.mse <- vitd_simNeg.focus %>%
  filter(method != "residual",
         method != "lmr",
         method != "rank(season)") %>%
  group_by(exposure, outcome, replicate, method) %>%
  summarise(
    mse = mean(lace^2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    method = factor(method, levels = c(
      "lmr", "rank", "rank(tdi)", "rank(age)", "rank(sex)",
      "rank(fi)", "rank(month)", "rank(all1)", "rank(all2)"
    )),
    Methods = recode(
      method,
      "lmr"         = "LMR",
      "rank"        = "Rank",
      "rank(tdi)"   = "Rank (TDI)",
      "rank(age)"   = "Rank (age)",
      "rank(sex)"   = "Rank (sex)",
      "rank(fi)"    = "Rank (FI)",
      "rank(month)" = "Rank (month)",
      "rank(all1)"  = "Rank (FO)",
      "rank(all2)"  = "Rank (HO)"
    ),
    Methods = factor(
      Methods,
      levels = c("LMR", "Rank", "Rank (TDI)", "Rank (age)", "Rank (sex)",
                 "Rank (FI)", "Rank (month)", "Rank (FO)", "Rank (HO)")
    )
  )

bmi.mse <- bmi_simNeg %>%
  filter(method != "residual",
         method != "lmr") %>%
  group_by(exposure, outcome, replicate, method) %>%
  summarise(
    mse = mean(lace^2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    method = factor(method, levels = c(
      "rank", "rank(tdi)", "rank(age)", "rank(sex)",
      "rank(fi)", "rank(met)", "rank(all1)", "rank(all2)"
    )),
    Methods = recode(
      method,
      "rank"       = "Rank",
      "rank(tdi)"  = "Rank (TDI)",
      "rank(age)"  = "Rank (age)",
      "rank(sex)"  = "Rank (sex)",
      "rank(fi)"   = "Rank (FI)",
      "rank(met)"  = "Rank (MET)",
      "rank(all1)" = "Rank (FO)",
      "rank(all2)" = "Rank (HO)"
    ),
    Methods = factor(
      Methods,
      levels = c("Rank", "Rank (TDI)", "Rank (age)", "Rank (sex)",
                 "Rank (FI)", "Rank (MET)", "Rank (FO)", "Rank (HO)")
    )
  )

ldl.mse <- ldl_simNeg %>%
  filter(method != "residual",
         method != "lmr",
         method != "rank(all1-tdi)",
         method != "rank(all2-tdi)") %>%
  group_by(exposure, outcome, replicate, method) %>%
  summarise(
    mse = mean(lace^2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    method = factor(method, levels = c(
      "rank", "rank(tdi)", "rank(age)", "rank(sex)",
      "rank(fi)", "rank(all1)", "rank(all2)"
    )),
    Methods = recode(
      method,
      "rank"       = "Rank",
      "rank(tdi)"  = "Rank (TDI)",
      "rank(age)"  = "Rank (age)",
      "rank(sex)"  = "Rank (sex)",
      "rank(fi)"   = "Rank (FI)",
      "rank(all1)" = "Rank (FO)",
      "rank(all2)" = "Rank (HO)"
    ),
    Methods = factor(
      Methods,
      levels = c("Rank", "Rank (TDI)", "Rank (age)", "Rank (sex)",
                 "Rank (FI)", "Rank (FO)", "Rank (HO)")
    )
  )

#========================================
# Combine and normalise
#========================================

all.mse <- bind_rows(vitd.mse, vitd.focus.mse, bmi.mse, ldl.mse) %>%
  mutate(
    Exposure = recode(
      exposure,
      "vitd"        = "25(OH)D [Replicated]",
      "vitd[focus]" = "25(OH)D [Focused]",
      "bmi"         = "BMI",
      "ldl"         = "LDL-C"
    ),
    Exposure = factor(
      Exposure,
      levels = c("25(OH)D [Replicated]", "25(OH)D [Focused]", "BMI", "LDL-C")
    ),
    Methods = factor(
      Methods,
      levels = c("LMR", "Rank", "Rank (TDI)", "Rank (age)", "Rank (sex)",
                 "Rank (FI)", "Rank (MET)", "Rank (month)", "Rank (FO)", "Rank (HO)")
    )
  )

baseline_rank <- all.mse %>%
  filter(Methods == "Rank") %>%
  select(exposure, outcome, replicate, mse_rank = mse)

all.mse <- all.mse %>%
  left_join(baseline_rank, by = c("exposure", "outcome", "replicate")) %>%
  mutate(mse_norm = mse / mse_rank) %>%
  select(-mse_rank)

#========================================
# Plot data
#========================================

plot_dat <- all.mse %>%
  filter(Methods != "Rank") %>%
  mutate(
    Exposure_spaced = factor(
      Exposure,
      levels = c(
        "25(OH)D [Replicated]",
        "gap0",
        "25(OH)D [Focused]",
        "gap1",
        "BMI",
        "gap2",
        "LDL-C"
      )
    )
  )

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

plot_mse_norm <- ggplot(
  plot_dat,
  aes(
    x = Exposure_spaced,
    y = mse_norm,
    fill = Methods,
    colour = Methods
  )
) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_boxplot(
    width = 0.45,
    outlier.size = 0.6,
    position = position_dodge(width = 0.90),
    show.legend = TRUE
  ) +
  stat_summary(
    fun = median,
    geom = "crossbar",
    width = 0.5,
    colour = "white",
    linewidth = 0.3,
    position = position_dodge(width = 0.90),
    show.legend = FALSE
  ) +
  labs(
    x = NULL,
    y = "Normalized MSE"
  ) +
  scale_x_discrete(
    drop = FALSE,
    labels = c(
      "25(OH)D [Replicated]",
      "",
      "25(OH)D [Focused]",
      "",
      "BMI",
      "",
      "LDL-C"
    )
  ) +
  scale_y_continuous(
    limits = c(0, 1.3),
    breaks = seq(0, 1.2, by = 0.2)
  ) +
  scale_fill_manual(values = method_colors) +
  scale_color_manual(values = method_colors) +
  theme_minimal() +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

plot_mse_norm

#========================================
# Save figure
#========================================

ggsave(
  plot = plot_mse_norm,
  filename = "manuscript/figure7.svg",
  device = "svg",
  scale = 0.9,
  width = 11,
  height = 5.5,
  units = "in"
)
