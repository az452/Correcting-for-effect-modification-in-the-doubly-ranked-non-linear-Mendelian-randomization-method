
#===============================================================================
#
# Figure 5: effect modification of the genetic score-exposure association
#
#===============================================================================

# Set your working directory to the project folder before running this script.
rm(list = ls())

library(dplyr)
library(forestplot)
library(grid)
library(ggpubr)
library(ggplotify)

#========================================
# Panel A: Vitamin D (replicated score)
#========================================
load("rdata_vitd_interaction_results_20250405_grsSD.RData")

tdi_results <- vitd.tdi$strata_results
tdi_results$mod <- "TDI"

age_results <- vitd.age$strata_results
age_results$mod <- "Age"

sex_results <- vitd.sex$strata_results
sex_results$mod <- "Sex"

fi_results <- vitd.fi$strata_results
fi_results$mod <- "FI"

season_results <- vitd.month$strata_results
season_results$mod <- "Month"

dat_strata <- bind_rows(
  tdi_results,
  age_results,
  sex_results,
  fi_results,
  season_results
) %>%
  mutate(
    blci = b - qnorm(0.975) * se,
    buci = b + qnorm(0.975) * se
  ) %>%
  mutate(
    mod = factor(mod, levels = c("Age", "FI", "Sex", "TDI", "Month"))
  )

dat_pval <- bind_rows(
  data.frame(mod = "TDI",   pval = vitd.tdi$pval_gv),
  data.frame(mod = "Age",   pval = vitd.age$pval_gv),
  data.frame(mod = "Sex",   pval = vitd.sex$pval_gv),
  data.frame(mod = "FI",    pval = vitd.fi$pval_gv),
  data.frame(mod = "Month", pval = vitd.month$pval_gv)
)

dat_strata <- merge(dat_strata, dat_pval, by = "mod", all.x = TRUE) %>%
  mutate(
    season_sort = case_when(
      strata == "Winter" ~ 1,
      strata == "Spring" ~ 2,
      strata == "Summer" ~ 3,
      strata == "Autumn" ~ 4,
      TRUE ~ 99
    )
  ) %>%
  arrange(mod, season_sort, strata) %>%
  group_by(mod) %>%
  mutate(P_inter = if_else(row_number() == n(), pval, NA_real_)) %>%
  ungroup() %>%
  select(mod, strata, v25, v50, v75, b, se, blci, buci, P_inter)

header_rows <- dat_strata %>%
  group_by(mod) %>%
  summarise(section = as.character(first(mod)), .groups = "drop") %>%
  mutate(
    subgroup = NA,
    v25 = NA, v50 = NA, v75 = NA,
    b = NA, se = NA, blci = NA, buci = NA, P_inter = NA
  )

subgroup_rows <- dat_strata %>%
  rename(section = mod, subgroup = strata)

df_vitd <- bind_rows(lapply(unique(subgroup_rows$section), function(sec) {
  header <- header_rows %>% filter(section == sec)
  subs   <- subgroup_rows %>% filter(section == sec)
  bind_rows(header, subs)
}))

df_vitd <- df_vitd %>%
  mutate(
    label1 = if_else(
      is.na(subgroup),
      section,
      case_when(
        section == "Age" & !is.na(v50) ~
          sprintf("%s (%.0f, %.0f to %.0f)", subgroup, v50, v25, v75),
        !is.na(v50) ~
          sprintf("%s (%.2f, %.2f to %.2f)", subgroup, v50, v25, v75),
        TRUE ~ subgroup
      )
    ),
    label2 = if_else(is.na(b), "", sprintf("%.2f (%.2f, %.2f)", b, blci, buci)),
    label3 = if_else(is.na(P_inter), "", sprintf("%.2e", P_inter))
  )

df_vitd$label3[df_vitd$label3 == ""] <- NA_character_

labeltext_vitd <- list(
  as.list(c("Effect modifier", df_vitd$label1)),
  as.list(c("Beta (95% CI)", df_vitd$label2)),
  c(list(expression(bold(P)[bold("interaction")])), as.list(df_vitd$label3))
)

is_summary_vitd <- c(TRUE, is.na(df_vitd$subgroup))
mean_vec_vitd   <- c(NA, df_vitd$b)
lower_vec_vitd  <- c(NA, df_vitd$blci)
upper_vec_vitd  <- c(NA, df_vitd$buci)

xticks_vitd <- seq(2.5, 6, by = 0.5)

vitd_plot <- forestplot(
  labeltext  = labeltext_vitd,
  mean       = mean_vec_vitd,
  lower      = lower_vec_vitd,
  upper      = upper_vec_vitd,
  is.summary = is_summary_vitd,
  zero       = NA,
  xlog       = FALSE,
  clip       = c(2.5, 6),
  xlim       = c(2.5, 6),
  xticks     = xticks_vitd,
  boxsize    = 0.2,
  graph.pos  = 2,
  align      = c("l", "l", "l"),
  colgap     = unit(0.7, "cm"),
  col        = fpColors(box = "black", line = "black", summary = "black"),
  txt_gp     = fpTxtGp(
    label = gpar(cex = 0.7),
    ticks = gpar(cex = 0.7),
    xlab  = gpar(cex = 0.7)
  ),
  xlab       = "VitaminD-GS [Replicated] - 25(OH)D\n(nmol/L per 1 SD increase in GS)",
  title      = NA
)

#========================================
# Panel B: Vitamin D (focused score)
#========================================
load("rdata_vitd_interaction_results_focusedScore_20250405_grsSD.RData")

tdi_results <- vitd.tdi.focus$strata_results
tdi_results$mod <- "TDI"

age_results <- vitd.age.focus$strata_results
age_results$mod <- "Age"

sex_results <- vitd.sex.focus$strata_results
sex_results$mod <- "Sex"

fi_results <- vitd.fi.focus$strata_results
fi_results$mod <- "FI"

season_results <- vitd.month.focus$strata_results
season_results$mod <- "Month"

dat_strata <- bind_rows(
  tdi_results,
  age_results,
  sex_results,
  fi_results,
  season_results
) %>%
  mutate(
    blci = b - qnorm(0.975) * se,
    buci = b + qnorm(0.975) * se
  ) %>%
  mutate(
    mod = factor(mod, levels = c("Age", "FI", "Sex", "TDI", "Month"))
  )

dat_pval <- bind_rows(
  data.frame(mod = "TDI",   pval = vitd.tdi.focus$pval_gv),
  data.frame(mod = "Age",   pval = vitd.age.focus$pval_gv),
  data.frame(mod = "Sex",   pval = vitd.sex.focus$pval_gv),
  data.frame(mod = "FI",    pval = vitd.fi.focus$pval_gv),
  data.frame(mod = "Month", pval = vitd.month.focus$pval_gv)
)

dat_strata <- merge(dat_strata, dat_pval, by = "mod", all.x = TRUE) %>%
  mutate(
    season_sort = case_when(
      strata == "Winter" ~ 1,
      strata == "Spring" ~ 2,
      strata == "Summer" ~ 3,
      strata == "Autumn" ~ 4,
      TRUE ~ 99
    )
  ) %>%
  arrange(mod, season_sort, strata) %>%
  group_by(mod) %>%
  mutate(P_inter = if_else(row_number() == n(), pval, NA_real_)) %>%
  ungroup() %>%
  select(mod, strata, v25, v50, v75, b, se, blci, buci, P_inter)

header_rows <- dat_strata %>%
  group_by(mod) %>%
  summarise(section = as.character(first(mod)), .groups = "drop") %>%
  mutate(
    subgroup = NA,
    v25 = NA, v50 = NA, v75 = NA,
    b = NA, se = NA, blci = NA, buci = NA, P_inter = NA
  )

subgroup_rows <- dat_strata %>%
  rename(section = mod, subgroup = strata)

df_vitd <- bind_rows(lapply(unique(subgroup_rows$section), function(sec) {
  header <- header_rows %>% filter(section == sec)
  subs   <- subgroup_rows %>% filter(section == sec)
  bind_rows(header, subs)
}))

df_vitd <- df_vitd %>%
  mutate(
    label1 = if_else(
      is.na(subgroup),
      section,
      case_when(
        section == "Age" & !is.na(v50) ~
          sprintf("%s (%.0f, %.0f to %.0f)", subgroup, v50, v25, v75),
        !is.na(v50) ~
          sprintf("%s (%.2f, %.2f to %.2f)", subgroup, v50, v25, v75),
        TRUE ~ subgroup
      )
    ),
    label2 = if_else(is.na(b), "", sprintf("%.2f (%.2f, %.2f)", b, blci, buci)),
    label3 = if_else(is.na(P_inter), "", sprintf("%.2e", P_inter))
  )

df_vitd$label3[df_vitd$label3 == ""] <- NA_character_

labeltext_vitd <- list(
  as.list(c("Effect modifier", df_vitd$label1)),
  as.list(c("Beta (95% CI)", df_vitd$label2)),
  c(list(expression(bold(P)[bold("interaction")])), as.list(df_vitd$label3))
)

is_summary_vitd <- c(TRUE, is.na(df_vitd$subgroup))
mean_vec_vitd   <- c(NA, df_vitd$b)
lower_vec_vitd  <- c(NA, df_vitd$blci)
upper_vec_vitd  <- c(NA, df_vitd$buci)

xticks_vitd <- seq(2.5, 6, by = 0.5)

vitd_plot.focus <- forestplot(
  labeltext  = labeltext_vitd,
  mean       = mean_vec_vitd,
  lower      = lower_vec_vitd,
  upper      = upper_vec_vitd,
  is.summary = is_summary_vitd,
  zero       = NA,
  xlog       = FALSE,
  clip       = c(2.5, 6),
  xlim       = c(2.5, 6),
  xticks     = xticks_vitd,
  boxsize    = 0.2,
  graph.pos  = 2,
  align      = c("l", "l", "l"),
  colgap     = unit(0.7, "cm"),
  col        = fpColors(box = "black", line = "black", summary = "black"),
  txt_gp     = fpTxtGp(
    label = gpar(cex = 0.7),
    ticks = gpar(cex = 0.7),
    xlab  = gpar(cex = 0.7)
  ),
  xlab       = "VitaminD-GS [Focused] - 25(OH)D\n(nmol/L per 1 SD increase in GS)",
  title      = NA
)

#========================================
# Panel C: BMI
#========================================
load("rdata_bmi_interaction_results_20250407.RData")

tdi_results <- bmi.tdi$strata_results
tdi_results$mod <- "TDI"

met_results <- bmi.met$strata_results
met_results$mod <- "MET"

age_results <- bmi.age$strata_results
age_results$mod <- "Age"

sex_results <- bmi.sex$strata_results
sex_results$mod <- "Sex"

fi_results <- bmi.fi$strata_results
fi_results$mod <- "FI"

dat_strata <- bind_rows(
  tdi_results,
  met_results,
  age_results,
  sex_results,
  fi_results
) %>%
  mutate(
    blci = b - qnorm(0.975) * se,
    buci = b + qnorm(0.975) * se
  ) %>%
  select(mod, strata, v25, v50, v75, b, se, blci, buci) %>%
  mutate(
    mod = factor(mod, levels = c("Age", "FI", "Sex", "TDI", "MET"))
  ) %>%
  arrange(mod, strata)

dat_pval <- bind_rows(
  data.frame(mod = "TDI", pval = bmi.tdi$pval_gv),
  data.frame(mod = "MET", pval = bmi.met$pval_gv),
  data.frame(mod = "Age", pval = bmi.age$pval_gv),
  data.frame(mod = "Sex", pval = bmi.sex$pval_gv),
  data.frame(mod = "FI",  pval = bmi.fi$pval_gv)
)

dat_strata <- merge(dat_strata, dat_pval, by = "mod", all.x = TRUE) %>%
  arrange(mod, strata) %>%
  group_by(mod) %>%
  mutate(P_inter = if_else(row_number() == n(), pval, NA_real_)) %>%
  ungroup()

header_rows <- dat_strata %>%
  group_by(mod) %>%
  summarise(section = first(mod), .groups = "drop") %>%
  mutate(
    subgroup = NA,
    v25 = NA, v50 = NA, v75 = NA,
    b = NA, se = NA, blci = NA, buci = NA, P_inter = NA
  )

subgroup_rows <- dat_strata %>%
  rename(section = mod, subgroup = strata)

df_bmi <- bind_rows(lapply(unique(subgroup_rows$section), function(sec) {
  header <- header_rows %>% filter(section == sec)
  subs   <- subgroup_rows %>% filter(section == sec)
  bind_rows(header, subs)
}))

df_bmi <- df_bmi %>%
  mutate(
    label1 = if_else(
      is.na(subgroup),
      section,
      case_when(
        section == "Age" & !is.na(v50) ~
          sprintf("%s (%.0f, %.0f to %.0f)", subgroup, v50, v25, v75),
        section == "MET" & !is.na(v50) ~
          sprintf("%s (%.0f, %.0f to %.0f)", subgroup, v50, v25, v75),
        !is.na(v50) ~
          sprintf("%s (%.2f, %.2f to %.2f)", subgroup, v50, v25, v75),
        TRUE ~ subgroup
      )
    ),
    label2 = if_else(is.na(b), "", sprintf("%.2f (%.2f, %.2f)", b, blci, buci)),
    label3 = if_else(is.na(P_inter), "", sprintf("%.2e", P_inter))
  )

df_bmi$label3[df_bmi$label3 == ""] <- NA_character_

labeltext_bmi <- list(
  as.list(c("Effect modifier", df_bmi$label1)),
  as.list(c("Beta (95% CI)", df_bmi$label2)),
  c(list(expression(bold(P)[bold("interaction")])), as.list(df_bmi$label3))
)

is_summary_bmi <- c(TRUE, is.na(df_bmi$subgroup))
mean_vec_bmi   <- c(NA, df_bmi$b)
lower_vec_bmi  <- c(NA, df_bmi$blci)
upper_vec_bmi  <- c(NA, df_bmi$buci)

xticks_bmi <- seq(2.5, 5.5, by = 0.5)

bmi_plot <- forestplot(
  labeltext  = labeltext_bmi,
  mean       = mean_vec_bmi,
  lower      = lower_vec_bmi,
  upper      = upper_vec_bmi,
  is.summary = is_summary_bmi,
  zero       = NA,
  xlog       = FALSE,
  clip       = c(2.5, 5.5),
  xlim       = c(2.5, 5.5),
  xticks     = xticks_bmi,
  boxsize    = 0.2,
  graph.pos  = 2,
  align      = c("l", "l", "l"),
  colgap     = unit(0.7, "cm"),
  col        = fpColors(box = "black", line = "black", summary = "black"),
  txt_gp     = fpTxtGp(
    label = gpar(cex = 0.7),
    ticks = gpar(cex = 0.7),
    xlab  = gpar(cex = 0.7)
  ),
  xlab       = "BMI-GS - BMI\n(kg/m² per 1 unit increase in GS)",
  title      = NA
)

#========================================
# Panel D: LDL
#========================================
load("rdata_ldl_interaction_results_20250428.RData")

tdi_results <- ldl.tdi$strata_results
tdi_results$mod <- "TDI"

age_results <- ldl.age$strata_results
age_results$mod <- "Age"

sex_results <- ldl.sex$strata_results
sex_results$mod <- "Sex"

fi_results <- ldl.fi$strata_results
fi_results$mod <- "FI"

dat_strata <- bind_rows(
  tdi_results,
  age_results,
  sex_results,
  fi_results
) %>%
  mutate(
    blci = b - qnorm(0.975) * se,
    buci = b + qnorm(0.975) * se
  ) %>%
  select(mod, strata, v25, v50, v75, b, se, blci, buci)

dat_pval <- bind_rows(
  data.frame(mod = "TDI", pval = ldl.tdi$pval_gv),
  data.frame(mod = "Age", pval = ldl.age$pval_gv),
  data.frame(mod = "Sex", pval = ldl.sex$pval_gv),
  data.frame(mod = "FI",  pval = ldl.fi$pval_gv)
)

dat_strata <- merge(dat_strata, dat_pval, by = "mod", all.x = TRUE) %>%
  group_by(mod) %>%
  mutate(P_inter = if_else(row_number() == n(), pval, NA_real_)) %>%
  ungroup()

header_rows <- dat_strata %>%
  group_by(mod) %>%
  summarise(section = first(mod), .groups = "drop") %>%
  mutate(
    subgroup = NA,
    v25 = NA, v50 = NA, v75 = NA,
    b = NA, se = NA, blci = NA, buci = NA, P_inter = NA
  )

subgroup_rows <- dat_strata %>%
  rename(section = mod, subgroup = strata)

df_ldl <- bind_rows(lapply(unique(subgroup_rows$section), function(sec) {
  header <- header_rows %>% filter(section == sec)
  subs   <- subgroup_rows %>% filter(section == sec)
  bind_rows(header, subs)
}))

df_ldl <- df_ldl %>%
  mutate(
    label1 = if_else(
      is.na(subgroup),
      section,
      case_when(
        section == "Age" & !is.na(v50) ~
          sprintf("%s (%.0f, %.0f to %.0f)", subgroup, v50, v25, v75),
        !is.na(v50) ~
          sprintf("%s (%.2f, %.2f to %.2f)", subgroup, v50, v25, v75),
        TRUE ~ subgroup
      )
    ),
    label2 = if_else(is.na(b), "", sprintf("%.2f (%.2f, %.2f)", b, blci, buci)),
    label3 = if_else(is.na(P_inter), "", sprintf("%.2e", P_inter))
  )

df_ldl$label3[df_ldl$label3 == ""] <- NA_character_

labeltext_ldl <- list(
  as.list(c("Effect modifier", df_ldl$label1)),
  as.list(c("Beta (95% CI)", df_ldl$label2)),
  c(list(expression(bold(P)[bold("interaction")])), as.list(df_ldl$label3))
)

is_summary_ldl <- c(TRUE, is.na(df_ldl$subgroup))
mean_vec_ldl   <- c(NA, df_ldl$b)
lower_vec_ldl  <- c(NA, df_ldl$blci)
upper_vec_ldl  <- c(NA, df_ldl$buci)

xticks_ldl <- seq(0.5, 0.9, by = 0.1)

ldl_plot <- forestplot(
  labeltext  = labeltext_ldl,
  mean       = mean_vec_ldl,
  lower      = lower_vec_ldl,
  upper      = upper_vec_ldl,
  is.summary = is_summary_ldl,
  zero       = NA,
  xlog       = FALSE,
  clip       = c(0.5, 0.9),
  xlim       = c(0.5, 0.9),
  xticks     = xticks_ldl,
  boxsize    = 0.2,
  graph.pos  = 2,
  align      = c("l", "l", "l"),
  colgap     = unit(0.7, "cm"),
  col        = fpColors(box = "black", line = "black", summary = "black"),
  txt_gp     = fpTxtGp(
    label = gpar(cex = 0.7),
    ticks = gpar(cex = 0.7),
    xlab  = gpar(cex = 0.7)
  ),
  xlab       = "LDL-GS - LDL\n(mmol/L per 1 unit increase in GS)",
  title      = NA
)

#========================================
# Combine panels
#========================================
vitd_grob       <- grid::grid.grabExpr({ print(vitd_plot) })
vitd_grob.focus <- grid::grid.grabExpr({ print(vitd_plot.focus) })
bmi_grob        <- grid::grid.grabExpr({ print(bmi_plot) })
ldl_grob        <- grid::grid.grabExpr({ print(ldl_plot) })

vitd_gg       <- as.ggplot(vitd_grob)
vitd_gg.focus <- as.ggplot(vitd_grob.focus)
bmi_gg        <- as.ggplot(bmi_grob)
ldl_gg        <- as.ggplot(ldl_grob)

panel_plot <- ggarrange(
  vitd_gg, vitd_gg.focus, bmi_gg, ldl_gg,
  ncol = 1,
  heights = c(10, 10, 10, 9),
  align = "v",
  labels = c("A", "B", "C", "D"),
  font.label = list(size = 11, face = "bold")
)

panel_plot

#========================================
# Save figure
#========================================
ggsave(
  plot = panel_plot,
  filename = "manuscript/figure5.svg",
  device = "svg",
  scale = 0.9,
  width = 8,
  height = 18,
  units = "in"
)
