
#===============================================================================
#
# Figure 2: simulation results under three correlation scenarios
# corr(E1, E2) = 0, 0.3, and 0.6
#
#===============================================================================


rm(list = ls())
library(ggplot2)
library(dplyr)

# Set your working directory to the project folder before running this script.
load("rdata_nlmr_gxe_sim_results_modifier_20251009.RData")


# rho = 0
dat1 <- sim_results_rho0 %>%
  mutate(
    strata_label = factor(strata, levels = 0:10, labels = c("Overall", 1:10)),
    method = factor(method, levels = c("lmr", "rank+X")),
    Methods = recode(method,
                     "lmr"    = "LMR (0)",
                     "rank+X" = "Rank (0)"
    )
  )

# rho = 0.3
dat2 <- sim_results_rho03 %>%
  mutate(
    strata_label = factor(strata, levels = 0:10, labels = c("Overall", 1:10)),
    method = factor(method, levels = c("lmr", "rank+X")),
    Methods = recode(method,
                     "lmr"    = "LMR (0.3)",
                     "rank+X" = "Rank (0.3)"
    )
  )

# rho = 0.6
dat3 <- sim_results_rho06 %>%
  mutate(
    strata_label = factor(strata, levels = 0:10, labels = c("Overall", 1:10)),
    method = factor(method, levels = c("lmr", "rank+X")),
    Methods = recode(method,
                     "lmr"    = "LMR (0.6)",
                     "rank+X" = "Rank (0.6)"
    )
  )

# Combine data and set legend order
dat <- bind_rows(dat1, dat2, dat3) %>%
  mutate(
    Methods = factor(
      Methods,
      levels = c(
        "LMR (0)", "LMR (0.3)", "LMR (0.6)",
        "Rank (0)", "Rank (0.3)", "Rank (0.6)"
      )
    )
  )


# Colours
method_pal <- c(
  "LMR (0)"    = "#A6CEE5FF",
  "LMR (0.3)"  = "#59A3CFFF",
  "LMR (0.6)"  = "#0072B5FF",
  "Rank (0)"   = "#F5D5B3FF",
  "Rank (0.3)" = "#ECB173FF",
  "Rank (0.6)" = "#E18727FF"
)

pos <- position_dodge(width = 0.65)

# Plot
plot <- ggplot(
  dat,
  aes(x = strata_label, y = lace, fill = Methods, colour = Methods)
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.3,
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
  theme(
    legend.position = "right"
  ) +
  scale_fill_manual(values = method_pal, drop = FALSE) +
  scale_color_manual(values = method_pal, drop = FALSE)

plot


ggsave(
  plot = plot,
  filename = "manuscript/figure2.svg",
  scale = 0.9,
  width = 15,
  height = 6,
  units = "in"
)
