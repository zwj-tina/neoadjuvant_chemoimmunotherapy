# R
library(ggplot2)
library(ggpubr)
library(dplyr)
library(rstatix)
library(tidyr)
library(lmerTest)
library(emmeans)
# Kenward-Roger
has_pbkrtest <- requireNamespace("pbkrtest", quietly = TRUE)

# 1. (Raw Data)
df_raw <- data.frame(
  Organoid = c(
    rep("Sample 1", 8),
    rep("Sample 2", 16),
    rep("Sample 3", 12)
  ),
  Treatment = c(
    rep(c("NC", "Pembro", "WX261-1", "Pembro+WX261-1"), each = 2),
    rep(c("NC", "Pembro", "WX261-1", "Pembro+WX261-1"), each = 4),
    rep(c("NC", "Pembro", "WX261-1", "Pembro+WX261-1"), each = 3)
  ),
  RawValue = c(
    5.076862, 4.532720,
    4.310866, 3.951289,
    2.476043, 2.841675,
    3.413037, 3.638607,
    
    0.909810, 1.011095, 1.139970, 1.049008,
    0.751899, 0.831850, 1.033690, 0.905589,
    0.690844, 0.747675, 0.741792, 0.628485,
    1.111004, 0.934249, 0.159663, 0.673583,
    
    2.195545, 2.778424, 3.818130,
    2.552180, 2.926712, 2.849991,
    2.154485, 2.720134, 3.315702,
    2.057212, 2.042949, 1.919889
  )
)

# 2. NC
# —— (technical replicate level)
df_tech <- df_raw %>%
  group_by(Organoid) %>%
  mutate(Viability = RawValue / mean(RawValue[Treatment == "NC"])) %>%
  ungroup()

df_tech$Treatment <- factor(
  df_tech$Treatment,
  levels = c("NC", "Pembro", "WX261-1", "Pembro+WX261-1")
)

# log2log-normal
df_tech <- df_tech %>%
  mutate(log2Viability = log2(Viability))

# 3.
# —— n = 3
# log2
# Viability
df_bio <- df_tech %>%
  group_by(Organoid, Treatment) %>%
  summarise(
    Viability     = mean(Viability),
    log2Viability = mean(log2Viability), # log =
    n_tech = n(),
    .groups = "drop"
  )

print(df_bio) # Organoid x Treatment 3 x 4 = 12

# 4. +
# rstatix::t_test(paired=TRUE)""
df_bio <- df_bio %>%
  arrange(Treatment, Organoid)

check_order <- df_bio %>%
  group_by(Treatment) %>%
  summarise(order_key = paste(Organoid, collapse = ","), .groups = "drop")

if (n_distinct(check_order$order_key) != 1) {
  stop("Fatal error: the Organoid order differs across Treatment groups; the paired t-test would be misaligned. Please check the data.")
} else {
  cat("Pairing order check passed: the Organoid order is identical across Treatment groups; paired tests can be run safely.\n")
}

# 4
df_wide_check <- df_bio %>%
  select(Organoid, Treatment, Viability) %>%
  pivot_wider(names_from = Treatment, values_from = Viability)
print(df_wide_check)

# 5. n = 3 t
# 5.1 tlog2log2 fold change
# 5BH
stat_test_paired <- df_bio %>%
  t_test(
    log2Viability ~ Treatment,
    paired = TRUE,
    comparisons = list(
      c("NC", "Pembro+WX261-1")
    )
  ) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance(p.col = "p.adj") %>% # p
  add_xy_position(x = "Treatment")

print(stat_test_paired)

# filter
# 5""
stat_test_plot <- stat_test_paired

# 5.2 /Organoid
# log2(Viability)t
model <- lmer(log2Viability ~ Treatment + (1 | Organoid), data = df_tech)
cat("\n=== Mixed-effects model ANOVA (Satterthwaite degrees of freedom) ===\n")
print(anova(model))

# n=3Satterthwaite
# Kenward-Rogerpbkrtest
if (has_pbkrtest) {
  cat("\n=== Mixed-effects model ANOVA (Kenward-Roger degrees of freedom, sensitivity analysis) ===\n")
  print(anova(model, ddf = "Kenward-Roger"))
} else {
  cat("\n[Note] pbkrtest is not installed; skipping the Kenward-Roger sensitivity analysis.\n")
}
cat("\n[Note] With only a few random-effect levels (Organoid = 3), the variance component estimates are unstable; ")
cat("the mixed-model results serve only as a supplementary validation of the paired t-test, not as the basis for the main conclusion.\n")

cat("\n=== Mixed-effects model post-hoc pairwise comparisons (Tukey adjustment, log2 scale) ===\n")
emm <- emmeans(model, pairwise ~ Treatment, adjust = "tukey")
print(emm)

# 6.
# (A) n=3 ± SD
# (B) +
# log2""

treat_colors <- c(
  "NC" = "#D1D5DB",
  "Pembro" = "#60A5FA",
  "WX261-1" = "#34D399",
  "Pembro+WX261-1" = "#F87171"
)

df_bio$Treatment <- factor(
  df_bio$Treatment,
  levels = c("NC", "Pembro", "WX261-1", "Pembro+WX261-1")
)

p <- ggbarplot(
  df_bio,
  x = "Treatment",
  y = "Viability",
  fill = "Treatment",
  color = "black",
  add = "mean_sd",
  palette = treat_colors,
  xlab = "Treatment",
  ylab = "Relative Viability",
  title = "Drug Sensitivity across Patient-Derived Organoids (n = 3)"
) +
  geom_point(
    data = df_bio,
    aes(x = Treatment, y = Viability, group = Organoid),
    position = position_dodge(width = 0),
    size = 2.5,
    shape = 21,
    fill = "white",
    color = "black",
    alpha = 0.9
  ) +
  geom_line(
    data = df_bio,
    aes(x = Treatment, y = Viability, group = Organoid),
    color = "grey40",
    linewidth = 0.5,
    linetype = "dashed"
  ) +
  stat_pvalue_manual(
    stat_test_plot,
    label = "p.adj = {p.adj}", # p
    size = 3.2,
    tip.length = 0.01,
    step.increase = 0.1
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 12, face = "bold", color = "black"),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.title.x = element_text(size = 13, face = "bold"),
    axis.title.y = element_text(size = 13, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 11)
  )

print(p)

ggsave(
  "MSD/figure/Organiod_combined_biological_n3.pdf",
  plot = p, width = 5.5, height = 5.5
)

