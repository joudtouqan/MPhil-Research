#!/usr/bin/env Rscript
# =============================================================================
# Publication-style figures: raw adiposity traits by site, stratified by sex
#   Fig 1  boxplots  : trait (rows) x sex (cols), boxes by study site + ANOVA F
#   Fig 2  histograms: trait distributions overlaid by sex
#
# Raw (pre-INT, pre-adjustment) traits, same source & filters as Table 2 /
# the phenotype prep: analysis_regeneron.csv -> !is.na(genid) -> ages>=18 (67,887).
# Palette: Okabe-Ito subset (colourblind-safe, validated). Box borders + x-axis
# site labels provide secondary encoding for the two low-contrast fills.
# =============================================================================

suppressMessages({
  library(data.table); library(dplyr); library(tidyr)
  library(ggplot2); library(purrr); library(broom)
})

RAW_PATH <- "/home/jt962/rds/rds-post_qc_data-pNR2rM6BWWA/believe/phenotype/data_freeze_3/analysis_regeneron.csv"
OUT_DIR  <- "/home/jt962/rds/hpc-work/believe_adiposity/01_phenotype_prep/publication_figures"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- load + filter (identical to prep) --------------------------------------
cols <- c("genid","study","ages","sex","bmi","waist","hip","whr","fatmass","fatperc")
d <- as.data.frame(fread(RAW_PATH, select = cols, showProgress = FALSE)) |>
  filter(!is.na(genid), ages >= 18)

# readable, ordered site + sex labels
site_lab <- c(BELIEVE = "BELIEVE", BELURBAN = "Urban",
              BELRURAL = "Rural",  BELSLUM  = "Slum")
d <- d |>
  mutate(site = factor(site_lab[study], levels = c("BELIEVE","Urban","Rural","Slum")),
         sex  = factor(sex, levels = c("Male","Female")))

# trait display order + units
trait_levels <- c(bmi="BMI (kg/m²)", waist="Waist (cm)", hip="Hip (cm)",
                  whr="WHR", fatmass="Fat mass (kg)", fatperc="Fat (%)")

long <- d |>
  select(site, sex, all_of(names(trait_levels))) |>
  pivot_longer(all_of(names(trait_levels)), names_to = "trait", values_to = "value") |>
  filter(!is.na(value), !is.na(sex)) |>
  mutate(trait = factor(trait_levels[trait], levels = trait_levels))

# ---- palettes (validated) ----------------------------------------------------
pal_site <- c(BELIEVE="#E69F00", Urban="#009E73", Rural="#0072B2", Slum="#CC79A7")
pal_sex  <- c(Male="#0072B2", Female="#D55E00")

# ---- shared publication theme ------------------------------------------------
theme_pub <- theme_bw(base_size = 12, base_family = "sans") +
  theme(
    panel.grid.minor  = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border      = element_rect(colour = "grey75", linewidth = 0.4),
    strip.background  = element_rect(fill = "grey94", colour = NA),
    strip.text        = element_text(face = "bold", size = 10.5),
    strip.text.y      = element_text(angle = 0),
    axis.text.x       = element_text(angle = 35, hjust = 1, size = 9),
    legend.position   = "top",
    legend.title      = element_text(face = "bold"),
    plot.title        = element_text(face = "bold", size = 13),
    plot.subtitle     = element_text(colour = "grey35", size = 10)
  )

# =============================================================================
# FIGURE 1 : boxplots, trait x sex, boxes by site, per-panel ANOVA F (site effect)
# =============================================================================
fstats <- long |>
  group_by(trait, sex) |>
  summarise(tidy_aov = list(tidy(aov(value ~ site))), .groups = "drop") |>
  mutate(row = map(tidy_aov, ~ filter(.x, term == "site"))) |>
  unnest(row) |>
  mutate(label = sprintf("F=%.0f, p=%.0e", statistic, p.value))

fig1 <- ggplot(long, aes(site, value, fill = site)) +
  geom_boxplot(outlier.size = 0.25, outlier.alpha = 0.25,
               linewidth = 0.35, colour = "grey25", width = 0.65) +
  geom_text(data = fstats, aes(x = -Inf, y = Inf, label = label),
            hjust = -0.06, vjust = 1.5, size = 2.7, colour = "grey30",
            inherit.aes = FALSE) +
  facet_grid(trait ~ sex, scales = "free_y", switch = "y") +
  scale_fill_manual(values = pal_site, name = "Study site") +
  labs(title = "Raw adiposity traits by study site, stratified by sex",
       subtitle = "BELIEVE cohort, post-QC N = 67,887 (age ≥ 18). Box = median & IQR, whiskers 1.5×IQR.\nF, p: one-way ANOVA of site effect within each sex stratum.",
       x = NULL, y = NULL) +
  theme_pub +
  theme(strip.placement = "outside", panel.spacing.y = unit(0.5, "lines"))

ggsave(file.path(OUT_DIR, "fig1_boxplots_traits_by_site_sex.pdf"), fig1,
       width = 8.5, height = 12, device = cairo_pdf)
ggsave(file.path(OUT_DIR, "fig1_boxplots_traits_by_site_sex.png"), fig1,
       width = 8.5, height = 12, dpi = 300)

# =============================================================================
# FIGURE 2 : histograms of each trait distribution, overlaid by sex
# =============================================================================
fig2 <- ggplot(long, aes(value, fill = sex, colour = sex)) +
  geom_histogram(aes(y = after_stat(density)), bins = 45,
                 position = "identity", alpha = 0.45, linewidth = 0.2) +
  facet_wrap(~ trait, scales = "free", ncol = 3) +
  scale_fill_manual(values = pal_sex, name = "Sex") +
  scale_colour_manual(values = pal_sex, name = "Sex") +
  labs(title = "Distribution of raw adiposity traits by sex",
       subtitle = "BELIEVE cohort, post-QC N = 67,887 (18+). Density-normalised histograms; Male vs Female overlaid.",
       x = "Trait value", y = "Density") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 9),
        panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3))

ggsave(file.path(OUT_DIR, "fig2_histograms_traits_by_sex.pdf"), fig2,
       width = 10, height = 6.5, device = cairo_pdf)
ggsave(file.path(OUT_DIR, "fig2_histograms_traits_by_sex.png"), fig2,
       width = 10, height = 6.5, dpi = 300)

# ---- console log -------------------------------------------------------------
cat("N used (18+, genid):", nrow(d), "\n")
cat("Per-trait non-missing rows in long form:\n")
print(long |> count(trait))
cat("\nANOVA (site effect) F / p by trait x sex:\n")
print(fstats |> select(trait, sex, statistic, p.value) |> as.data.frame(), digits = 4)
cat("\nFigures written to:", OUT_DIR, "\n")
cat("  fig1_boxplots_traits_by_site_sex.{pdf,png}\n  fig2_histograms_traits_by_sex.{pdf,png}\n")
