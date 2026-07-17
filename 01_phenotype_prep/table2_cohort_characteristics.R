#!/usr/bin/env Rscript
# =============================================================================
# Table 2 (cohort characteristics) regeneration
#   (a) BELIEVE + BELURBAN merged into a single "Urban" setting
#   (b) sex-stratified companion + setting x sex cross-tab (appendix candidate)
#
# Uses RAW (pre-INT, pre-adjustment) trait values from the analysis-ready
# phenotype source. Descriptive only: no residualised / INT values here.
#
# Source is the same file + filters the SAIGE Step 1 phenotype prep consumes:
#   read analysis_regeneron.csv -> filter !is.na(genid) -> filter ages >= 18
#   => 67,887 post-QC IIDs (52,009 / 6,685 / 4,041 / 5,152 by subgroup).
# NB: believe_phenotypes_saige.txt has only 67,881 (6 dropped for missing PCs)
#     and carries INT traits, so it is NOT used for descriptives.
# =============================================================================

suppressMessages({
  library(data.table)
  library(dplyr)
})

RAW_PATH <- "/home/jt962/rds/rds-post_qc_data-pNR2rM6BWWA/believe/phenotype/data_freeze_3/analysis_regeneron.csv"
OUT_DIR  <- "/home/jt962/rds/hpc-work/believe_adiposity/01_phenotype_prep/table2_cohort_characteristics"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- traits & display config -------------------------------------------------
# name in data, label, decimal places
TRAITS <- list(
  list(v = "bmi",     lab = "BMI (kg/m2)",   dp = 1),
  list(v = "waist",   lab = "Waist (cm)",    dp = 1),
  list(v = "hip",     lab = "Hip (cm)",      dp = 1),
  list(v = "whr",     lab = "WHR",           dp = 2),   # 2 dp: keeps Overall inside site range
  list(v = "fatmass", lab = "Fat mass (kg)", dp = 1),
  list(v = "fatperc", lab = "Fat % (%)",     dp = 1)
)
AGE_DP <- 1

# ============================================================================
# Step 0 : load + filter, print gate
# ============================================================================
cols <- c("genid","study","idno","ages","sex","bmi","waist","hip","whr","fatmass","fatperc")
raw  <- as.data.frame(fread(RAW_PATH, select = cols, showProgress = FALSE))

d <- raw %>% filter(!is.na(genid)) %>% filter(ages >= 18)

cat("================ STEP 0: SOURCE & GATE ================\n")
cat("Resolved path :", RAW_PATH, "\n")
cat("nrow (post-QC):", nrow(d), "\n")
cat("table(subgroup):\n"); print(table(d$study, useNA = "ifany"))

stopifnot(nrow(d) == 67887)
sg <- table(d$study)
stopifnot(sg["BELIEVE"] == 52009, sg["BELRURAL"] == 6685,
          sg["BELSLUM"] == 4041,  sg["BELURBAN"] == 5152)
cat("GATE 0 PASSED: nrow == 67,887 and subgroup counts match.\n\n")

# ============================================================================
# Step 1 : recode setting, gate
# ============================================================================
d <- d %>% mutate(
  rsetting = dplyr::case_when(
    study %in% c("BELIEVE", "BELURBAN") ~ "Urban",
    study == "BELRURAL" ~ "Rural",
    study == "BELSLUM"  ~ "Slum"
  ),
  rsetting = factor(rsetting, levels = c("Urban", "Rural", "Slum")),
  sex      = factor(sex, levels = c("Male", "Female"))
)

urban_n  <- sum(d$rsetting == "Urban")
urban_f  <- sum(d$rsetting == "Urban" & d$sex == "Female")
cat("================ STEP 1: SETTING RECODE ================\n")
cat("table(rsetting):\n"); print(table(d$rsetting, useNA = "ifany"))
cat("Urban N =", urban_n, " | Urban female =", urban_f, "\n")
stopifnot(urban_n == 57161, urban_f == 33258)
cat("GATE 1 PASSED: Urban N == 57,161 and Urban female == 33,258.\n\n")

# ============================================================================
# helpers
# ============================================================================
nfmt <- function(x) formatC(x, format = "d", big.mark = ",")
pfmt <- function(p) if (is.na(p)) "-" else if (p < 2.2e-16) "<2.2e-16" else sprintf("%.3g", p)

# mean (SD) with per-cell non-missing n, complete-case for that trait
cell_cont <- function(x, dp) {
  ok <- !is.na(x)
  n  <- sum(ok)
  if (n == 0) return(list(str = "-", n = 0L, mean = NA_real_))
  m  <- mean(x[ok]); s <- sd(x[ok])
  list(str = sprintf("%.*f (%.*f) [n=%s]", dp, m, dp, s, nfmt(n)),
       n = as.integer(n), mean = m)
}

# n (%) for a logical vector over a denominator of non-missing sex
cell_pct <- function(is_target, valid) {
  n <- sum(is_target & valid); den <- sum(valid)
  list(str = sprintf("%s (%.1f%%)", nfmt(n), 100 * n / den), n = as.integer(n))
}

# ============================================================================
# Step 2 : Table 2a -- by setting (Overall, Urban, Rural, Slum)
# ============================================================================
groups_2a <- list(
  Overall = d,
  Urban   = d %>% filter(rsetting == "Urban"),
  Rural   = d %>% filter(rsetting == "Rural"),
  Slum    = d %>% filter(rsetting == "Slum")
)
gnames_2a <- names(groups_2a)

t2a <- list()
# N row
t2a[["N (participants)"]] <- c(sapply(groups_2a, nrow) |> vapply(nfmt, ""), P = "")

# Age
age_cells <- lapply(groups_2a, function(g) cell_cont(g$ages, AGE_DP))
p_age <- kruskal.test(ages ~ rsetting, data = d)$p.value
t2a[["Age (years)"]] <- c(sapply(age_cells, `[[`, "str"), P = pfmt(p_age))

# Female n (%)
fem_cells <- lapply(groups_2a, function(g) cell_pct(g$sex == "Female", !is.na(g$sex)))
p_sex <- chisq.test(table(d$rsetting, d$sex))$p.value
t2a[["Female, n (%)"]] <- c(sapply(fem_cells, `[[`, "str"), P = pfmt(p_sex))

# Continuous traits
overall_in_range <- list()
for (tr in TRAITS) {
  cells <- lapply(groups_2a, function(g) cell_cont(g[[tr$v]], tr$dp))
  # Kruskal-Wallis across the 3 settings (complete-case)
  sub <- d[!is.na(d[[tr$v]]), c(tr$v, "rsetting")]
  p   <- kruskal.test(sub[[tr$v]], sub$rsetting)$p.value
  t2a[[tr$lab]] <- c(sapply(cells, `[[`, "str"), P = pfmt(p))
  # range check (unrounded)
  site_means <- c(cells$Urban$mean, cells$Rural$mean, cells$Slum$mean)
  ov <- cells$Overall$mean
  overall_in_range[[tr$lab]] <- (ov >= min(site_means) & ov <= max(site_means))
}

t2a_df <- data.frame(Characteristic = names(t2a),
                     do.call(rbind, t2a), row.names = NULL, check.names = FALSE)

# ============================================================================
# Step 3 : Table 2b -- by sex (Overall, Male, Female)
# ============================================================================
groups_2b <- list(
  Overall = d,
  Male    = d %>% filter(sex == "Male"),
  Female  = d %>% filter(sex == "Female")
)
t2b <- list()
t2b[["N (participants)"]] <- c(sapply(groups_2b, nrow) |> vapply(nfmt, ""), P = "")

age_b <- lapply(groups_2b, function(g) cell_cont(g$ages, AGE_DP))
p_age_sex <- wilcox.test(ages ~ sex, data = d)$p.value
t2b[["Age (years)"]] <- c(sapply(age_b, `[[`, "str"), P = pfmt(p_age_sex))

for (tr in TRAITS) {
  cells <- lapply(groups_2b, function(g) cell_cont(g[[tr$v]], tr$dp))
  sub <- d[!is.na(d[[tr$v]]), c(tr$v, "sex")]
  p   <- wilcox.test(sub[[tr$v]] ~ sub$sex)$p.value
  t2b[[tr$lab]] <- c(sapply(cells, `[[`, "str"), P = pfmt(p))
}
t2b_df <- data.frame(Characteristic = names(t2b),
                     do.call(rbind, t2b), row.names = NULL, check.names = FALSE)

# ============================================================================
# Step 3b : Table 2c -- setting x sex cross-tab (6 columns, appendix candidate)
# ============================================================================
combos <- list(
  "Urban Male"   = d %>% filter(rsetting == "Urban", sex == "Male"),
  "Urban Female" = d %>% filter(rsetting == "Urban", sex == "Female"),
  "Rural Male"   = d %>% filter(rsetting == "Rural", sex == "Male"),
  "Rural Female" = d %>% filter(rsetting == "Rural", sex == "Female"),
  "Slum Male"    = d %>% filter(rsetting == "Slum",  sex == "Male"),
  "Slum Female"  = d %>% filter(rsetting == "Slum",  sex == "Female")
)
t2c <- list()
t2c[["N (participants)"]] <- sapply(combos, nrow) |> vapply(nfmt, "")
t2c[["Age (years)"]] <- sapply(combos, function(g) cell_cont(g$ages, AGE_DP)$str)
for (tr in TRAITS)
  t2c[[tr$lab]] <- sapply(combos, function(g) cell_cont(g[[tr$v]], tr$dp)$str)
t2c_df <- data.frame(Characteristic = names(t2c),
                     do.call(rbind, t2c), row.names = NULL, check.names = FALSE)

# ============================================================================
# Step 4 : outputs (TSV + markdown)
# ============================================================================
write.table(t2a_df, file.path(OUT_DIR, "table2a_by_setting.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(t2b_df, file.path(OUT_DIR, "table2b_by_sex.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(t2c_df, file.path(OUT_DIR, "table2c_setting_by_sex.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

md_table <- function(df) {
  hdr <- paste0("| ", paste(colnames(df), collapse = " | "), " |")
  sep <- paste0("| ", paste(rep("---", ncol(df)), collapse = " | "), " |")
  bod <- apply(df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |"))
  paste(c(hdr, sep, bod), collapse = "\n")
}
md <- paste0(
  "### Table 2a. Cohort characteristics by setting (Urban = BELIEVE + BELURBAN)\n\n",
  "Continuous: mean (SD) [non-missing n]. P: Kruskal-Wallis (continuous) / chi-square (sex) across Urban/Rural/Slum.\n\n",
  md_table(t2a_df),
  "\n\n### Table 2b. Cohort characteristics by sex\n\n",
  "Continuous: mean (SD) [non-missing n]. P: Wilcoxon rank-sum (Male vs Female).\n\n",
  md_table(t2b_df),
  "\n\n### Table 2c (appendix candidate). Setting x sex cross-tabulation\n\n",
  "Continuous: mean (SD) [non-missing n].\n\n",
  md_table(t2c_df), "\n"
)
writeLines(md, file.path(OUT_DIR, "table2_markdown.md"))

# ============================================================================
# console log: every N and p-value
# ============================================================================
cat("================ STEP 2/3: COMPUTED N & P ================\n")
cat("-- Table 2a group Ns --\n")
for (g in gnames_2a) cat(sprintf("  %-8s N = %s\n", g, nfmt(nrow(groups_2a[[g]]))))
cat("-- Table 2a female Ns --\n")
for (g in gnames_2a) cat(sprintf("  %-8s female = %s\n", g, nfmt(fem_cells[[g]]$n)))
cat("-- Table 2a per-trait non-missing N by setting --\n")
for (tr in c(list(list(v="ages",lab="Age",dp=AGE_DP)), TRAITS)) {
  ns <- sapply(groups_2a, function(g) sum(!is.na(g[[tr$v]])))
  cat(sprintf("  %-14s %s\n", tr$lab,
              paste(sprintf("%s=%s", gnames_2a, sapply(ns, nfmt)), collapse = "  ")))
}
cat("-- Table 2a p-values (KW / chi-sq) --\n")
cat(sprintf("  Age  : %.4g (Kruskal-Wallis)\n", p_age))
cat(sprintf("  Sex  : %.4g (chi-square)\n", p_sex))
for (tr in TRAITS) {
  sub <- d[!is.na(d[[tr$v]]), c(tr$v, "rsetting")]
  cat(sprintf("  %-8s: %.4g (Kruskal-Wallis)\n", tr$lab,
              kruskal.test(sub[[tr$v]], sub$rsetting)$p.value))
}
cat("-- Table 2b group Ns & Wilcoxon p --\n")
for (g in names(groups_2b)) cat(sprintf("  %-8s N = %s\n", g, nfmt(nrow(groups_2b[[g]]))))
cat(sprintf("  Age  : %.4g\n", p_age_sex))
for (tr in TRAITS) {
  sub <- d[!is.na(d[[tr$v]]), c(tr$v, "sex")]
  cat(sprintf("  %-8s: %.4g\n", tr$lab, wilcox.test(sub[[tr$v]] ~ sub$sex)$p.value))
}

# ============================================================================
# Verification gates
# ============================================================================
cat("\n================ VERIFICATION GATES ================\n")
setting_ns <- sapply(groups_2a[c("Urban","Rural","Slum")], nrow)
cat(sprintf("Column Ns sum (Urban+Rural+Slum) = %s  [expect 67,887]\n", nfmt(sum(setting_ns))))
stopifnot(sum(setting_ns) == 67887)

fem_sum <- sum(sapply(fem_cells[c("Urban","Rural","Slum")], `[[`, "n"))
cat(sprintf("Female counts sum (Urban+Rural+Slum) = %s  [expect 40,096]\n", nfmt(fem_sum)))
stopifnot(fem_sum == 40096)

cat("Overall mean within [min,max] of the 3 setting means:\n")
for (nm in names(overall_in_range))
  cat(sprintf("  %-14s %s\n", nm, ifelse(overall_in_range[[nm]], "OK", "*** OUT OF RANGE ***")))
stopifnot(all(unlist(overall_in_range)))

cat("\nALL GATES PASSED. Outputs written to:\n  ", OUT_DIR, "\n", sep = "")
cat("Files: table2a_by_setting.tsv, table2b_by_sex.tsv, table2c_setting_by_sex.tsv, table2_markdown.md\n")
