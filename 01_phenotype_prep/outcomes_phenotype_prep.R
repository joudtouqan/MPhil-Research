#!/usr/bin/env Rscript
# ============================================================================
# In-sample OUTCOME phenotype prep for BELIEVE genetic-correlation matrix
#   HbA1c  -> hba1c_INT  (continuous: covariate-adjust + INT, mirrors adiposity)
#   T2D    -> t2d         (binary 0/1, from analysis var `diab2`)
#   CHD    -> chd         (binary 0/1, from analysis var `chd`)
# Output: believe_outcomes_saige.txt with covariates for SAIGE binary Step 1.
# ============================================================================
suppressMessages({ library(dplyr) })

RAW_PATH <- "/home/jt962/rds/rds-post_qc_data-pNR2rM6BWWA/believe/phenotype/data_freeze_3/analysis_regeneron.csv"
PCA_DIR  <- "/home/jt962/rds/rds-post_qc_data-pNR2rM6BWWA/believe/genotype/imputed/pca_flashpca_Freeze_Two"
OUTDIR   <- "/home/jt962/rds/hpc-work/believe_adiposity/01_phenotype_prep"
OUTFILE  <- file.path(OUTDIR, "believe_outcomes_saige.txt")

# ---- 1. Read raw phenotypes -------------------------------------------------
raw <- read.csv(RAW_PATH, stringsAsFactors = FALSE)
d <- raw %>%
  select(genid, idno, study, ages, sex, hba1c, diab2, chd) %>%
  filter(!is.na(genid))                       # must be genotyped
cat("Rows with genid:", nrow(d), "\n")

# ---- 2. Covariates ----------------------------------------------------------
d <- d %>% mutate(
  sex_numeric = case_when(sex == "Male" ~ 0, sex == "Female" ~ 1, TRUE ~ NA_real_),
  age_sq      = ages^2
) %>% filter(ages >= 18)                       # same adult filter as adiposity prep
cat("Rows after age>=18:", nrow(d), "\n")

# ---- 3. Binary outcomes: coerce to clean 0/1 --------------------------------
d$t2d <- suppressWarnings(as.integer(d$diab2))
d$chd <- suppressWarnings(as.integer(d$chd))
d$t2d[!d$t2d %in% c(0L, 1L)] <- NA
d$chd[!d$chd %in% c(0L, 1L)] <- NA
cat("T2D  cases/controls:", sum(d$t2d == 1, na.rm=TRUE), "/", sum(d$t2d == 0, na.rm=TRUE), "\n")
cat("CHD  cases/controls:", sum(d$chd == 1, na.rm=TRUE), "/", sum(d$chd == 0, na.rm=TRUE), "\n")

# ---- 4. HbA1c: numeric, plausibility filter --------------------------------
d$hba1c <- suppressWarnings(as.numeric(d$hba1c))
d$hba1c[d$hba1c < 3 | d$hba1c > 20] <- NA       # clinical plausibility (% units)
cat("HbA1c non-missing:", sum(!is.na(d$hba1c)), "\n")

# ---- 5. Merge PCs (unrelated + projected related) ---------------------------
pcs  <- read.table(file.path(PCA_DIR, "pcs.txt"),        header = TRUE) %>%
          rename(genid = IID) %>% select(genid, PC1:PC20)
proj <- read.table(file.path(PCA_DIR, "projections.txt"), header = TRUE) %>%
          rename(genid = IID) %>% select(genid, PC1:PC20)
pca_all <- bind_rows(pcs, proj) %>% distinct(genid, .keep_all = TRUE)
cat("Samples with PCs:", nrow(pca_all), "\n")

n_before <- nrow(d)
d <- inner_join(d, pca_all, by = "genid")
cat("N before PC merge:", n_before, " after:", nrow(d), "\n")

# ---- 6. HbA1c: covariate-adjust (age, age^2, sex, study, PC1-20) then INT ---
int_transform <- function(x) {
  r <- rank(x, na.last = "keep", ties.method = "average")
  qnorm((r - 0.5) / sum(!is.na(x)))
}
pc_cols <- paste0("PC", 1:20)
f_hba1c <- as.formula(paste("hba1c ~ ages + age_sq + sex_numeric + study +",
                            paste(pc_cols, collapse = " + ")))
d$hba1c_adj <- residuals(lm(f_hba1c, data = d, na.action = na.exclude))
d$hba1c_INT <- int_transform(d$hba1c_adj)
cat("HbA1c INT non-missing:", sum(!is.na(d$hba1c_INT)), "\n")

# ---- 7. Assemble SAIGE phenotype file --------------------------------------
out <- d %>%
  mutate(
    FID = paste0("CAMBRIDGE-BELIEVE_", genid, "_", genid),
    IID = paste0("CAMBRIDGE-BELIEVE_", genid, "_", genid)
  ) %>%
  select(FID, IID, hba1c_INT, t2d, chd,
         ages, age_sq, sex_numeric, study,
         PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10,
         PC11, PC12, PC13, PC14, PC15, PC16, PC17, PC18, PC19, PC20)

write.table(out, OUTFILE, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
cat("Saved:", OUTFILE, "  (", nrow(out), "rows x", ncol(out), "cols )\n")

# ---- 8. Report N with each outcome AND genotyped covariates -----------------
cat("\n=== analysis-ready N (non-missing outcome + full covariates) ===\n")
cc <- complete.cases(out[, c("ages","age_sq","sex_numeric","study", pc_cols)])
cat("HbA1c:", sum(!is.na(out$hba1c_INT) & cc), "\n")
cat("T2D  : cases", sum(out$t2d == 1 & cc, na.rm=TRUE),
    " controls", sum(out$t2d == 0 & cc, na.rm=TRUE), "\n")
cat("CHD  : cases", sum(out$chd == 1 & cc, na.rm=TRUE),
    " controls", sum(out$chd == 0 & cc, na.rm=TRUE), "\n")
