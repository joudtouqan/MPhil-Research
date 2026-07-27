#!/usr/bin/env Rscript
# ============================================================================
# EUR comparator steps 6-7: power analysis 
#
# Power: EUR Nikpay ncase=60801 / ncontrol=123504; R2 from EUR diagnostics.
# Comparison: for traits in both ancestries, Cochran's Q between the two IVW
# estimates (df=1) + SAS:EUR ratio. Trait map:
#   bmi(EUR)=BMI(SAS); fatperc(EUR)=Fat percentage(SAS);
#   whr(EUR)~WHRadjBMI(SAS)  [CAVEAT: WHR != WHRadjBMI; closest available match]
# ============================================================================
.libPaths("/rds/user/jt962/hpc-work/software/R/4.6.0_libs")
suppressMessages({library(data.table)})

EUR <- "/home/jt962/rds/hpc-work/believe_adiposity/09_MR/eur_comparator"
SAS <- "/home/jt962/rds/hpc-work/believe_adiposity/09_MR/Two_sample_MR"
RES <- file.path(EUR, "results"); POW <- file.path(EUR, "power")
dir.create(POW, showWarnings = FALSE)

## ---- power (Brion 2013, binary outcome), EUR Nikpay ----
NCASE <- 60801; NCTRL <- 123504; N <- NCASE + NCTRL; K <- NCASE/N; alpha <- 0.05
thres <- qchisq(1 - alpha, 1)
powfun  <- function(OR, R2){ b01 <- K*(OR/(1+K*(OR-1)) - 1); ncp <- N*R2*b01^2/(K*(1-K) - b01^2); pchisq(thres, 1, ncp = ncp, lower.tail = FALSE) }
minOR80 <- function(R2) uniroot(function(or) powfun(or, R2) - 0.80, c(1.0001, 5))$root

dg  <- fread(file.path(RES, "diagnostics_all_traits.tsv"))
res <- fread(file.path(RES, "mr_all_traits.tsv"))[method == "IVW (random)"]
pw  <- merge(dg[, .(trait, label, R2 = R2_exposure, nSNP = n_retained, N_exposure, F_overall)],
             res[, .(trait, obsOR = OR)], by = "trait")[order(-R2)]
pw[, `:=`(power_obsOR = mapply(powfun, obsOR, R2),
          power_OR1.2  = mapply(powfun, 1.2, R2),
          power_OR1.5  = mapply(powfun, 1.5, R2),
          minOR_80pct  = sapply(R2, minOR80))]
fwrite(pw, file.path(POW, "table_power_eur.tsv"), sep = "\t")
cat("================ EUR power (Nikpay 60,801/123,504) ================\n"); print(pw, digits = 4)

## ---- trans-ancestry comparison ----
eur <- fread(file.path(RES, "mr_all_traits.tsv"))[method == "IVW (random)"]
sas <- fread(file.path(SAS, "results/mr_all_traits.tsv"))[method == "IVW (random)"]
# SAS R2/power for the asymmetry caveat
saspow <- fread(file.path(SAS, "results/diagnostics/table_power.tsv"))

map <- data.table(
  trait     = c("bmi", "fatperc", "whr"),
  eur_label = c("BMI (Yengo)", "Body fat % (Pan-UKB EUR)", "WHR (Pan-UKB EUR)"),
  sas_label = c("BMI", "Fat percentage", "WHRadjBMI"),
  note      = c("direct match", "direct match", "WHR(EUR) vs WHRadjBMI(SAS) - not identical"))

cmp <- rbindlist(lapply(seq_len(nrow(map)), function(i){
  e <- eur[label == map$eur_label[i]]; s <- sas[label == map$sas_label[i]]
  sp <- saspow[label == map$sas_label[i]]
  Q  <- (e$b - s$b)^2 / (e$se^2 + s$se^2); Qp <- pchisq(Q, df = 1, lower.tail = FALSE)
  data.table(trait = map$trait[i], note = map$note[i],
    EUR_OR = e$OR, EUR_lci = e$OR_lci, EUR_uci = e$OR_uci, EUR_p = e$pval,
    SAS_OR = s$OR, SAS_lci = s$OR_lci, SAS_uci = s$OR_uci, SAS_p = s$pval,
    SAS_power_obsOR = sp$power_obsOR, SAS_minOR80 = sp$minOR_80pct,
    ratio_SAS_EUR_logOR = s$b / e$b,
    Q_between = Q, Q_between_p = Qp,
    CIs_overlap = !(e$OR_lci > s$OR_uci | s$OR_lci > e$OR_uci))
}))
fwrite(cmp, file.path(POW, "transancestry_comparison.tsv"), sep = "\t")
cat("\n================ Trans-ancestry SAS vs EUR (IVW) ================\n")
print(cmp[, .(trait, EUR_OR, EUR_p, SAS_OR, SAS_p, SAS_power_obsOR, Q_between, Q_between_p, CIs_overlap, note)], digits = 3)
cat("\nReadout: Q_between_p<0.05 => estimates differ beyond chance (candidate ancestry diff).\n")
cat("SAS_power_obsOR shows the SAS arm's power at its own observed OR -- low power means a\n")
cat("null SAS estimate is uninformative, so EUR-sig vs SAS-null is NOT alone an ancestry difference.\n")
cat("DONE\n")
