#!/usr/bin/env Rscript
# ============================================================================
# Steps 4-6: harmonise BELIEVE (exposure) vs Genes & Health CHD (outcome),
# run two-sample MR, diagnostics, and plots, for all six adiposity traits.
#
# Harmonisation is done explicitly
#   exposure effect allele = Allele2, other = Allele1, eaf = AF_Allele2  
#   outcome  effect allele = ALLELE1, other = ALLELE0, eaf = A1FREQ      
# Both datasets are TOPMed/GRCh38 -> same strand; join on chr:pos, match alleles,
# flip outcome beta + eaf when the effect alleles are swapped, drop ambiguous
# palindromes and allele/frequency-inconsistent SNPs.
#
# MR estimators from the MendelianRandomization package (v0.10.0):
#   IVW (random), MR-Egger (+intercept), weighted median, weighted mode (mr_mbe).
# Diagnostics: per-SNP F, Cochran's Q (IVW & Egger), Egger intercept, leave-one-out.
# Outcome BETA is log-odds -> OR = exp(b).
# ============================================================================
.libPaths("/rds/user/jt962/hpc-work/software/R/4.6.0_libs")
suppressMessages({library(data.table); library(MendelianRandomization); library(ggplot2); library(patchwork)})

MR   <- "/home/jt962/rds/hpc-work/believe_adiposity/09_MR"
INST <- file.path(MR, "instruments")
HARM <- file.path(MR, "harmonised")
RES  <- file.path(MR, "results")
PLT  <- file.path(MR, "plots")
dir.create(RES, showWarnings = FALSE); dir.create(PLT, showWarnings = FALSE)

traits <- c("bmi","waist","hip","whradjbmi","fatmass","fatperc")
labs   <- c(bmi="BMI", waist="Waist circumference", hip="Hip circumference",
            whradjbmi="WHRadjBMI", fatmass="Fat mass", fatperc="Fat percentage")

AF_TOL    <- 0.20   # gross allele-frequency inconsistency threshold
PAL_LO    <- 0.40   # palindrome ambiguous-frequency window
PAL_HI    <- 0.60

outcome <- fread(file.path(HARM, "outcome_chd_extract.tsv"))
# outcome cols: chr_pos CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ BETA SE LOG10P P N
out2 <- outcome[, .(key = chr_pos, o_ea = toupper(ALLELE1), o_oa = toupper(ALLELE0),
                    o_eaf = as.numeric(A1FREQ), o_beta = as.numeric(BETA),
                    o_se = as.numeric(SE), o_p = as.numeric(P), o_n = as.numeric(N))]

is_palin <- function(a, b) paste0(a, b) %in% c("AT","TA","CG","GC")

all_results <- list(); all_diag <- list()

for (tr in traits) {
  message("==== ", tr, " ====")
  exp <- fread(file.path(INST, paste0(tr, ".instruments.tsv")))
  exp2 <- exp[, .(SNP, chr = CHR, pos = POS,
                  ea = toupper(Allele2), oa = toupper(Allele1),
                  eaf = as.numeric(AF_Allele2), beta = as.numeric(BETA),
                  se = as.numeric(SE), pval = as.numeric(p.value), n = as.numeric(N),
                  key = paste0(CHR, ":", POS))]
  n_inst <- nrow(exp2)

  m <- merge(exp2, out2, by = "key", all.x = TRUE)
  n_missing_out <- length(setdiff(exp2$key, out2$key))      # not present in outcome at all
  mf <- m[!is.na(o_beta)]

  # allele matching (same orientation or swapped)
  mf[, same := ea == o_ea & oa == o_oa]
  mf[, flip := ea == o_oa & oa == o_ea]
  mf[, matched := same | flip]
  n_unmatched <- nrow(mf[matched == FALSE])

  mf <- mf[matched == TRUE]
  # harmonise outcome to exposure effect allele
  mf[, o_beta_h := ifelse(flip, -o_beta, o_beta)]
  mf[, o_eaf_h  := ifelse(flip, 1 - o_eaf, o_eaf)]
  # de-dup any multiallelic position collisions (keep one matched record per SNP)
  setorder(mf, SNP); mf <- mf[!duplicated(SNP)]

  # palindrome + frequency QC
  mf[, palin := is_palin(ea, oa)]
  mf[, af_diff := abs(eaf - o_eaf_h)]
  mf[, drop_palin_ambig := palin & eaf > PAL_LO & eaf < PAL_HI]
  mf[, drop_af := af_diff > AF_TOL]                          # catches strand-flipped palindromes too
  mf[, keep := !drop_palin_ambig & !drop_af]

  n_palin_ambig <- nrow(mf[drop_palin_ambig == TRUE])
  n_af_drop     <- nrow(mf[drop_af == TRUE & drop_palin_ambig == FALSE])
  dat <- mf[keep == TRUE]
  n_keep <- nrow(dat)

  fwrite(mf, file.path(HARM, paste0(tr, ".harmonised.tsv")), sep = "\t")

  # instrument strength
  dat[, Fstat := (beta / se)^2]
  meanF <- mean(dat$Fstat); minF <- min(dat$Fstat)

  # ---- MR ----
  mri <- mr_input(bx = dat$beta, bxse = dat$se, by = dat$o_beta_h, byse = dat$o_se,
                  snps = dat$SNP, exposure = labs[[tr]], outcome = "CHD (G&H)")
  ivw   <- mr_ivw(mri, model = "random")
  egg   <- mr_egger(mri)
  wmed  <- mr_median(mri, weighting = "weighted")
  wmode <- mr_mbe(mri, weighting = "weighted")

  mkrow <- function(method, b, se, p, nsnp) data.table(
    trait = tr, label = labs[[tr]], method = method, nSNP = nsnp,
    b = b, se = se, pval = p,
    OR = exp(b), OR_lci = exp(b - 1.96 * se), OR_uci = exp(b + 1.96 * se))

  res <- rbindlist(list(
    mkrow("IVW (random)",      ivw@Estimate,  ivw@StdError,      ivw@Pvalue,      ivw@SNPs),
    mkrow("Weighted median",   wmed@Estimate, wmed@StdError,     wmed@Pvalue,     wmed@SNPs),
    mkrow("Weighted mode",     wmode@Estimate,wmode@StdError,    wmode@Pvalue,    wmode@SNPs),
    mkrow("MR-Egger",          egg@Estimate,  egg@StdError.Est,  egg@Pvalue.Est,  egg@SNPs)
  ))
  fwrite(res, file.path(RES, paste0("mr_results_", tr, ".tsv")), sep = "\t")
  all_results[[tr]] <- res

  diag <- data.table(
    trait = tr, label = labs[[tr]],
    n_instruments = n_inst, n_missing_outcome = n_missing_out,
    n_unmatched_alleles = n_unmatched, n_palindromic_ambiguous = n_palin_ambig,
    n_af_inconsistent = n_af_drop, n_retained = n_keep,
    mean_F = meanF, min_F = minF,
    Q_ivw = ivw@Heter.Stat[1], Q_ivw_p = ivw@Heter.Stat[2],
    Q_egger = egg@Heter.Stat[1], Q_egger_p = egg@Heter.Stat[2],
    egger_intercept = egg@Intercept, egger_intercept_se = egg@StdError.Int,
    egger_intercept_p = egg@Pvalue.Int)
  fwrite(diag, file.path(RES, paste0("diagnostics_", tr, ".tsv")), sep = "\t")
  all_diag[[tr]] <- diag

  # ---- 2x2 grid, panels A) scatter B) forest C) leave-one-out D) funnel ----.
  big <- n_keep > 30
  base_theme <- theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
                      axis.title = element_text(size = 12), axis.text = element_text(size = 10),
                      plot.margin = margin(6, 10, 6, 6))
  pA <- mr_plot(mri, line = "ivw", orientate = TRUE, interactive = FALSE) +
        ggtitle("Genetic associations") + base_theme
  pB <- mr_forest(mri, methods = c("ivw","median","egger"), ordered = TRUE, snp_estimates = !big) +
        ggtitle(if (big) "Summary estimates by method" else "Per-variant & summary estimates") + base_theme
  if (!big) {
    pC <- mr_loo(mri) + ggtitle("Leave-one-out (IVW)") + base_theme
  } else {
    Sxy <- sum(dat$beta * dat$o_beta_h / dat$o_se^2); Sxx <- sum(dat$beta^2 / dat$o_se^2)
    loo_b <- (Sxy - dat$beta * dat$o_beta_h / dat$o_se^2) / (Sxx - dat$beta^2 / dat$o_se^2)
    loo <- data.frame(b = sort(loo_b), idx = seq_len(n_keep))
    pC <- ggplot(loo, aes(b, idx)) +
      geom_vline(xintercept = ivw@Estimate, linetype = 2, colour = "firebrick", linewidth = 0.7) +
      geom_point(size = 0.8, colour = "grey20") +
      labs(x = "Leave-one-out IVW estimate (log-OR)",
           y = "Variants (each omitted once, ordered)", title = "Leave-one-out (IVW)") +
      annotate("text", x = ivw@Estimate, y = n_keep, label = " full IVW", hjust = 0,
               colour = "firebrick", size = 3.6, fontface = "italic") +
      theme_bw() + base_theme +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
            panel.grid.minor = element_blank())
  }
  pD <- mr_funnel(mri) + ggtitle("Funnel plot") + base_theme
  fig <- (pA | pB) / (pC | pD) +
    plot_annotation(tag_levels = list(c("A)","B)","C)","D)")),
      title = sprintf("%s → CHD (G&H)   (%d genetic instruments)", labs[[tr]], n_keep),
      theme = theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5,
                                              margin = margin(b = 4)))) &
    theme(plot.tag = element_text(face = "bold", size = 16), plot.tag.position = c(0.01, 0.99))
  try(ggsave(file.path(PLT, paste0(tr, "_MR_2x2.png")), fig, width = 13, height = 11,
             dpi = 300, bg = "white"), silent = TRUE)

  message(sprintf("  retained %d/%d instruments; IVW OR=%.3f (%.3f-%.3f) p=%.2e; meanF=%.0f",
                  n_keep, n_inst, exp(ivw@Estimate), exp(ivw@Estimate-1.96*ivw@StdError),
                  exp(ivw@Estimate+1.96*ivw@StdError), ivw@Pvalue, meanF))
}

fwrite(rbindlist(all_results), file.path(RES, "mr_all_traits.tsv"), sep = "\t")
fwrite(rbindlist(all_diag),    file.path(RES, "diagnostics_all_traits.tsv"), sep = "\t")
message("DONE")
