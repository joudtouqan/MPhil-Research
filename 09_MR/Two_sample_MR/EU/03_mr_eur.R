#!/usr/bin/env Rscript
# ============================================================================
# EUR comparator step 3: harmonise EUR exposures vs Nikpay 2015 CAD and run
# two-sample MR. Same estimator set and allele-aware logic as the SAS arm
#   - harmonise on rsID (not chr:pos)
#   - outcome MAF merged from 1000G EUR (Nikpay has none)
#   - outcome beta_cad = log-odds -> OR = exp(b)
# Exposures: bmi (Yengo, per-SNP N), whr & fatperc (Pan-UKB EUR, irnt per-SD).
# Nikpay: ncase=60801, ncontrol=123504.
# ============================================================================
.libPaths("/rds/user/jt962/hpc-work/software/R/4.6.0_libs")
suppressMessages({library(data.table); library(MendelianRandomization); library(ggplot2); library(patchwork)})

EUR  <- "/home/jt962/rds/hpc-work/believe_adiposity/09_MR/eur_comparator"
EXPO <- file.path(EUR, "exposures"); INST <- file.path(EUR, "instruments")
HARM <- file.path(EUR, "harmonised"); RES <- file.path(EUR, "results"); PLT <- file.path(EUR, "plots")
dir.create(RES, showWarnings = FALSE); dir.create(PLT, showWarnings = FALSE)

traits <- c("bmi","whr","fatperc")
labs   <- c(bmi="BMI (Yengo)", whr="WHR (Pan-UKB EUR)", fatperc="Body fat % (Pan-UKB EUR)")
NCASE <- 60801; NCTRL <- 123504

AF_TOL <- 0.20; PAL_LO <- 0.42; PAL_HI <- 0.58            # brief: EUR palindrome window 0.42-0.58
is_palin <- function(a,b) paste0(a,b) %in% c("AT","TA","CG","GC")

outcome <- fread(file.path(HARM, "outcome_cad_extract.tsv"))   # rsid o_chr o_pos o_ea o_oa o_beta o_se o_p
frq     <- fread(file.path(HARM, "eur_freq_extract.tsv"))       # rsid ref_REF ref_ALT ref_ALT_FREQ
# outcome EAF from EUR ref: freq of the outcome effect allele (o_ea)
om <- merge(outcome, frq, by = "rsid", all.x = TRUE)
om[, o_eaf := fifelse(o_ea == ref_ALT, as.numeric(ref_ALT_FREQ),
                fifelse(o_ea == ref_REF, 1 - as.numeric(ref_ALT_FREQ), NA_real_))]
out2 <- om[, .(rsid, o_ea, o_oa, o_eaf, o_beta = as.numeric(o_beta),
               o_se = as.numeric(o_se), o_p = as.numeric(o_p))]

all_res <- list(); all_diag <- list()
for (tr in traits) {
  f <- file.path(INST, paste0(tr, ".instruments.tsv"))
  if (!file.exists(f)) { message("skip ", tr, " (no instruments file)"); next }
  exp <- fread(f)   # rsid chr pos ea oa eaf beta se p n
  exp2 <- exp[, .(rsid, chr, pos, ea = toupper(ea), oa = toupper(oa),
                  eaf = as.numeric(eaf), beta = as.numeric(beta), se = as.numeric(se),
                  pval = as.numeric(p), n = as.numeric(n))]
  n_inst <- nrow(exp2)

  m <- merge(exp2, out2, by = "rsid", all.x = TRUE)
  n_missing_out <- length(setdiff(exp2$rsid, out2$rsid))
  mf <- m[!is.na(o_beta)]

  mf[, same := ea == o_ea & oa == o_oa]
  mf[, flip := ea == o_oa & oa == o_ea]
  mf[, matched := same | flip]
  n_unmatched <- nrow(mf[matched == FALSE])
  mf <- mf[matched == TRUE]
  mf[, o_beta_h := ifelse(flip, -o_beta, o_beta)]
  mf[, o_eaf_h  := ifelse(flip, 1 - o_eaf, o_eaf)]
  setorder(mf, rsid); mf <- mf[!duplicated(rsid)]

  mf[, palin := is_palin(ea, oa)]
  mf[, af_diff := abs(eaf - o_eaf_h)]
  mf[, drop_palin_ambig := palin & eaf > PAL_LO & eaf < PAL_HI]
  mf[, drop_af := !is.na(af_diff) & af_diff > AF_TOL]
  mf[, keep := !drop_palin_ambig & !drop_af]
  n_palin_ambig <- nrow(mf[drop_palin_ambig == TRUE])
  n_af_drop     <- nrow(mf[drop_af == TRUE & drop_palin_ambig == FALSE])
  dat <- mf[keep == TRUE]
  n_keep <- nrow(dat)
  fwrite(mf, file.path(HARM, paste0(tr, ".harmonised.tsv")), sep = "\t")

  dat[, Fstat := (beta/se)^2]
  meanF <- mean(dat$Fstat); minF <- min(dat$Fstat)

  mri   <- mr_input(bx = dat$beta, bxse = dat$se, by = dat$o_beta_h, byse = dat$o_se,
                    snps = dat$rsid, exposure = labs[[tr]], outcome = "CAD (Nikpay)")
  ivw   <- mr_ivw(mri, model = "random")
  egg   <- mr_egger(mri)
  wmed  <- mr_median(mri, weighting = "weighted")
  wmode <- mr_mbe(mri, weighting = "weighted")

  mkrow <- function(method,b,se,p,nsnp) data.table(
    trait=tr,label=labs[[tr]],method=method,nSNP=nsnp,b=b,se=se,pval=p,
    OR=exp(b),OR_lci=exp(b-1.96*se),OR_uci=exp(b+1.96*se))
  res <- rbindlist(list(
    mkrow("IVW (random)",   ivw@Estimate,  ivw@StdError,    ivw@Pvalue,     ivw@SNPs),
    mkrow("Weighted median",wmed@Estimate, wmed@StdError,   wmed@Pvalue,    wmed@SNPs),
    mkrow("Weighted mode",  wmode@Estimate,wmode@StdError,  wmode@Pvalue,   wmode@SNPs),
    mkrow("MR-Egger",       egg@Estimate,  egg@StdError.Est,egg@Pvalue.Est, egg@SNPs)))
  fwrite(res, file.path(RES, paste0("mr_results_", tr, ".tsv")), sep = "\t")
  all_res[[tr]] <- res

  # Steiger (exposure SD scale R2; outcome binary liability R2)
  dat[, r2_exp := 2*eaf*(1-eaf)*beta^2]
  dat[, vg_out := 2*o_eaf_h*(1-o_eaf_h)*o_beta_h^2]
  dat[, r2_out := vg_out/(vg_out + pi^2/3)]
  R2 <- sum(dat$r2_exp); Nexp <- dat$n[1]
  Foverall <- (R2/(1-R2))*((Nexp - n_keep - 1)/n_keep)

  diag <- data.table(trait=tr,label=labs[[tr]],ncase=NCASE,ncontrol=NCTRL,
    n_instruments=n_inst,n_missing_outcome=n_missing_out,n_unmatched_alleles=n_unmatched,
    n_palindromic_ambiguous=n_palin_ambig,n_af_inconsistent=n_af_drop,n_retained=n_keep,
    mean_F=meanF,min_F=minF,R2_exposure=R2,F_overall=Foverall,N_exposure=Nexp,
    n_steiger_correct=sum(dat$r2_exp>dat$r2_out),
    Q_ivw=ivw@Heter.Stat[1],Q_ivw_p=ivw@Heter.Stat[2],
    egger_intercept=egg@Intercept,egger_intercept_se=egg@StdError.Int,egger_intercept_p=egg@Pvalue.Int)
  fwrite(diag, file.path(RES, paste0("diagnostics_", tr, ".tsv")), sep = "\t")
  all_diag[[tr]] <- diag

  # ----  2x2 grid, panels A) scatter B) forest C) leave-one-out D) funnel ----
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
      title = sprintf("%s → CAD (Nikpay)   (%d genetic instruments)", labs[[tr]], n_keep),
      theme = theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5,
                                              margin = margin(b = 4)))) &
    theme(plot.tag = element_text(face = "bold", size = 16), plot.tag.position = c(0.01, 0.99))
  try(ggsave(file.path(PLT, paste0(tr, "_MR_2x2.png")), fig, width = 13, height = 11,
             dpi = 300, bg = "white"), silent = TRUE)

  message(sprintf("  %s: kept %d/%d; IVW OR=%.3f (%.3f-%.3f) p=%.2e; meanF=%.0f; R2=%.4f Foverall=%.0f",
    tr,n_keep,n_inst,exp(ivw@Estimate),exp(ivw@Estimate-1.96*ivw@StdError),
    exp(ivw@Estimate+1.96*ivw@StdError),ivw@Pvalue,meanF,R2,Foverall))
}
fwrite(rbindlist(all_res),  file.path(RES,"mr_all_traits.tsv"), sep="\t")
fwrite(rbindlist(all_diag), file.path(RES,"diagnostics_all_traits.tsv"), sep="\t")
message("DONE")
