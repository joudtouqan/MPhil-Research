#!/usr/bin/env Rscript
# Align each BELIEVE independent signal to each comparator (matched trait only),
# orienting the comparator effect allele to effect allele (flip beta + 1-EAF
# when the comparator reports my NEA as its effect allele).
suppressMessages(library(data.table))

D   <- "/home/jt962/rds/hpc-work/believe_adiposity/08_cross_ancestry"
OUT <- file.path(D, "output")
EX  <- file.path(OUT, "extracts")

mas <- fread(file.path(OUT, "master_signals.tsv"), colClasses=list(character="rsID"))
h19 <- fread(file.path(D, "data/sigid_hg19.map"),
             header=FALSE, col.names=c("signal_id","h19chr","h19pos"))
mas <- merge(mas, h19[, .(signal_id, hg19_key=paste0(h19chr,":",h19pos))],
             by="signal_id", all.x=TRUE, sort=FALSE)

# hg38 key (native build of BELIEVE signals AND of the Genes & Health 51k GWAS)
mas[, hg38_key := paste0(chr, ":", pos)]
up <- function(x) toupper(x)
mas[, `:=`(EA=up(EA), NEA=up(NEA))]
mas[, palindromic := (paste0(EA,NEA) %in% c("AT","TA","CG","GC"))]

# comparator configs ---------------------------------------------------------
cfgs <- list(
  GIANT_BMI       = list(file="pulit_bmi.tsv",       trait="bmi_INT",         key="hg19_key", build="hg19", matched_by="pos+allele",
                          ea="Tested", nea="Other", freqEA="Freq", beta="BETA", se="SE", p="P", n="N", info="INFO"),
  GIANT_WHRadjBMI = list(file="pulit_whradjbmi.tsv", trait="whr_bmi_adj_INT", key="hg19_key", build="hg19", matched_by="pos+allele",
                          ea="Tested", nea="Other", freqEA="Freq", beta="BETA", se="SE", p="P", n="N", info="INFO"),
  GIANT_WC        = list(file="shungin_wc.tsv",      trait="waist_INT",       key="rsID",     build="rsID(GRCh37-era)", matched_by="rsID+allele",
                          ea="Allele1", nea="Allele2", freqEA="FreqAllele1", beta="b", se="se", p="p", n="N", info=NA),
  GIANT_HIP       = list(file="shungin_hip.tsv",     trait="hip_INT",         key="rsID",     build="rsID(GRCh37-era)", matched_by="rsID+allele",
                          ea="Allele1", nea="Allele2", freqEA="FreqAllele1", beta="b", se="se", p="p", n="N", info=NA),
  BBJ_BMI         = list(file="bbj_bmi.tsv",          trait="bmi_INT",         key="hg19_key", build="hg19", matched_by="pos+allele",
                          ea="ALT", nea="REF", freqEA="Frq", beta="BETA", se="SE", p="P", n=NA, info="Rsq", fixedN=158284L),
  GNH_BMI         = list(file="gnh_bmi.tsv",          trait="bmi_INT",         key="hg38_key", build="hg38", matched_by="pos+allele",
                          ea="ALLELE1", nea="ALLELE0", freqEA="A1FREQ", beta="BETA", se="SE", p="P", n="N", info=NA),
  GNH_WC          = list(file="gnh_wc.tsv",           trait="waist_INT",       key="hg38_key", build="hg38", matched_by="pos+allele",
                          ea="ALLELE1", nea="ALLELE0", freqEA="A1FREQ", beta="BETA", se="SE", p="P", n="N", info=NA)
)

align_one <- function(cn, cfg) {
  fp <- file.path(EX, cfg$file)
  if (!file.exists(fp) || file.size(fp) == 0L) {
    cat(sprintf("[skip] %s: %s not present yet - comparator omitted.\n", cn, cfg$file))
    return(NULL)
  }
  comp <- fread(fp, colClasses="character")
  for (col in c(cfg$freqEA,cfg$beta,cfg$se,cfg$p,if(!is.na(cfg$n))cfg$n,if(!is.na(cfg$info))cfg$info))
    comp[[col]] <- as.numeric(comp[[col]])
  comp[, cEA := up(get(cfg$ea))]; comp[, cNEA := up(get(cfg$nea))]
  keycol <- if (cfg$key %in% c("hg19_key","hg38_key")) "chr_pos" else "rsID"
  sig <- mas[trait==cfg$trait]
  res <- rbindlist(lapply(seq_len(nrow(sig)), function(i){
    s <- sig[i]; k <- s[[cfg$key]]
    cand <- comp[get(keycol)==k]
    base <- data.table(my_signal_id=s$signal_id, trait=s$trait, comparator=cn,
      comp_build=cfg$build, matched_by=cfg$matched_by, my_EA=s$EA, my_NEA=s$NEA,
      my_EAF=s$EAF, my_BETA=s$my_BETA, my_SE=s$my_SE, my_P=s$my_P,
      palindromic=s$palindromic, nearest_gene=s$nearest_gene)
    if (!nrow(cand)) return(cbind(base, present="N", allele_status="absent",
      comp_EA=NA, comp_NEA=NA, comp_EAF=NA_real_, comp_BETA_aligned=NA_real_,
      comp_SE=NA_real_, comp_P=NA_real_, comp_N=NA_real_, comp_INFO=NA_real_, flip=NA))
    # choose the candidate whose allele set matches mine
    myset <- sort(c(s$EA,s$NEA))
    cand[, setmatch := mapply(function(a,b) identical(sort(c(a,b)), myset), cEA, cNEA)]
    m <- cand[setmatch==TRUE]
    if (!nrow(m)) { r <- cand[1]   # present at locus but alleles differ
      return(cbind(base, present="Y", allele_status="allele_mismatch",
        comp_EA=r$cEA, comp_NEA=r$cNEA, comp_EAF=NA_real_, comp_BETA_aligned=NA_real_,
        comp_SE=NA_real_, comp_P=NA_real_, comp_N=NA_real_, comp_INFO=NA_real_, flip=NA)) }
    r <- m[1]
    flip <- (r$cEA == s$NEA)            # comparator effect allele is my NEA -> flip
    b  <- r[[cfg$beta]]; fr <- r[[cfg$freqEA]]
    bA <- if (flip) -b else b
    frA<- if (flip) 1-fr else fr
    Nv <- if (!is.na(cfg$n)) r[[cfg$n]] else cfg$fixedN
    Iv <- if (!is.na(cfg$info)) r[[cfg$info]] else NA_real_
    cbind(base, present="Y", allele_status=if(flip)"matched_flipped" else "matched_same",
      comp_EA=s$EA, comp_NEA=s$NEA, comp_EAF=round(frA,4), comp_BETA_aligned=bA,
      comp_SE=r[[cfg$se]], comp_P=r[[cfg$p]], comp_N=Nv, comp_INFO=Iv, flip=flip)
  }))
  res
}

all <- rbindlist(lapply(names(cfgs), function(cn) align_one(cn, cfgs[[cn]])), fill=TRUE)

# direction + replication tiers ----------------------------------------------
all[, same_direction := ifelse(present=="Y" & allele_status %in% c("matched_same","matched_flipped"),
                               sign(comp_BETA_aligned)==sign(my_BETA), NA)]
# Bonferroni per comparator over allele-matched signals tested
all[, matched := allele_status %in% c("matched_same","matched_flipped")]
all[, m_tested := sum(matched), by=comparator]
all[, bonf_alpha := 0.05/m_tested]
all[, repl_nominal     := matched & same_direction & comp_P < 0.05]
all[, repl_bonferroni  := matched & same_direction & comp_P < bonf_alpha]
all[, repl_genomewide  := matched & same_direction & comp_P < 5e-8]

setcolorder(all, c("my_signal_id","trait","comparator","comp_build","matched_by",
  "nearest_gene","my_EA","my_NEA","my_EAF","my_BETA","my_SE","my_P",
  "comp_EA","comp_EAF","comp_BETA_aligned","comp_SE","comp_P","comp_N","comp_INFO",
  "present","allele_status","flip","palindromic","same_direction",
  "m_tested","bonf_alpha","repl_nominal","repl_bonferroni","repl_genomewide"))
setorder(all, comparator, trait, my_signal_id)
fwrite(all, file.path(OUT,"per_signal_x_comparator.tsv"), sep="\t")
cat(sprintf("Wrote per_signal_x_comparator.tsv: %d rows\n", nrow(all)))

## ---- FTO validation gate ---------------------------------------------------
cat("\n================ FTO VALIDATION GATE (rs55872725 = 16:53775211 hg38) ===========\n")
fto <- all[my_signal_id=="16:53775211:C:T" & comparator %in% c("GIANT_BMI","BBJ_BMI")]
print(fto[, .(comparator, my_EA, my_BETA, comp_EA, comp_BETA_aligned, comp_P, flip, same_direction)])
ok <- nrow(fto)>0 && all(fto$same_direction==TRUE) && all(fto$comp_P < 5e-8)
cat(sprintf("FTO gate: %s\n", if(ok) "PASS (same direction + genome-wide sig in GIANT & BBJ)" else "*** FAIL - STOP, alignment suspect ***"))

## ---- per-comparator counts -------------------------------------------------
cat("\n================ PER-COMPARATOR COUNTS ===========\n")
cnt <- all[, .(
  signals_for_trait = .N,
  present           = sum(present=="Y"),
  allele_matched    = sum(matched),
  allele_mismatch   = sum(allele_status=="allele_mismatch"),
  same_direction    = sum(same_direction==TRUE, na.rm=TRUE),
  repl_nominal      = sum(repl_nominal, na.rm=TRUE),
  repl_bonferroni   = sum(repl_bonferroni, na.rm=TRUE),
  repl_genomewide   = sum(repl_genomewide, na.rm=TRUE),
  bonf_alpha        = signif(unique(bonf_alpha),3)
), by=.(comparator, trait)]
print(cnt)
fwrite(cnt, file.path(OUT,"comparator_counts.tsv"), sep="\t")

## ---- concordance table (matched only) --------------------------------------
conc <- all[matched==TRUE, .(my_signal_id, trait, comparator, nearest_gene,
  my_EA, my_BETA, my_SE, my_P, my_EAF,
  comp_BETA_aligned, comp_SE, comp_P, comp_EAF, same_direction,
  repl_nominal, repl_bonferroni)]
fwrite(conc, file.path(OUT,"concordance_betas.tsv"), sep="\t")
cat(sprintf("\nWrote concordance_betas.tsv: %d matched signal-comparator pairs\n", nrow(conc)))

## ---- SAS-specific candidates -----------------------------------------------
eur <- all[comparator %in% c("GIANT_BMI","GIANT_WHRadjBMI","GIANT_WC","GIANT_HIP")]
sas <- eur[my_P < 5e-8 & (present=="N" | allele_status=="allele_mismatch" |
                          (matched & (same_direction==FALSE | comp_P > 0.05)))]
sas_out <- sas[, .(my_signal_id, trait, comparator, nearest_gene, my_BETA, my_P,
  present, allele_status, comp_BETA_aligned, comp_P, same_direction)]
setorder(sas_out, trait, my_signal_id)
fwrite(sas_out, file.path(OUT,"sas_specific_candidates.tsv"), sep="\t")
cat(sprintf("SAS-specific CANDIDATES (EUR comparator, genome-wide in BELIEVE): %d (see file)\n", nrow(sas_out)))
cat("\nNOTE: traits with NO external comparator: fatmass_INT, fatperc_INT (not in GIANT/BBJ).\n")
