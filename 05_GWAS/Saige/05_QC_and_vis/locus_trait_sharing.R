#!/usr/bin/env Rscript
# =============================================================================
# Locus x trait sharing: how the consolidated physical loci distribute across
# the six adiposity traits (per Adam's Results comment).
#
# DENOMINATOR: 110 distance-pruned LEAD signals (NOT the 115 signals, which add
# 5 conditionally independent secondaries and would double-count FTO / chr5).
#
# LOCUS CONVENTION: +/-500 kb cross-trait single-linkage merge of the 110 leads
# => 45 loci. This is a TIGHTER consolidation than the manuscript's 75 "physical
# loci" (the annotation_table 'locus' column, which merges only near-identical
# positions). Both are stated so nothing silently contradicts the text.
#
# Heatmap: ggplot2 geom_tile (ComplexHeatmap/pheatmap unavailable in this env).
# UpSet:   UpSetR / ComplexUpset unavailable -> not hand-built; intersection
#          sizes reported as a table instead.
# =============================================================================

suppressMessages({ library(data.table); library(ggplot2) })

LEADS <- "/home/jt962/rds/hpc-work/believe_adiposity/05_GWAS/Saige/05_post_gwas/HDS/all_traits.leads.tsv"
ANN   <- "/home/jt962/rds/hpc-work/believe_adiposity/07_annotation/annotation_table.tsv"
OUT   <- "/home/jt962/rds/hpc-work/believe_adiposity/05_GWAS/Saige/05_post_gwas/HDS/locus_trait_sharing"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
WIN <- 500000L

# ---- trait display (UK spelling), fixed column order ------------------------
trait_disp <- c(bmi_INT="BMI", waist_INT="Waist", hip_INT="Hip",
                whr_bmi_adj_INT="WHRadjBMI", fatmass_INT="Fat mass",
                fatperc_INT="Fat %")
trait_order <- unname(trait_disp)

# ============================================================================
# Stage 0/1 : recompute 45 loci (+/-500kb cross-trait single-linkage), GATE
# ============================================================================
leads <- fread(LEADS)                     # columns from header, not memory
stopifnot(nrow(leads) == 110)
setorder(leads, CHR, POS)
leads[, locus_merge := {
  id <- integer(.N); cur <- 1L; id[1] <- 1L
  if (.N > 1) for (i in 2:.N)
    if (POS[i]-POS[i-1] <= WIN) id[i] <- cur else { cur <- cur+1L; id[i] <- cur }
  id
}, by = CHR]
leads[, locus_key := paste0("chr", CHR, "_", locus_merge)]

n_loci <- uniqueN(leads$locus_key)
cat("== GATE 0: locus count from +/-500kb merge of 110 leads =", n_loci, "==\n")
if (n_loci != 45) stop(sprintf("STOP: merge yielded %d loci, not 45. Not plotting.", n_loci))
cat("PASS: exactly 45 loci from 110 leads.\n\n")

# ---- nearest gene per lead (match to annotation by CHR+POS+trait) -----------
ann <- fread(ANN)
setnames(ann, "trait", "TRAIT")
lead_gene <- ann[signal_type=="lead", .(CHR=chr, POS=pos, TRAIT, nearest_gene)]
leads <- merge(leads, lead_gene, by=c("CHR","POS","TRAIT"), all.x=TRUE, sort=FALSE)

# ---- locus-level summary: representative = min-P lead -----------------------
loci <- leads[, {
  k <- which.min(P)
  .(CHR=CHR[1], repr_pos=POS[k], repr_gene=nearest_gene[k],
    n_traits=uniqueN(TRAIT), n_leads=.N, minP=min(P))
}, by=locus_key]
setorder(loci, CHR, repr_pos)

# ---- multi-signal loci (contain a conditional secondary) -> asterisk --------
sec <- ann[signal_type=="secondary", .(chr, pos)]
loci[, multi_signal := FALSE]
for (i in seq_len(nrow(sec)))
  loci[CHR==sec$chr[i] & abs(repr_pos-sec$pos[i])<=WIN, multi_signal := TRUE]

# row label: chr:pos . gene (* if multi-signal)
loci[, gene_lab := ifelse(is.na(repr_gene) | repr_gene=="", "", paste0(" · ", repr_gene))]
loci[, row_label := sprintf("chr%d:%d%s%s", CHR, repr_pos, gene_lab,
                            ifelse(multi_signal, " *", ""))]

# ============================================================================
# Stage 1 report
# ============================================================================
cat("== Stage 1 report ==\n")
cat("Loci:", nrow(loci), " | denominator: 110 leads (not 115 signals)\n")
cat("Traits per locus  min/median/max:",
    min(loci$n_traits), median(loci$n_traits), max(loci$n_traits), "\n")
cat("Single-trait loci:", sum(loci$n_traits==1), "\n")
top <- loci[order(-n_traits, minP)][1]
cat(sprintf("Most-trait locus : %s  (%d traits, minP=%.2g)\n",
            top$row_label, top$n_traits, top$minP))
cat("Multi-signal loci (asterisked):", sum(loci$multi_signal), "->",
    paste(loci[multi_signal==TRUE, row_label], collapse=" ; "), "\n")
mhc <- loci[CHR==6 & repr_pos>=29e6 & repr_pos<=34e6]
cat("MHC (chr6:29-34Mb) rows in matrix:", nrow(mhc),
    "-> convention = +/-500kb cross-trait merge (NOT collapsed to single HLA)\n\n")

# ============================================================================
# Stage 1 : build 45 x 6 binary matrix (+ -log10P matrix)
# ============================================================================
leads[, trait_lab := trait_disp[TRAIT]]
M <- matrix(0L, nrow=nrow(loci), ncol=length(trait_order),
            dimnames=list(loci$locus_key, trait_order))
L <- matrix(NA_real_, nrow=nrow(loci), ncol=length(trait_order),
            dimnames=list(loci$locus_key, trait_order))
for (i in seq_len(nrow(leads))) {
  r <- leads$locus_key[i]; cc <- leads$trait_lab[i]
  M[r, cc] <- 1L
  L[r, cc] <- max(L[r, cc], -log10(leads$P[i]), na.rm=TRUE)
}

# ---- row order by clustering on binary matrix -------------------------------
row_order_method <- "hierarchical clustering (binary/Jaccard distance, complete linkage)"
ord <- tryCatch({
  hc <- hclust(dist(M, method="binary"), method="complete"); hc$order
}, error=function(e) NULL)
if (is.null(ord) || anyNA(ord)) {
  ord <- order(loci$CHR, loci$repr_pos)
  row_order_method <- "chromosome and position (clustering unstable)"
}
loci_ord <- loci[ord]
M <- M[ord, , drop=FALSE]; L <- L[ord, , drop=FALSE]

# ============================================================================
# Stage 4 outputs: matrix TSV
# ============================================================================
mat_out <- data.table(locus=loci_ord$locus_key, chr=loci_ord$CHR,
                      pos=loci_ord$repr_pos, nearest_gene=loci_ord$repr_gene,
                      n_traits=loci_ord$n_traits, multi_signal=loci_ord$multi_signal)
mat_out <- cbind(mat_out, as.data.table(M))
fwrite(mat_out, file.path(OUT, "locus_trait_matrix_45x6.tsv"), sep="\t")

# ---- intersection-size table (trait combinations, desc) ---------------------
combo <- apply(M, 1, function(r) paste(trait_order[r==1], collapse=" + "))
isize <- as.data.table(table(combo))[order(-N)]
setnames(isize, c("trait_combination","n_loci"))
fwrite(isize, file.path(OUT, "intersection_sizes.tsv"), sep="\t")
cat("== Intersection sizes (trait combination -> n loci), desc ==\n")
print(isize)
cat("\nUpSet plot: UpSetR / ComplexUpset NOT installed -> plot not produced ",
    "(not hand-built, per brief). Table above is the checkable substitute.\n\n", sep="")

# ============================================================================
# Stage 2 : heatmap (ggplot geom_tile) -- binary (primary) + -log10P version
# ============================================================================
lab_levels <- loci_ord$row_label            # bottom->top = plotting order
df <- data.table(
  locus = factor(rep(loci_ord$row_label, times=length(trait_order)),
                 levels=lab_levels),
  trait = factor(rep(trait_order, each=nrow(M)), levels=trait_order),
  present = as.vector(M),
  logp    = as.vector(L))

cap_common <- paste0(
  "45 consolidated loci x 6 adiposity traits. Denominator: 110 distance-pruned lead ",
  "signals (not 115 independent signals). Loci = +/-500 kb cross-trait single-linkage ",
  "merge of the 110 leads; this is tighter than the manuscript's 75 'physical loci'. ",
  "Rows ordered by ", row_order_method, ". * = locus carries a conditionally ",
  "independent secondary signal. MHC (chr6:30-32.6 Mb) kept as 3 loci (not collapsed).")

p_bin <- ggplot(df, aes(trait, locus, fill=factor(present))) +
  geom_tile(colour="white", linewidth=0.6) +
  scale_fill_manual(values=c("0"="grey92","1"="#0072B2"),
                    labels=c("absent","lead present"), name=NULL) +
  scale_x_discrete(position="top") +
  labs(title="Locus x trait sharing across the six adiposity traits",
       subtitle="Binary: a genome-wide lead for the trait at the locus",
       caption=strwrap(cap_common, width=110) |> paste(collapse="\n"),
       x=NULL, y=NULL) +
  theme_minimal(base_size=10) +
  theme(axis.text.y=element_text(size=6.6),
        axis.text.x.top=element_text(face="bold", angle=0),
        panel.grid=element_blank(), legend.position="bottom",
        plot.caption=element_text(hjust=0, size=6.5, colour="grey35"),
        plot.title=element_text(face="bold"))
ggsave(file.path(OUT,"heatmap_locus_trait_binary.pdf"), p_bin,
       width=7.2, height=10, device=cairo_pdf)
ggsave(file.path(OUT,"heatmap_locus_trait_binary.png"), p_bin,
       width=7.2, height=10, dpi=200)

p_lp <- ggplot(df[present==1], aes(trait, locus, fill=logp)) +
  geom_tile(data=df, aes(trait, locus), fill="grey92", colour="white", linewidth=0.6) +
  geom_tile(colour="white", linewidth=0.6) +
  scale_fill_viridis_c(option="magma", direction=-1,
                       name=expression(-log[10]*italic(P)), end=0.92) +
  scale_x_discrete(position="top") +
  labs(title="Locus x trait sharing (shaded by lead significance)",
       subtitle=expression("Fill = "*-log[10]*italic(P)*" of the trait's lead at the locus"),
       caption=strwrap(cap_common, width=110) |> paste(collapse="\n"),
       x=NULL, y=NULL) +
  theme_minimal(base_size=10) +
  theme(axis.text.y=element_text(size=6.6),
        axis.text.x.top=element_text(face="bold"),
        panel.grid=element_blank(), legend.position="bottom",
        plot.caption=element_text(hjust=0, size=6.5, colour="grey35"),
        plot.title=element_text(face="bold"))
ggsave(file.path(OUT,"heatmap_locus_trait_logP.pdf"), p_lp,
       width=7.2, height=10, device=cairo_pdf)
ggsave(file.path(OUT,"heatmap_locus_trait_logP.png"), p_lp,
       width=7.2, height=10, dpi=200)

cat("== outputs written to", OUT, "==\n")
cat("  heatmap_locus_trait_binary.{pdf,png}\n  heatmap_locus_trait_logP.{pdf,png}\n")
cat("  locus_trait_matrix_45x6.tsv\n  intersection_sizes.tsv\n")
