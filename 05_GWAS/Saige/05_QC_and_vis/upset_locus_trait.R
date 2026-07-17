#!/usr/bin/env Rscript
# =============================================================================
# UpSet plot of locus x trait sharing (45 consolidated loci, 6 adiposity traits)
# + 6x6 pairwise co-occurrence matrix (the direct answer to Adam's question).
#
# Reuses the existing 45x6 binary matrix (locus_trait_matrix_45x6.tsv); the
# locus assignment is NOT re-derived here.
#
# Route C (hand-built, ggplot2 + patchwork): ComplexHeatmap / ggupset / UpSetR /
# ComplexUpset are all unavailable in this environment.
#
# IMPORTANT distinction (stated in captions):
#   * UpSet bars = EXCLUSIVE intersections: "BMI+Hip = 8" means loci in BMI and
#     hip and NO other trait.
#   * Pairwise matrix cell(i,j) = loci with a lead in both i and j regardless of
#     other traits (INCLUSIVE overlap). These are different quantities.
# =============================================================================

suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })

DIR <- "/home/jt962/rds/hpc-work/believe_adiposity/05_GWAS/Saige/05_post_gwas/HDS/locus_trait_sharing"
MAT <- file.path(DIR, "locus_trait_matrix_45x6.tsv")

# ============================================================================
# Stage 0 : read matrix, GATE 0
# ============================================================================
dt <- fread(MAT)                                   # columns from header
trait_cols <- c("BMI","Waist","Hip","WHRadjBMI","Fat mass","Fat %")
stopifnot(all(trait_cols %in% names(dt)))
M <- as.matrix(dt[, ..trait_cols]); storage.mode(M) <- "integer"
rownames(M) <- dt$locus

if (nrow(M) != 45L) stop(sprintf("STOP: matrix has %d rows, not 45.", nrow(M)))
if (sum(M) != 110L) stop(sprintf("STOP: total cells = %d, not 110 (leads).", sum(M)))
cat("== GATE 0 == rows =", nrow(M), " trait cols =", ncol(M),
    " total cells =", sum(M), " (reconciles to 110 leads)\n")
cat("set sizes (loci per trait):\n"); print(colSums(M))
cat("\n")

# ============================================================================
# Stage 2 : exclusive intersections (UpSet semantics)
# ============================================================================
combo_str <- apply(M, 1, function(r) paste(trait_cols[r == 1], collapse=" + "))
combo_n   <- rowSums(M)
int_tab <- as.data.table(table(combo_str))[, .(trait_combination = combo_str, n_loci = N)]
# stable order: size desc, then #traits asc, then alphabetical
int_tab[, ntr := lengths(strsplit(trait_combination, " \\+ "))]
setorder(int_tab, -n_loci, ntr, trait_combination)
int_tab[, ntr := NULL]

fwrite(int_tab, file.path(DIR, "intersection_sizes.tsv"), sep="\t")
cat("== Stage 2: exclusive intersection sizes (desc) ==\n"); print(int_tab)

# ---- GATE 2 : sums + known values -------------------------------------------
stopifnot(sum(int_tab$n_loci) == 45L)
known <- c("BMI + Hip"=8, "BMI"=6,
           "BMI + Waist + Hip + Fat mass + Fat %"=6,
           "WHRadjBMI"=4, "Fat %"=3)
for (k in names(known)) {
  got <- int_tab[trait_combination==k, n_loci]
  if (length(got)==0 || got != known[[k]])
    stop(sprintf("STOP GATE 2: '%s' = %s, expected %d.",
                 k, ifelse(length(got)==0,"absent",got), known[[k]]))
}
cat("\n== GATE 2 PASSED == intersections sum to 45; known values reproduced.\n\n")

# ============================================================================
# Stage 3 : pairwise co-occurrence (INCLUSIVE) 6x6
# ============================================================================
cooc <- t(M) %*% M                                  # cell(i,j)= loci in both i & j
co_dt <- data.table(trait = rownames(cooc), as.data.table(cooc))
fwrite(co_dt, file.path(DIR, "pairwise_cooccurrence_6x6.tsv"), sep="\t")
cat("== Stage 3: pairwise co-occurrence (diag = per-trait loci) ==\n")
print(cooc)
cat("\n")

# ============================================================================
# Build the hand-built UpSet (Route C)
# ============================================================================
ACC <- "#0072B2"; GREY <- "grey85"
set_sizes <- colSums(M)
# sets ordered by size; largest at TOP of the matrix
trait_lv  <- names(sort(set_sizes, decreasing = FALSE))
int_lv    <- int_tab$trait_combination                # x order (size desc)
K <- length(int_lv)

# membership long table (intersection x trait)
memb <- CJ(intersection = int_lv, trait = trait_lv)
memb[, intersection := factor(intersection, levels = int_lv)]
memb[, trait := factor(trait, levels = trait_lv)]
memb[, member := as.integer(mapply(function(ic, tr)
        tr %in% strsplit(ic, " \\+ ")[[1]], as.character(intersection), as.character(trait)))]
# vertical connectors: min/max member row per intersection
conn <- memb[member==1, .(ymin = min(as.integer(trait)),
                          ymax = max(as.integer(trait))), by = intersection]

int_tab[, intersection := factor(trait_combination, levels = int_lv)]

# --- (1) top: intersection size bars ---
p_top <- ggplot(int_tab, aes(intersection, n_loci)) +
  geom_col(fill = ACC, width = 0.65) +
  geom_text(aes(label = n_loci), vjust = -0.3, size = 2.8) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(y = "Loci in\nintersection", x = NULL) +
  theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 8))

# --- (2) bottom-right: membership dot matrix ---
p_mat <- ggplot(memb, aes(intersection, trait)) +
  geom_point(aes(colour = factor(member)), size = 3) +
  geom_segment(data = conn, aes(x = intersection, xend = intersection,
               y = ymin, yend = ymax), colour = ACC, linewidth = 0.7,
               inherit.aes = FALSE) +
  geom_point(data = memb[member==1], colour = ACC, size = 3) +
  scale_colour_manual(values = c("0"=GREY, "1"=ACC), guide = "none") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.major.y = element_line(colour = "grey93"),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 8))

# --- (3) bottom-left: set-size bars (loci per trait) ---
ss <- data.table(trait = factor(names(set_sizes), levels = trait_lv),
                 n = as.integer(set_sizes))
p_set <- ggplot(ss, aes(n, trait)) +
  geom_col(fill = "grey55", width = 0.6) +
  geom_text(aes(label = n), hjust = 1.2, colour = "white", size = 2.7) +
  scale_x_reverse(expand = expansion(mult = c(0.12, 0))) +
  labs(x = "Loci per\ntrait", y = NULL) +
  theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), axis.title.x = element_text(size = 8))

cap <- paste0(
  "45 consolidated loci x 6 adiposity traits. Denominator: 110 distance-pruned lead signals ",
  "(not 115 independent signals). Loci = +/-500 kb cross-trait single-linkage merge of the 110 ",
  "leads. UpSet bars are EXCLUSIVE intersections (e.g. 'BMI + Hip' = loci with leads for BMI and ",
  "hip and no other trait); for pairwise overlap regardless of other traits see the 6x6 ",
  "co-occurrence matrix. All 19 observed intersections shown (no truncation). MHC (chr6:30-32.6 Mb) ",
  "kept as 3 loci (not collapsed to a single HLA region).")

design <- "
#A
BC
"
upset <- p_top + p_set + p_mat +
  plot_layout(design = design, widths = c(1, 4), heights = c(1.6, 3)) +
  plot_annotation(
    title = "UpSet: locus x trait sharing across the six adiposity traits",
    caption = paste(strwrap(cap, width = 118), collapse = "\n"),
    theme = theme(plot.title = element_text(face = "bold", size = 12),
                  plot.caption = element_text(hjust = 0, size = 6.5, colour = "grey35")))

ggsave(file.path(DIR, "upset_locus_trait.pdf"), upset, width = 11, height = 7, device = cairo_pdf)
ggsave(file.path(DIR, "upset_locus_trait.png"), upset, width = 11, height = 7, dpi = 200)

# ---- pairwise tile plot (cheap) ---------------------------------------------
co_long <- as.data.table(as.table(cooc)); setnames(co_long, c("t1","t2","n"))
co_long[, t1 := factor(t1, levels = trait_cols)]
co_long[, t2 := factor(t2, levels = rev(trait_cols))]
p_co <- ggplot(co_long, aes(t1, t2, fill = n)) +
  geom_tile(colour = "white", linewidth = 1) +
  geom_text(aes(label = n), size = 3.4,
            colour = ifelse(co_long$n > max(cooc)*0.6, "white", "grey15")) +
  scale_fill_viridis_c(option = "mako", direction = -1, name = "Loci") +
  labs(title = "Pairwise locus co-occurrence between traits",
       subtitle = "Cell = loci with a lead in BOTH traits (inclusive of other traits); diagonal = loci per trait",
       caption = "Distinct from the UpSet exclusive intersections. 45 loci, 110 leads. MHC = 3 loci.",
       x = NULL, y = NULL) +
  coord_equal() +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 8.5, colour = "grey35"),
        plot.caption = element_text(hjust = 0, size = 7, colour = "grey35"))
ggsave(file.path(DIR, "pairwise_cooccurrence.pdf"), p_co, width = 6.5, height = 6, device = cairo_pdf)
ggsave(file.path(DIR, "pairwise_cooccurrence.png"), p_co, width = 6.5, height = 6, dpi = 200)

cat("== outputs written to", DIR, "==\n")
cat("  upset_locus_trait.{pdf,png}\n  pairwise_cooccurrence.{pdf,png}\n",
    "  intersection_sizes.tsv\n  pairwise_cooccurrence_6x6.tsv\n", sep="")
