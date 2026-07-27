# MPhil-Research

All analysis **scripts** for the BELIEVE adiposity MPhil project. This tree mirrors the
**results** tree at `~/rds/hpc-work/believe_adiposity/`, stage for stage, so each script sits at
the same `NN_stage/...` path as the outputs it produces. Each results folder there carries a
`README.md` describing the test/method behind it.

## Core analysis pipeline (script ↔ Methods section)

These  scripts are the primary analysis chain — the ones the Methods describe. Everything
else in the tree is diagnostics, figure regeneration, sensitivity re-runs, or exploratory work.
They run in the order listed; each consumes the previous stage's output.

### 1. Phenotype preparation

| Script | Methods step |
|---|---|
| [`01_phenotype_prep/outcomes_phenotype_prep.R`](01_phenotype_prep/outcomes_phenotype_prep.R) | Builds the SAIGE phenotype/covariate file. Adults (age ≥ 18) with genotypes; continuous traits regressed on age, age², sex, study site and PC1–20, then rank-based inverse-normal transformed; binary outcomes (T2D, CHD) coerced to clean 0/1; PCs merged from the FlashPCA unrelated set plus projections for relateds. Reports analysis-ready N per outcome. |

### 2. Variant QC

| Script | Methods step |
|---|---|
| [`03_variant_QC/3_make_filtered_pgen_hds.sh`](03_variant_QC/3_make_filtered_pgen_hds.sh) | Builds the analysis-ready **HDS (dosage) PGEN** by subsetting the imputed dosage data to the variant whitelist already passed by GT-based QC (MAF, HWE, missingness, biallelic SNVs). This dosage set is the substrate for GWAS, conditional analysis, MR clumping and fine-mapping. |

### 3. GWAS — SAIGE mixed model

| Script | Methods step |
|---|---|
| [`05_GWAS/Saige/01_LD_pruning/ld_pruning.sh`](05_GWAS/Saige/01_LD_pruning/ld_pruning.sh) | Relatedness-aware LD pruning for the GRM. Unrelated subset by KING coefficient < 0.022; variants filtered to MAF > 5%, HWE > 1e-10, missingness < 5%, biallelic SNVs, with MHC and the chr8/chr17 inversions excluded (hg38); pruned at `--indep-pairwise 1000 80 0.1`; the pruned variant list is then re-extracted into the **full** cohort. |
| [`05_GWAS/Saige/02_GRM/build_sparse_grm.sh`](05_GWAS/Saige/02_GRM/build_sparse_grm.sh) | Sparse genetic relationship matrix (`createSparseGRM.R`), autosomes only, relatedness cutoff 0.05, 2,000 random markers for kinship estimation. |
| [`05_GWAS/Saige/03_step1_saige/saige_step1.sh`](05_GWAS/Saige/03_step1_saige/saige_step1.sh) | **SAIGE Step 1 — null model.** Fits a null logistic/linear mixed model per trait using the sparse GRM, for the six adiposity traits (BMI, WHRadjBMI, waist, hip, fat mass, fat %). `--invNormalize=FALSE` because the INT is already applied at the phenotype-prep stage. |
| [`05_GWAS/Saige/04_step2_saige/saige_step2_HDS_perchrom.sh`](05_GWAS/Saige/04_step2_saige/saige_step2_HDS_perchrom.sh) | **SAIGE Step 2 — the primary association scan.** Single-variant score tests with saddlepoint approximation on imputed dosages; MAF ≥ 0.005, MAC ≥ 20, INFO ≥ 0.3, LOCO off. One SLURM array task per chromosome, each extracting its chromosome from the genome-wide BGEN with `bgenix` first (avoids re-streaming 307M variants 22 times) and deleting it afterwards; six phenotypes run three-at-a-time per chromosome with `.done` sentinels for restartability. |
| [`05_GWAS/Saige/06_locus_def/task1_get_leads.py`](05_GWAS/Saige/06_locus_def/task1_get_leads.py) | **Locus definition.** Distance-based clumping in gwaslab (`get_lead`, p < 5e-8, 500 kb window): a new locus begins where the gap to the previous significant variant exceeds 500 kb, and the lead is the lowest-p variant in the locus. Defines boundaries and leads only — *not* independence. Effect allele = `Allele2` (SAIGE's β/AF reference). Includes an FTO positive-control check. |

### 4. Conditional analysis — independent signals

| Script | Methods step |
|---|---|
| [`06_Conditional_analysis/02_run_locus_iterate.sh`](06_Conditional_analysis/02_run_locus_iterate.sh) | **Iterative (stepwise) conditional analysis** for one trait × locus. Extracts a ±1 Mb window from the HDS dosage PGEN, then repeatedly re-runs SAIGE Step 2 conditioning on the current signal set: any variant surviving at p_conditional < 5e-8 is a candidate, the strongest non-degenerate one is promoted into the condition set, and the loop repeats (max 10 rounds) until nothing survives. The final condition set = lead + promoted secondaries = the independent signals at that locus. Conditioning is run on PGEN because SAIGE's `--condition` is silently unsupported for BGEN input. |

### 5. Cross-ancestry comparison

| Script | Methods step |
|---|---|
| [`08_cross_ancestry/03_align.R`](08_cross_ancestry/03_align.R) | **Allele-aware alignment and replication testing** of every BELIEVE independent signal against each matched-trait comparator (GIANT BMI/WHRadjBMI/WC/hip, Biobank Japan BMI, Genes & Health BMI/WC). Handles per-comparator build and key differences (hg19 position, hg38 position, or rsID), flips comparator β and EAF when its effect allele is the BELIEVE non-effect allele, and flags palindromes. Reports direction concordance and replication at nominal, per-comparator Bonferroni, and genome-wide thresholds; emits SAS-specific candidates (genome-wide in BELIEVE but absent/non-replicating in European comparators). Includes an FTO alignment gate that must pass before results are used. |

### 6. Mendelian randomisation — adiposity → coronary disease

**South Asian arm (BELIEVE exposure → Genes & Health CHD outcome):**

| Script | Methods step |
|---|---|
| [`09_MR/Two_sample_MR/SAS/01_build_instruments.sh`](09_MR/Two_sample_MR/SAS/01_build_instruments.sh) | **Instrument selection.** Per trait: p < 5e-8, MAF ≥ 0.01, imputation INFO ≥ 0.7, extended MHC (chr6:25.7–33.4 Mb, GRCh38) excluded; then LD-clumped in PLINK2 against the **in-sample** BELIEVE reference at r² < 0.001 over 10 Mb. |
| [`09_MR/Two_sample_MR/SAS/03_mr_analysis.R`](09_MR/Two_sample_MR/SAS/03_mr_analysis.R) | **Harmonisation, MR estimation and diagnostics.** Explicit harmonisation on chr:pos (both datasets TOPMed/GRCh38, same strand) with outcome β and EAF flipped on allele swap, ambiguous palindromes (EAF 0.40–0.60) dropped, and allele-frequency-inconsistent SNPs (|ΔEAF| > 0.20) dropped. Estimators from `MendelianRandomization`: random-effects IVW (primary), weighted median, weighted mode, MR-Egger. Diagnostics: per-SNP F, Cochran's Q for IVW and Egger, Egger intercept, leave-one-out. Outcome β is log-odds, so OR = exp(β). Produces the four-panel scatter/forest/LOO/funnel figure per trait. |

**European comparator arm (Yengo BMI, Pan-UKB WHR & body fat % → Nikpay CAD):**

| Script | Methods step |
|---|---|
| [`09_MR/Two_sample_MR/EU/03_mr_eur.R`](09_MR/Two_sample_MR/EU/03_mr_eur.R) | Same estimator set and allele-aware logic as the South Asian arm, with the differences the data force: harmonisation on rsID rather than chr:pos, outcome allele frequencies merged from a 1000 Genomes EUR reference (Nikpay reports none), and a narrower palindrome exclusion window (0.42–0.58). Nikpay: 60,801 cases / 123,504 controls. |
| [`09_MR/Two_sample_MR/EU/04_power_and_transancestry.R`](09_MR/Two_sample_MR/EU/04_power_and_transancestry.R) | **Power and trans-ancestry heterogeneity.** Power by the Brion (2013) binary-outcome approximation, plus the minimum OR detectable at 80% power. For traits present in both ancestries, tests the South Asian vs European IVW estimates against each other with a one-degree-of-freedom Cochran's Q and reports the log-OR ratio and CI overlap — this is what distinguishes a genuine ancestry difference from a power difference. WHR (EUR) vs WHRadjBMI (SAS) is flagged as the closest available, non-identical match. |

### 7. Fine-mapping

| Script | Methods step |
|---|---|
| [`11_finemap/sweep_hds/prep_locus_hds.sh`](11_finemap/sweep_hds/prep_locus_hds.sh) | **FINEMAP input preparation** for one locus using **in-sample** BELIEVE LD. Takes variants within ±250 kb of the index with MAF ≥ 0.01, computes the signed r correlation matrix in PLINK with `--keep-allele-order`, and writes the z-file in exactly the LD matrix's variant order with effect allele = ALT = SAIGE `Allele2`. Row-order and allele-orientation agreement between the z-file and the LD matrix is the correctness condition for FINEMAP. |
| [`11_finemap/sweep_hds/sweep_parse_hds.py`](11_finemap/sweep_hds/sweep_parse_hds.py) | **Fine-mapping summary.** Iterates every (trait, lead) window in the manifest and extracts credible-set sizes, posterior inclusion probabilities, the lead variant's PIP (matched on exact variant ID, not position), and the log10 Bayes factor for ≥ one causal variant; flags leads with low INFO or absent from the LD reference. |


## Code/output split
**Every** script lives here; **every** output lives under `~/rds/hpc-work/believe_adiposity/`.
Scripts address their outputs by absolute path into the results tree, 
