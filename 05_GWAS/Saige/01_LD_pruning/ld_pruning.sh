#!/bin/bash
#SBATCH --job-name=ld_prune
#SBATCH --time=4:00:00
#SBATCH --partition=icelake
#SBATCH --account=BUTTERWORTH-SL2-CPU
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jt962@cam.ac.uk
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=/home/jt962/rds/hpc-work/believe_adiposity/05_GWAS/Saige/01_LD_pruning/logs/ld_prune_%j.log
#SBATCH --error=/home/jt962/rds/hpc-work/believe_adiposity/05_GWAS/Saige/01_LD_pruning/logs/ld_prune_%j.err

set -euo pipefail
module load plink/2.00-alpha

# ── Paths ─────────────────────────────────────────────────────────────────────
DATA=/home/jt962/rds/rds-post_qc_data-pNR2rM6BWWA/believe/genotype/genomewide/plink/aug_2023
PREFIX=CAMBRIDGE-BELIEVE_Freeze_Two.GxS
WORK=/home/jt962/rds/hpc-work/believe_adiposity/05_GWAS/Saige/01_LD_pruning

# ── Build the long-range LD exclusion file (hg38) ─────────────────────────────

cat > $WORK/hg38_high_LD.txt <<'EOF'
6  28510120  33480577  MHC
8  7000000   13000000  chr8_inversion
17 40000000  47000000  chr17_inversion
EOF

# ── Step 1: identify unrelated participants (Henry: KING coef < 0.022) ────────
plink2 \
  --bfile $DATA/$PREFIX \
  --king-cutoff 0.022 \
  --threads 8 \
  --out $WORK/unrelated

# ── Step 2: filter variants then prune ─────────────────
#   MAF > 5%, HWE > 1e-10, missingness < 5%, biallelic SNVs only,
#   high-LD regions removed, unrelated subset only.
#   Pruning parameters: window=1000, step=80, r² threshold=0.1.
plink2 \
  --bfile $DATA/$PREFIX \
  --keep $WORK/unrelated.king.cutoff.in.id \
  --maf 0.05 \
  --hwe 1e-10 \
  --geno 0.05 \
  --snps-only just-acgt \
  --max-alleles 2 --min-alleles 2 \
  --exclude range $WORK/hg38_high_LD.txt \
  --indep-pairwise 1000 80 0.1 \
  --threads 8 \
  --out $WORK/ld_pruned

echo "LD-pruned variants retained: $(wc -l < $WORK/ld_pruned.prune.in)"

# ── Step 3: extract the pruned variants back into the FULL cohort ────────────
# do the pruning on unrelateds, then apply the variant list to the full cohort for sparse GRM construction.
plink2 \
  --bfile $DATA/$PREFIX \
  --extract $WORK/ld_pruned.prune.in \
  --make-pgen \
  --threads 8 \
  --out $WORK/ld_pruned_full_cohort