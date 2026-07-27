
# ============================================================
# Rebuild the BELIEVE Freeze Two analysis-ready PGEN using HDS (dosage) input
# instead of GT (hard calls). Variant QC was already applied when building the
# GT pgen, so here I only subset the HDS pgen to the variants that passed.
#
# Output: freeze_two_common_qc_dosages.{pgen,pvar,psam} 
# ============================================================

set -euo pipefail
module load plink/2.00-alpha

# ── Paths ─────────────────────────────────────────────────────────────────────
export DATA=~/rds/rds-post_qc_data-pNR2rM6BWWA/believe/genotype/imputed
export WORK=~/rds/hpc-work/believe_adiposity/03_variant_QC
export PREFIX=CAMBRIDGE-BELIEVE_Freeze_Two.GxS.TOPMED_dosages

mkdir -p $WORK/pgen $WORK/logs

TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOG=$WORK/logs/${TIMESTAMP}_make_pgen_hds.log

echo "Run started:  $(date)"  | tee -a $LOG

# ── Build HDS-based filtered PGEN ─────────────────────────────────────────────
# Use the snplist from the GT QC run as the variant whitelist. All thresholds
# (--maf, --hwe, --geno, --max-alleles, etc.) were already enforced when that
# snplist was built, so they're not re-applied here.

plink2 \
  --pgen   $DATA/${PREFIX}.HDS.pgen \
  --pvar   $DATA/${PREFIX}.HDS.pvar \
  --psam   $DATA/${PREFIX}.HDS.COLLAB.psam \
  --extract $WORK/pgen/freeze_two_common_qc.snplist \
  --make-pgen \
  --out    $WORK/pgen/freeze_two_common_qc_dosages \
  --threads 8 \
  --memory 28000 \
  2>&1 | tee -a $LOG

echo "Run finished: $(date)"  | tee -a $LOG

# ── Sanity checks ────────────────────────────────────────────────────────────
echo "" | tee -a $LOG
echo "=== Output summary ===" | tee -a $LOG
echo "Variant count:" | tee -a $LOG
grep -cv '^#' $WORK/pgen/freeze_two_common_qc_dosages.pvar | tee -a $LOG
echo "Sample count:" | tee -a $LOG
grep -cv '^#' $WORK/pgen/freeze_two_common_qc_dosages.psam | tee -a $LOG
echo "Match to GT pgen (should be identical):" | tee -a $LOG
echo "  GT  variants: $(grep -cv '^#' $WORK/pgen/freeze_two_common_qc.pvar)" | tee -a $LOG
echo "  HDS variants: $(grep -cv '^#' $WORK/pgen/freeze_two_common_qc_dosages.pvar)" | tee -a $LOG
echo "  GT  samples : $(grep -cv '^#' $WORK/pgen/freeze_two_common_qc.psam)" | tee -a $LOG
echo "  HDS samples : $(grep -cv '^#' $WORK/pgen/freeze_two_common_qc_dosages.psam)" | tee -a $LOG
