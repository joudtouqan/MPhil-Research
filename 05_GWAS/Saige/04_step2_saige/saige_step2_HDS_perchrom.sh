
# ============================================================
# SAIGE Step 2 - single-variant tests on HDS dosages (BGEN)
# One array task per chromosome (1..22).
# Each task:
#   1. bgenix-extracts that chromosome into a small per-chrom BGEN
#   2. runs SAIGE step2 for all 6 phenotypes on that BGEN
#   3. deletes the per-chrom BGEN to stay within disk quota
# This avoids streaming the whole 307M-variant genome-wide BGEN for every chromosome
# ============================================================

set -euo pipefail

CHR=${SLURM_ARRAY_TASK_ID}

PHENOS=(bmi_INT whr_bmi_adj_INT waist_INT hip_INT fatmass_INT fatperc_INT)

BGENIX=/rds/user/jt962/hpc-work/micromamba_root/envs/bgenix_env/bin/bgenix
SAIGE_SIF=/home/jt962/rds/hpc-work/software/saige/saige_1.5.0.sif
BGEN=/home/jt962/rds/rds-post_qc_data-pNR2rM6BWWA/believe/genotype/imputed/CAMBRIDGE-BELIEVE_Freeze_Two.GxS.TOPMED_dosages.bgen
BGI=/home/jt962/rds/hpc-work/believe_adiposity/03_variant_QC/CAMBRIDGE-BELIEVE_Freeze_Two.GxS.TOPMED_dosages.bgen.bgi
SAMPLE=/home/jt962/rds/hpc-work/believe_adiposity/03_variant_QC/believe_freeze_two_saige_samples.txt
STEP1DIR=/home/jt962/rds/hpc-work/believe_adiposity/05_GWAS/Saige/03_step1_saige/output
BASE=/home/jt962/rds/hpc-work/believe_adiposity/05_GWAS/Saige/04_step2_saige/HDS
OUTDIR=${BASE}/output
LOGDIR=${BASE}/logs
TMPBGENDIR=${BASE}/perchrom_bgen

mkdir -p "$OUTDIR" "$LOGDIR" "$TMPBGENDIR"

CHRBGEN=${TMPBGENDIR}/chr${CHR}.bgen
cleanup() { rm -f "$CHRBGEN" "${CHRBGEN}.bgi"; }
trap cleanup EXIT

echo "=================================================="
echo "SAIGE Step 2 HDS per-chrom : chr${CHR}"
echo "Host: $(hostname)   Start: $(date)"
echo "=================================================="

# ---- 1. Extract this chromosome from the genome-wide BGEN ----
echo "[$(date)] Extracting chr${CHR} via bgenix ..."
"$BGENIX" -g "$BGEN" -i "$BGI" -incl-range "${CHR}:0-999999999" > "$CHRBGEN"
echo "[$(date)] Extraction done. Size: $(du -h "$CHRBGEN" | cut -f1)"

echo "[$(date)] Indexing per-chrom BGEN ..."
"$BGENIX" -g "$CHRBGEN" -index -clobber
echo "[$(date)] Index built."

module load singularity 2>/dev/null || true
unset R_LIBS R_LIBS_USER R_LIBS_SITE LD_PRELOAD

# ---- 2. Run SAIGE step2 for each phenotype ----
# Phenotypes are run 3-at-a-time in parallel (memory-light: sparse GRM).
# Each (pheno,chr) writes its own SAIGE log; a sentinel ".done" marks success
# so reruns skip completed work.
run_pheno() {
    local PHENO=$1
    local OUT=${OUTDIR}/${PHENO}_chr${CHR}.txt
    local DONE=${OUT}.done
    local PLOG=${LOGDIR}/saige_${PHENO}_chr${CHR}.log
    if [ -f "$DONE" ]; then echo "[$(date)] SKIP ${PHENO} chr${CHR} (already done)"; return 0; fi
    echo "[$(date)] START SAIGE step2: ${PHENO} chr${CHR}"
    singularity exec \
        --no-home \
        --bind /home/jt962/rds:/home/jt962/rds \
        --bind /home/jt962/MPhil-Research:/home/jt962/MPhil-Research \
        --env R_LIBS_USER= --env R_LIBS= --env R_LIBS_SITE= --env LD_PRELOAD= \
        "$SAIGE_SIF" \
        step2_SPAtests.R \
        --bgenFile="$CHRBGEN" \
        --bgenFileIndex="${CHRBGEN}.bgi" \
        --sampleFile="$SAMPLE" \
        --AlleleOrder=ref-first \
        --chrom="${CHR}" \
        --minMAF=0.005 \
        --minMAC=20 \
        --minInfo=0.3 \
        --LOCO=FALSE \
        --is_fastTest=TRUE \
        --is_imputed_data=TRUE \
        --GMMATmodelFile="${STEP1DIR}/${PHENO}_step1.rda" \
        --varianceRatioFile="${STEP1DIR}/${PHENO}_step1.varianceRatio.txt" \
        --SAIGEOutputFile="$OUT" > "$PLOG" 2>&1
    if grep -q "Analysis done" "$PLOG"; then
        touch "$DONE"
        echo "[$(date)] DONE ${PHENO} chr${CHR}: $(wc -l < "$OUT") lines"
    else
        echo "[$(date)] FAILED ${PHENO} chr${CHR} (see $PLOG)" >&2
        return 1
    fi
}
export -f run_pheno
export CHRBGEN SAMPLE CHR OUTDIR LOGDIR STEP1DIR SAIGE_SIF

# Launch all phenotypes, but cap concurrency at 3 simultaneous SAIGE runs.
RC=0
printf '%s\n' "${PHENOS[@]}" | xargs -I{} -P 3 bash -c 'run_pheno "$@"' _ {} || RC=$?
echo "[$(date)] All phenotype runs finished for chr${CHR} (xargs rc=$RC)"

# Report per-pheno completion for this chromosome
echo "[$(date)] Completion status chr${CHR}:"
for PHENO in "${PHENOS[@]}"; do
    if [ -f "${OUTDIR}/${PHENO}_chr${CHR}.txt.done" ]; then
        echo "  OK   ${PHENO}_chr${CHR}.txt ($(wc -l < "${OUTDIR}/${PHENO}_chr${CHR}.txt") lines)"
    else
        echo "  MISS ${PHENO}_chr${CHR}.txt"
    fi
done

echo "=================================================="
echo "[$(date)] chr${CHR} complete (cleanup on exit)."
sacct -j "$SLURM_JOB_ID" --format=JobID,JobName%18,MaxRSS,TotalCPU,Elapsed,State
echo "=================================================="
