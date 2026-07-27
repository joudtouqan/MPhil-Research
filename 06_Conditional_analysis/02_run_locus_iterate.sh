#!/bin/bash
#SBATCH --job-name=cond_iter
#SBATCH --output=/home/jt962/rds/hpc-work/believe_adiposity/06_Conditional_analysis/logs/iter_%x_%j.out
#SBATCH --error=/home/jt962/rds/hpc-work/believe_adiposity/06_Conditional_analysis/logs/iter_%x_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --partition=icelake-himem
#SBATCH --account=BUTTERWORTH-SL2-CPU
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jt962@cam.ac.uk

# ============================================================
# STEP 2 : iterative conditional analysis for ONE trait x locus.
# ============================================================

set -euo pipefail

TRAIT=$1
CHR=$2
LEAD_POS=$3
LEAD_ID=$4
LOCUS=${5:-locus${CHR}_${LEAD_POS}}

WIN_HALF=1000000
WIN_START=$(( LEAD_POS - WIN_HALF )); [ "$WIN_START" -lt 1 ] && WIN_START=1
WIN_END=$(( LEAD_POS + WIN_HALF ))

# ---- source HDS dosage PGEN trio ----
IMPDIR=/home/jt962/rds/rds-post_qc_data-pNR2rM6BWWA/believe/genotype/imputed
SRC_PGEN=${IMPDIR}/CAMBRIDGE-BELIEVE_Freeze_Two.GxS.TOPMED_dosages.HDS.pgen
SRC_PVAR=${IMPDIR}/CAMBRIDGE-BELIEVE_Freeze_Two.GxS.TOPMED_dosages.HDS.pvar
SRC_PSAM=${IMPDIR}/CAMBRIDGE-BELIEVE_Freeze_Two.GxS.TOPMED_dosages.HDS.COLLAB.psam

# ---- SAIGE / paths ----
SAIGE_SIF=/home/jt962/rds/hpc-work/software/saige/saige_1.5.0.sif
SAMPLE=/home/jt962/rds/hpc-work/believe_adiposity/03_variant_QC/believe_freeze_two_saige_samples.txt
STEP1DIR=/home/jt962/rds/hpc-work/believe_adiposity/05_GWAS/Saige/03_step1_saige/output
GWASDIR=/home/jt962/rds/hpc-work/believe_adiposity/05_GWAS/Saige/04_step2_saige/HDS/output
BASE=/home/jt962/rds/hpc-work/believe_adiposity/06_Conditional_analysis
OUTDIR=${BASE}/output/step2/${LOCUS}/${TRAIT}
LOGDIR=${BASE}/logs
PGENDIR=${BASE}/window_pgen
mkdir -p "$OUTDIR" "$LOGDIR" "$PGENDIR"

TAG=${LOCUS}_${TRAIT}_chr${CHR}_${LEAD_POS}
WINPREFIX=${PGENDIR}/${TAG}_pm1Mb

echo "=================================================="
echo "STEP2 iterate : ${LOCUS} / ${TRAIT}"
echo "Host: $(hostname)   Start: $(date)"
echo "Lead: ${LEAD_ID}   window ${CHR}:${WIN_START}-${WIN_END} (+/-1Mb)"
echo "Outdir: ${OUTDIR}"
echo "=================================================="

# ---- 1. Extract window from HDS PGEN ----
module load rhel8/default-icl 2>/dev/null || true
module load ceuadmin/plink/2.0_20240105 2>/dev/null || true
echo "[$(date)] plink2: $(command -v plink2)"
if [ ! -f "${WINPREFIX}.pgen" ]; then
    plink2 --pgen "$SRC_PGEN" --pvar "$SRC_PVAR" --psam "$SRC_PSAM" \
        --chr "$CHR" --from-bp "$WIN_START" --to-bp "$WIN_END" \
        --make-pgen --out "$WINPREFIX"
else
    echo "[$(date)] window pgen already present, reusing."
fi
NWIN=$(grep -vc '^#' "${WINPREFIX}.pvar")
echo "[$(date)] window variants: ${NWIN}"

# ---- 2. Confirm lead present (rebuild CHROM:POS:REF:ALT from window pvar) ----
LEADCHK=$(awk -F'\t' -v p="$LEAD_POS" '$2==p{print $1":"$2":"$4":"$5}' "${WINPREFIX}.pvar" | grep -Fx "$LEAD_ID" || true)
if [ -z "$LEADCHK" ]; then
    echo "[$(date)] ERROR: lead ${LEAD_ID} not found as CHROM:POS:REF:ALT in window pvar. Rows at pos:" >&2
    awk -F'\t' -v p="$LEAD_POS" '$2==p{print "  "$1":"$2":"$4":"$5}' "${WINPREFIX}.pvar" >&2
    exit 1
fi
echo "[$(date)] lead confirmed in window: ${LEAD_ID}"

module load singularity 2>/dev/null || true
unset R_LIBS R_LIBS_USER R_LIBS_SITE LD_PRELOAD
cd "$OUTDIR"

run_saige_cond() {   # $1=outfile  $2=comma-joined condition set
    local OUT=$1 COND=$2
    local PLOG=${OUT%.txt}.log
    singularity exec --no-home \
        --bind /home/jt962/rds:/home/jt962/rds \
        --bind /home/jt962/MPhil-Research:/home/jt962/MPhil-Research \
        --env R_LIBS_USER= --env R_LIBS= --env R_LIBS_SITE= --env LD_PRELOAD= \
        --env R_ENVIRON_USER=/dev/null --env R_PROFILE_USER=/dev/null \
        "$SAIGE_SIF" step2_SPAtests.R \
        --pgenPrefix="$WINPREFIX" \
        --sampleFile="$SAMPLE" \
        --AlleleOrder=ref-first \
        --chrom="${CHR}" \
        --minMAF=0.005 --minMAC=20 --minInfo=0.3 \
        --LOCO=FALSE --is_imputed_data=TRUE \
        --condition="$COND" \
        --GMMATmodelFile="${STEP1DIR}/${TRAIT}_step1.rda" \
        --varianceRatioFile="${STEP1DIR}/${TRAIT}_step1.varianceRatio.txt" \
        --SAIGEOutputFile="$OUT" > "$PLOG" 2>&1
    grep -q "Analysis done" "$PLOG"
}

# ---- 3. Iterate ----
COND_SET="$LEAD_ID"
SIGNALS_FILE=${OUTDIR}/${TAG}.signals.tsv
SURV_FILE=${OUTDIR}/${TAG}.survivors_byround.tsv
: > "$SURV_FILE"
echo -e "round\tcond_id\tpromoted_secondary\tsecondary_p.value_c\tn_surv_5e8\tn_suspect" > "$SIGNALS_FILE"

MAXROUNDS=10
round=0
while [ "$round" -lt "$MAXROUNDS" ]; do
    round=$(( round + 1 ))
    OUT=${OUTDIR}/${TAG}.cond_round${round}.txt
    echo "[$(date)] ROUND ${round}: conditioning on => ${COND_SET}"
    if ! run_saige_cond "$OUT" "$COND_SET"; then
        echo "[$(date)] SAIGE FAILED round ${round} (see ${OUT%.txt}.log)" >&2
        tail -20 "${OUT%.txt}.log" >&2
        exit 1
    fi

    # Find survivors: p.value_c (col18) < 5e-8, excluding conditioned markers.
    # canonical id = col1:col2:col4:col5 (CHROM:POS:REF:ALT).
    # Emit: id, p.value (col13), p.value_c (col18), AF (col7), suspect flag.
    awk -F'\t' -v cset=",${COND_SET}," -v rnd="$round" '
        NR==1 {next}
        {
            id=$1":"$2":"$4":"$5
            if (index(cset, ","id",") > 0) next       # skip already-conditioned
            pc=$18+0
            if (pc < 5e-8) {
                susp = (pc==0) ? "SUSPECT" : "ok"
                printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n", rnd, id, $13, $18, $7, $8, susp
            }
        }' "$OUT" >> "$SURV_FILE"

    # Among this round survivors (rnd==round): pick top non-suspect by p.value_c.
    PROMO=$(awk -F'\t' -v r="$round" '$1==r && $7=="ok"{print $2"\t"$4}' "$SURV_FILE" \
            | sort -t$'\t' -k2,2g | head -1)
    NSURV=$(awk -F'\t' -v r="$round" '$1==r{c++} END{print c+0}' "$SURV_FILE")
    NSUSP=$(awk -F'\t' -v r="$round" '$1==r && $7=="SUSPECT"{c++} END{print c+0}' "$SURV_FILE")

    if [ -n "$PROMO" ]; then
        PID=$(echo "$PROMO" | cut -f1); PPC=$(echo "$PROMO" | cut -f2)
        echo -e "${round}\t${COND_SET}\t${PID}\t${PPC}\t${NSURV}\t${NSUSP}" >> "$SIGNALS_FILE"
        echo "[$(date)] round ${round}: ${NSURV} survivors (<5e-8), ${NSUSP} suspect. Promote ${PID} (p.value_c=${PPC})"
        COND_SET="${COND_SET},${PID}"
    else
        echo -e "${round}\t${COND_SET}\t<none>\tNA\t${NSURV}\t${NSUSP}" >> "$SIGNALS_FILE"
        echo "[$(date)] round ${round}: ${NSURV} survivors (<5e-8), ${NSUSP} suspect, 0 promotable. STOP."
        break
    fi
done

echo "[$(date)] Iteration complete. Final condition set: ${COND_SET}"
echo "[$(date)] signals: ${SIGNALS_FILE}"
echo "[$(date)] survivors-by-round: ${SURV_FILE}"
echo "[$(date)] independent signals at this locus = lead + promoted secondaries:"
echo "    ${COND_SET}" | tr ',' '\n' | sed 's/^/      /'
echo "=================================================="
sacct -j "$SLURM_JOB_ID" --format=JobID,JobName%14,MaxRSS,Elapsed,State 2>/dev/null || true
echo "=================================================="
