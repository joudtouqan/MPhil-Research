

PHENOS=(bmi_INT whr_bmi_adj_INT waist_INT hip_INT fatmass_INT fatperc_INT)
PHENO_NAME=${PHENOS[$((SLURM_ARRAY_TASK_ID - 1))]}

SAIGE_SIF=/home/jt962/rds/hpc-work/software/saige/saige_1.5.0.sif
LDP=/home/jt962/rds/hpc-work/believe_adiposity/05_GWAS/Saige/02_GRM/bfile/ld_pruned_autosomes
PHENO=/home/jt962/rds/hpc-work/believe_adiposity/01_phenotype_prep/believe_phenotypes_saige.txt
GRM=/home/jt962/rds/hpc-work/believe_adiposity/05_GWAS/Saige/02_GRM/grm/sparse_grm_relatednessCutoff_0.05_2000_randomMarkersUsed.sparseGRM.mtx
GRM_IDS=${GRM}.sampleIDs.txt
OUTDIR=/home/jt962/rds/hpc-work/believe_adiposity/05_GWAS/Saige/03_step1_saige/output

mkdir -p "$OUTDIR" logs

echo "=== SAIGE Step 1: $PHENO_NAME ==="
echo "Start: $(date)"

export APPTAINERENV_R_PROFILE_USER=/dev/null
export APPTAINERENV_R_ENVIRON_USER=/dev/null

module load singularity 2>/dev/null || true
unset R_LIBS R_LIBS_USER R_LIBS_SITE LD_PRELOAD

singularity exec \
    --no-home \
    --bind /home/jt962/rds:/home/jt962/rds \
    --env R_LIBS_USER= \
    --env R_LIBS= \
    --env R_LIBS_SITE= \
    --env LD_PRELOAD= \
    "$SAIGE_SIF" \
    step1_fitNULLGLMM.R \
    --plinkFile="$LDP" \
    --sparseGRMFile="$GRM" \
    --sparseGRMSampleIDFile="$GRM_IDS" \
    --useSparseGRMtoFitNULL=TRUE \
    --phenoFile="$PHENO" \
    --phenoCol="$PHENO_NAME" \
    --sampleIDColinphenoFile=IID \
    --traitType=quantitative \
    --invNormalize=FALSE \
    --outputPrefix="${OUTDIR}/${PHENO_NAME}_step1" \
    --nThreads=8 \
    --IsOverwriteVarianceRatioFile=TRUE

echo "End: $(date)"
